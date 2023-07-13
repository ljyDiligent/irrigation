 subroutine CalcIrrigationNeeded(this, bounds, num_exposedvegp, filter_exposedvegp, &
       elai, t_soisno, eff_porosity, h2osoi_liq, volr, rof_prognostic)
    !
    ! !DESCRIPTION:
    ! Calculate whether and how much irrigation is needed for each column. However, this
    ! does NOT actually set the irrigation flux.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds

    ! number of points in filter_exposedvegp
    integer, intent(in) :: num_exposedvegp

    ! patch filter for non-snow-covered veg
    integer, intent(in) :: filter_exposedvegp(:)

    ! one-sided leaf area index with burying by snow [patch]
    real(r8), intent(in) :: elai( bounds%begp: )

    ! col soil temperature (K) [col, nlevgrnd] (note that this does NOT contain the snow levels)
    real(r8), intent(in) :: t_soisno( bounds%begc: , 1: )

    ! effective porosity (0 to 1) [col, nlevgrnd]
    real(r8), intent(in) :: eff_porosity( bounds%begc: , 1: )

    ! column liquid water (kg/m2) [col, nlevgrnd] (note that this does NOT contain the snow levels)
    real(r8), intent(in) :: h2osoi_liq( bounds%begc: , 1: )
    ! river water volume (m3) (ignored if rof_prognostic is .false.)
    real(r8), intent(in) :: volr( bounds%begg: )

    ! whether we're running with a prognostic ROF component; this is needed to determine
    ! whether we can limit irrigation based on river volume.
    logical, intent(in) :: rof_prognostic

    !
    ! !LOCAL VARIABLES:
    integer :: fp   ! patch filter index
    integer :: fc   ! column filter index
    integer :: p    ! patch index
    integer :: c    ! column index
    integer :: g    ! gridcell index
    integer :: j    ! level

    ! Filter for columns where we need to check for irrigation
    type(filter_col_type) :: check_for_irrig_col_filter

    ! target soil moisture for this layer [kg/m2]
    real(r8) :: h2osoi_liq_target
	real(r8) :: h2osoi_liq_target_satu

    ! soil moisture at wilting point for this layer [kg/m2]
    real(r8) :: h2osoi_liq_wilting_point

    ! Total of h2osoi down to the depth of irrigation in each column [kg/m2]
    real(r8) :: h2osoi_liq_tot(bounds%begc:bounds%endc)

    ! Total of h2osoi_liq_target down to the depth of irrigation in each column [kg/m2]
    real(r8) :: h2osoi_liq_target_tot(bounds%begc:bounds%endc)
	real(r8) :: h2osoi_liq_target_satu_tot(bounds%begc:bounds%endc)

    ! Total of h2osoi_liq at wilting point down to the depth of irrigation in each column
    ! [kg/m2]
    real(r8) :: h2osoi_liq_wilting_point_tot(bounds%begc:bounds%endc)

    ! h2osoi_liq at the threshold for irrigation in this column [kg/m2]
    real(r8) :: h2osoi_liq_at_threshold

    ! difference between desired soil moisture level for each column and current soil
    ! moisture level [kg/m2] [i.e., mm]
    real(r8) :: deficit(bounds%begc:bounds%endc)
	real(r8) :: deficit_satu(bounds%begc:bounds%endc)
    	real(r8) :: deficit_pool(bounds%begc:bounds%endc)


    ! deficit limited by river volume [kg/m2] [i.e., mm]
    real(r8) :: deficit_volr_limited(bounds%begc:bounds%endc)
    real(r8) :: deficit_satu_volr_limited(bounds%begc:bounds%endc)
    real(r8) :: deficit_pool_volr_limited(bounds%begc:bounds%endc)

    ! where do we need to check soil moisture to see if we need to irrigate?
    logical  :: check_for_irrig_patch(bounds%begp:bounds%endp)
    logical  :: check_for_irrig_col(bounds%begc:bounds%endc)

    ! set to true once we have reached the max allowable depth for irrigation in a given
    ! column
    logical  :: reached_max_depth(bounds%begc:bounds%endc)

    ! Whether we should limit deficits by available volr
    logical :: limit_irrigation

    character(len=*), parameter :: subname = 'CalcIrrigationNeeded'
    !-----------------------------------------------------------------------
    
    ! Enforce expected array sizes
    SHR_ASSERT_ALL_FL((ubound(elai) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(t_soisno) == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(eff_porosity) == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(h2osoi_liq) == (/bounds%endc, nlevgrnd/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(volr) == (/bounds%endg/)), sourcefile, __LINE__)


    ! Determine if irrigation is needed (over irrigated soil columns)
    
    ! First, determine in what grid cells we need to bother 'measuring' soil water, to see
    ! if we need irrigation
    check_for_irrig_col(bounds%begc:bounds%endc) = .false.
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       g = patch%gridcell(p)

       check_for_irrig_patch(p) = this%PointNeedsCheckForIrrig( &
            pft_type=patch%itype(p), elai=elai(p), &
            londeg=grc%londeg(g))
       if (check_for_irrig_patch(p)) then
          c = patch%column(p)
          check_for_irrig_col(c) = .true.
       end if
    end do

    check_for_irrig_col_filter = col_filter_from_logical_array(bounds, &
         check_for_irrig_col(bounds%begc:bounds%endc))

    ! Initialize some variables
    do fc = 1, check_for_irrig_col_filter%num
       c = check_for_irrig_col_filter%indices(fc)

       reached_max_depth(c) = .false.
       h2osoi_liq_tot(c) = 0._r8
       h2osoi_liq_target_tot(c) = 0._r8
       h2osoi_liq_target_satu_tot(c) = 0._r8
       h2osoi_liq_wilting_point_tot(c) = 0._r8
    end do

    ! Now 'measure' soil water for the grid cells identified above and see if the soil is
    ! dry enough to warrant irrigation
    do j = 1,nlevsoi
       do fc = 1, check_for_irrig_col_filter%num
          c = check_for_irrig_col_filter%indices(fc)

          if (.not. reached_max_depth(c)) then
             if (col%z(c,j) > this%params%irrig_depth) then
                reached_max_depth(c) = .true.
             else if (j > col%nbedrock(c)) then
                reached_max_depth(c) = .true.
             else if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                ! if level L was frozen, then we don't look at any levels below L
                reached_max_depth(c) = .true.
             else
                h2osoi_liq_tot(c) = h2osoi_liq_tot(c) + h2osoi_liq(c,j)

                h2osoi_liq_target = this%RelsatToH2osoi( &
                     relsat = this%relsat_target_col(c,j), &
                     eff_porosity = eff_porosity(c,j), &
                     dz = col%dz(c,j))
                h2osoi_liq_target_tot(c) = h2osoi_liq_target_tot(c) + &
                     h2osoi_liq_target 
                h2osoi_liq_target_satu = this%RelsatToH2osoi( &
                     relsat = 1._r8, &
                     eff_porosity = eff_porosity(c,j), &
                     dz = col%dz(c,j))
                h2osoi_liq_target_satu_tot(c) = h2osoi_liq_target_satu_tot(c) + &
                     h2osoi_liq_target_satu

                h2osoi_liq_wilting_point = this%RelsatToH2osoi( &
                     relsat = this%relsat_wilting_point_col(c,j), &
                     eff_porosity = eff_porosity(c,j), &
                     dz = col%dz(c,j))
                h2osoi_liq_wilting_point_tot(c) = h2osoi_liq_wilting_point_tot(c) + &
                     h2osoi_liq_wilting_point
             end if
          end if     ! if (.not. reached_max_depth(c))
       end do        ! do fc
    end do           ! do j

    ! Compute deficits
    ! First initialize deficits to 0 everywhere; this is needed for later averaging up to gridcell
    deficit(bounds%begc:bounds%endc) = 0._r8
    deficit_satu(bounds%begc:bounds%endc) = 0._r8
    deficit_pool(bounds%begc:bounds%endc) = 0._r8

    do fc = 1, check_for_irrig_col_filter%num
       c = check_for_irrig_col_filter%indices(fc)
       h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot(c) + &
            this%params%irrig_threshold_fraction * &
            (h2osoi_liq_target_tot(c) - h2osoi_liq_wilting_point_tot(c))
       if (h2osoi_liq_tot(c) < h2osoi_liq_at_threshold) then 
         deficit(c) = (h2osoi_liq_target_tot(c) - h2osoi_liq_tot(c))
         deficit_satu(c) = h2osoi_liq_target_satu_tot(c) - h2osoi_liq_tot(c)   
         deficit_pool(c) = h2osoi_liq_target_satu_tot(c) - h2osoi_liq_tot(c) 
         if (deficit(c) < 0._r8) then
             write(iulog,*) subname//' ERROR: deficit < 0'
             write(iulog,*) 'This implies that irrigation target is less than irrigatio&
                  &n threshold, which should never happen'
             call endrun(decomp_index=c, clmlevel=namec, msg='deficit < 0 '// &
                  errMsg(sourcefile, __LINE__))
         end if
      else if (h2osoi_liq_tot(c) < h2osoi_liq_target_satu_tot(c)) then 
         deficit_pool(c) = h2osoi_liq_target_satu_tot(c) - h2osoi_liq_tot(c) 

      else
         ! We're above the threshold - so don't irrigate
         deficit(c) = 0._r8
         deficit_satu(c) = 0._r8
         deficit_pool(c) = 0._r8
      end if
    end do

    ! Limit deficits by available volr, if desired. Note that we cannot do this limiting
    ! if running without a prognostic river model, since we need river volume for this
    ! limiting.
    !
    ! NOTE(wjs, 2016-11-22) In principle we could base this on rof_present rather than
    ! rof_prognostic, but that would depend on the data runoff (drof) model sending river
    ! volume, which it currently does not.
    limit_irrigation = (this%params%limit_irrigation_if_rof_enabled .and. rof_prognostic)
    if (limit_irrigation) then
       call this%CalcDeficitVolrLimited( &
            bounds = bounds, &
            check_for_irrig_col_filter = check_for_irrig_col_filter, &
            deficit = deficit(bounds%begc:bounds%endc), &
            volr = volr(bounds%begg:bounds%endg), &
            deficit_volr_limited = deficit_volr_limited(bounds%begc:bounds%endc))
    else
       deficit_volr_limited(bounds%begc:bounds%endc) = deficit(bounds%begc:bounds%endc)
    end if

    if (limit_irrigation) then
       call this%CalcDeficitVolrLimited( &
            bounds = bounds, &
            check_for_irrig_col_filter = check_for_irrig_col_filter, &
            deficit = deficit_satu(bounds%begc:bounds%endc), &
            volr = volr(bounds%begg:bounds%endg), &
            deficit_volr_limited = deficit_satu_volr_limited(bounds%begc:bounds%endc))
    else
       deficit_satu_volr_limited(bounds%begc:bounds%endc) = deficit_satu(bounds%begc:bounds%endc)
    end if
    if (limit_irrigation) then
       call this%CalcDeficitVolrLimited( &
            bounds = bounds, &
            check_for_irrig_col_filter = check_for_irrig_col_filter, &
            deficit = deficit_pool(bounds%begc:bounds%endc), &
            volr = volr(bounds%begg:bounds%endg), &
            deficit_volr_limited = deficit_pool_volr_limited(bounds%begc:bounds%endc))
    else
       deficit_pool_volr_limited(bounds%begc:bounds%endc) = deficit_pool(bounds%begc:bounds%endc)
    end if

    ! Convert deficits to irrigation rate
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       c = patch%column(p)

       if (check_for_irrig_patch(p)) then

       ! Convert units from mm to mm/sec
       if (col%itype(c)==(200+cft_lb+47)) then
         this%sfc_irrig_rate_patch(p) = deficit_pool_volr_limited(c) / &
             (this%dtime*this%irrig_nsteps_per_day)
         this%irrig_rate_demand_patch(p) = deficit_pool(c) / &
                 (this%dtime*this%irrig_nsteps_per_day)

       else if(this%irrig_method_patch(p) == irrig_method_drip  .or.  this%irrig_method_patch(p) == irrig_method_sprinkler) then
         this%sfc_irrig_rate_patch(p) = deficit_volr_limited(c) / &
             (this%dtime*this%irrig_nsteps_per_day)
         this%irrig_rate_demand_patch(p) = deficit(c) / &
                 (this%dtime*this%irrig_nsteps_per_day)
       else if(this%irrig_method_patch(p) == irrig_method_flood) then
         this%sfc_irrig_rate_patch(p) = deficit_satu_volr_limited(c) / &
             (this%dtime*this%irrig_nsteps_per_day)
         this%irrig_rate_demand_patch(p) = deficit_satu(c) / &
                 (this%dtime*this%irrig_nsteps_per_day)
       else
         call endrun(msg=' ERROR: irrig_method_patch set to invalid value ' // &
            errMsg(sourcefile, __LINE__))
        endif
          ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
          ! in this case, we'll irrigate by 0 for the given number of time steps
          this%n_irrig_steps_left_patch(p) = this%irrig_nsteps_per_day
       end if
    end do

  end subroutine CalcIrrigationNeeded
