  subroutine CalcIrrigationNeeded(this, bounds, num_exposedvegp, filter_exposedvegp, &

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
    do fc = 1, check_for_irrig_col_filter%num
       c = check_for_irrig_col_filter%indices(fc)

       h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot(c) + &
            this%params%irrig_threshold_fraction * &
            (h2osoi_liq_target_tot(c) - h2osoi_liq_wilting_point_tot(c))
       if (h2osoi_liq_tot(c) < h2osoi_liq_at_threshold) then
          deficit(c) = h2osoi_liq_target_tot(c) - h2osoi_liq_tot(c)
          ! deficit shouldn't be less than 0: if it is, that implies that the
          ! irrigation target is less than the irrigation threshold, which is not
          ! supposed to happen
          if (deficit(c) < 0._r8) then
             write(iulog,*) subname//' ERROR: deficit < 0'
             write(iulog,*) 'This implies that irrigation target is less than irrigatio&
                  &n threshold, which should never happen'
             call endrun(subgrid_index=c, subgrid_level=subgrid_level_column, msg='deficit < 0 '// &
                  errMsg(sourcefile, __LINE__))
          end if
       else
          ! We're above the threshold - so don't irrigate
          deficit(c) = 0._r8
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

    ! Convert deficits to irrigation rate
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       c = patch%column(p)

       if (check_for_irrig_patch(p)) then

          ! Convert units from mm to mm/sec
          this%sfc_irrig_rate_patch(p) = deficit_volr_limited(c) / &
               (this%dtime*this%irrig_nsteps_per_day)
          this%irrig_rate_demand_patch(p) = deficit(c) / &
               (this%dtime*this%irrig_nsteps_per_day)

          ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
          ! in this case, we'll irrigate by 0 for the given number of time steps
          this%n_irrig_steps_left_patch(p) = this%irrig_nsteps_per_day
       end if
    end do

  end subroutine CalcIrrigationNeeded
