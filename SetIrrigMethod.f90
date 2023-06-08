subroutine SetIrrigMethod(this, bounds)
    !
    ! !DESCRIPTION:
    ! Set this%irrig_method_patch based on surface dataset
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p ! patch index
    integer :: g ! gridcell index
    integer :: m ! patch itype

    character(len=*), parameter :: subname = 'SetIrrigMethod'
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       g = patch%gridcell(p)
       m = patch%itype(p)
       if (m >= lbound(irrig_method, 2) .and. m <= ubound(irrig_method, 2) &
            .and. pftcon%irrigated(m) == 1._r8) then
          this%irrig_method_patch(p) = irrig_method(g,m)
          ! ensure irrig_method is valid; if not set, use drip irrigation
          if(irrig_method(g,m) == irrig_method_unset) then
             this%irrig_method_patch(p) = this%params%irrig_method_default
          else if (irrig_method(g,m) /= irrig_method_drip .and. irrig_method(g,m) /= irrig_method_sprinkler) then
             write(iulog,*) subname //' invalid irrigation method specified'
             call endrun(subgrid_index=g, subgrid_level=subgrid_level_gridcell, msg='bad irrig_method '// &
                  errMsg(sourcefile, __LINE__))
          end if
       else
          this%irrig_method_patch(p) = this%params%irrig_method_default
       end if
    end do

  end subroutine SetIrrigMethod
