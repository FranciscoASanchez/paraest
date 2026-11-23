!--------------------------------------------------------------------------------------------------------------------------
!
!  Module model_output: writes the attractive and associative parameters of
!  the model (initial and final), in the layout of PARGR of the GCA package
!  (v1.14), with groups identified by their database number and name. The
!  precision is higher than PARGR's, so that the final values can be copied
!  to an input file.
!
!--------------------------------------------------------------------------------------------------------------------------
module model_output
  use gca_parameters, only: NG, NGA, NU, MASSOC, ONUASOC, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, &
                            ENASS, RKASS, groupID, assocID, groupName, assocName
  implicit none
  private
  public :: write_model_params

contains

subroutine write_model_params(u, NC)
  integer, intent(in) :: u, NC
  integer :: i, j, k, l, NST
  double precision :: Tspl, epx
  character(300) :: f

  ! ---- Pure groups
  write(u,'(/,X,A,/,106X,A)') 'Pure group attractive parameters:', 'Groups per component'
  write(f,'(A,I0,A)') '(3X,"No",3X,"Name",8X,"T*(K)",6X,"q",7X,"g*(atm cm6/mol2)",7X,"g''",11X,"g""",10X,"Ts(K)",5X,"exp",', &
                      NC, '(I5))'
  write(u,f) (j, j = 1, NC)
  write(f,'(A,I0,A)') '(X,107("-"),', NC, '("-----"))'
  write(u,f)
  do i = 1, NG
    call spline_temperature(TSTR(i), GDOT(i), GDDOT(i), Tspl, epx)
    write(u,'(X,I4,2X,A10,F8.2,F9.4,3X,F14.2,2X,2F13.6,F12.2,F7.3,*(I5))') groupID(i), groupName(i), TSTR(i), &
      QAREA(i), GSTR(i), GDOT(i), GDDOT(i), Tspl, epx, (NU(i,j), j = 1, NC)
  enddo

  ! ---- Binary (lower triangle; k* and k' are symmetric)
  if (NG > 1) then
    write(u,'(/,X,A,/)') 'Binary interaction parameter at T* matrix, [k*]'
    call matrix_header()
    do i = 1, NG
      write(u,'(X,I4,2X,A10,"|",*(F10.6))') groupID(i), groupName(i), RKIJ(i,:i)
    enddo
    write(u,'(/,X,A,/)') "Binary interaction parameter temperature dependence matrix, [k']"
    call matrix_header()
    do i = 1, NG
      write(u,'(X,I4,2X,A10,"|",*(F10.6))') groupID(i), groupName(i), RKTIJ(i,:i)
    enddo
    write(u,'(/,X,A,/)') 'Non-randomness parameters matrix, [alpha] (row i, column j: alpha_ij)'
    call matrix_header()
    do i = 1, NG
      write(f,'(A)') '(X,I4,2X,A10,"|"'
      do j = 1, NG
        if (j == i) then
          f = trim(f)//',"   -------"'
        else
          f = trim(f)//',F10.4'
        endif
      enddo
      f = trim(f)//')'
      write(u,f) groupID(i), groupName(i), (ALFAIJ(i,j), j = 1, i - 1), (ALFAIJ(i,j), j = i + 1, NG)
    enddo
  endif

  ! ---- Association
  write(u,'(/,X,A,I3)') 'Number of associating groups =', NGA
  if (NGA == 0) return
  write(u,'(/,X,A,/,32X,A)') 'Associating group configuration', 'Groups per component'
  write(u,'(3X,"No",3X,"Name",7X,"Sites",4X,*(I5))') (j, j = 1, NC)
  write(f,'(A,I0,A)') '(X,30("-"),', NC, '("-----"))'
  write(u,f)
  do i = 1, NGA
    write(u,'(X,I4,2X,A10,I4,6X,*(I5))') assocID(i), assocName(i), MASSOC(i), (nint(ONUASOC(i,j)), j = 1, NC)
  enddo
  NST = sum(MASSOC(:NGA))
  write(u,'(/,X,A,/)') 'Association energy matrix, [epsilon/R] (K)'
  call site_header()
  do i = 1, NGA
    do k = 1, MASSOC(i)
      write(u,'(X,I4,2X,A10,I2," |",*(F12.4))') assocID(i), assocName(i), k, &
        ((ENASS(k,i,l,j), l = 1, MASSOC(j)), j = 1, i - 1), (ENASS(k,i,l,i), l = 1, k)
    enddo
  enddo
  write(u,'(/,X,A,/)') 'Association volume matrix, [kappa] (cm3/mol)'
  call site_header()
  do i = 1, NGA
    do k = 1, MASSOC(i)
      write(u,'(X,I4,2X,A10,I2," |",*(ES12.5))') assocID(i), assocName(i), k, &
        ((RKASS(k,i,l,j), l = 1, MASSOC(j)), j = 1, i - 1), (RKASS(k,i,l,i), l = 1, k)
    enddo
  enddo

contains

  subroutine matrix_header()
    write(u,'(18X,*(I10))') groupID(:NG)
    write(u,'(18X,*(A10))') (adjustr(groupName(j)), j = 1, NG)
  endsubroutine matrix_header

  subroutine site_header()
    integer :: jj, kk
    write(u,'(21X,*(I12))') ((assocID(jj), kk = 1, MASSOC(jj)), jj = 1, NGA)
    write(u,'(21X,*(A12))') ((adjustr(assocName(jj)), kk = 1, MASSOC(jj)), jj = 1, NGA)
    write(u,'(21X,*(I12))') ((kk, kk = 1, MASSOC(jj)), jj = 1, NGA)
  endsubroutine site_header

endsubroutine write_model_params

endmodule model_output
