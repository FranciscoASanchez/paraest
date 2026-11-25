!--------------------------------------------------------------------------------------------------------------------------
!
!  Module model_output: writes the attractive and associative parameters of
!  the model (initial and final), in the layout of PARGR of the GCA package
!  (v1.14), with groups identified by their position in the system (Pos, used in
!  the fitted-parameter specification), database number and name. The
!  precision is higher than PARGR's, so that the final values can be copied
!  to an input file.
!
!--------------------------------------------------------------------------------------------------------------------------
module model_output
  use gca_parameters, only: NG, NGA, NU, MASSOC, ONUASOC, QAREA, TSTR, GSTR, GDOT, GDDOT, RKIJ, RKTIJ, ALFAIJ, &
                            ENASS, RKASS, groupID, assocID, groupName, assocName
  implicit none
  private
  public :: write_model_params, keep_pargr_warnings, write_pargr_warnings

  ! Warnings that PARGR wrote (binary parameters not in the database)
  character(200), allocatable :: pargrWarning(:)

contains

subroutine write_model_params(u, NC)
  integer, intent(in) :: u, NC
  integer :: i, j, k, l, NST
  double precision :: Tspl, epx
  character(300) :: f

  ! ---- Pure groups
  write(u,'(/,X,A,/,111X,A)') 'Pure group attractive parameters:', 'Groups per component'
  write(f,'(A,I0,A)') '(2X,"Pos",3X,"ID",2X,"Name",8X,"T*(K)",6X,"q",7X,"g*(atm cm6/mol2)",7X,"g''",11X,"g""",10X,"Ts(K)",5X,"exp",', &
                      NC, '(I5))'
  write(u,f) (j, j = 1, NC)
  write(f,'(A,I0,A)') '(X,112("-"),', NC, '("-----"))'
  write(u,f)
  do i = 1, NG
    call spline_temperature(TSTR(i), GDOT(i), GDDOT(i), Tspl, epx)
    write(u,'(X,I4,I5,2X,A10,F8.2,F9.4,3X,F14.2,2X,2F13.6,F12.2,F7.3,*(I5))') i, groupID(i), groupName(i), TSTR(i), &
      QAREA(i), GSTR(i), GDOT(i), GDDOT(i), Tspl, epx, (NU(i,j), j = 1, NC)
  enddo

  ! ---- Binary (lower triangle; k* and k' are symmetric)
  if (NG > 1) then
    write(u,'(/,X,A,/)') 'Binary interaction parameter at T* matrix, [k*]'
    call matrix_header()
    do i = 1, NG
      write(u,'(X,I4,I5,2X,A10,"|",*(F10.6))') i, groupID(i), groupName(i), RKIJ(i,:i)
    enddo
    write(u,'(/,X,A,/)') "Binary interaction parameter temperature dependence matrix, [k']"
    call matrix_header()
    do i = 1, NG
      write(u,'(X,I4,I5,2X,A10,"|",*(F10.6))') i, groupID(i), groupName(i), RKTIJ(i,:i)
    enddo
    write(u,'(/,X,A,/)') 'Non-randomness parameters matrix, [alpha] (row i, column j: alpha_ij)'
    call matrix_header()
    do i = 1, NG
      write(f,'(A)') '(X,I4,I5,2X,A10,"|"'
      do j = 1, NG
        if (j == i) then
          f = trim(f)//',"   -------"'
        else
          f = trim(f)//',F10.4'
        endif
      enddo
      f = trim(f)//')'
      write(u,f) i, groupID(i), groupName(i), (ALFAIJ(i,j), j = 1, i - 1), (ALFAIJ(i,j), j = i + 1, NG)
    enddo
  endif

  ! ---- Association
  write(u,'(/,X,A,I3)') 'Number of associating groups =', NGA
  if (NGA == 0) return
  write(u,'(/,X,A,/,37X,A)') 'Associating group configuration', 'Groups per component'
  write(u,'(2X,"Pos",3X,"ID",2X,"Name",7X,"Sites",4X,*(I5))') (j, j = 1, NC)
  write(f,'(A,I0,A)') '(X,35("-"),', NC, '("-----"))'
  write(u,f)
  do i = 1, NGA
    write(u,'(X,I4,I5,2X,A10,I4,6X,*(I5))') i, assocID(i), assocName(i), MASSOC(i), (nint(ONUASOC(i,j)), j = 1, NC)
  enddo
  NST = sum(MASSOC(:NGA))
  write(u,'(/,X,A,/)') 'Association energy matrix, [epsilon/R] (K)'
  call site_header()
  do i = 1, NGA
    do k = 1, MASSOC(i)
      write(u,'(X,I4,I5,2X,A10,I2," |",*(F12.4))') i, assocID(i), assocName(i), k, &
        ((ENASS(k,i,l,j), l = 1, MASSOC(j)), j = 1, i - 1), (ENASS(k,i,l,i), l = 1, k)
    enddo
  enddo
  write(u,'(/,X,A,/)') 'Association volume matrix, [kappa] (cm3/mol)'
  call site_header()
  do i = 1, NGA
    do k = 1, MASSOC(i)
      write(u,'(X,I4,I5,2X,A10,I2," |",*(ES12.5))') i, assocID(i), assocName(i), k, &
        ((RKASS(k,i,l,j), l = 1, MASSOC(j)), j = 1, i - 1), (RKASS(k,i,l,i), l = 1, k)
    enddo
  enddo

contains

  subroutine matrix_header()
    integer :: jj
    write(u,'(20X,"Pos",*(I10))') (jj, jj = 1, NG)
    write(u,'(20X,"ID ",*(I10))') groupID(:NG)
    write(u,'(23X,*(A10))') (adjustr(groupName(j)), j = 1, NG)
  endsubroutine matrix_header

  subroutine site_header()
    integer :: jj, kk
    write(u,'(23X,"Pos",*(I12))') ((jj, kk = 1, MASSOC(jj)), jj = 1, NGA)
    write(u,'(23X,"ID ",*(I12))') ((assocID(jj), kk = 1, MASSOC(jj)), jj = 1, NGA)
    write(u,'(26X,*(A12))') ((adjustr(assocName(jj)), kk = 1, MASSOC(jj)), jj = 1, NGA)
    write(u,'(22X,"Site",*(I12))') ((kk, kk = 1, MASSOC(jj)), jj = 1, NGA)
  endsubroutine site_header

endsubroutine write_model_params

!------------------------------------------------------------------------------------------------------------------------------------
!
!  Keeps the warnings that PARGR wrote to unit u (binary parameters not found
!  in the database: the two lines of the notice and the list of pairs of
!  groups) and closes u, so
!  that the scratch file does not stay behind if the program stops later.
!
subroutine keep_pargr_warnings(u)
  integer, intent(in) :: u
  character(200) :: line
  integer :: ios, n, pass
  do pass = 1, 2
    rewind(u)
    n = 0
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'not available in the database') > 0 .or. index(line, 'Standard (predictive)') > 0 &
          .or. index(line, ')  and ') > 0) then
        n = n + 1
        if (pass == 2) pargrWarning(n) = line
      endif
    enddo
    if (pass == 1) then
      if (allocated(pargrWarning)) deallocate(pargrWarning)
      allocate(pargrWarning(n))
    endif
  enddo
  close(u)
endsubroutine keep_pargr_warnings

subroutine write_pargr_warnings(u)
  integer, intent(in) :: u
  integer :: i
  if (.not. allocated(pargrWarning)) return
  do i = 1, size(pargrWarning)
    write(u,'(A)') trim(pargrWarning(i))
  enddo
endsubroutine write_pargr_warnings

endmodule model_output
