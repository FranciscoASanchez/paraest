!> General utilities of paraest.
module utilities
  implicit none

contains

  !> Converts an integer to text, without blanks: 12 -> '12'.
  !>
  !> Used to build formats whose repeat count is known only at run time
  !> (what Intel wrote as <n>):
  !>
  !>     write(*,'(/,'//i2s(ng)//'(F7.4,2X),"  !g*")') (g(i), i = 1, ng)
  !>
  !> The result is a deferred-length string (character(:), allocatable):
  !> its length adjusts to the number of digits.
  function i2s(n) result(text)
    integer, intent(in) :: n
    character(:), allocatable :: text
    character(20) :: aux

    write(aux, '(I0)') n
    text = trim(aux)
  end function i2s

  !> Builds a repeated group of a format, followed by a comma:
  !>
  !>     rep(3, 'F7.4,2X')  ->  '3(F7.4,2X),'
  !>     rep(0, 'F7.4,2X')  ->  ''            (the group disappears)
  !>
  !> Typical use:
  !>
  !>     write(*,'(/,'//rep(ng,'F7.4,2X')//'"  !g*")') (g(i), i = 1, ng)
  !>
  !> With n = 0, '0(...)' cannot be written in a format; Intel accepted it
  !> in its <n> form, so the empty string is returned and the rest of the
  !> format is written anyway.
  function rep(n, descriptor) result(text)
    integer, intent(in) :: n
    character(*), intent(in) :: descriptor
    character(:), allocatable :: text

    if (n > 0) then
      text = i2s(n)//'('//descriptor//'),'
    else
      text = ''
    end if
  end function rep

  !> Stops the program with a clear message if a quantity read from the
  !> data file exceeds the maximum size of the arrays:
  !>
  !>     call check_limit(nnp, NPTIPOM, 'interaction parameters (NNP)')
  !>
  !> Without this check the program would write outside the arrays, and the
  !> results (or the titles) would come out with garbage, without warning.
  subroutine check_limit(n, nmax, what)
    integer, intent(in) :: n, nmax
    character(*), intent(in) :: what

    if (n > nmax) then
      write(*, '(/, "Error: there are ", i0, " ", a, "; the maximum is ", i0, ".")') n, what, nmax
      write(*, '("To increase it, change src/main/dimensions.f90 and rebuild.")')
      stop 1
    end if
  end subroutine check_limit

end module utilities
