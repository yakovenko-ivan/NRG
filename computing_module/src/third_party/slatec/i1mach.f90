function i1mach ( i )

	use, intrinsic	:: iso_fortran_env, only :	stdin	=>	input_unit, &
												stdout	=>	output_unit, &
												stderr	=>	error_unit
!*****************************************************************************80
!
!! I1MACH returns integer machine constants.
!
!  Discussion:
!
!    Input/output unit numbers.
!
!      I1MACH(1) = the standard input unit.
!      I1MACH(2) = the standard output unit.
!      I1MACH(3) = the standard punch unit.
!      I1MACH(4) = the standard error message unit.
!
!    Words.
!
!      I1MACH(5) = the number of bits per integer storage unit.
!      I1MACH(6) = the number of characters per integer storage unit.
!
!    Integers.
!
!    Assume integers are represented in the S digit base A form:
!
!      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
!
!    where 0 <= X(1:S-1) < A.
!
!      I1MACH(7) = A, the base.
!      I1MACH(8) = S, the number of base A digits.
!      I1MACH(9) = A^S-1, the largest integer.
!
!    Floating point numbers
!
!    Assume floating point numbers are represented in the T digit 
!    base B form:
!
!      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
!
!    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
!
!      I1MACH(10) = B, the base.
!
!    Single precision
!
!      I1MACH(11) = T, the number of base B digits.
!      I1MACH(12) = EMIN, the smallest exponent E.
!      I1MACH(13) = EMAX, the largest exponent E.
!
!    Double precision
!
!      I1MACH(14) = T, the number of base B digits.
!      I1MACH(15) = EMIN, the smallest exponent E.
!      I1MACH(16) = EMAX, the largest exponent E.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 16.
!
!    Output, integer ( kind = 4 ) I1MACH, the value of the chosen parameter.
!
  implicit none
  
  integer ( kind = 4 )				::	i1mach
  integer ( kind = 4 )	,intent(in)	::	i

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i1mach = 0
   ! stop
  else if ( i == 1 ) then
    i1mach = stdin
  else if ( i == 2 ) then
    i1mach = stdout
  else if ( i == 3 ) then
    i1mach = 7
  else if ( i == 4 ) then
    i1mach = stderr
  else if ( i == 5 ) then
    i1mach = 32
  else if ( i == 6 ) then
    i1mach = 4
  else if ( i == 7 ) then
    i1mach = RADIX(1)
  else if ( i == 8 ) then
    i1mach = DIGITS(1)
  else if ( i == 9 ) then
    i1mach = huge(1)
  else if ( i == 10 ) then
    i1mach = RADIX(1.0D+0)
  else if ( i == 11 ) then
    i1mach = DIGITS(1.0E+0)
  else if ( i == 12 ) then
    i1mach = MINEXPONENT(1.0E+0)
  else if ( i == 13 ) then
    i1mach = MAXEXPONENT(1.0E+0)
  else if ( i == 14 ) then
    i1mach = DIGITS(1.0D+0)
  else if ( i == 15 ) then
    i1mach = MINEXPONENT(1.0D+0)
  else if ( i == 16 ) then
    i1mach = MAXEXPONENT(1.0D+0)
  else if ( 16 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i1mach = 0
    stop
  end if

  return
end