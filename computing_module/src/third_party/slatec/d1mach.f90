function d1mach ( i )

!*****************************************************************************80
!
!! D1MACH returns double precision real machine constants.
!
!  Discussion:
!
!    Assuming that the internal representation of a double precision real
!    number is in base B, with T the number of base-B digits in the mantissa,
!    and EMIN the smallest possible exponent and EMAX the largest possible 
!    exponent, then
!
!      D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
!      D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
!      D1MACH(3) = B^(-T), the smallest relative spacing.
!      D1MACH(4) = B^(1-T), the largest relative spacing.
!      D1MACH(5) = log10(B).
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
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 8 ) D1MACH, the value of the chosen parameter.
!
  implicit none

  real		( kind = 8 )				::	d1mach
  integer	( kind = 4 )	,intent(in)	::	i
  real		( kind = 8 )				:: base

  base = radix(1.0D+0)
  
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  else if ( i == 1 ) then
    !d1mach = 4.450147717014403D-308
	d1mach = tiny(1.0D+0)	! RADIX(1.0D+0)**(MINEXPONENT(1.0D+0)-1.0D+0)
	continue
  else if ( i == 2 ) then
    !d1mach = 8.988465674311579D+307
	d1mach = huge(1.0D+0)	! RADIX(1.0D+0)**(MAXEXPONENT(1.0D+0)*(1.0D+0 - (RADIX(1.0D+0))**(-DIGITS(1.0D+0))))
	continue
  else if ( i == 3 ) then
    !d1mach = 1.110223024625157D-016
	d1mach = epsilon(1.0D+0)	/ base ! RADIX(1.0D+0)**(-DIGITS(1.0D+0))
  else if ( i == 4 ) then
    !d1mach = 2.220446049250313D-016
	d1mach = epsilon(1.0D+0)	! RADIX(1.0D+0)**(1.0D+0-DIGITS(1.0D+0))
	continue
  else if ( i == 5 ) then
    !d1mach = 0.301029995663981D+000
	d1mach = dlog10(base)
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  end if

  return
end