module supplementary_routines
    
    use kind_parameters
    
    implicit none
    
    contains 
   
    function lagrange_poly(x,x_0,y,dx)
        real(dp),    dimension(:)	,intent(in)		:: y
		real(dp),    intent(in)						:: x_0   
        real(dp),    intent(in)						:: x
        real(dp),    intent(in)						:: dx
        
        real(dp)		:: lagrange_poly
        real(dp)		:: k
        
        real(dp)	:: a, b, c, d, e, f, g
        real(dp),    dimension(:), allocatable :: P
        real(dp),    dimension(:), allocatable :: x_arr
        real(dp) :: factorial, l
        integer     :: i, j, order
        
        allocate(P(0:size(y)-1))
        allocate(x_arr(0:size(y)-1))
        
        P = y
        factorial = 1
        
        order = size(y)
        
        if (order == 2) then
            x_arr(0) = x_0 - dx
            x_arr(1) = x_0
        else
            do i = 0, order-1
                x_arr(i) = x_0 - dx * ((order-1)/2 - i)
            end do
        end if
        
        do i = 0, order-1
            factorial   = factorial*(i+1)
            do j = 0, order-1
                if(j /= i) then
                    l = (x - x_arr(j))/(x_arr(i) - x_arr(j))
                    P(i)    = P(i) * l
                end if
            end do
        end do
        
        lagrange_poly = 0.0
        do i = 0, order-1
            lagrange_poly = lagrange_poly + P(i)
        end do
        
		!a =  (1.0/24.0)*(y(0)-4.0*y(1)+6.0*y(2)-4.0*y(3)+y(4))/dx**4
		!b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		!c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6.0*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
		!d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
		!e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
  !!  
  !      k = a*x**4 + b*x**3+c*x**2+d*x+e
  !!      
		!lagrange_poly = 4.0*a*x**3+3.0*b*x**2+2.0*c*x+d
        
    end function

    function dlagrange_poly(x,x_0,y,dx)
        real(dp),    dimension(0:6)	,intent(in)		:: y
		real(dp),    intent(in)						:: x_0   
        real(dp),    intent(in)						:: x
        real(dp),    intent(in)						:: dx
        
        real(dp)		:: dlagrange_poly
        real(dp)		:: k
        
        real(dp)	:: a, b, c, d, e, f, g
        real(dp),    dimension(0:6) :: P
        real(dp),    dimension(0:6) :: x_arr, summ
        real(dp) :: factorial, prod
        integer     :: i, j, l, order
        
        P = y
        factorial = 1
        
        order = size(y)
        summ = 0.0
        
        do i = 0, order-1
            x_arr(i) = x_0 - dx * ((order-1)/2 - i)
        end do

        do i = 0, order-1
            do j = 0, order-1
                prod = 1.0_dp
                do l = 0, order-1
                    if ((l /= i).and.(l /= j)) then
                        prod = prod * (x - x_arr(l))
                    end if
                end do
                if ( j /= i ) then
                    summ(i) = summ(i)  + prod
                end if
            end do
        end do        

        do i = 0, order-1
            P(i) = summ(i)
            factorial   = factorial*(i+1)
            prod = 1.0_dp
            do j = 0, order-1
                if(j /= i) then
                    P(i) = P(i) / (x_arr(i) - x_arr(j))
                end if
            end do
        end do
        
        dlagrange_poly = 0.0
        do i = 0, order-1
            dlagrange_poly = dlagrange_poly + P(i) * y(i)
        end do
        
		a =  (1.0/24.0)*(y(0)-4.0*y(1)+6*y(2)-4.0*y(3)+y(4))/dx**4
		b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
!		d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
!		e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
    
!		dlagrange_poly = 12.0*a*x**2+6.0*b*x+2.0*c
        
    end function
    
    
    function d2lagrange_poly(x,x_0,y,dx)
        real(dp),    dimension(0:6)	,intent(in)		:: y
		real(dp),    intent(in)						:: x_0   
        real(dp),    intent(in)						:: x
        real(dp),    intent(in)						:: dx
        
        real(dp)		:: d2lagrange_poly
        real(dp)		:: k
        
        real(dp)	:: a, b, c, d, e, f, g
        real(dp),    dimension(0:6) :: P
        real(dp),    dimension(0:6) :: x_arr, summ2
        real(dp) :: factorial, prod, summ1
        integer     :: i, j, l, m,  order
        
        P = y
        factorial = 1
        
        order = size(y)
        summ1 = 0.0
        summ2 = 0.0
        
        do i = 0, order-1
            x_arr(i) = x_0 - dx * ((order-1)/2 - i)
        end do

        do i = 0, order-1
            summ2(i) = 0.0_dp
            do j = 0, order-1
                summ1 = 0.0_dp
                do l = 0, order-1
                    prod = 1.0_dp
                    do m = 0, order-1 
                        if ((m /= i).and.(m /= j).and.(m /= l)) then
                            prod = prod * (x - x_arr(m)) / (x_arr(i) - x_arr(m))
                        end if
                    end do
                    if ((l /= i).and.(l /= j)) then
                        summ1 = summ1 + 1.0_dp / (x_arr(i) - x_arr(l)) * prod
                    end if     
                end do
                if ( j /= i ) then
                    summ2(i) = summ2(i)  + 1.0_dp / (x_arr(i) - x_arr(j)) * summ1
                end if
            end do
        end do        

        do i = 0, order-1
            P(i) = summ2(i)
        end do
        
        d2lagrange_poly = 0.0
        do i = 0, order-1
            d2lagrange_poly = d2lagrange_poly + P(i) * y(i)
        end do
        
		a =  (1.0/24.0)*(y(0)-4.0*y(1)+6*y(2)-4.0*y(3)+y(4))/dx**4
		b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
!		d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
!		e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
    
!		d2lagrange_poly = 12.0*a*x**2+6.0*b*x+2.0*c
        
    end function    
    
    function five_point_poly(x,x_0,y,dx)
        real(dp),    dimension(0:4)	,intent(in)		:: y
		real(dp),    intent(in)						:: x_0   
        real(dp),    intent(in)						:: x
        real(dp),    intent(in)						:: dx
        
        real(dp)		:: five_point_poly
        real(dp)		:: k
        
        real(dp)	:: a, b, c, d, e
        
		a =  (1.0/24.0)*(y(0)-4.0*y(1)+6.0*y(2)-4.0*y(3)+y(4))/dx**4
		b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6.0*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
		d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
		e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
    
        k = a*x**4 + b*x**3+c*x**2+d*x+e
        
		five_point_poly = 4.0*a*x**3+3.0*b*x**2+2.0*c*x+d
        
    end function
        
	function dfive_point_poly(x,x_0,y,dx)
        real(dp),    dimension(0:4)	,intent(in)		:: y
		real(dp),    intent(in)						:: x_0   
        real(dp),    intent(in)						:: x
        real(dp),    intent(in)						:: dx
        
        real(dp)		:: dfive_point_poly
        
        real(dp)	:: a, b, c, d, e
        
		a =  (1.0/24.0)*(y(0)-4.0*y(1)+6*y(2)-4.0*y(3)+y(4))/dx**4
		b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
!		d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
!		e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
    
		dfive_point_poly = 12.0*a*x**2+6.0*b*x+2.0*c
    end function   
    
    function Newton(x_0,y,dx,tol)
        real(dp),    dimension(0:6)	        		:: y
		real(dp),    intent(in)						:: x_0   
        real(dp),    intent(in)						:: dx
        real(dp),	intent(in)	                    :: tol
        
		real(dp)		:: Newton
        
        real(dp)		:: f, df, x, xn
        integer			:: n
        
        xn = x_0
        x = x_0
        
        do n = 1, 200
            f	= dlagrange_poly(xn,x_0,y,dx)
            if (abs(f)<tol) exit
            df	= d2lagrange_poly(xn,x_0,y,dx)
            x	= xn - f/df
            xn	= x
        end do
        
        Newton = x
!        print *, 'Lagrange' , x
        
        !xn = x_0
        !x = x_0        
        !do n = 1, 200
        !    f	= five_point_poly(xn,x_0,y,dx)
        !    if (abs(f)<tol) exit
        !    df	= dfive_point_poly(xn,x_0,y,dx)
        !    x	= xn - f/df
        !    xn	= x
        !end do
        !
        !Newton = x
!        print *, 'Old' , x
        
    end function
    

    character(len=20) function str_e(k)
	    use kind_parameters

    !   "Convert an real to string."
        real(dp), intent(in) :: k
    
        write (str_e, '(e10.3)') k    
        str_e = adjustl(str_e)
    
    end function str_e
    
    character(len=20) function str_r(k)
	    use kind_parameters

    !   "Convert an real to string."
        real(dp), intent(in) :: k
    
        write (str_r, '(f10.2)') k    
        str_r = adjustl(str_r)
    
    end function str_r
    
    character(len=20) function str_i(k)
	    use kind_parameters

    !   "Convert an integer to string."
        integer, intent(in) :: k
    
        write (str_i, '(I2)') k
        str_i = adjustl(str_i)
    
    end function str_i 
    
end module