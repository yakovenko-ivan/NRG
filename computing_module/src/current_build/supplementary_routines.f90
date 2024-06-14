module supplementary_routines
    
    use kind_parameters
    
    implicit none
    
    contains 
   
    function five_point_poly(x,x_0,y,dx)
        real(dkind),    dimension(0:4)	,intent(in)		:: y
		real(dkind),    intent(in)						:: x_0   
        real(dkind),    intent(in)						:: x
        real(dkind),    intent(in)						:: dx
        
        real(dkind)		:: five_point_poly
        real(dkind)		:: k
        
        real(dkind)	:: a, b, c, d, e
        
		a =  (1.0/24.0)*(y(0)-4.0*y(1)+6.0*y(2)-4.0*y(3)+y(4))/dx**4
		b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6.0*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
		d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
		e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
    
        k = a*x**4 + b*x**3+c*x**2+d*x+e
        
		five_point_poly = 4.0*a*x**3+3.0*b*x**2+2.0*c*x+d
        
    end function
        
	function dfive_point_poly(x,x_0,y,dx)
        real(dkind),    dimension(0:4)	,intent(in)		:: y
		real(dkind),    intent(in)						:: x_0   
        real(dkind),    intent(in)						:: x
        real(dkind),    intent(in)						:: dx
        
        real(dkind)		:: dfive_point_poly
        
        real(dkind)	:: a, b, c, d, e
        
		a =  (1.0/24.0)*(y(0)-4.0*y(1)+6*y(2)-4.0*y(3)+y(4))/dx**4
		b = -(1.0/12.0)*(dx*y(0)-2.0*dx*y(1)+2.0*dx*y(3)-dx*y(4)+2.0*x_0*y(0)-8.0*x_0*y(1)+12.0*x_0*y(2)-8.0*x_0*y(3)+2.0*x_0*y(4))/dx**4
		c = -(1.0/24.0)*(dx**2*y(0)-16.0*dx**2.0*y(1)+30.0*dx**2*y(2)-16.0*dx**2*y(3)+dx**2*y(4)-6.0*dx*x_0*y(0)+12.0*dx*x_0*y(1)-12.0*dx*x_0*y(3)+6*dx*x_0*y(4)-6.0*x_0**2*y(0)+24.0*x_0**2*y(1)-36.0*x_0**2*y(2)+24.0*x_0**2*y(3)-6.0*x_0**2*y(4))/dx**4
!		d = (1.0/12.0)*(dx**3*y(0)-8.0*dx**3*y(1)+8.0*dx**3*y(3)-dx**3*y(4)+dx**2*x_0*y(0)-16*dx**2.0*x_0*y(1)+30.0*dx**2*x_0*y(2)-16.0*dx**2*x_0*y(3)+dx**2*x_0*y(4)-3*dx*x_0**2*y(0)+6.0*dx*x_0**2*y(1)-6.0*dx*x_0**2*y(3)+3.0*dx*x_0**2*y(4)-2.0*x_0**3*y(0)+8.0*x_0**3*y(1)-12.0*x_0**3*y(2)+8*x_0**3*y(3)-2.0*x_0**3*y(4))/dx**4
!		e = (1.0/24.0)*(24.0*dx**4*y(2)-2.0*dx**3*x_0*y(0)+16.0*dx**3*x_0*y(1)-16*dx**3*x_0*y(3)+2.0*dx**3*x_0*y(4)-dx**2*x_0**2*y(0)+16.0*dx**2*x_0**2*y(1)-30.0*dx**2*x_0**2*y(2)+16.0*dx**2*x_0**2*y(3)-dx**2*x_0**2*y(4)+2.0*dx*x_0**3*y(0)-4.0*dx*x_0**3*y(1)+4.0*dx*x_0**3*y(3)-2.0*dx*x_0**3*y(4)+x_0**4*y(0)-4.0*x_0**4*y(1)+6.0*x_0**4*y(2)-4.0*x_0**4*y(3)+x_0**4*y(4))/dx**4
    
		dfive_point_poly = 12.0*a*x**2+6.0*b*x+2.0*c
    end function   
    
    function Newton(x_0,y,dx,tol)
        real(dkind),    dimension(0:4)	,intent(in)		:: y
		real(dkind),    intent(in)						:: x_0   
        real(dkind),    intent(in)						:: dx
        real(dkind),	intent(in)	                    :: tol
        
		real(dkind)		:: Newton
        
        real(dkind)		:: f, df, x, xn
        integer			:: n
        
        xn = x_0
        
        do n = 1, 20
            f	= five_point_poly(xn,x_0,y,dx)
            if (abs(f)<tol) exit
            df	= dfive_point_poly(xn,x_0,y,dx)
            x	= xn - f/df
            xn	= x
        end do
    
        Newton = x
        
    end function
end module