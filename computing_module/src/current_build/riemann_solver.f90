module riemann_solver_class

	use kind_parameters
	use global_data

	use data_manager_class

	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	private
	public  riemann_solver, riemann_solver_c

	type 	:: riemann_solver
        
		real(dkind)					:: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r
        real(dkind)					:: rho, p, gamma, u
        integer                     :: contact_direction
        
        real(dkind)                 :: mu_l, mu_r, h_l, h_r, s_l, s_r, G_l, G_r
        
        real(dkind)                 :: conv_rate
        
        integer                     :: counter
     
	contains
		procedure				:: set_parameters
        procedure				:: clear
        procedure				:: solve
        procedure, private		:: F_lr
        procedure, private		:: F_l
        procedure, private		:: F_r
        procedure, private		:: a_l
        procedure, private		:: a_r
        procedure, private		:: F_lr_prime
        procedure, private		:: F_l_prime
        procedure, private		:: F_r_prime
        procedure, private		:: a_l_prime
        procedure, private		:: a_r_prime
        
        procedure               :: get_density
        procedure               :: get_pressure
        procedure               :: get_gamma
        procedure               :: get_velocity
        procedure               :: rightward_contact
        procedure               :: get_activation_count
	end type

	interface riemann_solver_c
		module procedure constructor
	end interface

contains

	type(riemann_solver) function constructor (manager)

		type(data_manager)						,intent(inout)	:: manager

		constructor%rho_l = 0.0_dkind
        constructor%rho_r = 0.0_dkind
        constructor%p_l = 0.0_dkind
        constructor%p_r = 0.0_dkind
        constructor%gamma_l = 0.0_dkind
        constructor%gamma_r = 0.0_dkind
        constructor%u_l = 0.0_dkind
        constructor%u_r = 0.0_dkind
        
        constructor%mu_l = 0.0_dkind
        constructor%mu_r = 0.0_dkind
        constructor%h_l = 0.0_dkind
        constructor%h_r = 0.0_dkind
        constructor%s_l = 0.0_dkind
        constructor%s_r = 0.0_dkind
        constructor%G_l = 0.0_dkind
        constructor%G_r = 0.0_dkind
        
        constructor%rho = 0.0_dkind
        constructor%p = 0.0_dkind
        constructor%gamma = 0.0_dkind
        constructor%u = 0.0_dkind
        constructor%contact_direction = 0
        
        constructor%conv_rate = 0.1_dkind
        
        constructor%counter = 0

    end function

    
	subroutine set_parameters(this, rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r)

		class(riemann_solver)   ,intent(inout)  :: this
		real(dkind)             ,intent(in)     :: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r
        
        this%rho_l = rho_l
        this%rho_r = rho_r
        this%p_l = p_l
        this%p_r = p_r
        this%gamma_l = gamma_l
        this%gamma_r = gamma_r
		this%u_l = u_l
        this%u_r = u_r

        this%mu_l = (gamma_l - 1.0_dkind) / (2.0_dkind * gamma_l)
        this%mu_r = (gamma_r - 1.0_dkind) / (2.0_dkind * gamma_r)
        this%h_l = (gamma_l - 1.0_dkind) / (gamma_l + 1.0_dkind)
        this%h_r = (gamma_r - 1.0_dkind) / (gamma_r + 1.0_dkind)
        this%s_l = p_l * rho_l**(-gamma_l)
        this%s_r = p_r * rho_r**(-gamma_r)
        this%G_l = 2.0 * sqrt(gamma_l) * this%s_l**(0.5_dkind/gamma_l) / (gamma_l - 1.0_dkind)
        this%G_r = 2.0 * sqrt(gamma_r) * this%s_r**(0.5_dkind/gamma_r) / (gamma_r - 1.0_dkind)
	
    end subroutine
    
    
    subroutine clear(this)

		class(riemann_solver) ,intent(inout) :: this
        
        this%rho_l = 0.0_dkind
        this%rho_r = 0.0_dkind
        this%p_l = 0.0_dkind
        this%p_r = 0.0_dkind
        this%gamma_l = 0.0_dkind
        this%gamma_r = 0.0_dkind
		this%u_l = 0.0_dkind
        this%u_r = 0.0_dkind

        this%mu_l = 0.0_dkind
        this%mu_r = 0.0_dkind
        this%h_l = 0.0_dkind
        this%h_r = 0.0_dkind
        this%s_l = 0.0_dkind
        this%s_r = 0.0_dkind
        this%G_l = 0.0_dkind
        this%G_r = 0.0_dkind
        
        this%rho = 0.0_dkind
        this%p = 0.0_dkind
        this%gamma = 0.0_dkind
        this%u = 0.0_dkind
	
    end subroutine
    
    
    subroutine solve(this)

		class(riemann_solver) ,intent(inout) :: this
        
        real(dkind)     :: delta_p = 1.0e-2_dkind
        real(dkind)     :: p_0
        
        real(dkind)     :: p_c, u_c
        real(dkind)     :: func_value, func_deriv_value
        
        real(dkind)     :: max_l, max_r, max
        real(dkind)     :: relative_shock_speed, rho_s, beta
        
        real(dkind)     :: c_l, c_r, c_m, rho_m, s_vac_l, s_vac_r
        
        integer         :: iter
        
        this%counter = this%counter + 1
        
        if (this%p_l > 0.0_dkind .and. this%p_r > 0.0_dkind) then           ! Standart case
            
        ! Contact discontinuity
            
        p_0 = (this%p_l + this%p_r) / 2.0_dkind * min1(this%p_l, this%p_r) / max1(this%p_l, this%p_r)
        
        iter = 1
        do while (iter < 1000)
            func_value = F_lr(this, p_0)
            func_deriv_value = F_lr_prime(this, p_0)
            
            p_c = p_0 - this%conv_rate * func_value/func_deriv_value
            if ((p_c < 0.0_dkind) .or. p_c /= abs(p_c)) p_c = 0.0_dkind
                
            if ((abs(p_0 - p_c) < 1.0_dkind)) exit
                
            p_0 = p_c
            iter = iter + 1
        end do
        if (isnan(p_c)) then
            print *, "Riemann Solver: p_c is nan"
            p_c = (this%p_l + this%p_r) / 2.0_dkind * min1(this%p_l, this%p_r) / max1(this%p_l, this%p_r)
            !pause
        end if  
        
        u_c = (this%F_l(p_c) + this%F_r(p_c)) * 0.5_dkind
            
        max_l = 0.0_dkind
        max_r = 0.0_dkind
        max = 0.0_dkind
        
        if (u_c > 0.0_dkind) then                                                           ! CD to the right
            this%gamma = this%gamma_l
            this%contact_direction = 1
            if (p_c > this%p_l) then                                                        ! Shock wave to the left  (2 domains)
                relative_shock_speed = this%u_l + (p_c - this%p_l) / (u_c - this%u_l) / this%rho_l
                if (relative_shock_speed > 0.0_dkind) then                                  ! At x = 0.0 - left unpertrubed domain
                    this%rho = this%rho_l
                    this%u = this%u_l
                    this%p = this%p_l
                else
                    this%rho = this%rho_l/this%a_l(p_c)                                     ! At x = 0.0 - domain between SW and CD
                    this%u = u_c
                    this%p = p_c
                end if
                max_l = abs(relative_shock_speed)
            else                                                                            ! Rarefication wave to the left (3 domains)
                c_l = sqrt(this%gamma_l * this%p_l / this%rho_l)
                rho_m = (p_c / this%s_l)**(1.0_dkind / this%gamma_l)
                c_m = sqrt(this%gamma_l * p_c / rho_m)
                if (this%u_l - c_l > 0.0_dkind) then                                        ! At x = 0.0 - left unpertrubed domain
                    this%rho = this%rho_l
                    this%u = this%u_l
                    this%p = this%p_l
                elseif (this%u_l - c_l <= 0.0_dkind .and. u_c - c_m > 0.0_dkind) then       ! At x = 0.0 - rarefication fan
                    this%u = ((this%gamma_l - 1.0_dkind) * this%u_l + 2.0_dkind * c_l) / (this%gamma_l + 1.0_dkind)
                    this%p = ((2.0_dkind * c_l / (this%gamma_l - 1.0_dkind) + this%u_l - this%u) / this%G_l)**(1.0_dkind/this%mu_l)
                    this%rho = (this%p / this%s_l)**(1.0_dkind / this%gamma_l)
                else                                                                        ! At x = 0.0 - steady flow domain
                    this%rho = rho_m
                    this%u = u_c
                    this%p = p_c
                end if
                max_l = max1(abs(this%u_l - c_l), abs(u_c - c_m))
                max_l = max1(abs(u_c), max_l)
            end if
        else                                                                                ! CD to the left
            this%gamma = this%gamma_r
            this%contact_direction = -1
            if (p_c > this%p_r) then                                                        ! Shock wave to the right  (2 domains)
                relative_shock_speed = this%u_r + (p_c - this%p_r) / (u_c - this%u_r) / this%rho_r
                if (relative_shock_speed < 0.0_dkind) then                                  ! At x = 0.0 - right unpertrubed domain
                    this%rho = this%rho_r
                    this%u = this%u_r
                    this%p = this%p_r
                else
                    this%rho = this%rho_r/this%a_r(p_c)                                     ! At x = 0.0 - domain between SW and CD
                    this%u = u_c
                    this%p = p_c
                end if
                max_r = abs(relative_shock_speed)
            else                                                                            ! Rarefication wave to the right (3 domains)
                c_r = sqrt(this%gamma_r * this%p_r / this%rho_r)
                rho_m = (p_c / this%s_r)**(1.0_dkind / this%gamma_r)
                c_m = sqrt(this%gamma_r * p_c / rho_m)
                if (this%u_r + c_r < 0.0_dkind) then                                        ! At x = 0.0 - right unpertrubed domain
                    this%rho = this%rho_r
                    this%u = this%u_r
                    this%p = this%p_r
                elseif (this%u_r + c_r >= 0.0_dkind .and. u_c + c_m < 0.0_dkind) then       ! At x = 0.0 - rarefication fan
                    this%u = ((this%gamma_r - 1.0_dkind) * this%u_r - 2.0_dkind * c_r) / (this%gamma_r + 1.0_dkind)
                    this%p = ((2.0_dkind * c_r / (this%gamma_r - 1.0_dkind) - this%u_r + this%u) / this%G_r)**(1.0_dkind/this%mu_r)
                    this%rho = (this%p / this%s_r)**(1.0_dkind / this%gamma_r)
                else                                                                        ! At x = 0.0 - steady flow domain
                    this%rho = rho_m
                    this%u = u_c
                    this%p = p_c
                end if
                max_r = max1(abs(this%u_r + c_r), abs(u_c + c_m))
                max_r = max1(abs(u_c), max_r)
            end if
        end if
        max = max1(max_l, max_r)
            
        elseif (this%p_l <= 0.0_dkind .and. this%p_r > 0.0_dkind) then      ! Vacuum on the left side
            
            this%contact_direction = 1
            c_r         = sqrt(this%gamma_r * this%p_r / this%rho_r)
            s_vac_r     = this%u_r - this%G_r * this%p_r**this%mu_r
            
            if (this%u_r + c_r <= 0.0_dkind) then                           ! At x = 0.0 - right unpertrubed domain
                this%rho = this%rho_r
                this%u = this%u_r
                this%p = this%p_r
            elseif (s_vac_r < 0.0_dkind) then                               ! At x = 0.0 - rarefication fan
                this%u = ((this%gamma_r - 1.0_dkind) * this%u_r - 2.0_dkind * c_r) / (this%gamma_r + 1.0_dkind)
                this%p = ((2.0_dkind * c_r / (this%gamma_r - 1.0_dkind) - this%u_r + this%u) / this%G_r)**(1.0_dkind/this%mu_r)
                this%rho = (this%p / this%s_r)**(1.0_dkind / this%gamma_r)
            else                                                            ! At x = 0.0 - vacuum
                this%rho = 0.0_dkind
                this%u = s_vac_r
                this%p = 0.0_dkind
            end if
            
        elseif (this%p_l > 0.0_dkind .and. this%p_r <= 0.0_dkind) then      ! Vacuum on the right side
            
            this%contact_direction = -1
            c_l         = sqrt(this%gamma_l * this%p_l / this%rho_l)
            s_vac_l     = this%u_l + this%G_l * this%p_l**this%mu_l
            
            if (this%u_l - c_l >= 0.0_dkind) then                           ! At x = 0.0 - left unpertrubed domain
                this%rho = this%rho_l
                this%u = this%u_l
                this%p = this%p_l
            elseif (s_vac_l > 0.0_dkind) then                               ! At x = 0.0 - rarefication fan
                this%u = ((this%gamma_l - 1.0_dkind) * this%u_l + 2.0_dkind * c_l) / (this%gamma_l + 1.0_dkind)
                this%p = ((2.0_dkind * c_l / (this%gamma_l - 1.0_dkind) + this%u_l - this%u) / this%G_l)**(1.0_dkind/this%mu_l)
                this%rho = (this%p / this%s_l)**(1.0_dkind / this%gamma_l)
            else                                                            ! At x = 0.0 - vacuum
                this%rho = 0.0_dkind
                this%u = s_vac_l
                this%p = 0.0_dkind
            end if
            
        else                                                                ! Vacuum on the both sides
            this%rho = 0.0_dkind
            this%u = 0.0_dkind
            this%p = 0.0_dkind
        end if
        
    end subroutine
    
    
    function F_lr(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: F_lr
        real(dkind) ,intent(in)     :: p
        
        F_lr = this%F_l(p) - this%F_r(p)
    end function
    
    function F_l(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: F_l
        real(dkind) ,intent(in)     :: p

        if (abs(p - this%p_l) < 1.0_dkind) then
            F_l = this%u_l
        elseif (p < this%p_l) then
            F_l = this%u_l - this%G_l * (p**this%mu_l - this%p_l**this%mu_l)
        else
            F_l = this%u_l - sqrt((p - this%p_l) * (1.0_dkind - this%a_l(p)) / this%rho_l)
        end if
    end function
    
    function F_r(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: F_r
        real(dkind) ,intent(in)     :: p
        
        if (abs(p - this%p_r) < 1.0_dkind) then
            F_r = this%u_r
        elseif (p < this%p_r) then
            F_r = this%u_r + this%G_r * (p**this%mu_r - this%p_r**this%mu_r)
        else
            F_r = this%u_r + sqrt((p - this%p_r) * (1.0_dkind - this%a_r(p)) / this%rho_r)
        end if
    end function
	
    function a_l(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: a_l
        real(dkind) ,intent(in)     :: p
        
        a_l = (this%p_l + this%h_l * p) / (p + this%h_l * this%p_l)
    end function
    
    function a_r(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: a_r
        real(dkind) ,intent(in)     :: p
        
        a_r = (this%p_r + this%h_r * p) / (p + this%h_r * this%p_r)
    end function
    
    function F_lr_prime(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: F_lr_prime
        real(dkind) ,intent(in)     :: p
        
        F_lr_prime = this%F_l_prime(p) - this%F_r_prime(p)
    end function
    
    function F_l_prime(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: F_l_prime
        real(dkind) ,intent(in)     :: p
        
        if (abs(p - this%p_r) < 1.0_dkind) then
            F_l_prime = - this%G_l * this%mu_l * p**(this%mu_l - 1.0_dkind)
        elseif (p < this%p_l) then
            F_l_prime = - this%G_l * this%mu_l * p**(this%mu_l - 1.0_dkind)
        else
            F_l_prime = - 0.5_dkind * (1.0_dkind - this%a_l(p) - (p - this%p_l) * this%a_l_prime(p)) / sqrt((p - this%p_l) * (1.0_dkind - this%a_l(p)) / this%rho_l) / this%rho_l
        end if
    end function
    
    function F_r_prime(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: F_r_prime
        real(dkind) ,intent(in)     :: p

        if (abs(p - this%p_r) < 1.0) then
            F_r_prime = this%G_r * this%mu_r * p**(this%mu_r - 1.0_dkind)
        elseif (p < this%p_r) then
            F_r_prime = this%G_r * this%mu_r * p**(this%mu_r - 1.0_dkind)
        else
            F_r_prime = 0.5_dkind * (1.0_dkind - this%a_r(p) - (p - this%p_r) * this%a_r_prime(p)) / sqrt((p - this%p_r) * (1.0_dkind - this%a_r(p)) / this%rho_r) / this%rho_r
        end if
    end function
    
    function a_l_prime(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: a_l_prime
        real(dkind) ,intent(in)     :: p
        
        a_l_prime = this%p_l * (this%h_l**2 - 1.0_dkind) / (p + this%h_l * this%p_l)**2
    end function
    
    function a_r_prime(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: a_r_prime
        real(dkind) ,intent(in)     :: p
        
        a_r_prime = this%p_r * (this%h_r**2 - 1.0_dkind) / (p + this%h_r * this%p_r)**2
    end function
    
    
    function get_density(this)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: get_density
        
        get_density = this%rho
    end function
    
    function get_pressure(this)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: get_pressure
        
        get_pressure = this%p
    end function
    
    function get_gamma(this)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: get_gamma
        
        get_gamma = this%gamma
    end function
    
    function get_velocity(this)
        class(riemann_solver) ,intent(in) :: this
        real(dkind)                 :: get_velocity
        
        get_velocity = this%u
    end function
    
    function rightward_contact(this)
        class(riemann_solver) ,intent(in) :: this
        logical                 :: rightward_contact
        
        if (this%contact_direction == 1) then
            rightward_contact = .true.
        else
            rightward_contact = .false.
        end if
    end function
    
    function get_activation_count(this)
        class(riemann_solver) ,intent(in) :: this
        integer                 :: get_activation_count
        
        get_activation_count = this%counter
    end function
    
end module