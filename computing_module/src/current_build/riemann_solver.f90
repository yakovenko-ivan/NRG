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
        
		real(dp)					:: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r
        real(dp)					:: rho, p, gamma, u
        
        real(dp)                 :: mu_l, mu_r, h_l, h_r, s_l, s_r, G_l, G_r
        
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
        
        procedure               :: get_density
        procedure               :: get_pressure
        procedure               :: get_gamma
        procedure               :: get_velocity
        procedure               :: get_activation_count
	end type

	interface riemann_solver_c
		module procedure constructor
	end interface

contains

	type(riemann_solver) function constructor (manager)

		type(data_manager)						,intent(inout)	:: manager

		constructor%rho_l = 0.0_dp
        constructor%rho_r = 0.0_dp
        constructor%p_l = 0.0_dp
        constructor%p_r = 0.0_dp
        constructor%gamma_l = 0.0_dp
        constructor%gamma_r = 0.0_dp
        constructor%u_l = 0.0_dp
        constructor%u_r = 0.0_dp
        
        constructor%mu_l = 0.0_dp
        constructor%mu_r = 0.0_dp
        constructor%h_l = 0.0_dp
        constructor%h_r = 0.0_dp
        constructor%s_l = 0.0_dp
        constructor%s_r = 0.0_dp
        constructor%G_l = 0.0_dp
        constructor%G_r = 0.0_dp
        
        constructor%rho = 0.0_dp
        constructor%p = 0.0_dp
        constructor%gamma = 0.0_dp
        constructor%u = 0.0_dp
        
        constructor%counter = 0

    end function

    
	subroutine set_parameters(this, rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r)

		class(riemann_solver)   ,intent(inout)  :: this
		real(dp)             ,intent(in)     :: rho_l, rho_r, p_l, p_r, gamma_l, gamma_r, u_l, u_r
        
        this%rho_l = rho_l
        this%rho_r = rho_r
        this%p_l = p_l
        this%p_r = p_r
        this%gamma_l = gamma_l
        this%gamma_r = gamma_r
		this%u_l = u_l
        this%u_r = u_r

        this%mu_l = (gamma_l - 1.0_dp) / (2.0_dp * gamma_l)
        this%mu_r = (gamma_r - 1.0_dp) / (2.0_dp * gamma_r)
        this%h_l = (gamma_l - 1.0_dp) / (gamma_l + 1.0_dp)
        this%h_r = (gamma_r - 1.0_dp) / (gamma_r + 1.0_dp)
        this%s_l = p_l * rho_l**(-gamma_l)
        this%s_r = p_r * rho_r**(-gamma_r)
        this%G_l = 2.0 * sqrt(gamma_l) * this%s_l**(0.5_dp/gamma_l) / (gamma_l - 1.0_dp)
        this%G_r = 2.0 * sqrt(gamma_r) * this%s_r**(0.5_dp/gamma_r) / (gamma_r - 1.0_dp)
	
    end subroutine
    
    
    subroutine clear(this)

		class(riemann_solver) ,intent(inout) :: this
        
        this%rho_l = 0.0_dp
        this%rho_r = 0.0_dp
        this%p_l = 0.0_dp
        this%p_r = 0.0_dp
        this%gamma_l = 0.0_dp
        this%gamma_r = 0.0_dp
		this%u_l = 0.0_dp
        this%u_r = 0.0_dp

        this%mu_l = 0.0_dp
        this%mu_r = 0.0_dp
        this%h_l = 0.0_dp
        this%h_r = 0.0_dp
        this%s_l = 0.0_dp
        this%s_r = 0.0_dp
        this%G_l = 0.0_dp
        this%G_r = 0.0_dp
        
        this%rho = 0.0_dp
        this%p = 0.0_dp
        this%gamma = 0.0_dp
        this%u = 0.0_dp
	
    end subroutine
    
    
    subroutine solve(this)

		class(riemann_solver) ,intent(inout) :: this
        
        real(dp)     :: delta_p = 1.0e-2_dp
        real(dp)     :: p_0 = 100.0_dp
        
        real(dp)     :: p_c, u_c
        real(dp)     :: func_value, func_deriv_value
        
        real(dp)     :: max_l, max_r, max
        real(dp)     :: relative_shock_speed, rho_s, beta
        
        real(dp)     :: c_l, c_r, c_m, rho_m
        
        integer         :: iter
        
        this%counter = this%counter + 1
        
        ! ���������� ������
        do iter = 1,20
            func_value = F_lr(this, p_0)
            func_deriv_value = (this%F_lr(p_0 + delta_p) - this%F_lr(p_0 - delta_p)) / (2.0 * delta_p)
            
            p_c = p_0 - func_value/func_deriv_value
            p_0 = p_c
        end do
        
        u_c = this%F_l(p_c)
        
        ! ����������� ���������� �� ������
        max_l = 0.0_dp
        max_r = 0.0_dp
        max = 0.0_dp
        
        if (u_c > 0.0_dp) then                           ! �� ������
            this%gamma = this%gamma_l
            if (p_c > this%p_l) then                       ! �� �����  (2 �������)
                relative_shock_speed = this%u_l + (p_c - this%p_l) / (u_c - this%u_l) / this%rho_l      ! �������� ��
                if (relative_shock_speed > 0.0_dp) then  ! ����� ������������� �������
                    this%rho = this%rho_l
                    this%u = this%u_l
                    this%p = this%p_l
                else
                    this%rho = this%rho_l/this%a_l(p_c)     ! ������� ����� �� � ��
                    this%u = u_c
                    this%p = p_c
                end if
                max_l = abs(relative_shock_speed)
            else                                            ! ����� ���������� �����  (3 �������)
                c_l = sqrt(this%gamma_l * this%p_l / this%rho_l)        ! �������� ����� �����
                rho_m = (p_c / this%s_l)**(1.0_dp / this%gamma_l)    ! �������������� ���������
                c_m = sqrt(this%gamma_l * p_c / rho_m)                  ! �������� ����� � �������������� �������
                if (this%u_l - c_l > 0.0_dp) then        ! ����� ������������� �������
                    this%rho = this%rho_l
                    this%u = this%u_l
                    this%p = this%p_l
                elseif (this%u_l - c_l <= 0.0_dp .and. u_c - c_m > 0.0_dp) then   ! ������� ����� ����������
                    this%u = ((this%gamma_l - 1.0_dp) * this%u_l + 2.0_dp * c_l) / (this%gamma_l + 1.0_dp)
                    this%p = ((2.0_dp * c_l / (this%gamma_l - 1.0_dp) + this%u_l - this%u) / this%G_l)**(1.0_dp/this%mu_l)
                    this%rho = (this%p / this%s_l)**(1.0_dp / this%gamma_l)
                else                                        ! ������� ��������������� �������
                    this%rho = rho_m
                    this%u = u_c
                    this%p = p_c
                end if
                max_l = max1(abs(this%u_l - c_l), abs(u_c - c_m))
                max_l = max1(abs(u_c), max_l)
            end if
        else                                                ! �� �����
            this%gamma = this%gamma_r
            if (p_c > this%p_r) then                       ! �� ������  (2 �������)
                beta = (this%u_r - u_c) * (this%u_r - u_c) / (p_c - this%p_r)
                rho_s = this%rho_r / (1.0_dp - beta * this%rho_r)            ! ��������� �� ��
                relative_shock_speed = (u_c * rho_s - this%u_r * this%rho_r) / (rho_s - this%rho_r)      ! �������� ��
                !relative_shock_speed = this%u_r + (p_c - this%p_r) / (u_c - this%u_r) / this%rho_r
                if (relative_shock_speed < 0.0_dp) then  ! ������ ������������� �������
                    this%rho = this%rho_r
                    this%u = this%u_r
                    this%p = this%p_r
                else
                    this%rho = this%rho_r/this%a_r(p_c)     ! ������� ����� �� � ��
                    this%u = u_c
                    this%p = p_c
                end if
                max_r = abs(relative_shock_speed)
            else                                            ! ����� ���������� ������  (3 �������)
                c_r = sqrt(this%gamma_r * this%p_r / this%rho_r)        ! �������� ����� �����
                rho_m = (p_c / this%s_r)**(1.0_dp / this%gamma_r)     ! �������������� ���������
                c_m = sqrt(this%gamma_r * p_c / rho_m)                  ! �������� ����� � �������������� �������
                if (this%u_r + c_r < 0.0_dp) then        ! ������ ������������� �������
                    this%rho = this%rho_r
                    this%u = this%u_r
                    this%p = this%p_r
                elseif (this%u_r + c_r >= 0.0_dp .and. u_c + c_m < 0.0_dp) then   ! ������� ����� ����������
                    this%u = ((this%gamma_r - 1.0_dp) * this%u_r - 2.0_dp * c_r) / (this%gamma_r + 1.0_dp)
                    this%p = ((2.0_dp * c_r / (this%gamma_r - 1.0_dp) - this%u_r + this%u) / this%G_r)**(1.0_dp/this%mu_r)
                    this%rho = (this%p / this%s_r)**(1.0_dp / this%gamma_r)
                else                                        ! ������� ��������������� �������
                    this%rho = rho_m
                    this%u = u_c
                    this%p = p_c
                end if
                max_r = max1(abs(this%u_r + c_r), abs(u_c + c_m))
                max_r = max1(abs(u_c), max_r)
            end if
        end if
        max = max1(max_l, max_r)
        
    end subroutine
    
    
    function F_lr(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: F_lr
        real(dp) ,intent(in)     :: p
        
        F_lr = this%F_l(p) - this%F_r(p)
    end function
    
    function F_l(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: F_l
        real(dp) ,intent(in)     :: p

        if (p < this%p_l) then
            F_l = this%u_l - this%G_l * (p**this%mu_l - this%p_l**this%mu_l)
        else
            F_l = this%u_l - sqrt((p - this%p_l) * (1.0_dp - this%a_l(p)) / this%rho_l)
        end if
    end function
    
    function F_r(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: F_r
        real(dp) ,intent(in)     :: p
        
        if (p < this%p_r) then
            F_r = this%u_r + this%G_r * (p**this%mu_r - this%p_r**this%mu_r)
        else
            F_r = this%u_r + sqrt((p - this%p_r) * (1.0_dp - this%a_r(p)) / this%rho_r)
        end if
    end function
	
    function a_l(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: a_l
        real(dp) ,intent(in)     :: p
        
        a_l = (this%p_l + this%h_l * p) / (p + this%h_l * this%p_l)
    end function
    
    function a_r(this, p)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: a_r
        real(dp) ,intent(in)     :: p
        
        a_r = (this%p_r + this%h_r * p) / (p + this%h_r * this%p_r)
    end function
    
    
    function get_density(this)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: get_density
        
        get_density = this%rho
    end function
    
    function get_pressure(this)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: get_pressure
        
        get_pressure = this%p
    end function
    
    function get_gamma(this)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: get_gamma
        
        get_gamma = this%gamma
    end function
    
    function get_velocity(this)
        class(riemann_solver) ,intent(in) :: this
        real(dp)                 :: get_velocity
        
        get_velocity = this%u
    end function
    
    function get_activation_count(this)
        class(riemann_solver) ,intent(in) :: this
        integer                 :: get_activation_count
        
        get_activation_count = this%counter
    end function
    
end module