program computing_module

!	use IFPORT , only: SETENVQQ

	use global_data
	use kind_parameters
	use computational_domain_class
	use chemical_properties_class
	use thermophysical_properties_class
	use solver_options_class
	use computational_mesh_class
	use mpi_communications_class
	use data_manager_class
	use boundary_conditions_class
	use field_scalar_class
	use field_vector_class
	use data_save_class
	use data_io_class
	use post_processor_manager_class

	use cpm_solver_class	
	use cabaret_solver_class
	use fds_low_mach_solver_class
	
    use benchmarking
	
#ifdef mpi
	use MPI
#endif

	implicit none

#ifdef OMP
	include "omp_lib.h"
#endif

	type(computational_domain)					:: problem_domain
	type(mpi_communications)					:: problem_mpi_support
	type(computational_mesh)			,target	:: problem_mesh
	type(data_manager)							:: problem_manager
	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics
	type(cabaret_solver)						:: problem_cabaret_solver
	type(cpm_solver)							:: problem_cpm_solver
	type(fds_solver)							:: problem_fds_solver

	type(solver_options)						:: problem_solver_options
	
	type(boundary_conditions)	,target			:: problem_boundaries
	type(field_scalar_cons)		,target			:: p, T, rho, e_i, E_f, mol_mix_conc
	type(field_vector_cons)		,target			:: v, Y

	type(data_io)								:: problem_data_io
	type(data_save)								:: problem_data_save

	type(post_processor_manager)				:: problem_post_proc_manager

    type(timer)     :: main_clock
 
	integer	:: c(8)
	integer	:: day, h, m , s

	integer		:: log_unit, mpi_io_unit

	integer		:: num_threads, max_threads, processor_rank, num, nth

	integer		:: iter, num_iterations, test
	real(dp)	:: calculation_time, time_step
	logical		:: stop_flag, precision_flag, benchmarking, there
	
	integer		:: error

    character(len=*), parameter :: VERSION = '1.0'
    character(len=32)           :: arg, sub_arg
    integer                     :: eq_index, arg_val
    integer                     :: i
    
    integer                     :: bench_io
    character(len=100)          :: bechmark_result_file
    logical                     :: new_bench_io
    
    num_threads     = 1
    benchmarking    = .false.
    
    do i = 1, command_argument_count()
        call get_command_argument(i, arg)

        sub_arg     = arg
        eq_index    = index(sub_arg,'=')
        if (eq_index /= 0) then
            sub_arg = arg(:eq_index-1)
            read (arg(eq_index+1:),'(I6)') arg_val
        end if

        select case (sub_arg)
            case ('-v', '--version')
                print '(A)', 'Joint institute for High Temperature of RAS'
                print '(A)', 'Numerical Reactive Gasdynamics (NRG) package'
                print '(2A)', 'Version ', VERSION
                stop

            case ('-h', '--help')
                call print_help()
                stop

            case ('--num_threads')
                num_threads = arg_val
                
            case ('--benchmark')
                num_iterations = arg_val
                benchmarking    = .true.
                
                write(bechmark_result_file,'(A,I6.6,A)') 'bench_ni_', num_iterations, '.dat'
                inquire(file = bechmark_result_file, exist = there)
                
                if (.not.there) then
                    open(newunit = bench_io, file = bechmark_result_file, status = 'new', form = 'formatted')
                    new_bench_io = .true.
                else
                    open(newunit = bench_io, file = bechmark_result_file, status = 'old', position = 'append', form = 'formatted')
                end if
                
                !do test = 1, 10
                !    write(bechmark_result_file,'(A,I2.2,A,I6.6,A,I2.2,A)') 'bench_nt',num_threads,'_ni_',num_iterations,'_test_',test,'.dat'
                !    inquire(file = bechmark_result_file, exist = there)
                !    if (.not.there) then
                !        open(newunit = bench_io, file = bechmark_result_file, status = 'replace', form = 'formatted')
                !        exit
                !    end if
                !end do
                
            case default
                print '(2A, /)', 'unrecognised command-line option: ', arg
                call print_help()
                stop
        end select
    end do
    
#ifdef mpi
	call mpi_init(error)
#endif

#ifdef OMP
    if(num_threads /= 0) then
        print *, 'OpenMP initialization. Max threads:', num_threads
    else
        print *, 'OpenMP initialization. Threads number undefined, please set --num_threads option.'
        print *, 'Using default threads number', num_threads
    end if

    call omp_set_num_threads(num_threads)
    max_threads = omp_get_max_threads()
	
	!$omp parallel private(nth, num)
	nth		=	omp_get_num_threads()
	num		=	omp_get_thread_num()
	if (num == 0) then
		print *, "Total number of threads: ", nth
	else
		print *, "Thread number: ", num
	end if
	!$omp end parallel
#endif

	open(newunit = log_unit, file = problem_setup_log_file, status = 'old', form = 'formatted', position = 'append')

	problem_domain 			= computational_domain_c()

	problem_chemistry		= chemical_properties_c()
	problem_thermophysics	= thermophysical_properties_c(problem_chemistry)
	
	problem_mpi_support		= mpi_communications_c(problem_domain)

	problem_solver_options	= solver_options_c()

	problem_manager			= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics,problem_solver_options)
    
	call problem_manager%create_boundary_conditions(problem_boundaries)
	call problem_manager%create_computational_mesh(problem_mesh)
	call problem_manager%create_scalar_field(p				,'pressure'						,'P')
	call problem_manager%create_scalar_field(T 				,'temperature'					,'T')
	call problem_manager%create_scalar_field(rho			,'density'						,'rho')
	call problem_manager%create_scalar_field(e_i			,'internal_energy'				,'e_i')
	call problem_manager%create_scalar_field(E_f			,'full_energy'					,'E_f')
	call problem_manager%create_scalar_field(mol_mix_conc	,'mixture_molar_concentration'	,'mix_mol_conc')
	call problem_manager%create_vector_field(v				,'velocity'						,'v',	'spatial')
	call problem_manager%create_vector_field(Y				,'specie_molar_concentration '	,'Y',	'chemical')

!	problem_data_io				= data_io_c(problem_manager,calculation_time)

	select case(problem_solver_options%get_solver_name())
		case('cpm')
			problem_cpm_solver = cpm_solver_c(	problem_manager,  &
												problem_data_io			= problem_data_io)
		case('CABARET')
			problem_cabaret_solver = cabaret_solver_c(	problem_manager,  &
														problem_data_io			= problem_data_io)
		case('fds_low_mach')											
			problem_fds_solver = fds_solver_c(	problem_manager,  &
												problem_data_io			= problem_data_io)																		
	end select

	problem_data_save			= data_save_c(problem_manager)
	problem_post_proc_manager	= post_processor_manager_c(problem_manager)	

	processor_rank = problem_domain%get_processor_rank()
	if (processor_rank == 0) then
		call problem_domain				%write_log(log_unit)													
		call problem_chemistry			%write_log(log_unit)		
		call problem_thermophysics		%write_log(log_unit)
		call problem_post_proc_manager	%write_log(log_unit)
		call problem_data_save			%write_log(log_unit)
		call problem_data_io			%write_log(log_unit)
		call problem_boundaries			%write_log(log_unit)
		call problem_solver_options		%write_log(log_unit)
	end if
	
    call problem_manager%create_timer(main_clock,'Main cycle time', 'MAIN')

	iter = 0
	stop_flag		= .false.
	precision_flag	= .false.
	do while(.not.stop_flag)
        
        call main_clock%tic()
        
		iter = iter + 1

		select case(problem_solver_options%get_solver_name())
			case('cpm')
				call problem_cpm_solver%solve_problem()		
				calculation_time	= problem_cpm_solver%get_time()		
				time_step			= problem_cpm_solver%get_time_step()
			case('CABARET')
				call problem_cabaret_solver%solve_problem()
				calculation_time	= problem_cabaret_solver%get_time()
				time_step			= problem_cabaret_solver%get_time_step()
			case('fds_low_mach')											
				call problem_fds_solver%solve_problem(iter,stop_flag)		
				calculation_time	= problem_fds_solver%get_time()		
				time_step			= problem_fds_solver%get_time_step()					
		end select
		
		if ((mod(iter-1,1) == 0).and.(processor_rank == 0)) then
			print *, ' Calculation time : ', calculation_time
			print *, ' Current time step : ', time_step
			print *, ' Amount of iterations : ', iter
  
			call date_and_time(values=c)
			day=c(3)
			h=c(5)
			m=c(6)
			s=c(7)
			print *, ' Current time = ', day,'  ',h,':',m,':',s
        end if
		
!		if ((precision_flag).and.(calculation_time > 165.0e-09_dp)) then
!		if (time_step < 3.66e-08_dp * problem_cpm_solver%get_CFL_coefficient()) then
!			call problem_data_save%set_save_time(1.0_dp)
!			call problem_cpm_solver%set_CFL_coefficient(0.1_dp)
!		else
!			call problem_cpm_solver%set_CFL_coefficient(0.75_dp)
!		end if
	
!		if (calculation_time > 1.0e-05_dp) then
!			call problem_data_save%set_save_time(200.0_dp)
!            stop
!        end if    
            
        call problem_data_io%output_all_data(calculation_time			,stop_flag)	
		call problem_post_proc_manager%process_data(calculation_time	,stop_flag)
		call problem_data_save%save_all_data(calculation_time			,stop_flag)

        call main_clock%toc(new_iter=.true.)

        if ((benchmarking).and.(iter == num_iterations)) then
            exit
        end if
	end do

    if (benchmarking) then
        call problem_manager%print_all_clocks(num_threads,lun=bench_io,new_io=new_bench_io)
    end if
    
#ifdef mpi
	call mpi_finalize(error)
#endif

contains
    
    subroutine print_help()
        print '(a, /)', 'command-line options:'
        print '(a)',    '  -v, --version     print version information and exit'
        print '(a, /)', '  -h, --help        print usage information and exit'
    end subroutine print_help

end program
