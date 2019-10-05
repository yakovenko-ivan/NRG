program computing_module

	use IFPORT

	use global_data
	use kind_parameters
	use computational_domain_class
	use mpi_communications_class	
	use chemical_properties_class
	use thermophysical_properties_class	
	use data_manager_class
	use data_io_class

	implicit none

	type(computational_domain)					:: problem_domain
	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics	
	type(data_manager)							:: problem_manager
	type(mpi_communications)					:: problem_mpi_support
	type(data_io)								:: problem_data_io

	integer	:: processors_number

	character(len=100)	:: system_command, file_name, data_output_folder, file_path, header_file_name, data_file_name
	character(len=100)	:: field_name

	real(dkind)		,dimension(:,:,:)	,allocatable	:: read_buffer
	real(dkind)		,dimension(:,:,:)	,allocatable	:: write_buffer

	real(dkind)					:: calculation_time

	integer		,dimension(3)	:: processor_grid

	integer		:: sep_index, argument_index, flow_index
	integer		:: dir_io, write_io 
	integer	,dimension(:)	,allocatable	:: read_io

	character(len=10)	:: full_command_line
	character(len=10)	:: a_single_argument

	integer	,dimension(3,2)	:: global_utter_cells_bounds ,local_utter_cells_bounds
	integer	,dimension(3,2)	:: global_inner_faces_bounds ,local_inner_faces_bounds	
	integer		:: load_counter
	integer		:: processor_rank, dimensions
	integer		:: i, j, k, dim
	integer	,dimension(3,2)	:: ind_r, ind_w
	integer		:: proc_i, proc_j, proc_k, proc
	integer		:: var = 0
	integer		:: io_stat
	integer		:: error = 0

	if (command_argument_count() == 0) then
		stop "Error: please, specify number of processors used with proc_num argument"
	end if

	call get_command_argument(2, value = a_single_argument, status = io_stat)

	print *, trim(a_single_argument)
	read(a_single_argument,*) processors_number

	problem_domain 			= computational_domain_c()
	problem_mpi_support		= mpi_communications_c(problem_domain)
	
	problem_chemistry		= chemical_properties_c()
	problem_thermophysics	= thermophysical_properties_c(problem_chemistry)

	problem_manager			= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics)

	problem_data_io			= data_io_c(problem_manager,calculation_time)

	call problem_domain%set_mpi_communicator_size(processors_number)
	call problem_domain%generate_processors_grid()

	dimensions 		= problem_domain%get_domain_dimensions()
	processor_grid 	= problem_domain%get_processor_number()

	allocate(read_io(0:processor_grid(1)*processor_grid(2)*processor_grid(3)-1))

	data_output_folder	= trim(problem_data_io%get_data_output_folder())
	load_counter		= problem_data_io%get_load_counter()

	write(system_command,'(A,A,A,I3.3,A)') 'ls . -d ', trim(data_output_folder) ,trim(fold_sep), load_counter, '/*/ > dir.txt'
	call execute_command_line(system_command) 

	open(newunit = dir_io, file = 'dir.txt', status = 'old', form = 'formatted')
	read(dir_io,'(A)') file_path

	do 
		read(dir_io,'(A)',iostat = error) file_path
		if (error /= 0) exit

		sep_index = index(file_path,trim(fold_sep))

		field_name = file_path(sep_index+1:)
		sep_index = index(field_name,trim(fold_sep))
		field_name = field_name(sep_index+1:)
		sep_index = index(field_name,trim(fold_sep))		
		field_name = field_name(:sep_index-1)

		write(file_name,'(A,A,I3.3,A,A,A)') trim(data_output_folder) ,trim(fold_sep), load_counter, trim(fold_sep), trim(field_name), '.dat'
		open(newunit = write_io, file = file_name, status = 'replace', form = 'binary', position = 'append') 

		do proc_i = 0, processor_grid(1) -1 
		do proc_j = 0, processor_grid(2) -1
		do proc_k = 0, processor_grid(3) -1
			processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
			write(data_file_name,'(A,A,I3.3,A,A,A,A,A,I4.4,A)') trim(data_output_folder) ,trim(fold_sep), load_counter, trim(fold_sep), trim(field_name) ,trim(fold_sep), trim(field_name) ,'_proc_', processor_rank, '.dat'
			open(newunit = read_io(processor_rank), file = data_file_name, status = 'old', form = 'binary') 
		end do
		end do
		end do

		flow_index = index(field_name,'flow')

		ind_r = 1
		ind_w = 1
		ind_w(:dimensions,:) = 0

		if (flow_index == 0) then
			global_utter_cells_bounds		= problem_domain%get_global_utter_cells_bounds()

			allocate(write_buffer(	global_utter_cells_bounds(1,1):global_utter_cells_bounds(1,2), &
									global_utter_cells_bounds(2,1):global_utter_cells_bounds(2,2), &
									global_utter_cells_bounds(3,1):global_utter_cells_bounds(3,2)))

			write_buffer = 0.0_dkind

			do 

				ind_w = 1
				ind_w(:dimensions,:) = 0

				do proc_i = 0, processor_grid(1)-1
				do proc_j = 0, processor_grid(2)-1
				do proc_k = 0, processor_grid(3)-1
					processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
					
					call problem_domain%set_processor_grid_coord((/proc_i,proc_j,proc_k/))
					call problem_domain%decompose_domain()
					local_utter_cells_bounds	= problem_domain%get_local_utter_cells_bounds()

					allocate(read_buffer(	local_utter_cells_bounds(1,1):local_utter_cells_bounds(1,2), &
											local_utter_cells_bounds(2,1):local_utter_cells_bounds(2,2), &
											local_utter_cells_bounds(3,1):local_utter_cells_bounds(3,2)))

					read(read_io(processor_rank),iostat = error) read_buffer
					
					!print *, read_buffer(:,1,1)

					if (error == 0) then
						
						do dim = 1,dimensions
							proc = proc_i*I_m(1,dim) + proc_j*I_m(2,dim) + proc_k*I_m(3,dim)
							ind_r(dim,1) = local_utter_cells_bounds(dim,1) + ceiling(proc/(proc+1.0))
							ind_r(dim,2) = local_utter_cells_bounds(dim,2) - ceiling((processor_grid(dim)-1-proc)/(processor_grid(dim)-1-proc+1.0))

							ind_w(dim,2) = ind_w(dim,1) + ind_r(dim,2) - ind_r(dim,1) 
						end do

						write_buffer(ind_w(1,1):ind_w(1,2), ind_w(2,1):ind_w(2,2), ind_w(3,1):ind_w(3,2)) = read_buffer(ind_r(1,1):ind_r(1,2), ind_r(2,1):ind_r(2,2), ind_r(3,1):ind_r(3,2))
					end if
					!print *, '     '
					!print *, write_buffer(:,1,1)

					deallocate(read_buffer)
					
					ind_w(3,1) = ind_w(3,2) + I_m(3,dimensions)
				end do
				ind_w(3,1) = I_m(1,dimensions) + I_m(2,dimensions) 
				ind_w(2,1) = ind_w(2,2) + I_m(2,dimensions) + I_m(3,dimensions) 
				end do
				ind_w(2,1) = I_m(1,dimensions)
				ind_w(1,1) = ind_w(1,2) + I_m(1,dimensions) + I_m(2,dimensions) + I_m(3,dimensions) 
				end do

				if (error /= 0) exit
				write(write_io) (((write_buffer(i,j,k),	i = global_utter_cells_bounds(1,1),global_utter_cells_bounds(1,2)), &
														j = global_utter_cells_bounds(2,1),global_utter_cells_bounds(2,2)), &
														k = global_utter_cells_bounds(3,1),global_utter_cells_bounds(3,2))
			!	print *, '     '
			!	print *, write_buffer(:,1,1)
				continue
			end do


		else
			global_inner_faces_bounds		= problem_domain%get_global_inner_faces_bounds()

			allocate(write_buffer(	global_inner_faces_bounds(1,1):global_inner_faces_bounds(1,2), &
									global_inner_faces_bounds(2,1):global_inner_faces_bounds(2,2), &
									global_inner_faces_bounds(3,1):global_inner_faces_bounds(3,2)))

			do 

				ind_w = 1
			!	ind_w(:dimensions,:) = 0

				do proc_i = 0, processor_grid(1)-1
				do proc_j = 0, processor_grid(2)-1
				do proc_k = 0, processor_grid(3)-1
					processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
					
					call problem_domain%set_processor_grid_coord((/proc_i,proc_j,proc_k/))
					call problem_domain%decompose_domain()
					local_inner_faces_bounds	= problem_domain%get_local_inner_faces_bounds()

					allocate(read_buffer(	local_inner_faces_bounds(1,1):local_inner_faces_bounds(1,2), &
											local_inner_faces_bounds(2,1):local_inner_faces_bounds(2,2), &
											local_inner_faces_bounds(3,1):local_inner_faces_bounds(3,2)))

					read(read_io(processor_rank),iostat = error) read_buffer
					
					if (error == 0) then

						do dim = 1,dimensions
							proc = proc_i*I_m(1,dim) + proc_j*I_m(2,dim) + proc_k*I_m(3,dim)
							ind_r(dim,1) = local_inner_faces_bounds(dim,1) !+ ceiling(proc*I_m(1,dim)/(proc+1.0))
							ind_r(dim,2) = local_inner_faces_bounds(dim,2) !- ceiling((processor_grid(dim)-1-proc)/(processor_grid(dim)-1-proc+1.0))

							ind_w(dim,2) = ind_w(dim,1) + ind_r(dim,2) - ind_r(dim,1) 
						end do

						write_buffer(ind_w(1,1):ind_w(1,2), ind_w(2,1):ind_w(2,2), ind_w(3,1):ind_w(3,2)) = read_buffer(ind_r(1,1):ind_r(1,2), ind_r(2,1):ind_r(2,2), ind_r(3,1):ind_r(3,2))
					end if

					deallocate(read_buffer)
				ind_w(3,1) = ind_w(3,2) 
				end do
				ind_w(3,1) = 1 !I_m(1,dimensions) + I_m(2,dimensions) 
				ind_w(2,1) = ind_w(2,2) 
				end do
				ind_w(2,1) = 1 !I_m(1,dimensions)
				ind_w(1,1) = ind_w(1,2) 
				end do

				if (error /= 0) exit
				write(write_io) (((write_buffer(i,j,k),	i = global_inner_faces_bounds(1,1),global_inner_faces_bounds(1,2)), &
														j = global_inner_faces_bounds(2,1),global_inner_faces_bounds(2,2)), &
														k = global_inner_faces_bounds(3,1),global_inner_faces_bounds(3,2))

				
			end do

		end if

		do proc_i = 0, processor_grid(1) -1 
		do proc_j = 0, processor_grid(2) -1
		do proc_k = 0, processor_grid(3) -1
			processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
			close(unit = read_io(processor_rank)) 
		end do
		end do
		end do

		deallocate(write_buffer)
	end do

end program
