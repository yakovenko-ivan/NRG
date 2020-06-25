program computing_module

!	use IFPORT

	use global_data
	use kind_parameters
	use computational_domain_class
	use mpi_communications_class	
	use chemical_properties_class
	use thermophysical_properties_class	
	use data_manager_class
	use data_save_class

	implicit none

	type(computational_domain)					:: problem_domain
	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics	
	type(data_manager)							:: problem_manager
	type(mpi_communications)					:: problem_mpi_support
	type(data_save)								:: problem_data_save

	integer	:: processors_number

	character(len=100)	:: system_command, file_name, data_save_folder, file_path, header_file_name, data_file_name
	character(len=100)	:: save_time

	real(rkind)		,dimension(:,:,:)	,allocatable	:: read_buffer
	real(rkind)		,dimension(:,:,:)	,allocatable	:: write_buffer

	integer		,dimension(3)	:: processor_grid

	integer		:: variables_number

	integer		:: sep_index, argument_index
	integer		:: dir_io, write_io 
	integer	,dimension(:)	,allocatable	:: read_io

	character(len=10)	:: full_command_line
	character(len=10)	:: a_single_argument

	integer	,dimension(3,2)	:: global_inner_cells_bounds ,local_inner_cells_bounds
	integer		:: processor_rank
	integer		:: i, j, k
	integer		:: i_ind, j_ind, k_ind
	integer		:: proc_i, proc_j, proc_k, proc
	integer		:: var = 0
	integer		:: sta
	integer		:: error = 0

	if (command_argument_count() == 0) then
		stop "Error: please, specify number of processors used with proc_num argument"
	end if

	call get_command_argument(2, value = a_single_argument, status = sta)

	print *, trim(a_single_argument)
	read(a_single_argument,*) processors_number

	problem_domain 			= computational_domain_c()
	problem_mpi_support		= mpi_communications_c(problem_domain)
	
	problem_chemistry		= chemical_properties_c()
	problem_thermophysics	= thermophysical_properties_c(problem_chemistry)

	problem_manager			= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics)

	problem_data_save		= data_save_c(problem_manager)

	variables_number = 0
	variables_number = variables_number + problem_domain%get_domain_dimensions()				! Mesh 
	variables_number = variables_number + problem_data_save%get_visible_fields_number()

	call problem_domain%set_mpi_communicator_size(processors_number)
	call problem_domain%generate_processors_grid()

	processor_grid = problem_domain%get_processor_number()

	allocate(read_io(0:processor_grid(1)*processor_grid(2)*processor_grid(3)-1))

	global_inner_cells_bounds		= problem_domain%get_global_inner_cells_bounds()

	allocate(write_buffer(	global_inner_cells_bounds(1,1):global_inner_cells_bounds(1,2), &
							global_inner_cells_bounds(2,1):global_inner_cells_bounds(2,2), &
							global_inner_cells_bounds(3,1):global_inner_cells_bounds(3,2)))

	data_save_folder = trim(problem_data_save%get_data_save_folder())
	write(system_command,'(A,A,A)') 'ls . -d ', trim(data_save_folder) , '/*/ > dir.txt'
	call system(system_command) 

	open(newunit = dir_io, file = 'dir.txt', status = 'old', form = 'formatted')
	read(dir_io,'(A)') file_path

	do 
		read(dir_io,'(A)',iostat = error) file_path
		if (error /= 0) exit

		sep_index = index(file_path,trim(fold_sep))

		save_time = file_path(sep_index+1:)
		sep_index = index(save_time,trim(fold_sep))
		save_time = save_time(:sep_index-1)

		write(header_file_name,'(A,A,A,A,A,A)') trim(data_save_folder) , trim(fold_sep), trim(save_time) ,trim(fold_sep), trim(save_time) , '_header.plt'
		write(file_name,'(A,A,A,A)') trim(data_save_folder) , trim(fold_sep), trim(save_time), '.plt'
		write(system_command,'(A,A,A,A)') 'cp ',trim(header_file_name), '  ', trim(file_name)
		call system(system_command) 
	
		open(newunit = write_io, file = file_name, status = 'old', form = 'binary', position = 'append') 

		do proc_i = 0, processor_grid(1) -1 
		do proc_j = 0, processor_grid(2) -1
		do proc_k = 0, processor_grid(3) -1
			processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
			write(data_file_name,'(A,A,A,A,A,A,I4.4,A)') trim(data_save_folder) , trim(fold_sep), trim(save_time) ,trim(fold_sep), trim(save_time) ,'_proc_', processor_rank, '_data.plt'
			open(newunit = read_io(processor_rank), file = data_file_name, status = 'old', form = 'binary') 
		end do
		end do
		end do

		do 

			i_ind = 1
			j_ind = 1
			k_ind = 1

			do proc_i = 0, processor_grid(1)-1
				do proc_j = 0, processor_grid(2)-1
					do proc_k = 0, processor_grid(3)-1
						processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
						
						call problem_domain%set_processor_grid_coord((/proc_i,proc_j,proc_k/))
						call problem_domain%decompose_domain()
						local_inner_cells_bounds	= problem_domain%get_local_inner_cells_bounds()

						allocate(read_buffer(	local_inner_cells_bounds(1,1):local_inner_cells_bounds(1,2), &
												local_inner_cells_bounds(2,1):local_inner_cells_bounds(2,2), &
												local_inner_cells_bounds(3,1):local_inner_cells_bounds(3,2)))

						read(read_io(processor_rank),iostat = error) read_buffer
						
						!print *, read_buffer

						if (error == 0) then
							write_buffer(i_ind :i_ind + local_inner_cells_bounds(1,2) - local_inner_cells_bounds(1,1), &
										 j_ind :j_ind + local_inner_cells_bounds(2,2) - local_inner_cells_bounds(2,1), &
										 k_ind :k_ind + local_inner_cells_bounds(3,2) - local_inner_cells_bounds(3,1)) = read_buffer
						end if

						!print *, '                '
						!print *, read_buffer(:,1,1)

						deallocate(read_buffer)
						k_ind = k_ind + local_inner_cells_bounds(3,2) - local_inner_cells_bounds(3,1) + 1
					end do
				k_ind = 1
				j_ind = j_ind + local_inner_cells_bounds(2,2) - local_inner_cells_bounds(2,1) + 1
				end do
			j_ind = 1
			i_ind = i_ind + local_inner_cells_bounds(1,2) - local_inner_cells_bounds(1,1) + 1
			end do

			write(write_io) (((write_buffer(i,j,k),	i = global_inner_cells_bounds(1,1),global_inner_cells_bounds(1,2)), &
													j = global_inner_cells_bounds(2,1),global_inner_cells_bounds(2,2)), &
													k = global_inner_cells_bounds(3,1),global_inner_cells_bounds(3,2))

			if (error /= 0) exit
		end do

		close(write_io)
		do proc_i = 0, processor_grid(1) -1 
		do proc_j = 0, processor_grid(2) -1
		do proc_k = 0, processor_grid(3) -1
			processor_rank = (proc_k)   +  processor_grid(3) * (proc_j)   + processor_grid(3) * processor_grid(2) * (proc_i)
			close(read_io(processor_rank))
		end do
		end do
		end do
	end do

end program
