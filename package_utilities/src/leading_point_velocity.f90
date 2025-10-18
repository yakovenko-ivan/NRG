program Leading_point
    
	use IFPORT

	use global_data
	use kind_parameters
	use computational_domain_class
	use mpi_communications_class	
	use chemical_properties_class
	use thermophysical_properties_class	
    use solver_options_class
	use data_manager_class
	use data_save_class

	implicit none
    
	type(computational_domain)					:: problem_domain
	type(chemical_properties)			,target	:: problem_chemistry
	type(thermophysical_properties)		,target	:: problem_thermophysics
	type(solver_options)						:: problem_solver_options
	type(data_manager)							:: problem_data_manager
	type(mpi_communications)					:: problem_mpi_support
	type(data_save)								:: problem_data_save
    
    integer					:: CONVERTSTR3

	character(len=200)	:: system_command, file_name, data_save_folder, file_path, header_file_name, data_file_name
	integer	            :: save_time
    character(len=20)	:: cstr
    real(dp)         :: data_save_time

    real				:: header(3000)		!1399 !976!1493	 !1399 !1587				!# Header служит для считывания заголовков файлов с данными. Размер заголовка может меняться в зависимости от постановки, такчто 
	integer				:: num_of_strings
	real				:: string														!# величину header приходится каждый раз подбирать заново. Для этого можно выставить заведомо большое количество элементов массива,
																						!# например header(1000), и сделать предваритьельный прогон считывания заголовка. Далее внизу помечен блок программы (стр. 171), 
																						!# где нужно сделать точку останова и посмотреть значения в массиве, чтобы определить сколько элементов должно быть в действительности 
																						!# (см. комментарии к стр. 171). После этого, завышенное значение кол-ва элементов здесь нужно изменить на реальное и перезапустить программу.  

	real(sp)		,dimension(:,:,:)	,allocatable	:: read_buffer, read_buffer2, read_buffer3

    integer ,parameter		                        :: variables_number_in_file = 19!29!20!31	!# Сколько переменных было в файлах с расчетами
	integer ,parameter		                        :: variables_number = 11			!# Сколько переменных будет в файле с x-t разверткой, названия переменных см. стр. 88. Первая переменная - всегда время! 

    character(len=20) ,dimension(variables_number)	:: variables_names

    integer	,dimension(40)	:: istr
    
	integer     							:: file_number
	integer									:: write_unit, write_unit2
	integer	,dimension(:)	,allocatable	:: read_units
	integer									:: dir_io
	
	integer	,dimension(3,2)	:: global_utter_cells_bounds, global_inner_cells_bounds
	integer	:: i, j, k, n, variable_counter, file_counter, lp_index,lp_index2,lp_index3
	real	:: r

	integer					        :: sta	
	
	real					        :: min_grad, max_grad_x, max_grad_y, grad_x, grad_y
	integer					        :: lp_number, lp_neighbour, lp_copies, lp_tip, lp_bound
	real(dp)	,dimension(2000,2)	:: lp_coord
	real(dp)	,dimension(2)		:: lp_copy
	real(dp)				        :: tip_coord_y
    real(dp)				        :: tip_coord_x
	real(dp)				        :: cell_size
	
	
	problem_domain 			= computational_domain_c()
	problem_mpi_support		= mpi_communications_c(problem_domain)
	
	problem_chemistry		= chemical_properties_c()
	problem_thermophysics	= thermophysical_properties_c(problem_chemistry)

    problem_solver_options  = solver_options_c()	
    
	problem_data_manager	= data_manager_c(problem_domain,problem_mpi_support,problem_chemistry,problem_thermophysics,problem_solver_options)

	problem_data_save		= data_save_c(problem_data_manager)
	
	global_inner_cells_bounds	= problem_domain%get_global_inner_cells_bounds()	
	
	cell_size = 2e-05_dp
	
	data_save_folder = trim(problem_data_save%get_data_save_folder())
!	write(system_command,'(A,A,A)') 'ls ', trim(data_save_folder) , ' | wc -l > dir.txt'				!# Linux
	write(system_command,'(A,A,A)') trim(data_save_folder),'\script.bat'								!# Windows. Необходимо наличие скрипт файла script.bat в папке data_save !!
	call execute_command_line(system_command)															!#          script.bat считает количество *.plt файлов в директории
	open(newunit = dir_io, file = 'dir.txt', status = 'old', form = 'formatted')
	read(dir_io,'(I5)') file_number	
	close(dir_io)
	
!	write(system_command,'(A,A,A)') 'ls ', trim(data_save_folder) , ' > dir.txt'						!# Linux
	write(system_command,'(A,A,A)') 'del dir.txt && dir /B ', trim(data_save_folder) , ' >> dir.txt'	!# Windows
	call execute_command_line(system_command) 
	open(newunit = dir_io, file = 'dir.txt', status = 'old', form = 'formatted')
	
	allocate(read_buffer(	global_inner_cells_bounds(1,1):global_inner_cells_bounds(1,2), &
						    global_inner_cells_bounds(2,1):global_inner_cells_bounds(2,2), &
						    global_inner_cells_bounds(3,1):global_inner_cells_bounds(3,2)))
						
!	allocate(tip_coord_x(	global_inner_cells_bounds(2,1):global_inner_cells_bounds(2,2)))
                            
	allocate(read_units(file_number))
	
	write(file_name,'(A)') 'lp_coord(t).dat' 
	open(newunit=write_unit, file=file_name, status='replace')
	
	write(write_unit,'(A)') 'VARIABLES="time" "lp_x"' 

	num_of_strings = 0
	
	do file_counter = 1, file_number 

		read(dir_io,'(A)') file_path
		write(file_name,'(A,A,A)') trim(data_save_folder),trim(fold_sep),trim(file_path)
		open(newunit=read_units(file_counter),file=file_name,status='old',form='binary')
		
		if(file_counter == 1) then 
			do while ((abs(string - cell_size/2) > cell_size/100))
				read(read_units(file_counter)) string
				num_of_strings = num_of_strings + 1
                if (isnan(string) == .true.) string = 0.0
			end do
			rewind(read_units(file_counter))
		end if
		
		read(read_units(file_counter)) header(1:num_of_strings-1)                           !# Считывание заголовка. После него на предварительном пуске нужно ставить точку останова и смотреть где в считанном массиве начнут идти данные расчетной сетки.
																							!# Первыми в файлах с данными идут координаты x, так что массив header должен заканчиваться до считывания координаты x первой ячейки (1,1,1). 
																							!# Именно номер элемента перед элементом с координатой x ячейки (1,1,1) должен стать общим числом элементов при объявлении массива header (стр. 30). 
		
		read(file_path(:6),*)  save_time
	
		data_save_time			= save_time/problem_data_save%get_save_time_coefficient()
		
		
		do variable_counter = 1, 4
			read(read_units(file_counter),iostat = sta) (((read_buffer(i,j,k),  i = global_inner_cells_bounds(1,1),global_inner_cells_bounds(1,2)), &
																				j = global_inner_cells_bounds(2,1),global_inner_cells_bounds(2,2)), &
																				k = global_inner_cells_bounds(3,1),global_inner_cells_bounds(3,2))			
		end do
																			
        max_grad_x = 0.0
        max_grad_y = 0.0

        tip_coord_y	= 0.0
        tip_coord_x = 0.0
        
		do i = global_inner_cells_bounds(1,1)+1, global_inner_cells_bounds(1,2)-1
		do j = global_inner_cells_bounds(2,1)+1, global_inner_cells_bounds(2,2)-1

            if(j==25) then
				grad_x = (read_buffer(i-1,j,1) - read_buffer(i+1,j,1))
				if ((grad_x > max_grad_x))	then
					max_grad_x = grad_x
					tip_coord_x = i*cell_size - 0.5*cell_size
				end if
            end if
            
		end do
		end do


		write(write_unit,'(3E14.7)') data_save_time, tip_coord_x

	end do

end program 
    
    
INTEGER FUNCTION CONVERTSTR3(pistr,pcstr) RESULT (Num)
INTEGER*4 pistr(20)
CHARACTER*20 pcstr

i=1
DO WHILE (ICHAR(pcstr(i:i)).NE.0)
  pistr(i)=ICHAR(pcstr(i:i))
  i=i+1
ENDDO
pistr(i)=0
Num=i
END