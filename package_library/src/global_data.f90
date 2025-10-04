module global_data
    use kind_parameters

    implicit none

#ifdef WIN
	character(len=100)	,parameter	:: fold_sep						= '\'
#else
	character(len=100)	,parameter	:: fold_sep						= '/'
#endif
	
	character(len=100)	,parameter	:: task_setup_folder						= 'task_setup'
	character(len=100)	,parameter	:: chemical_mechanisms_folder				= 'chemical_mechanisms'
	character(len=100)	,parameter	:: thermophysical_data_folder				= 'thermophysical_data'	
    character(len=100)  ,parameter  :: variables_data_file_name					= trim(task_setup_folder) // trim(fold_sep) // 'field_data.inf'
    character(len=100)  ,parameter  :: domain_data_file_name					= trim(task_setup_folder) // trim(fold_sep) // 'domain_data.inf'
	character(len=100)	,parameter	:: chemical_data_file_name					= trim(task_setup_folder) // trim(fold_sep) // 'chemical_data.inf'
	character(len=100)	,parameter	:: thermophysical_data_file_name			= trim(task_setup_folder) // trim(fold_sep) // 'thermophysical_setup.inf'
	character(len=100)	,parameter	:: solver_data_file_name					= trim(task_setup_folder) // trim(fold_sep)	// 'solver_setup.inf'
    character(len=100)  ,parameter  :: post_processor_manager_data_file_name	= trim(task_setup_folder) // trim(fold_sep) // 'post_processors_manager_setup.inf'
	character(len=100)	,parameter	:: post_processor_data_file_name			= trim(task_setup_folder) // trim(fold_sep) // 'post_processor#.inf'
	character(len=100)	,parameter	:: post_processor_data_format				= '.dat'
	character(len=100)	,parameter	:: data_save_data_file_name					= trim(task_setup_folder) // trim(fold_sep) // 'data_save_setup.inf'
	character(len=100)	,parameter	:: data_io_data_file_name					= trim(task_setup_folder) // trim(fold_sep) // 'data_io_setup.inf'
	character(len=100)	,parameter	:: data_io_data_format						= '.dat'
	character(len=100)	,parameter	:: boundary_data_file_name					= trim(task_setup_folder) // trim(fold_sep) // 'boundary_conditions_setup.inf'
	character(len=100)	,parameter	:: problem_setup_log_file					= 'problem_setup.log'

    character(len=10)   ,dimension(3)   ,parameter  :: spatial_axis_names_cartesian   = [character(len=10) :: 'x','y','z']
    character(len=10)   ,dimension(3)   ,parameter  :: spatial_axis_names_cylindrical = [character(len=10) :: 'r','phi','z']
    character(len=10)   ,dimension(3)   ,parameter  :: spatial_axis_names_spherical   = [character(len=10) :: 'r','theta','phi']

    integer             ,dimension(3,3) ,parameter  :: I_m = reshape((/1,0,0,0,1,0,0,0,1/), (/3,3/))

    real(dp)	,parameter  :: pi			= 3.141592653589793115997963468544185161590576171875_dp
    real(dp)	,parameter	:: r_gase_J		= 8.31446261815324_dp
    real(dp)	,parameter	:: r_gase_cal	= 1.98720425864083_dp
    real(dp)	,parameter	:: T_ref        = 298.15_dp
    real(dp)	,parameter	:: P_atm		= 101325.0_dp
	
end module
