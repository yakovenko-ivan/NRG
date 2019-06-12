module chemical_properties_class

	use kind_parameters
	use global_data

	implicit none

	private
	public  chemical_properties, chemical_properties_pointer, chemical_properties_c

	type 	:: chemical_properties_pointer
		type(chemical_properties)	,pointer	:: chem_ptr
	end type

	type 	:: chemical_properties
		character(len=50)		:: chemical_mechanism_file_name
		real(dkind)				:: default_enhanced_efficiencies
		character(len=15)		:: E_act_units	
	
		integer                                             ::  species_number
		character(len=20)   ,dimension(:)   ,allocatable    ::  species_names
		integer                                             ::  reactions_number
		integer     ,dimension(:)           ,allocatable    ::  reactions_type
		integer     ,dimension(:,:,:)       ,allocatable    ::  chemical_coeffs
		real(dkind) ,dimension(:)           ,allocatable    ::  E_act, A, beta, E_act_low, A_low, beta_low
		real(dkind) ,dimension(:,:)         ,allocatable    ::  enhanced_efficiencies, Troe_coeffs

		character(len=20)                                   ::  activation_energy_units
	contains
		!   procedure   :: get_chemical_species_number
		!   procedure   :: get_chemical_species_names
		!   procedure   :: get_chemical_reactions_number

		! IO Getters/Setters
		procedure	,private	:: read_properties
		procedure	,private	:: write_properties
		procedure	,private	:: set_properties
		
		procedure	,private	:: read_chemical_species_number
		procedure	,private	:: read_chemical_species_names
		procedure   ,private	:: read_chemical_reactions_number
		procedure   ,private	:: read_chemical_reactions_properties
		procedure   ,private	:: chemical_reaction_parser
		
		procedure	:: get_chemical_specie_index
		procedure	:: get_chemical_specie_name

		! Logger
		procedure	:: write_log
		
	end type

	interface   chemical_properties_c
		module procedure   constructor
		module procedure   constructor_file
	end interface

contains

    type(chemical_properties)   function constructor(chemical_mechanism_file_name,default_enhanced_efficiencies,E_act_units)

 		character(len=*)	,intent(in)	:: chemical_mechanism_file_name
		real(dkind)			,intent(in)	:: default_enhanced_efficiencies
		character(len=*)	,intent(in)	:: E_act_units

		integer	:: io_unit
		
		call constructor%set_properties(chemical_mechanism_file_name,default_enhanced_efficiencies,E_act_units)
		
		open(newunit = io_unit, file = chemical_data_file_name, status = 'replace', form = 'formatted')
		call constructor%write_properties(io_unit)
		close(io_unit)	
		
    end function

	type(chemical_properties)   function constructor_file()
		integer	:: io_unit
		
		open(newunit = io_unit, file = chemical_data_file_name, status = 'old', form = 'formatted')
		call constructor_file%read_properties(io_unit)
		close(io_unit)	
    end function
	
	subroutine read_properties(this,chemical_data_unit)
		class(chemical_properties)	,intent(inout)	:: this
		integer						,intent(in)		:: chemical_data_unit
	
		character(len=50)	:: chemical_mechanism_file_name
		real(dkind)			:: default_enhanced_efficiencies
		character(len=15)	:: E_act_units
		
		namelist /chemical_properties/ chemical_mechanism_file_name, default_enhanced_efficiencies, E_act_units
		
		read(unit = chemical_data_unit, nml = chemical_properties)
		call this%set_properties(chemical_mechanism_file_name, default_enhanced_efficiencies, E_act_units)
	end subroutine
	
	subroutine write_properties(this,chemical_data_unit)
		class(chemical_properties)	,intent(in)	:: this
		integer						,intent(in)	:: chemical_data_unit
	
		character(len=50)	:: chemical_mechanism_file_name
		real(dkind)			:: default_enhanced_efficiencies
		character(len=15)	:: E_act_units
		
		namelist /chemical_properties/ chemical_mechanism_file_name, default_enhanced_efficiencies, E_act_units
		
		chemical_mechanism_file_name	= this%chemical_mechanism_file_name
		default_enhanced_efficiencies	= this%default_enhanced_efficiencies
		E_act_units						= this%E_act_units
		
		write(unit = chemical_data_unit, nml = chemical_properties)
	end subroutine
	
	
	subroutine set_properties(this,chemical_mechanism_file_name,default_enhanced_efficiencies,E_act_units)
	
		class(chemical_properties)  ,intent(inout) :: this
		
 		character(len=*)	,intent(in)	:: chemical_mechanism_file_name
		real(dkind)			,intent(in)	:: default_enhanced_efficiencies
		character(len=*)	,intent(in)	:: E_act_units
		
		real(dkind) :: activation_energy_coefficient	
		
		integer :: chemical_mechanism_data_file_unit
		integer :: reaction_number
		integer :: error
		
		this%chemical_mechanism_file_name		=	chemical_mechanism_file_name
		this%default_enhanced_efficiencies		=	default_enhanced_efficiencies
		this%E_act_units						=	E_act_units
		
		open(newunit = chemical_mechanism_data_file_unit, file = trim(task_setup_folder) // trim(fold_sep) // trim(chemical_mechanisms_folder) // trim(fold_sep) // this%chemical_mechanism_file_name, status = 'old', iostat = error)

		this%species_number = this%read_chemical_species_number(chemical_mechanism_data_file_unit)
		allocate(this%species_names(this%species_number+1))
		call this%read_chemical_species_names(chemical_mechanism_data_file_unit)
		this%reactions_number = this%read_chemical_reactions_number(chemical_mechanism_data_file_unit)

		allocate(this%reactions_type(this%reactions_number), &
		this%E_act(this%reactions_number)     , this%A(this%reactions_number)       , this%beta(this%reactions_number)        ,&
		this%E_act_low(this%reactions_number) , this%A_low(this%reactions_number)   , this%beta_low(this%reactions_number)    ,&
		this%enhanced_efficiencies(this%reactions_number,this%species_number), this%Troe_coeffs(this%reactions_number,4)                                      ,&
		this%chemical_coeffs(4,this%reactions_number,2))

		this%reactions_type  = 0
		this%chemical_coeffs = 0
		this%E_act           = 0.0_dkind
		this%A               = 0.0_dkind
		this%beta            = 0.0_dkind
		this%E_act_low       = 0.0_dkind
		this%A_low           = 0.0_dkind
		this%beta_low        = 0.0_dkind
		this%Troe_coeffs     = 0.0_dkind
		this%enhanced_efficiencies   = this%default_enhanced_efficiencies

		call this%read_chemical_reactions_properties(chemical_mechanism_data_file_unit)
		close(chemical_mechanism_data_file_unit)

		select case (this%E_act_units)
			case('J.mol')
				activation_energy_coefficient = 1.0      			! [Eact] = J/mol
			case('kJ.mol')
			    activation_energy_coefficient = 1000.0      		! [Eact] = kJ/mol -> J/mol
			case('cal.mol')
			    activation_energy_coefficient = r_gase_J/r_gase_cal	! [Eact] = cal/mol -> J/mol
		end select

		do reaction_number = 1,this%reactions_number
			if(this%chemical_coeffs(1,reaction_number,1) == 2)	this%A		(reaction_number) = this%A		(reaction_number) * 1.0e-06_dkind		! Second order reaction 10^-6*[sm^3]-> [m^3]
			if(this%chemical_coeffs(1,reaction_number,1) == 3)	this%A		(reaction_number) = this%A		(reaction_number) * 1.0e-12_dkind		! Third order reaction 10^-12*[sm^6]-> [m^6]
			if(this%chemical_coeffs(1,reaction_number,1)+1 == 2)	this%A_low	(reaction_number) = this%A_low	(reaction_number) * 1.0e-06_dkind	
			if(this%chemical_coeffs(1,reaction_number,1)+1 == 3)	this%A_low	(reaction_number) = this%A_low	(reaction_number) * 1.0e-12_dkind
			this%E_act		(reaction_number) = this%E_act		(reaction_number) * activation_energy_coefficient
			this%E_act_low	(reaction_number) = this%E_act_low	(reaction_number) * activation_energy_coefficient
		end do
	end subroutine
	
	subroutine write_log(this,log_unit)
		class(chemical_properties)	,intent(in)	:: this
		integer						,intent(in)	:: log_unit
		
		write(log_unit,'(A)')		'************************************************************************************* '
		write(log_unit,'(A)')		' Chemical data setup : '
		write(log_unit,'(A,A)')		' Chemical mechanism file              : ',	this%chemical_mechanism_file_name
		write(log_unit,'(A,E8.1)')	' Enhanced efficiencies default value  : ',	this%default_enhanced_efficiencies
		write(log_unit,'(A,A)')     ' Activation energy units              : ',	this%E_act_units
		write(log_unit,'(A,I3)')	' Species number                       : ',	this%species_number
		write(log_unit,'(A)')		' Species_names                        : '
		write(log_unit,*)		this%species_names(1:this%species_number)
		write(log_unit,'(A,I4)')	' Reactions_number                     : ',	this%reactions_number
		write(log_unit,'(A)')		'************************************************************************************* '
		
	end subroutine

	
   !  pure integer function get_chemical_species_number(this)
   !      class(chemical_properties)    ,intent(in)           :: this
   !
   !      get_chemical_species_number = this%species_number
   !  end function
   !
   !  pure integer function get_chemical_reactions_number(this)
   !      class(chemical_properties)    ,intent(in)           :: this
   !
   !      get_chemical_reactions_number = this%reactions_number
   !  end function
   !
   !  pure function get_chemical_species_names(this)
   !      class(chemical_properties)    ,intent(in)          :: this
   !      character(len=20)    ,dimension(:)   ,allocatable   :: get_chemical_species_names
   !
   !      allocate(get_chemical_species_names(this%species_number))
   !
   !      get_chemical_species_names = this%species_names
   !  end function

    integer function read_chemical_reactions_number(this,file_unit)
        class(chemical_properties)  ,intent(inout) :: this
		  integer							,intent(in)		:: file_unit

        character(len=150)  :: string

        integer     :: reactions_number
        integer     :: sign_position, third_body_position

        integer     :: r_type
        integer     :: trim_length
        integer     :: cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs
        logical     :: us
        integer     :: error

        reactions_number = 0

        do while (string /= 'REACTIONS')
			read (file_unit,*)	string
        end do

        do
            read(file_unit,*,iostat=error)  string
            if(error /= 0) exit
            string = trim(string)

            third_body_position = index(string,'(+M)')

            if (third_body_position > 0) then

                trim_length = len(trim(string))
                string      = string(1:third_body_position-1) // string(third_body_position+4:trim_length-4)

            end if

            r_type          = 0
            sign_position   = index(string,'=')
            if(sign_position > 0) r_type = 1
            sign_position   = index(string,'=>')
            if(sign_position > 0) r_type = 2
            sign_position   = index(string,'<=')
            if(sign_position > 0) r_type = 3

            select case (r_type)
                case(1)
                    call this%chemical_reaction_parser(string,'=',  cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs, us)
                case(2)
                    call this%chemical_reaction_parser(string,'=>', cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs, us)
                case(3)
                    call this%chemical_reaction_parser(string,'<=', cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs, us)
                case default
                    cycle
            end select

            if (us) then
                cycle
            else
                reactions_number = reactions_number + 1
            end if

        end do

        rewind(file_unit)

        read_chemical_reactions_number = reactions_number
    end function

    subroutine read_chemical_reactions_properties(this, file_unit)
        class(chemical_properties)  ,intent(inout) :: this
		  integer							,intent(in)		:: file_unit

        character(len=150)  :: string

        integer     :: reactions_number
        integer     :: sign_position, third_body_position, fall_off_third_body_position, delimiter_position1, delimiter_position2

        integer     :: r_type, specie_number
        integer     :: trim_length, adjtrim_length
        integer     :: cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs
        logical     :: us, enhanced_efficiency
        integer     :: error

        reactions_number    = 0

        enhanced_efficiency = .false.

        do while (string /= 'REACTIONS')
			read (file_unit,*)	string
        end do

        do
            read(file_unit,*,iostat=error)  string
            if(error /= 0) exit
            string = trim(string)

            third_body_position             = index(string,'+M')

            fall_off_third_body_position    = index(string,'(+M)')

            if (fall_off_third_body_position > 0) then
                trim_length         = len(trim(string))
                string              = string(1:fall_off_third_body_position-1) // string(fall_off_third_body_position+4:trim_length-4)
                enhanced_efficiency = .true.
            end if

            r_type          = 0
            sign_position   = index(string,'=')
            if(sign_position > 0) r_type = 1
            sign_position   = index(string,'=>')
            if(sign_position > 0) r_type = 2
            sign_position   = index(string,'<=')
            if(sign_position > 0) r_type = 3

            select case (r_type)
                case(1)
                    call this%chemical_reaction_parser(string,'=',  cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs, us)
                    if (us) then
                        cycle
                    else
                        reactions_number = reactions_number + 1
                        if( third_body_position > 0) then
                            this%reactions_type(reactions_number) = 1
                            enhanced_efficiency = .true.
                        else
                            this%reactions_type(reactions_number) = 0
                        end if
                    end if
                case(2)
                    call this%chemical_reaction_parser(string,'=>', cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs, us)
                    if (us) then
                        cycle
                    else
                        reactions_number = reactions_number + 1
                        if( third_body_position > 0) then
                            this%reactions_type(reactions_number) = 6
                            enhanced_efficiency = .true.
                        else
                            this%reactions_type(reactions_number) = 5
                        end if
                    end if
                case(3)
                    call this%chemical_reaction_parser(string,'<=', cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3, cc_ls, cc_rs, us)
                    if (us) then
                        cycle
                    else
                        reactions_number = reactions_number + 1
                        if( third_body_position > 0) then
                            this%reactions_type(reactions_number) = 1
                            enhanced_efficiency = .true.
                        else
                            this%reactions_type(reactions_number) = 0
                        end if
                    end if
                case default
                    cycle
            end select

            this%chemical_coeffs(1,reactions_number,1) = cc_ls
            this%chemical_coeffs(2,reactions_number,1) = cc_l1
            this%chemical_coeffs(3,reactions_number,1) = cc_l2
            this%chemical_coeffs(4,reactions_number,1) = cc_l3

            this%chemical_coeffs(1,reactions_number,2) = cc_rs
            this%chemical_coeffs(2,reactions_number,2) = cc_r1
            this%chemical_coeffs(3,reactions_number,2) = cc_r2
            this%chemical_coeffs(4,reactions_number,2) = cc_r3

            backspace(file_unit)

            read(file_unit,*,iostat=error) string, this%A(reactions_number), this%beta(reactions_number), this%E_act(reactions_number)

            if(fall_off_third_body_position > 0) then
                read(file_unit,*,iostat=error) string
                if (string == 'LOW') then
                    this%reactions_type(reactions_number) = 2
                    if (r_type == 2) this%reactions_type(reactions_number) = 7

                    backspace(file_unit)

                    read(file_unit,*,iostat=error) string, this%A_low(reactions_number), this%beta_low(reactions_number), this%E_act_low(reactions_number)
                end if
                read(file_unit,*,iostat=error) string
                if (string == 'TROE') then
                    this%reactions_type(reactions_number) = 3
                    if (r_type == 2) this%reactions_type(reactions_number) = 7

                    backspace(file_unit)

                    read(file_unit,*,iostat=error) string , &
						  											this%Troe_coeffs(reactions_number,1), &
						  											this%Troe_coeffs(reactions_number,2), &
																	this%Troe_coeffs(reactions_number,3), &
																	this%Troe_coeffs(reactions_number,4)

                else
                   backspace(file_unit)
                end if
            end if

            if (error == 0) then

                read(file_unit,*,iostat=error) string
                delimiter_position1 = index(string , '/')
                if ((delimiter_position1 > 0).and.(enhanced_efficiency)) then
                    trim_length     = len(trim(string))
                    adjtrim_length  = len(adjustl(trim(string)))

                    do while (adjtrim_length > 0)

                        delimiter_position1 = index(string,'/')

                        specie_number       = this%get_chemical_specie_index(trim(adjustl(string(1:delimiter_position1-1))))
                        string              = string(delimiter_position1+1:trim_length)
                        delimiter_position2 = index(string,'/')

                        if (specie_number > 0) read(string(1:delimiter_position2-1),*) this%enhanced_efficiencies(reactions_number,specie_number)

                        string          = string(delimiter_position2+1:trim_length)
                        adjtrim_length  = len(adjustl(trim(string)))
                    end do
                else
                    backspace(file_unit)
                end if
            end if

        end do

        rewind(file_unit)

    end subroutine

    subroutine chemical_reaction_parser(this,reaction,sign,spec_num_left1,spec_num_left2,spec_num_left3,&
        spec_num_right1,spec_num_right2,spec_num_right3,sum_chem_coef_left,sum_chem_coef_right,unknown_specie_flag)

        class(chemical_properties)      ,intent(inout)  ::  this
        character(len=*)                ,intent(in)     ::  reaction, sign

        integer             ,intent(out)    ::  spec_num_left1,spec_num_left2,spec_num_left3
        integer             ,intent(out)    ::  spec_num_right1,spec_num_right2,spec_num_right3
        integer             ,intent(out)    ::  sum_chem_coef_left,sum_chem_coef_right

        logical             ,intent(out)    ::  unknown_specie_flag

        character(len=20)   ::  r_left1,r_left2,r_right1,r_right2
        character(len=20)   ::  r_left_spec1,r_left_spec2,r_left_spec3
        character(len=20)   ::  r_right_spec1,r_right_spec2,r_right_spec3

        integer :: pos, ia1, ia2, ia3, ib1, ib2, ib3, ixz, iyz

        ia1=0; ia2=0; ia3=0; ib1=0; ib2=0; ib3=0; ixz=0; iyz=0


	    spec_num_left1  = 0
	    spec_num_left2  = 0
	    spec_num_left3  = 0


	    spec_num_right1 = 0
	    spec_num_right2 = 0
	    spec_num_right3 = 0

	    sum_chem_coef_left  = 0
	    sum_chem_coef_right = 0

	    unknown_specie_flag      = .false.

        pos = index(reaction,sign)

 	    select case(sign)
		    case ('=')
			    r_left1  = reaction(1:pos-1)
			    r_right1 = reaction(pos+1:len(reaction))
		    case ('<=')
			    r_left1  = reaction(1:pos-1)
			    r_right1 = reaction(pos+3:len(reaction))
		    case ('=>')
			    r_left1  = reaction(1:pos-1)
			    r_right1 = reaction(pos+2:len(reaction))
	    end select

	    r_left1  = trim(r_left1)
	    r_right1 = trim(r_right1)

        pos = index(r_left1,'+')

	    if (pos.eq.0) then
		    if (r_left1(1:1).eq.'2') then
			    spec_num_left1 = this%get_chemical_specie_index(r_left1(2:len(r_left1)))     ! 2a=...
			    spec_num_left2 = spec_num_left1
			    sum_chem_coef_left = 2
		    else
			    spec_num_left1 = this%get_chemical_specie_index(r_left1)                     ! a=...
			    sum_chem_coef_left = 1
		    endif
		    if (spec_num_left1.eq.0) unknown_specie_flag = .true.

	    else

		    r_left_spec1 = r_left1(1:pos-1)

		    if (r_left_spec1(1:1).eq.'2') then
			    spec_num_left1 = this%get_chemical_specie_index(r_left_spec1(2:len(r_left_spec1)))   !2a(+b)=... !!or  2a(+b+c)=...  !!m.k. 3
			    spec_num_left2 = spec_num_left1
			    sum_chem_coef_left = 2
			    ixz=1
		    else
			    spec_num_left1 = this%get_chemical_specie_index(r_left_spec1)                        !a(+b)=... or  a(+b+c)=...
			    sum_chem_coef_left = 1
		    endif
		    if (spec_num_left1.eq.0) unknown_specie_flag=.true.

		    r_left2=r_left1(pos+1:)

		    pos = index(r_left2,'+')

		    if (pos.eq.0) then
			    if (r_left2(1:1).eq.'2') then
				    spec_num_left2 = this%get_chemical_specie_index(r_left2(2:))         !(a+)2b=...
				    spec_num_left3 = spec_num_left2
				    sum_chem_coef_left = 3
				    if (spec_num_left2.eq.0) unknown_specie_flag = .true.
			    else
				    if (ixz>0) then
					    spec_num_left3 = this%get_chemical_specie_index(r_left2)         !(2a+)b=...
					    sum_chem_coef_left = 3
					    if (spec_num_left3.eq.0) unknown_specie_flag=.true.
				    else
					    spec_num_left2 = this%get_chemical_specie_index(r_left2)         !(a+)b=...
					    sum_chem_coef_left = 2
					    if (spec_num_left2.eq.0) unknown_specie_flag=.true.
				    endif
			    endif
	        else

			    r_left_spec2=r_left2(1:(pos-1))

			    spec_num_left2 = this%get_chemical_specie_index(r_left_spec2) 	        !(a+)b+c=...
			    if (spec_num_left2.eq.0) unknown_specie_flag=.true.

			    r_left_spec3=r_left2(pos+1:)

			    spec_num_left3 = this%get_chemical_specie_index(r_left_spec3)

			    if (spec_num_left3.eq.0) unknown_specie_flag=.true.
			    sum_chem_coef_left = 3
		    endif

        endif

!!!	REACTION RIGHT SIDE !!!

	    pos = index(r_right1,'+')
	    if (pos.eq.0) then
		    if (r_right1(1:1).eq.'2') then
			    spec_num_right1 = this%get_chemical_specie_index(r_right1(2:len(r_right1)))  ! 2a=...
			    spec_num_right2 = spec_num_right1
			    sum_chem_coef_right = 2
		    else
			    spec_num_right1 = this%get_chemical_specie_index(r_right1)                   ! a=...
			    sum_chem_coef_right = 1
		    endif
		    if (spec_num_right1.eq.0) unknown_specie_flag=.true.
	    else

		    r_right_spec1=r_right1(1:pos-1)

		    if (r_right_spec1(1:1).eq.'2') then
			    spec_num_right1 = this%get_chemical_specie_index(r_right_spec1(2:len(r_right_spec1)))    !2a(+b)=... !!or  2a(+b+c)=...  !!m.k. 3
			    spec_num_right2 = spec_num_right1
			    sum_chem_coef_right = 2
			    iyz=1
		    else
			    spec_num_right1 = this%get_chemical_specie_index(r_right_spec1)                          !a(+b)=... or  a(+b+c)=...
			    sum_chem_coef_right = 1
		    endif
		    if (spec_num_right1.eq.0) unknown_specie_flag=.true.

		    r_right2=r_right1(pos+1:)

		    pos = index(r_right2,'+')

		    if (pos.eq.0) then
			    if (r_right2(1:1).eq.'2') then
				    spec_num_right2 = this%get_chemical_specie_index(r_right2(2:))       !(a+)2b=...
				    spec_num_right3 = spec_num_right2
				    sum_chem_coef_right = 3
				    if (spec_num_right2.eq.0) unknown_specie_flag=.true.
			    else
				    if (iyz>0) then
					    spec_num_right3 = this%get_chemical_specie_index(r_right2)      !(2a+)b=...
					    sum_chem_coef_right = 3
					    if (spec_num_right3.eq.0) unknown_specie_flag=.true.
				    else
					    spec_num_right2 = this%get_chemical_specie_index(r_right2)       !(a+)b=...
					    sum_chem_coef_right = 2
					    if (spec_num_right2.eq.0) unknown_specie_flag=.true.
				    endif
			    endif
	        else

			    r_right_spec2=r_right2(1:(pos-1))

			    spec_num_right2 = this%get_chemical_specie_index(r_right_spec2) 	        !(a+)b+c=...
			    if (spec_num_right2.eq.0) unknown_specie_flag=.true.

			    r_right_spec3=r_right2(pos+1:)

			    spec_num_right3 = this%get_chemical_specie_index(r_right_spec3)

			    if (spec_num_right3.eq.0) unknown_specie_flag=.true.
			    sum_chem_coef_right = 3
		    endif
	    endif
    end subroutine

    integer function get_chemical_specie_index(this,string)

        class(chemical_properties)          ,intent(inout)  ::  this
        character(len=*)                    ,intent(in)     ::  string

        integer ::  species_number
        integer ::  specie_number, i

        specie_number   = 0

	    do i = 1, this%species_number+1

		    if (string == this%species_names(i))	then
			    specie_number = i
			    exit
            endif

        enddo

        get_chemical_specie_index = specie_number

    end function
	
    character function get_chemical_specie_name(this,index)

        class(chemical_properties)		,intent(inout)  ::  this
        integer			                ,intent(in)		::  index

        integer ::  species_number
        integer ::  specie_number, i

        get_chemical_specie_name = this%species_names(index)

    end function	


    integer function read_chemical_species_number(this,file_unit)

		class(chemical_properties)          ,intent(inout)	:: this
		integer										,intent(in)		:: file_unit

		character(len=20)			:: string
		character(len=1)			:: letter, prev_letter

		integer						:: stat
		integer						:: species_number, i, j

        species_number = 0

		do while (string /= 'SPECIES')
			read (file_unit,*)	string
		end do

		reading_species_number : do
			i		= 3
			string	= ' '
			letter	= 'x'

			do while ((letter /= ' ').or.(string(i-2:i-2) == ' '))
				read(file_unit,'(A)',advance = 'no',iostat = stat) letter
				string(i:i) = letter
				i = i + 1
			end do

			string = trim(adjustl(string))

			if (string /= 'END') then
				species_number = species_number + 1
			else
				exit reading_species_number
			end if
        end do reading_species_number

        read_chemical_species_number = species_number

		rewind(file_unit)

    end function

	 subroutine read_chemical_species_names(this,file_unit)

		class(chemical_properties)	,intent(inout)	:: this
		integer							,intent(in)		:: file_unit

		character(len=20)			:: string
		character(len=1)			:: letter ,prev_letter

		integer						:: stat
		integer						:: i ,j

		do while (string /= 'SPECIES')
			read(file_unit,*)	string
		end do

		j = 1
		reading_species_names: do
			i		= 3
			string	= ' '
			letter	= 'x'
			do while ((letter /= ' ').or.(string(i-2:i-2) == ' '))
				read(file_unit,'(A)',advance = 'no',iostat=stat) letter
				string(i:i) = letter
				i = i + 1
			end do

			string = trim(adjustl(string))

			if (string /= 'END') then
				this%species_names(j) = string
				j = j + 1
			else
				exit reading_species_names
			end if
        end do reading_species_names

        rewind(file_unit)

		this%species_names(size(this%species_names)) = 'M'

	end subroutine

end module
