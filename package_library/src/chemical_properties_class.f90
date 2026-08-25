module chemical_properties_class

    use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end

    use kind_parameters
    use global_data

    implicit none

    private
    public :: chemical_properties, chemical_properties_pointer, chemical_properties_c

    integer, parameter :: mechanism_line_length = 1024
    integer, parameter :: parser_message_length = 512

    type :: chemical_properties_pointer
        type(chemical_properties), pointer :: chem_ptr
    end type chemical_properties_pointer

    type :: chemical_properties
        character(len=50) :: chemical_mechanism_file_name
        real(dp) :: default_enhanced_efficiencies
        character(len=15) :: E_act_units

        integer :: species_number
        character(len=20), dimension(:), allocatable :: species_names
        integer :: reactions_number
        integer, dimension(:), allocatable :: reactions_type
        integer, dimension(:,:,:), allocatable :: chemical_coeffs
        real(dp), dimension(:), allocatable :: E_act, A, beta
        real(dp), dimension(:), allocatable :: E_act_low, A_low, beta_low
        real(dp), dimension(:,:), allocatable :: enhanced_efficiencies
        real(dp), dimension(:,:), allocatable :: Troe_coeffs

        character(len=20) :: activation_energy_units

    contains

        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties

        procedure, private :: read_chemical_species_number
        procedure, private :: read_chemical_species_names
        procedure, private :: read_chemical_reactions_number
        procedure, private :: read_chemical_reactions_properties
        procedure, private :: chemical_reaction_parser

        procedure :: get_chemical_specie_index
        procedure :: get_chemical_specie_name
        procedure :: get_chemical_mechanism

        procedure :: write_log

    end type chemical_properties

    interface chemical_properties_c
        module procedure constructor
        module procedure constructor_file
    end interface chemical_properties_c

contains

    type(chemical_properties) function constructor( &
            chemical_mechanism_file_name, default_enhanced_efficiencies, &
            E_act_units)

        character(len=*), intent(in) :: chemical_mechanism_file_name
        real(dp), intent(in) :: default_enhanced_efficiencies
        character(len=*), intent(in) :: E_act_units

        integer :: io_unit, io_status
        character(len=parser_message_length) :: io_message

        call constructor%set_properties(chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units)

        open(newunit=io_unit, file=chemical_data_file_name, status='replace', &
            form='formatted', action='write', delim='quote', &
            iostat=io_status, iomsg=io_message)

        if (io_status /= 0) then
            write(error_unit,'(A)') &
                'Chemical properties: unable to create setup file: '// &
                trim(chemical_data_file_name)
            write(error_unit,'(A)') 'I/O message: '//trim(io_message)
            error stop 'Chemical properties setup-file write failure'
        end if

        call constructor%write_properties(io_unit)
        close(io_unit)

    end function constructor


    type(chemical_properties) function constructor_file()

        integer :: io_unit, io_status
        character(len=parser_message_length) :: io_message

        open(newunit=io_unit, file=chemical_data_file_name, status='old', &
            form='formatted', action='read', iostat=io_status, iomsg=io_message)

        if (io_status /= 0) then
            write(error_unit,'(A)') &
                'Chemical properties: unable to open setup file: '// &
                trim(chemical_data_file_name)
            write(error_unit,'(A)') 'I/O message: '//trim(io_message)
            error stop 'Chemical properties setup-file read failure'
        end if

        call constructor_file%read_properties(io_unit)
        close(io_unit)

    end function constructor_file


    subroutine read_properties(this, chemical_data_unit)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: chemical_data_unit

        character(len=50) :: chemical_mechanism_file_name
        real(dp) :: default_enhanced_efficiencies
        character(len=15) :: E_act_units
        integer :: io_status
        character(len=parser_message_length) :: io_message

        namelist /chemical_properties/ chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units

        read(unit=chemical_data_unit, nml=chemical_properties, &
            iostat=io_status, iomsg=io_message)

        if (io_status /= 0) then
            write(error_unit,'(A)') &
                'Chemical properties: unable to read /chemical_properties/ namelist.'
            write(error_unit,'(A)') 'I/O message: '//trim(io_message)
            error stop 'Chemical properties namelist read failure'
        end if

        call this%set_properties(chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units)

    end subroutine read_properties


    subroutine write_properties(this, chemical_data_unit)

        class(chemical_properties), intent(in) :: this
        integer, intent(in) :: chemical_data_unit

        character(len=50) :: chemical_mechanism_file_name
        real(dp) :: default_enhanced_efficiencies
        character(len=15) :: E_act_units

        namelist /chemical_properties/ chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units

        chemical_mechanism_file_name = this%chemical_mechanism_file_name
        default_enhanced_efficiencies = this%default_enhanced_efficiencies
        E_act_units = this%E_act_units

        write(unit=chemical_data_unit, nml=chemical_properties)

    end subroutine write_properties


    subroutine set_properties(this, chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units)

        class(chemical_properties), intent(inout) :: this

        character(len=*), intent(in) :: chemical_mechanism_file_name
        real(dp), intent(in) :: default_enhanced_efficiencies
        character(len=*), intent(in) :: E_act_units

        real(dp) :: activation_energy_coefficient

        integer :: chemical_mechanism_data_file_unit
        integer :: reaction_number
        integer :: io_status
        logical :: file_exists

        character(len=mechanism_line_length) :: mechanism_path
        character(len=parser_message_length) :: io_message

        this%chemical_mechanism_file_name = trim(chemical_mechanism_file_name)
        this%default_enhanced_efficiencies = default_enhanced_efficiencies
        this%E_act_units = trim(E_act_units)

        mechanism_path = trim(task_setup_folder)//trim(fold_sep)// &
            trim(chemical_mechanisms_folder)//trim(fold_sep)// &
            trim(this%chemical_mechanism_file_name)

        ! Do not allow a failed OPEN to propagate into the parser.  Previously
        ! the IOSTAT value from OPEN was ignored, so an incorrect terminal
        ! working directory could later appear as a misleading "end of record"
        ! error inside the mechanism parser.
        inquire(file=trim(mechanism_path), exist=file_exists)

        if (.not. file_exists) then
            write(error_unit,'(A)') &
                'Chemical properties: mechanism file does not exist.'
            write(error_unit,'(A)') &
                'Resolved relative path: '//trim(mechanism_path)
            write(error_unit,'(A)') &
                'Check the process working directory and task_setup folder.'
            error stop 'Chemical mechanism file not found'
        end if

        open(newunit=chemical_mechanism_data_file_unit, &
            file=trim(mechanism_path), status='old', form='formatted', &
            action='read', position='rewind', iostat=io_status, iomsg=io_message)

        if (io_status /= 0) then
            write(error_unit,'(A)') &
                'Chemical properties: unable to open mechanism file: '// &
                trim(mechanism_path)
            write(error_unit,'(A)') 'I/O message: '//trim(io_message)
            error stop 'Chemical mechanism open failure'
        end if

        this%species_number = &
            this%read_chemical_species_number(chemical_mechanism_data_file_unit)

        if (this%species_number <= 0) then
            error stop 'Chemical mechanism contains no species'
        end if

        allocate(this%species_names(this%species_number+1))
        this%species_names = ''

        call this%read_chemical_species_names(chemical_mechanism_data_file_unit)

        this%reactions_number = &
            this%read_chemical_reactions_number(chemical_mechanism_data_file_unit)

        if (this%reactions_number <= 0) then
            error stop 'Chemical mechanism contains no supported reactions'
        end if

        allocate(this%reactions_type(this%reactions_number), &
            this%E_act(this%reactions_number), &
            this%A(this%reactions_number), &
            this%beta(this%reactions_number), &
            this%E_act_low(this%reactions_number), &
            this%A_low(this%reactions_number), &
            this%beta_low(this%reactions_number), &
            this%enhanced_efficiencies(this%reactions_number, &
                this%species_number), &
            this%Troe_coeffs(this%reactions_number,4), &
            this%chemical_coeffs(4,this%reactions_number,2))

        this%reactions_type = 0
        this%chemical_coeffs = 0
        this%E_act = 0.0_dp
        this%A = 0.0_dp
        this%beta = 0.0_dp
        this%E_act_low = 0.0_dp
        this%A_low = 0.0_dp
        this%beta_low = 0.0_dp
        this%Troe_coeffs = 0.0_dp
        this%enhanced_efficiencies = this%default_enhanced_efficiencies

        call this%read_chemical_reactions_properties( &
            chemical_mechanism_data_file_unit)

        close(chemical_mechanism_data_file_unit)

        select case (trim(this%E_act_units))
        case ('J.mol')
            activation_energy_coefficient = 1.0_dp
        case ('kJ.mol')
            activation_energy_coefficient = 1000.0_dp
        case ('cal.mol')
            activation_energy_coefficient = r_gase_J/r_gase_cal
        case ('kcal.mol')
            activation_energy_coefficient = &
                r_gase_J/r_gase_cal*1000.0_dp
        case default
            write(error_unit,'(A,A)') &
                'Chemical properties: unsupported activation-energy units: ', &
                trim(this%E_act_units)
            error stop 'Unsupported activation-energy units'
        end select

        do reaction_number = 1, this%reactions_number

            if (this%chemical_coeffs(1,reaction_number,1) == 2) then
                this%A(reaction_number) = &
                    this%A(reaction_number)*1.0e-06_dp
            end if

            if (this%chemical_coeffs(1,reaction_number,1) == 3) then
                this%A(reaction_number) = &
                    this%A(reaction_number)*1.0e-12_dp
            end if

            if (this%chemical_coeffs(1,reaction_number,1)+1 == 2) then
                this%A_low(reaction_number) = &
                    this%A_low(reaction_number)*1.0e-06_dp
            end if

            if (this%chemical_coeffs(1,reaction_number,1)+1 == 3) then
                this%A_low(reaction_number) = &
                    this%A_low(reaction_number)*1.0e-12_dp
            end if

            this%E_act(reaction_number) = &
                this%E_act(reaction_number)*activation_energy_coefficient

            this%E_act_low(reaction_number) = &
                this%E_act_low(reaction_number)*activation_energy_coefficient

        end do

    end subroutine set_properties


    subroutine write_log(this, log_unit)

        class(chemical_properties), intent(in) :: this
        integer, intent(in) :: log_unit

        write(log_unit,'(A)') &
            '*************************************************************************************'
        write(log_unit,'(A)') ' Chemical data setup : '
        write(log_unit,'(A,A)') &
            ' Chemical mechanism file              : ', &
            trim(this%chemical_mechanism_file_name)
        write(log_unit,'(A,E8.1)') &
            ' Enhanced efficiencies default value  : ', &
            this%default_enhanced_efficiencies
        write(log_unit,'(A,A)') &
            ' Activation energy units              : ', &
            trim(this%E_act_units)
        write(log_unit,'(A,I3)') &
            ' Species number                       : ', this%species_number
        write(log_unit,'(A)') ' Species_names                        : '
        write(log_unit,*) this%species_names(1:this%species_number)
        write(log_unit,'(A,I4)') &
            ' Reactions_number                     : ', this%reactions_number
        write(log_unit,'(A)') &
            '*************************************************************************************'

    end subroutine write_log


    integer function read_chemical_reactions_number(this, file_unit)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: file_unit

        character(len=150) :: string

        integer :: reactions_number
        integer :: sign_position, third_body_position

        integer :: r_type
        integer :: trim_length
        integer :: cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3
        integer :: cc_ls, cc_rs
        logical :: us
        integer :: io_status

        reactions_number = 0

        call seek_section(file_unit, 'REACTIONS')

        do
            read(file_unit,*,iostat=io_status) string

            if (io_status == iostat_end) exit
            if (io_status /= 0) then
                write(error_unit,'(A)') &
                    'Chemical properties: error while counting reactions.'
                error stop 'Chemical reaction-count parser failure'
            end if

            string = trim(string)

            if (string == 'END') exit

            third_body_position = index(string,'(+M)')

            if (third_body_position > 0) then
                trim_length = len_trim(string)
                if (trim_length > 8) then
                    string = string(1:third_body_position-1)// &
                        string(third_body_position+4:trim_length-4)
                end if
            end if

            r_type = 0

            sign_position = index(string,'=')
            if (sign_position > 0) r_type = 1

            sign_position = index(string,'=>')
            if (sign_position > 0) r_type = 2

            sign_position = index(string,'<=')
            if (sign_position > 0) r_type = 3

            select case (r_type)
            case (1)
                call this%chemical_reaction_parser(string,'=', &
                    cc_l1,cc_l2,cc_l3,cc_r1,cc_r2,cc_r3, &
                    cc_ls,cc_rs,us)
            case (2)
                call this%chemical_reaction_parser(string,'=>', &
                    cc_l1,cc_l2,cc_l3,cc_r1,cc_r2,cc_r3, &
                    cc_ls,cc_rs,us)
            case (3)
                call this%chemical_reaction_parser(string,'<=', &
                    cc_l1,cc_l2,cc_l3,cc_r1,cc_r2,cc_r3, &
                    cc_ls,cc_rs,us)
            case default
                cycle
            end select

            if (.not. us) reactions_number = reactions_number + 1
        end do

        rewind(file_unit)

        read_chemical_reactions_number = reactions_number

    end function read_chemical_reactions_number


    subroutine read_chemical_reactions_properties(this, file_unit)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: file_unit

        character(len=150) :: string

        integer :: reactions_number
        integer :: sign_position
        integer :: third_body_position
        integer :: fall_off_third_body_position
        integer :: delimiter_position1, delimiter_position2

        integer :: r_type, specie_number
        integer :: trim_length, adjtrim_length
        integer :: cc_l1, cc_l2, cc_l3, cc_r1, cc_r2, cc_r3
        integer :: cc_ls, cc_rs
        logical :: us, enhanced_efficiency
        integer :: io_status

        reactions_number = 0
        enhanced_efficiency = .false.

        call seek_section(file_unit, 'REACTIONS')

        do
            read(file_unit,*,iostat=io_status) string

            if (io_status == iostat_end) exit
            if (io_status /= 0) then
                write(error_unit,'(A)') &
                    'Chemical properties: error while reading reactions.'
                error stop 'Chemical reaction parser failure'
            end if

            string = trim(string)

            if (string == 'END') exit

            third_body_position = index(string,'+M')
            fall_off_third_body_position = index(string,'(+M)')

            if (fall_off_third_body_position > 0) then
                trim_length = len_trim(string)
                if (trim_length > 8) then
                    string = string(1:fall_off_third_body_position-1)// &
                        string(fall_off_third_body_position+4:trim_length-4)
                end if
                enhanced_efficiency = .true.
            end if

            r_type = 0

            sign_position = index(string,'=')
            if (sign_position > 0) r_type = 1

            sign_position = index(string,'=>')
            if (sign_position > 0) r_type = 2

            sign_position = index(string,'<=')
            if (sign_position > 0) r_type = 3

            select case (r_type)

            case (1)
                call this%chemical_reaction_parser(string,'=', &
                    cc_l1,cc_l2,cc_l3,cc_r1,cc_r2,cc_r3, &
                    cc_ls,cc_rs,us)

                if (us) cycle

                reactions_number = reactions_number + 1

                if (third_body_position > 0) then
                    this%reactions_type(reactions_number) = 1
                    enhanced_efficiency = .true.
                else
                    this%reactions_type(reactions_number) = 0
                end if

            case (2)
                call this%chemical_reaction_parser(string,'=>', &
                    cc_l1,cc_l2,cc_l3,cc_r1,cc_r2,cc_r3, &
                    cc_ls,cc_rs,us)

                if (us) cycle

                reactions_number = reactions_number + 1

                if (third_body_position > 0) then
                    this%reactions_type(reactions_number) = 6
                    enhanced_efficiency = .true.
                else
                    this%reactions_type(reactions_number) = 5
                end if

            case (3)
                call this%chemical_reaction_parser(string,'<=', &
                    cc_l1,cc_l2,cc_l3,cc_r1,cc_r2,cc_r3, &
                    cc_ls,cc_rs,us)

                if (us) cycle

                reactions_number = reactions_number + 1

                if (third_body_position > 0) then
                    this%reactions_type(reactions_number) = 1
                    enhanced_efficiency = .true.
                else
                    this%reactions_type(reactions_number) = 0
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

            read(file_unit,*,iostat=io_status) string, &
                this%A(reactions_number), &
                this%beta(reactions_number), &
                this%E_act(reactions_number)

            if (io_status /= 0) then
                write(error_unit,'(A,I0)') &
                    'Chemical properties: cannot read Arrhenius data for reaction ', &
                    reactions_number
                error stop 'Chemical Arrhenius parser failure'
            end if

            if (fall_off_third_body_position > 0) then

                read(file_unit,*,iostat=io_status) string

                if (io_status /= 0) then
                    error stop 'Chemical LOW/TROE parser failure'
                end if

                if (string == 'LOW') then
                    this%reactions_type(reactions_number) = 2

                    if (r_type == 2) then
                        this%reactions_type(reactions_number) = 7
                    end if

                    backspace(file_unit)

                    read(file_unit,*,iostat=io_status) string, &
                        this%A_low(reactions_number), &
                        this%beta_low(reactions_number), &
                        this%E_act_low(reactions_number)

                    if (io_status /= 0) then
                        error stop 'Chemical LOW parser failure'
                    end if
                end if

                read(file_unit,*,iostat=io_status) string

                if (io_status /= 0) then
                    error stop 'Chemical TROE parser failure'
                end if

                if (string == 'TROE') then
                    this%reactions_type(reactions_number) = 3

                    if (r_type == 2) then
                        this%reactions_type(reactions_number) = 7
                    end if

                    backspace(file_unit)

                    read(file_unit,*,iostat=io_status) string, &
                        this%Troe_coeffs(reactions_number,1), &
                        this%Troe_coeffs(reactions_number,2), &
                        this%Troe_coeffs(reactions_number,3), &
                        this%Troe_coeffs(reactions_number,4)

                    if (io_status /= 0) then
                        error stop 'Chemical TROE coefficients parser failure'
                    end if
                else
                    backspace(file_unit)
                end if
            end if

            if (io_status == 0) then

                read(file_unit,*,iostat=io_status) string

                if (io_status == iostat_end) exit

                if (io_status /= 0) then
                    error stop 'Chemical enhanced-efficiency parser failure'
                end if

                delimiter_position1 = index(string,'/')

                if (delimiter_position1 > 0 .and. enhanced_efficiency) then

                    trim_length = len_trim(string)
                    adjtrim_length = len_trim(adjustl(string))

                    do while (adjtrim_length > 0)

                        delimiter_position1 = index(string,'/')

                        if (delimiter_position1 <= 1) exit

                        specie_number = this%get_chemical_specie_index( &
                            trim(adjustl(string(1:delimiter_position1-1))))

                        string = string(delimiter_position1+1:trim_length)
                        delimiter_position2 = index(string,'/')

                        if (delimiter_position2 <= 1) then
                            error stop &
                                'Malformed third-body efficiency specification'
                        end if

                        if (specie_number > 0) then
                            read(string(1:delimiter_position2-1),*) &
                                this%enhanced_efficiencies( &
                                    reactions_number,specie_number)
                        end if

                        string = string(delimiter_position2+1:trim_length)
                        adjtrim_length = len_trim(adjustl(string))

                    end do

                else
                    backspace(file_unit)
                end if
            end if

        end do

        rewind(file_unit)

    end subroutine read_chemical_reactions_properties


    subroutine chemical_reaction_parser(this, reaction, sign, &
            spec_num_left1, spec_num_left2, spec_num_left3, &
            spec_num_right1, spec_num_right2, spec_num_right3, &
            sum_chem_coef_left, sum_chem_coef_right, unknown_specie_flag)

        class(chemical_properties), intent(inout) :: this
        character(len=*), intent(in) :: reaction, sign

        integer, intent(out) :: spec_num_left1, spec_num_left2, spec_num_left3
        integer, intent(out) :: spec_num_right1, spec_num_right2, spec_num_right3
        integer, intent(out) :: sum_chem_coef_left, sum_chem_coef_right

        logical, intent(out) :: unknown_specie_flag

        character(len=20) :: r_left1, r_left2, r_right1, r_right2
        character(len=20) :: r_left_spec1, r_left_spec2, r_left_spec3
        character(len=20) :: r_right_spec1, r_right_spec2, r_right_spec3

        integer :: pos, ixz, iyz

        ixz = 0
        iyz = 0

        spec_num_left1 = 0
        spec_num_left2 = 0
        spec_num_left3 = 0

        spec_num_right1 = 0
        spec_num_right2 = 0
        spec_num_right3 = 0

        sum_chem_coef_left = 0
        sum_chem_coef_right = 0

        unknown_specie_flag = .false.

        pos = index(reaction,sign)

        select case (sign)
        case ('=')
            r_left1 = reaction(1:pos-1)
            r_right1 = reaction(pos+1:len(reaction))
        case ('<=')
            r_left1 = reaction(1:pos-1)
            r_right1 = reaction(pos+3:len(reaction))
        case ('=>')
            r_left1 = reaction(1:pos-1)
            r_right1 = reaction(pos+2:len(reaction))
        end select

        r_left1 = trim(r_left1)
        r_right1 = trim(r_right1)

        pos = index(r_left1,'+')

        if (pos == 0) then

            if (r_left1(1:1) == '2') then
                spec_num_left1 = this%get_chemical_specie_index( &
                    r_left1(2:len(r_left1)))
                spec_num_left2 = spec_num_left1
                sum_chem_coef_left = 2
            else
                spec_num_left1 = this%get_chemical_specie_index(r_left1)
                sum_chem_coef_left = 1
            end if

            if (spec_num_left1 == 0) unknown_specie_flag = .true.

        else

            r_left_spec1 = r_left1(1:pos-1)

            if (r_left_spec1(1:1) == '2') then
                spec_num_left1 = this%get_chemical_specie_index( &
                    r_left_spec1(2:len(r_left_spec1)))
                spec_num_left2 = spec_num_left1
                sum_chem_coef_left = 2
                ixz = 1
            else
                spec_num_left1 = &
                    this%get_chemical_specie_index(r_left_spec1)
                sum_chem_coef_left = 1
            end if

            if (spec_num_left1 == 0) unknown_specie_flag = .true.

            r_left2 = r_left1(pos+1:)
            pos = index(r_left2,'+')

            if (pos == 0) then

                if (r_left2(1:1) == '2') then
                    spec_num_left2 = this%get_chemical_specie_index(r_left2(2:))
                    spec_num_left3 = spec_num_left2
                    sum_chem_coef_left = 3
                    if (spec_num_left2 == 0) unknown_specie_flag = .true.
                else
                    if (ixz > 0) then
                        spec_num_left3 = &
                            this%get_chemical_specie_index(r_left2)
                        sum_chem_coef_left = 3
                        if (spec_num_left3 == 0) unknown_specie_flag = .true.
                    else
                        spec_num_left2 = &
                            this%get_chemical_specie_index(r_left2)
                        sum_chem_coef_left = 2
                        if (spec_num_left2 == 0) unknown_specie_flag = .true.
                    end if
                end if

            else

                r_left_spec2 = r_left2(1:pos-1)
                spec_num_left2 = &
                    this%get_chemical_specie_index(r_left_spec2)

                if (spec_num_left2 == 0) unknown_specie_flag = .true.

                r_left_spec3 = r_left2(pos+1:)
                spec_num_left3 = &
                    this%get_chemical_specie_index(r_left_spec3)

                if (spec_num_left3 == 0) unknown_specie_flag = .true.

                sum_chem_coef_left = 3

            end if

        end if

        ! Reaction right side.
        pos = index(r_right1,'+')

        if (pos == 0) then

            if (r_right1(1:1) == '2') then
                spec_num_right1 = this%get_chemical_specie_index( &
                    r_right1(2:len(r_right1)))
                spec_num_right2 = spec_num_right1
                sum_chem_coef_right = 2
            else
                spec_num_right1 = this%get_chemical_specie_index(r_right1)
                sum_chem_coef_right = 1
            end if

            if (spec_num_right1 == 0) unknown_specie_flag = .true.

        else

            r_right_spec1 = r_right1(1:pos-1)

            if (r_right_spec1(1:1) == '2') then
                spec_num_right1 = this%get_chemical_specie_index( &
                    r_right_spec1(2:len(r_right_spec1)))
                spec_num_right2 = spec_num_right1
                sum_chem_coef_right = 2
                iyz = 1
            else
                spec_num_right1 = &
                    this%get_chemical_specie_index(r_right_spec1)
                sum_chem_coef_right = 1
            end if

            if (spec_num_right1 == 0) unknown_specie_flag = .true.

            r_right2 = r_right1(pos+1:)
            pos = index(r_right2,'+')

            if (pos == 0) then

                if (r_right2(1:1) == '2') then
                    spec_num_right2 = &
                        this%get_chemical_specie_index(r_right2(2:))
                    spec_num_right3 = spec_num_right2
                    sum_chem_coef_right = 3
                    if (spec_num_right2 == 0) unknown_specie_flag = .true.
                else
                    if (iyz > 0) then
                        spec_num_right3 = &
                            this%get_chemical_specie_index(r_right2)
                        sum_chem_coef_right = 3
                        if (spec_num_right3 == 0) unknown_specie_flag = .true.
                    else
                        spec_num_right2 = &
                            this%get_chemical_specie_index(r_right2)
                        sum_chem_coef_right = 2
                        if (spec_num_right2 == 0) unknown_specie_flag = .true.
                    end if
                end if

            else

                r_right_spec2 = r_right2(1:pos-1)
                spec_num_right2 = &
                    this%get_chemical_specie_index(r_right_spec2)

                if (spec_num_right2 == 0) unknown_specie_flag = .true.

                r_right_spec3 = r_right2(pos+1:)
                spec_num_right3 = &
                    this%get_chemical_specie_index(r_right_spec3)

                if (spec_num_right3 == 0) unknown_specie_flag = .true.

                sum_chem_coef_right = 3

            end if

        end if

    end subroutine chemical_reaction_parser


    integer function get_chemical_specie_index(this, string)

        class(chemical_properties), intent(inout) :: this
        character(len=*), intent(in) :: string

        integer :: specie_number, i

        specie_number = 0

        do i = 1, this%species_number+1
            if (trim(string) == trim(this%species_names(i))) then
                specie_number = i
                exit
            end if
        end do

        get_chemical_specie_index = specie_number

    end function get_chemical_specie_index


    character function get_chemical_specie_name(this, index)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: index

        get_chemical_specie_name = this%species_names(index)

    end function get_chemical_specie_name


    character(len=20) function get_chemical_mechanism(this)

        class(chemical_properties), intent(inout) :: this

        character(len=20) :: mechanism_name
        integer :: index_type

        index_type = index(this%chemical_mechanism_file_name,'.txt')

        if (index_type > 1) then
            mechanism_name = this%chemical_mechanism_file_name(:index_type-1)
        else
            mechanism_name = trim(this%chemical_mechanism_file_name)
        end if

        get_chemical_mechanism = trim(mechanism_name)

    end function get_chemical_mechanism


    integer function read_chemical_species_number(this, file_unit)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: file_unit

        integer :: parsed_species_number

        call scan_species_section(file_unit, parsed_species_number)

        read_chemical_species_number = parsed_species_number

    end function read_chemical_species_number


    subroutine read_chemical_species_names(this, file_unit)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: file_unit

        integer :: parsed_species_number

        call scan_species_section(file_unit, parsed_species_number, &
            this%species_names(1:this%species_number))

        if (parsed_species_number /= this%species_number) then
            write(error_unit,'(A,I0,A,I0)') &
                'Chemical properties: species count changed between scans: ', &
                this%species_number, ' -> ', parsed_species_number
            error stop 'Inconsistent SPECIES section'
        end if

        this%species_names(this%species_number+1) = 'M'

    end subroutine read_chemical_species_names


    subroutine scan_species_section(file_unit, species_count, species_names)

        integer, intent(in) :: file_unit
        integer, intent(out) :: species_count
        character(len=20), intent(out), optional :: species_names(:)

        character(len=mechanism_line_length) :: line
        character(len=64) :: token
        character(len=parser_message_length) :: io_message

        integer :: io_status, cursor
        logical :: in_species_section, token_found, section_found

        species_count = 0
        in_species_section = .false.
        section_found = .false.

        rewind(file_unit)

        do

            read(file_unit,'(A)',iostat=io_status,iomsg=io_message) line

            if (io_status == iostat_end) exit

            if (io_status /= 0) then
                write(error_unit,'(A)') &
                    'Chemical properties: error while scanning SPECIES section.'
                write(error_unit,'(A)') 'I/O message: '//trim(io_message)
                error stop 'SPECIES parser read failure'
            end if

            call remove_chemkin_comment(line)

            cursor = 1

            do
                call next_whitespace_token(line, cursor, token, token_found)

                if (.not. token_found) exit

                if (.not. in_species_section) then
                    if (trim(token) == 'SPECIES') then
                        in_species_section = .true.
                        section_found = .true.
                    end if
                    cycle
                end if

                if (trim(token) == 'END') then
                    rewind(file_unit)
                    return
                end if

                species_count = species_count + 1

                if (present(species_names)) then

                    if (species_count > size(species_names)) then
                        write(error_unit,'(A)') &
                            'Chemical properties: too many species while rereading section.'
                        error stop 'SPECIES array overflow'
                    end if

                    if (len_trim(token) > len(species_names(species_count))) then
                        write(error_unit,'(A,A)') &
                            'Chemical properties: species name is too long: ', &
                            trim(token)
                        error stop 'Chemical species name too long'
                    end if

                    species_names(species_count) = trim(token)

                end if

            end do

        end do

        rewind(file_unit)

        if (.not. section_found) then
            error stop 'Chemical properties: SPECIES section was not found'
        end if

        error stop 'Chemical properties: SPECIES section has no END marker'

    end subroutine scan_species_section


    subroutine seek_section(file_unit, section_name)

        integer, intent(in) :: file_unit
        character(len=*), intent(in) :: section_name

        character(len=mechanism_line_length) :: line
        character(len=64) :: token
        character(len=parser_message_length) :: io_message

        integer :: io_status, cursor
        logical :: token_found

        rewind(file_unit)

        do
            read(file_unit,'(A)',iostat=io_status,iomsg=io_message) line

            if (io_status == iostat_end) exit

            if (io_status /= 0) then
                write(error_unit,'(A,A)') &
                    'Chemical properties: error while searching section ', &
                    trim(section_name)
                write(error_unit,'(A)') 'I/O message: '//trim(io_message)
                error stop 'Chemical mechanism section search failure'
            end if

            call remove_chemkin_comment(line)

            cursor = 1
            call next_whitespace_token(line, cursor, token, token_found)

            if (token_found) then
                if (trim(token) == trim(section_name)) return
            end if
        end do

        write(error_unit,'(A,A)') &
            'Chemical properties: section not found: ', trim(section_name)
        error stop 'Chemical mechanism section not found'

    end subroutine seek_section


    subroutine remove_chemkin_comment(line)

        character(len=*), intent(inout) :: line
        integer :: comment_position

        comment_position = index(line,'!')

        if (comment_position == 1) then
            line = ''
        else if (comment_position > 1) then
            line(comment_position:) = ' '
        end if

    end subroutine remove_chemkin_comment


    subroutine next_whitespace_token(line, cursor, token, found)

        character(len=*), intent(in) :: line
        integer, intent(inout) :: cursor
        character(len=*), intent(out) :: token
        logical, intent(out) :: found

        integer :: start_position, end_position, usable_length

        token = ''
        found = .false.
        usable_length = len_trim(line)

        do while (cursor <= usable_length)
            if (.not. is_parser_whitespace(line(cursor:cursor))) exit
            cursor = cursor + 1
        end do

        if (cursor > usable_length) return

        start_position = cursor

        do while (cursor <= usable_length)
            if (is_parser_whitespace(line(cursor:cursor))) exit
            cursor = cursor + 1
        end do

        end_position = cursor-1

        if (end_position < start_position) return

        if (end_position-start_position+1 > len(token)) then
            write(error_unit,'(A)') &
                'Chemical properties: token in mechanism file is too long.'
            error stop 'Chemical mechanism token too long'
        end if

        token = line(start_position:end_position)
        found = .true.

    end subroutine next_whitespace_token


    pure logical function is_parser_whitespace(character_value)

        character(len=1), intent(in) :: character_value

        is_parser_whitespace = &
            character_value == ' ' .or. &
            character_value == achar(9) .or. &
            character_value == achar(13)

    end function is_parser_whitespace

end module chemical_properties_class
