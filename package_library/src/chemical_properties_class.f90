module chemical_properties_class

    use, intrinsic :: iso_fortran_env, only: error_unit, iostat_end

    use kind_parameters
    use global_data

    implicit none

    private
    public :: chemical_properties, chemical_properties_pointer, chemical_properties_c

    integer, parameter :: mechanism_line_length = 4096
    integer, parameter :: parser_message_length = 512
    integer, parameter :: maximum_excluded_species = 64

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
        ! A reduced in-memory mechanism; the source CHEMKIN file is unchanged.
        character(len=20), allocatable :: excluded_species(:)
        integer :: original_species_number = 0
        integer :: original_reactions_number = 0
        integer, allocatable :: source_reaction_indices(:)

    contains

        procedure, private :: read_properties
        procedure, private :: write_properties
        procedure, private :: set_properties

        procedure, private :: read_chemical_species_number
        procedure, private :: read_chemical_species_names
        procedure, private :: read_chemical_reactions_number
        procedure, private :: read_chemical_reactions_properties
        procedure, private :: chemical_reaction_parser
        procedure, private :: count_standard_chemkin_reactions
        procedure, private :: read_standard_chemkin_reactions
        procedure, private :: reduce_mechanism

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
            E_act_units, excluded_species)

        character(len=*), intent(in) :: chemical_mechanism_file_name
        real(dp), intent(in) :: default_enhanced_efficiencies
        character(len=*), intent(in) :: E_act_units
        character(len=*), intent(in), optional :: excluded_species(:)

        integer :: io_unit, io_status
        character(len=parser_message_length) :: io_message

        call constructor%set_properties(chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units, excluded_species)

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
        integer :: excluded_species_count
        character(len=20) :: excluded_species(maximum_excluded_species)
        integer :: io_status
        character(len=parser_message_length) :: io_message

        namelist /chemical_properties/ chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units, &
            excluded_species_count, excluded_species

        excluded_species_count = 0
        excluded_species = ''
        read(unit=chemical_data_unit, nml=chemical_properties, &
            iostat=io_status, iomsg=io_message)

        if (io_status /= 0) then
            write(error_unit,'(A)') &
                'Chemical properties: unable to read /chemical_properties/ namelist.'
            write(error_unit,'(A)') 'I/O message: '//trim(io_message)
            error stop 'Chemical properties namelist read failure'
        end if

        if (excluded_species_count < 0 .or. &
            excluded_species_count > maximum_excluded_species) then
            error stop 'CHEMKIN: invalid excluded_species_count in setup file'
        end if
        call this%set_properties(chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units, &
            excluded_species(1:excluded_species_count))

    end subroutine read_properties


    subroutine write_properties(this, chemical_data_unit)

        class(chemical_properties), intent(in) :: this
        integer, intent(in) :: chemical_data_unit

        character(len=50) :: chemical_mechanism_file_name
        real(dp) :: default_enhanced_efficiencies
        character(len=15) :: E_act_units
        integer :: excluded_species_count
        character(len=20) :: excluded_species(maximum_excluded_species)

        namelist /chemical_properties/ chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units, &
            excluded_species_count, excluded_species

        chemical_mechanism_file_name = this%chemical_mechanism_file_name
        default_enhanced_efficiencies = this%default_enhanced_efficiencies
        E_act_units = this%E_act_units
        excluded_species_count = size(this%excluded_species)
        excluded_species = ''
        if (excluded_species_count > 0) &
            excluded_species(1:excluded_species_count) = this%excluded_species

        write(unit=chemical_data_unit, nml=chemical_properties)

    end subroutine write_properties


    subroutine set_properties(this, chemical_mechanism_file_name, &
            default_enhanced_efficiencies, E_act_units, excluded_species)

        class(chemical_properties), intent(inout) :: this

        character(len=*), intent(in) :: chemical_mechanism_file_name
        real(dp), intent(in) :: default_enhanced_efficiencies
        character(len=*), intent(in) :: E_act_units
        character(len=*), intent(in), optional :: excluded_species(:)

        real(dp) :: activation_energy_coefficient

        integer :: chemical_mechanism_data_file_unit
        integer :: reaction_number
        integer :: io_status
        logical :: file_exists
        character(len=20) :: amount_units
        real(dp) :: amount_factor

        character(len=mechanism_line_length) :: mechanism_path
        character(len=parser_message_length) :: io_message

        this%chemical_mechanism_file_name = trim(chemical_mechanism_file_name)
        this%default_enhanced_efficiencies = default_enhanced_efficiencies
        this%E_act_units = trim(E_act_units)
        if (present(excluded_species)) then
            if (size(excluded_species) > maximum_excluded_species) &
                error stop 'CHEMKIN: too many excluded species'
            this%excluded_species = excluded_species
        else
            allocate(this%excluded_species(0))
        end if

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

        ! Only the standard CHEMKIN grammar is supported.  File extensions
        ! do not select an alternative reader.
        amount_units='MOLES'

        this%species_number = &
            this%read_chemical_species_number(chemical_mechanism_data_file_unit)

        if (this%species_number <= 0) then
            error stop 'Chemical mechanism contains no species'
        end if

        allocate(this%species_names(this%species_number+1))
        this%species_names = ''

        call this%read_chemical_species_names(chemical_mechanism_data_file_unit)

        this%reactions_number = &
            this%count_standard_chemkin_reactions(chemical_mechanism_data_file_unit)

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

        call this%read_standard_chemkin_reactions( &
            chemical_mechanism_data_file_unit,amount_units)

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
        case ('kelvins')
            activation_energy_coefficient = r_gase_J
        case default
            write(error_unit,'(A,A)') &
                'Chemical properties: unsupported activation-energy units: ', &
                trim(this%E_act_units)
            error stop 'Unsupported activation-energy units'
        end select

        amount_factor=1.0_dp
        if (trim(amount_units)=='MOLECULES') &
            amount_factor=N_A

        do reaction_number = 1, this%reactions_number

            ! CHEMKIN MOLECULES: rates in particles/cm3/s; native NRG uses
            ! mol/cm3/s before the existing cm3 -> m3 conversion below.
            if (amount_factor /= 1.0_dp) then
                this%A(reaction_number)=this%A(reaction_number)* &
                    amount_factor**(this%chemical_coeffs(1,reaction_number,1)-1)
                if (this%reactions_type(reaction_number)==2 .or. &
                    this%reactions_type(reaction_number)==3 .or. &
                    this%reactions_type(reaction_number)==7) then
                    this%A_low(reaction_number)=this%A_low(reaction_number)* &
                        amount_factor**this%chemical_coeffs(1,reaction_number,1)
                end if
            end if

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

        this%original_species_number = this%species_number
        this%original_reactions_number = this%reactions_number
        allocate(this%source_reaction_indices(this%reactions_number))
        do reaction_number = 1, this%reactions_number
            this%source_reaction_indices(reaction_number) = reaction_number
        end do
        if (size(this%excluded_species) > 0) call this%reduce_mechanism()

    end subroutine set_properties


    ! Operate only on fully validated, fully parsed and SI-converted arrays.
    ! The source file remains standard CHEMKIN and is never modified.
    subroutine reduce_mechanism(this)
        class(chemical_properties), intent(inout) :: this
        integer :: ns, nr, i, j, r, side, part, idx, nactive, nkept, new_r
        integer, allocatable :: mapping(:), reaction_ids(:)
        integer, allocatable :: new_type(:), new_coeffs(:,:,:)
        logical, allocatable :: keep(:)
        character(len=20), allocatable :: new_names(:)
        real(dp), allocatable :: new_a(:), new_beta(:), new_e(:)
        real(dp), allocatable :: new_alow(:), new_blow(:), new_elow(:)
        real(dp), allocatable :: new_eff(:,:), new_troe(:,:)

        ns = this%species_number
        nr = this%reactions_number
        allocate(mapping(ns+1), keep(nr))
        mapping = 0
        do i = 1, size(this%excluded_species)
            if (len_trim(this%excluded_species(i)) == 0 .or. &
                trim(this%excluded_species(i)) == 'M') &
                error stop 'CHEMKIN: blank or pseudo-species exclusion'
            idx = this%get_chemical_specie_index(trim(this%excluded_species(i)))
            if (idx < 1 .or. idx > ns) then
                write(error_unit,'(A,A)') &
                    'CHEMKIN: excluded species not found in canonical mechanism: ', &
                    trim(this%excluded_species(i))
                error stop 'CHEMKIN: invalid excluded species'
            end if
            do j = 1, i-1
                if (trim(this%excluded_species(i)) == &
                    trim(this%excluded_species(j))) &
                    error stop 'CHEMKIN: duplicate excluded species'
            end do
            mapping(idx) = -1
        end do

        nactive = 0
        do i = 1, ns
            if (mapping(i) < 0) cycle
            nactive = nactive + 1
            mapping(i) = nactive
        end do
        if (nactive == 0) error stop 'CHEMKIN: no active species after exclusion'
        ! The generic collider M is a virtual species at ns+1.
        mapping(ns+1) = nactive+1

        keep = .true.
        do r = 1, nr
            do side = 1, 2
                do part = 2, this%chemical_coeffs(1,r,side)+1
                    idx = this%chemical_coeffs(part,r,side)
                    if (idx < 1 .or. idx > ns+1) &
                        error stop 'CHEMKIN: invalid parsed reaction participant'
                    if (mapping(idx) < 1) keep(r) = .false.
                end do
            end do
        end do
        nkept = count(keep)
        if (nkept == 0) error stop 'CHEMKIN: exclusion removes all reactions'

        allocate(new_names(nactive+1), new_type(nkept), &
            new_coeffs(4,nkept,2), reaction_ids(nkept), &
            new_a(nkept), new_beta(nkept), new_e(nkept), &
            new_alow(nkept), new_blow(nkept), new_elow(nkept), &
            new_eff(nkept,nactive), new_troe(nkept,4))
        do i = 1, ns
            if (mapping(i) > 0) &
                new_names(mapping(i)) = this%species_names(i)
        end do
        new_names(nactive+1) = 'M'
        new_r = 0
        do r = 1, nr
            if (.not. keep(r)) cycle
            new_r = new_r+1
            reaction_ids(new_r) = this%source_reaction_indices(r)
            new_type(new_r) = this%reactions_type(r)
            new_coeffs(:,new_r,:) = this%chemical_coeffs(:,r,:)
            do side = 1, 2
                do part = 2, new_coeffs(1,new_r,side)+1
                    idx = new_coeffs(part,new_r,side)
                    new_coeffs(part,new_r,side) = mapping(idx)
                end do
            end do
            new_a(new_r) = this%A(r)
            new_beta(new_r) = this%beta(r)
            new_e(new_r) = this%E_act(r)
            new_alow(new_r) = this%A_low(r)
            new_blow(new_r) = this%beta_low(r)
            new_elow(new_r) = this%E_act_low(r)
            new_troe(new_r,:) = this%Troe_coeffs(r,:)
            do i = 1, ns
                if (mapping(i) > 0) &
                    new_eff(new_r,mapping(i)) = &
                        this%enhanced_efficiencies(r,i)
            end do
        end do

        call move_alloc(new_names,this%species_names)
        call move_alloc(new_type,this%reactions_type)
        call move_alloc(new_coeffs,this%chemical_coeffs)
        call move_alloc(new_a,this%A)
        call move_alloc(new_beta,this%beta)
        call move_alloc(new_e,this%E_act)
        call move_alloc(new_alow,this%A_low)
        call move_alloc(new_blow,this%beta_low)
        call move_alloc(new_elow,this%E_act_low)
        call move_alloc(new_eff,this%enhanced_efficiencies)
        call move_alloc(new_troe,this%Troe_coeffs)
        call move_alloc(reaction_ids,this%source_reaction_indices)
        this%species_number = nactive
        this%reactions_number = nkept
    end subroutine reduce_mechanism

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
        write(log_unit,'(A,I0,A,I0)') &
            ' Source/reduced species               : ', &
            this%original_species_number, ' / ', this%species_number
        write(log_unit,'(A,I0,A,I0)') &
            ' Source/reduced reactions             : ', &
            this%original_reactions_number, ' / ', this%reactions_number
        if (size(this%excluded_species) > 0) then
            write(log_unit,'(A)',advance='no') ' Excluded species                      : '
            write(log_unit,*) this%excluded_species
            write(log_unit,'(A)') &
                ' Retained source reaction indices:'
            write(log_unit,*) this%source_reaction_indices
        end if
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


    character(len=20) function get_chemical_specie_name(this, index)

        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: index

        get_chemical_specie_name = this%species_names(index)

    end function get_chemical_specie_name


    character(len=20) function get_chemical_mechanism(this)

        class(chemical_properties), intent(inout) :: this

        character(len=20) :: mechanism_name
        integer :: index_type

        index_type = scan(trim(this%chemical_mechanism_file_name),'.',back=.true.)

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

    ! Read a standard CHEMKIN REACTIONS section without an intermediate format.
    ! The mechanism filename extension does not control parser selection.
    integer function count_standard_chemkin_reactions(this, file_unit) result(nreactions)
        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: file_unit
        character(len=mechanism_line_length) :: line
        integer :: ios
        logical :: terminated
        nreactions = 0
        terminated = .false.
        call seek_section(file_unit, 'REACTIONS')
        do
            read(file_unit,'(A)',iostat=ios) line
            if (ios == iostat_end) exit
            if (ios /= 0) error stop 'CHEMKIN: reading REACTIONS failed'
            call strip_chemkin_comment(line)
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (upper_chemkin(trim(line)) == 'END') then
                terminated = .true.
                exit
            end if
            if (index(line,'=') > 0) nreactions = nreactions + 1
        end do
        if (.not. terminated) error stop 'CHEMKIN: missing REACTIONS END'
        rewind(file_unit)
    end function count_standard_chemkin_reactions

    subroutine read_standard_chemkin_reactions(this, file_unit, amount_units)
        class(chemical_properties), intent(inout) :: this
        integer, intent(in) :: file_unit
        character(len=*), intent(out) :: amount_units
        character(len=mechanism_line_length) :: line, equation
        character(len=128) :: parts(128), keyword
        character(len=20) :: energy_units
        integer :: ios, n, r, i
        logical :: is_falloff, third_body, has_low, has_troe, terminated
        logical :: seen_reaction
        amount_units = 'MOLES'
        ! CHEMKIN defaults: CAL/MOLE and MOLES, independent of the legacy
        ! E_act_units constructor argument.  Header declarations override.
        energy_units = 'cal.mol'
        call seek_section(file_unit, 'REACTIONS')
        ! seek_section has consumed the header: read its options separately.
        call read_standard_chemkin_units(file_unit,energy_units,amount_units)
        this%E_act_units = energy_units
        r = 0
        has_low = .false.
        has_troe = .false.
        is_falloff = .false.
        third_body = .false.
        terminated = .false.
        seen_reaction = .false.
        do
            read(file_unit,'(A)',iostat=ios) line
            if (ios == iostat_end) exit
            if (ios /= 0) error stop 'CHEMKIN: reaction file read failure'
            call strip_chemkin_comment(line)
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            call split_chemkin_words(line,parts,n)
            keyword = upper_chemkin(trim(parts(1)))
            if (keyword == 'END') then
                if (is_falloff .and. .not. has_low) &
                    error stop 'CHEMKIN: missing LOW for falloff reaction'
                terminated = .true.
                exit
            end if
            if (index(line,'=') > 0) then
                if (is_falloff .and. .not. has_low) &
                    error stop 'CHEMKIN: missing LOW for falloff reaction'
                r = r + 1
                if (r > this%reactions_number) error stop 'CHEMKIN: inconsistent reaction count'
                if (n < 4) error stop 'CHEMKIN: missing Arrhenius parameters'
                equation = ''
                do i=1,n-3
                    equation = trim(equation)//trim(parts(i))
                end do
                call parse_standard_chemkin_equation(this,trim(equation),r,is_falloff,third_body)
                call parse_chemkin_real(parts(n-2),this%A(r))
                call parse_chemkin_real(parts(n-1),this%beta(r))
                call parse_chemkin_real(parts(n),this%E_act(r))
                if (this%A(r) < 0.0_dp) error stop 'CHEMKIN: negative Arrhenius A unsupported'
                has_low = .false.
                has_troe = .false.
                seen_reaction = .true.
                cycle
            end if
            if (.not. seen_reaction) error stop 'CHEMKIN: modifier before first reaction'
            if (keyword == 'DUP' .or. keyword == 'DUPLICATE') cycle
            select case (trim(keyword))
            case ('PLOG','CHEB','TCHEB','PCHEB','SRI','REV','FORD','RORD', &
                  'HIGH','LT','RLT','UNITS','STICK','COV','MWON','MWOFF')
                write(error_unit,'(A,A)') &
                    'CHEMKIN: unsupported kinetics directive: ',trim(keyword)
                error stop 'CHEMKIN: unsupported kinetics model'
            end select
            if (chemkin_starts_with(line,'LOW')) then
                if (.not. is_falloff .or. has_low) &
                    error stop 'CHEMKIN: unexpected or repeated LOW'
                call parse_chemkin_parameters(line,'LOW',3,3,parts,n)
                call parse_chemkin_real(parts(1),this%A_low(r))
                call parse_chemkin_real(parts(2),this%beta_low(r))
                call parse_chemkin_real(parts(3),this%E_act_low(r))
                if (this%A_low(r) <= 0.0_dp) error stop 'CHEMKIN: invalid LOW A'
                has_low = .true.
            else if (chemkin_starts_with(line,'TROE')) then
                if (.not. is_falloff .or. .not. has_low .or. has_troe) &
                    error stop 'CHEMKIN: unexpected or repeated TROE'
                call parse_chemkin_parameters(line,'TROE',3,4,parts,n)
                do i=1,n
                    call parse_chemkin_real(parts(i),this%Troe_coeffs(r,i))
                end do
                has_troe = .true.
                if (this%reactions_type(r) == 2) this%reactions_type(r) = 3
                ! Irreversible falloff remains type 7; presence of TROE_coeffs
                ! tells chemical_kinetics_core_class to apply Troe broadening.
            else if (index(line,'/') > 0) then
                if (.not. third_body) error stop 'CHEMKIN: efficiencies without third body'
                call parse_standard_chemkin_efficiencies(this,line,r)
            else
                write(error_unit,'(A,A)') 'CHEMKIN: unsupported directive: ',trim(line)
                error stop 'CHEMKIN: unsupported reaction model'
            end if
        end do
        if (.not. terminated .or. r /= this%reactions_number) &
            error stop 'CHEMKIN: incorrect reaction count or missing END'
        rewind(file_unit)
    end subroutine read_standard_chemkin_reactions

    subroutine read_standard_chemkin_units(file_unit,energy_units,amount_units)
        integer, intent(in) :: file_unit
        character(len=*), intent(inout) :: energy_units, amount_units
        character(len=mechanism_line_length) :: line
        character(len=128) :: words(128), option
        integer :: n, i, ios
        rewind(file_unit)
        do
            read(file_unit,'(A)',iostat=ios) line
            if (ios /= 0) error stop 'CHEMKIN: missing REACTIONS header'
            call strip_chemkin_comment(line)
            line = adjustl(line)
            call split_chemkin_words(line,words,n)
            if (n == 0) cycle
            if (upper_chemkin(trim(words(1))) == 'REACTIONS') exit
        end do
        do i=2,n
            option=upper_chemkin(trim(words(i)))
            select case(trim(option))
            case ('MOLES','MOLE')
                amount_units='MOLES'
            case ('MOLECULES','MOLECULE')
                amount_units='MOLECULES'
            case ('CAL/MOLE','CAL/MOL','CALORIES/MOLE')
                energy_units='cal.mol'
            case ('KCAL/MOLE','KCAL/MOL','KCALORIES/MOLE')
                energy_units='kcal.mol'
            case ('JOULES/MOLE','JOULE/MOLE','JOULES/MOL')
                energy_units='J.mol'
            case ('KJOULES/MOLE','KJOULES/MOL','KJOULE/MOLE')
                energy_units='kJ.mol'
            case ('KELVINS','KELVIN')
                energy_units='kelvins'
            case default
                write(error_unit,'(A,A)') 'CHEMKIN: unknown REACTIONS option ', trim(option)
                error stop 'CHEMKIN: unsupported units'
            end select
        end do
    end subroutine read_standard_chemkin_units

    subroutine parse_standard_chemkin_equation(this,equation,reaction,falloff,thirdbody)
        class(chemical_properties), intent(inout) :: this
        character(len=*), intent(in) :: equation
        integer, intent(in) :: reaction
        logical, intent(out) :: falloff,thirdbody
        character(len=mechanism_line_length) :: left,right,compact
        integer :: arrow, arrow_len, i, nleft, nright, nmleft, nmright
        logical :: reversible
        compact = ''
        do i=1,len_trim(equation)
            if (equation(i:i) /= ' ' .and. equation(i:i) /= achar(9)) &
                compact = trim(compact)//equation(i:i)
        end do
        arrow = index(compact,'<=>')
        reversible = arrow > 0
        arrow_len = 3
        if (.not. reversible) then
            arrow = index(compact,'=>')
            arrow_len = 2
        end if
        if (arrow == 0) then
            arrow = index(compact,'=')
            arrow_len = 1
            reversible = .true.
        end if
        if (arrow == 0 .or. (index(compact,'<=') > 0 .and. &
            index(compact,'<=>') == 0)) &
            error stop 'CHEMKIN: unsupported reaction arrow'
        left = compact(:arrow-1)
        right = compact(arrow+arrow_len:)
        falloff = index(left,'(+M)') > 0 .or. index(right,'(+M)') > 0
        if (falloff) then
            if (index(left,'(+M)') == 0 .or. index(right,'(+M)') == 0) &
                error stop 'CHEMKIN: unmatched (+M) collider'
            call erase_chemkin_falloff(left)
            call erase_chemkin_falloff(right)
        end if
        call parse_standard_chemkin_side(this,trim(left), &
            this%chemical_coeffs(:,reaction,1),nleft,nmleft)
        call parse_standard_chemkin_side(this,trim(right), &
            this%chemical_coeffs(:,reaction,2),nright,nmright)
        if (nmleft /= nmright) error stop 'CHEMKIN: unmatched explicit third body'
        if (nmleft > 1 .or. nmright > 1) error stop 'CHEMKIN: repeated collider M'
        if (falloff .and. nmleft > 0) error stop 'CHEMKIN: mixed M and (+M)'
        thirdbody = falloff .or. nmleft > 0
        if (falloff) then
            if (reversible) then
                this%reactions_type(reaction)=2
            else
                this%reactions_type(reaction)=7
            end if
        else if (thirdbody) then
            if (reversible) then
                this%reactions_type(reaction)=1
            else
                this%reactions_type(reaction)=6
            end if
        else
            if (reversible) then
                this%reactions_type(reaction)=0
            else
                this%reactions_type(reaction)=5
            end if
        end if
    end subroutine parse_standard_chemkin_equation

    subroutine erase_chemkin_falloff(side)
        character(len=*), intent(inout) :: side
        integer :: pos
        pos=index(side,'(+M)')
        if (pos > 0) side=side(:pos-1)//side(pos+4:)
    end subroutine erase_chemkin_falloff

    subroutine parse_standard_chemkin_side(this,side,output,nparticipants,ncolliders)
        class(chemical_properties), intent(inout) :: this
        character(len=*), intent(in) :: side
        integer, intent(out) :: output(4),nparticipants,ncolliders
        character(len=mechanism_line_length) :: remaining,token
        integer :: next, coefficient, i, spec, ios, total
        output=0
        total=0
        ncolliders=0
        remaining=trim(side)
        do
            if (len_trim(remaining)==0) error stop 'CHEMKIN: empty reaction participant'
            next=index(remaining,'+')
            if (next > 0) then
                token=remaining(:next-1)
                remaining=remaining(next+1:)
            else
                token=remaining
                remaining=''
            end if
            token=adjustl(token)
            coefficient=1
            ! An exact species name takes precedence over a numeric prefix: some
            ! Chemkin species identifiers begin with a digit.
            spec=this%get_chemical_specie_index(trim(token))
            i=1
            if (spec==0) then
            do while (i<=len_trim(token))
                if (index('0123456789',token(i:i))==0) exit
                i=i+1
            end do
            if (i>1) then
                read(token(:i-1),*,iostat=ios) coefficient
                if (ios/=0 .or. coefficient<1) error stop 'CHEMKIN: invalid stoichiometry'
                token=adjustl(token(i:))
            end if
            end if
            if (len_trim(token)==0) error stop 'CHEMKIN: missing species name'
            if (index(token,'.')/=0) error stop 'CHEMKIN: fractional stoichiometry unsupported'
            spec=this%get_chemical_specie_index(trim(token))
            if (spec==0) then
                write(error_unit,'(A,A)') 'CHEMKIN: unknown species: ',trim(token)
                error stop 'CHEMKIN: unknown reaction species'
            end if
            if (spec==this%species_number+1) ncolliders=ncolliders+coefficient
            if (total+coefficient>3) error stop 'CHEMKIN: >3 reactant/product occurrences unsupported'
            output(total+2:total+coefficient+1)=spec
            total=total+coefficient
            if (next==0) exit
        end do
        nparticipants=total
        output(1)=total
    end subroutine parse_standard_chemkin_side

    subroutine parse_standard_chemkin_efficiencies(this,record,reaction)
        class(chemical_properties), intent(inout) :: this
        character(len=*), intent(in) :: record
        integer, intent(in) :: reaction
        character(len=mechanism_line_length) :: rest,name,value
        integer :: slash1, slash2, spec
        rest=adjustl(record)
        if (rest(1:1)=="'" .or. rest(1:1)=='"') &
            error stop 'CHEMKIN: quoted efficiency lists are nonstandard'
        do
            rest=adjustl(rest)
            if (len_trim(rest)==0) exit
            if (len_trim(rest)==0) exit
            slash1=index(rest,'/')
            if (slash1<2) error stop 'CHEMKIN: malformed third-body efficiency'
            name=trim(adjustl(rest(:slash1-1)))
            spec=this%get_chemical_specie_index(trim(name))
            if (spec < 1 .or. spec > this%species_number) then
                write(error_unit,'(A,A)') 'CHEMKIN: unknown efficiency species: ',trim(name)
                error stop 'CHEMKIN: unknown third-body species'
            end if
            rest=adjustl(rest(slash1+1:))
            slash2=index(rest,'/')
            if (slash2<2) error stop 'CHEMKIN: malformed efficiency value'
            value=trim(adjustl(rest(:slash2-1)))
            call parse_chemkin_real(trim(value),this%enhanced_efficiencies(reaction,spec))
            if (this%enhanced_efficiencies(reaction,spec)<0.0_dp) &
                error stop 'CHEMKIN: negative third-body efficiency'
            rest=rest(slash2+1:)
        end do
    end subroutine parse_standard_chemkin_efficiencies

    subroutine parse_chemkin_parameters(line,keyword,min_n,max_n,parts,n)
        character(len=*), intent(in) :: line,keyword
        integer, intent(in) :: min_n,max_n
        character(len=*), intent(out) :: parts(:)
        integer, intent(out) :: n
        character(len=mechanism_line_length) :: rest
        integer :: first,last
        rest=adjustl(line(len_trim(keyword)+1:))
        if (len_trim(rest)==0) error stop 'CHEMKIN: missing reaction parameters'
        first=index(rest,'/')
        if (first==0) error stop 'CHEMKIN: standard LOW/TROE require slash delimiters'
        if (len_trim(rest(:first-1))>0) error stop 'CHEMKIN: malformed modifier'
        rest=adjustl(rest(first+1:))
        last=index(rest,'/')
        if (last==0) error stop 'CHEMKIN: missing closing delimiter'
        if (len_trim(rest(last+1:))>0) error stop 'CHEMKIN: trailing data after modifier'
        rest=rest(:last-1)
        call split_chemkin_words(rest,parts,n)
        if (n<min_n .or. n>max_n) error stop 'CHEMKIN: invalid parameter count'
    end subroutine parse_chemkin_parameters

    subroutine parse_chemkin_real(token,value)
        character(len=*), intent(in) :: token
        real(dp), intent(out) :: value
        integer :: ios
        read(token,*,iostat=ios) value
        if (ios/=0) then
            write(error_unit,'(A,A)') 'CHEMKIN: invalid number: ',trim(token)
            error stop 'CHEMKIN: non-numeric reaction parameter'
        end if
    end subroutine parse_chemkin_real

    subroutine strip_chemkin_comment(line)
        character(len=*), intent(inout) :: line
        integer :: i
        i=index(line,'!')
        if (i>0) line(i:)=' '
    end subroutine strip_chemkin_comment

    subroutine split_chemkin_words(line,words,n)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: words(:)
        integer, intent(out) :: n
        integer :: i,j,L
        L=len_trim(line)
        words=''
        n=0
        i=1
        do while (i<=L)
            if (line(i:i)==' ' .or. line(i:i)==achar(9) .or. line(i:i)==achar(13)) then
                i=i+1
                cycle
            end if
            j=i
            do while (j<=L)
                if (line(j:j)==' ' .or. line(j:j)==achar(9) .or. line(j:j)==achar(13)) exit
                j=j+1
            end do
            n=n+1
            if (n>size(words)) error stop 'CHEMKIN: too many words on one record'
            if (j-i>len(words(1))) error stop 'CHEMKIN: oversized input token'
            words(n)=line(i:j-1)
            i=j
        end do
    end subroutine split_chemkin_words

    pure function upper_chemkin(src) result(dst)
        character(len=*), intent(in) :: src
        character(len=len(src)) :: dst
        integer :: i,j
        dst=src
        do i=1,len(dst)
            j=iachar(dst(i:i))
            if (j>=iachar('a') .and. j<=iachar('z')) dst(i:i)=achar(j-32)
        end do
    end function upper_chemkin

    logical function chemkin_starts_with(line,key) result(found)
        character(len=*), intent(in) :: line,key
        character(len=mechanism_line_length) :: tmp
        integer :: n
        tmp=upper_chemkin(adjustl(line))
        n=len_trim(key)
        found=.false.
        if (len_trim(tmp)<n) return
        if (tmp(:n)/=key) return
        if (len_trim(tmp)==n) then
            found=.true.
        else
            found=tmp(n+1:n+1)==' ' .or. tmp(n+1:n+1)=='/'
        end if
    end function chemkin_starts_with

end module chemical_properties_class
