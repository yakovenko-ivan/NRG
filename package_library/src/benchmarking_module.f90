module benchmarking
    use :: kind_parameters
    use,intrinsic :: iso_fortran_env, only : stdout=>OUTPUT_UNIT

implicit none
private

type timer
       real(dp)             :: cpu_start
       real(dp)             :: cpu_end
       integer(i8)          :: clock_start
       integer(i8)          :: clock_end
       integer              :: wall_start(8)
       integer              :: wall_end(8)
       real(dp)             :: cum_elapsed
       real(dp)             :: cum_cpu
       integer              :: iter
       character(len=100)   :: description
       character(len=100)   :: short_description
   contains
       procedure  ::  tic               =>  clock_tic
       procedure  ::  toc               =>  clock_toc
       procedure  ::  print             =>  clock_print
       procedure  ::  walltime          =>  clock_walltime
       procedure  ::  cputime           =>  clock_cputime
       procedure  ::  get_header        =>  generate_header 
end type
   
type :: timer_pointer
    type(timer) ,   pointer    :: timer_ptr
end type timer_pointer

interface timer
     procedure :: clock_new
end interface timer

public :: timer, timer_pointer

character(len=*),parameter :: gen='(*(g0))'
character(len=*),parameter  :: all='(*(g0,1x))'

contains

! initialization constructor
    type(timer) function clock_new(description, short_description)
        character(len=*)  ,intent(in)     :: description
        character(len=*)  ,intent(in)     :: short_description
        
        call cpu_time(clock_new%cpu_start)
        call system_clock(clock_new%clock_start)
        call date_and_time(values=clock_new%wall_start)

        clock_new%cpu_end   = clock_new%cpu_start
        clock_new%clock_end = clock_new%clock_start
        clock_new%wall_end  = clock_new%wall_start

        clock_new%description       = description
        clock_new%short_description = short_description
    
        clock_new%cum_elapsed   = 0.0
        clock_new%cum_cpu       = 0.0
    
        clock_new%iter          = 0
    
    end function clock_new


    subroutine clock_tic(this)
        class(timer) :: this

       call cpu_time(this%cpu_start)
       call system_clock(this%clock_start)
       call date_and_time(values=this%wall_start)

       this%cpu_end   = this%cpu_start
       this%clock_end = this%clock_start
       this%wall_end  = this%wall_start

    end subroutine clock_tic

    subroutine clock_toc(this,new_iter)
        class(timer)                            :: this
        logical     ,intent(in),    optional    :: new_iter
    
        call cpu_time(this%cpu_end)
        call system_clock(this%clock_end)
        call date_and_time(values=this%wall_end)

        if(present(new_iter)) then
            if(new_iter)    this%iter = this%iter + 1
        end if
    
        this%cum_elapsed   = this%cum_elapsed + this%walltime()
        this%cum_cpu       = this%cum_cpu + this%cputime()
    
    end subroutine clock_toc

    subroutine clock_print(this,lun,num_threads)
        class(timer)    ,intent(in)                 :: this
        integer(i4)     ,intent(in),optional        :: lun
        integer(i2)     ,intent(in),optional        :: num_threads
        integer(i4)                                 :: lun_
        real(dp)                                    :: elapsed_time
        real(dp)                                    :: cpu_time
        real(dp)                                    :: cum_elapsed
        real(dp)                                    :: cum_cpu
        character(len=105)                          :: biggest
        integer(i8)                                 :: count_rate

        if(present(lun))then
            lun_=lun
        else
            lun_=stdout
        endif

        cum_elapsed = this%cum_elapsed
        cum_cpu     = this%cum_cpu
    
        elapsed_time           =  this%cum_elapsed / this%iter !this%walltime()
        cpu_time               =  this%cum_cpu / this%iter     !this%cputime()
    
        if(present(lun)) then
            if(this%cum_elapsed > 0.1) then
                write( lun_,'(2f10.5)', advance='no') cum_elapsed, cum_cpu!, elapsed_time, cpu_time
            end if
        else
            if(this%cum_elapsed > 0.1) then
                write( lun_,gen ) "###############"
                write( lun_,gen ) this%description
                write( lun_,gen ) "###############"
    
                ! try to make a reasonable format for the number of digits of precision
                call system_clock(count_rate=count_rate) ! Find the time rate
                write(biggest,'("(a,f0.",i0,")")')ceiling(log10(real(count_rate,kind=dp)))

                write( lun_,'(A,f10.5)')    'Total elapsed time (sec) ::', cum_elapsed
                write( lun_,'(A,f10.5)')    'Total CPU time     (sec) ::', cum_cpu    
                write( lun_,'(A,f10.5)')    'Average elapsed time (sec/iter)  ::', elapsed_time
                write( lun_,'(A,f10.5)')    'Average CPU time     (sec/iter)  ::', cpu_time
                write( lun_,'(A,I6)')       'Amount of iterations (iter)      ::', this%iter
                write( lun_,'(A,1x,f0.2)')  'Percentage         ::',(this%cum_cpu/this%cum_elapsed)*100
            end if
        end if
        
    end subroutine clock_print

    function clock_walltime(this) result(elapsed_time)
        class(timer)        :: this
        integer(kind=i8)    :: count_rate
        real(kind=dp)       :: elapsed_time
        real(kind=dp)       :: cpu_time
        
           call system_clock(count_rate=count_rate)
           elapsed_time = real(this%clock_end-this%clock_start,kind=dp)/real(count_rate,kind=dp)
           
    end function clock_walltime

    function  clock_cputime(this)  result(cpu_time)
        class(timer)         :: this
        real(kind=dp)        :: cpu_time
        
           cpu_time = real(this%cpu_end-this%cpu_start,kind=dp)
           
    end function clock_cputime
    
    function generate_header(this)  result(header)
        class(timer)        ,intent(in) :: this
        character(len=100)               :: header
        
        if(this%cum_elapsed > 0.1) then
            write(header,'(2A)') trim(this%short_description) // '_TELA    ', trim(this%short_description) // '_TCPU    ' !, trim(this%short_description) // '_ELA/i    ', trim(this%short_description) // '_CPU/i    '
        else
            header=''
        end if
    
    end function generate_header
end module benchmarking