module interface_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    type :: databloc_params_t
        ! Particle-related variables
        real(kind=dp)       :: mass
        real(kind=dp)       :: lj_epsilon
        real(kind=dp)       :: lj_sigma 
        ! Simulation-related variables
        real(kind=dp)       :: timestep      ! Reduced time
        real(kind=dp)       :: cutoff_set
        integer(kind=i64)   :: n_particles   
        real(kind=dp)       :: density       ! Reduced density
        integer(kind=i64)   :: n_steps       ! Equilibration steps
        integer(kind=i64)   :: n_steps_prod  ! Production steps  
        ! I/O variables
        integer(kind=i64)   :: write_file    ! How often do we write log information
        integer(kind=i64)   :: write_stats   ! How often do we write statistics
        character(len=2048) :: sim_name      ! Name of the simulation we are doing
        ! Thermostat variables
        real(kind=dp)       :: ref_temp      ! Reference temperature for thermostats and velocity initialization 
        real(kind=dp)       :: andersen_nu   ! Andersen parameter
        ! Simulation-dependent variables
        real(kind=dp)       :: box           ! Size of the box
    end type



    !abstract interface
    !    subroutine integrator_func(vel, pos, parambox, log_unit, therm_ptr)
    !        import :: databloc_params_t
    !        real(kind=8), intent(inout), dimension(:,:)     :: pos, vel
    !        integer(kind=8), intent(in)                     :: log_unit
    !        procedure(thermostat_func), pointer             :: therm_ptr
    !        type(databloc_params_t), intent(in)             :: parambox
    !    end subroutine
    !end interface

    !abstract interface 
    !    subroutine thermostat_func(vel, parambox)
    !        import :: databloc_params_t
    !        real(kind=8), dimension(:, :), intent(inout)    :: vel
    !        type(databloc_params_t), intent(in)             :: parambox
    !    end subroutine
    !end interface
end module interface_m