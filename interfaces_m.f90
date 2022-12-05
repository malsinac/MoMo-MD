module interface_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    type :: databloc_params_t
        ! Particle-related variables
        real(kind=dp)     :: mass
        ! Simulation-related variables
        real(kind=dp)     :: timestep
        real(kind=dp)     :: cutoff_set
        integer(kind=i64) :: n_particles
        real(kind=dp)     :: density
        integer(kind=i64) :: n_steps
        ! I/O variables
        integer(kind=i64) :: write_file
        ! Thermostat variables
        real(kind=dp)     :: ref_temp
        real(kind=dp)     :: andersen_nu
        ! Simulation-dependent variables
        real(kind=dp)     :: box
    end type

    abstract interface
        subroutine integrator_func(vel, pos, parambox, log_unit, therm_ptr)
            import :: databloc_params_t
            real(kind=8), intent(inout), dimension(:,:)     :: pos, vel
            integer(kind=8), intent(in)                     :: log_unit
            procedure(thermostat_func), pointer             :: therm_ptr
            type(databloc_params_t), intent(in)             :: parambox
        end subroutine
    end interface

    abstract interface 
        subroutine thermostat_func(vel, parambox)
            import :: databloc_params_t
            real(kind=8), dimension(:, :), intent(inout)    :: vel
            type(databloc_params_t), intent(in)             :: parambox
        end subroutine
    end interface
end module interface_m