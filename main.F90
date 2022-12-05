program main
    use, intrinsic :: iso_fortran_env, only: DP => REAL64, I64 => INT64
    use            :: therm_m,         only: calc_KE, calc_vdw_pbc, pbc, calc_vdw_force, compute_com_momenta
    use            :: init_cond_m,     only: init_positions_sc, bimodal_dist_velocities
    use            :: math_utils_m,    only: r8_normal_ab
    use            :: thermostats_m,   only: andersen_thermostat, null_thermostat
    use            :: integrtors_m,    only: velocity_verlet
    use            :: writers_m,       only: write_frame, write_system_information, write_velocities
    use            :: interface_m,     only: thermostat_func, integrator_func, databloc_params_t
    implicit none

    ! Derived types
    type(databloc_params_t) :: datablock
    ! Scalar values 
    real(kind=DP)       :: time0, time1
    integer(kind=I64)   :: stp, log_unit, vel_unit

    ! Arrays
    real(kind=DP), allocatable, dimension(:, :) :: positions, velocities

    ! Pointers
    procedure(thermostat_func), pointer :: thermostat_ptr => null()
    procedure(integrator_func), pointer :: integrator_ptr => null()

    ! Files
    open(newunit=log_unit, file="part1.log", access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=vel_unit, file="part1_vel.log", access='sequential', action='write',&
    status='replace', form='formatted')

    ! ~ Definim parametres ~
    datablock%mass = 1.0_dp
    datablock%timestep = 0.0001_dp
    datablock%cutoff_set = 2.5_dp
    datablock%write_file = 4_I64
    datablock%n_particles = 125_I64
    datablock%density = 0.7_DP
    datablock%n_steps = 100000_I64
    datablock%box = (real(datablock%n_particles, kind=DP)/datablock%density) ** (1.0_DP / 3.0_DP)
    datablock%ref_temp = 100_dp
    thermostat_ptr => null_thermostat
    integrator_ptr => velocity_verlet

    ! ~ Inicialitzem les posicions i velocitats~
    allocate(positions(datablock%n_particles, 3), velocities(datablock%n_particles, 3))
    
    call init_positions_sc(rho=datablock%density, pos=positions)
    call bimodal_dist_velocities(vel=velocities, temp=datablock%ref_temp, mass=datablock%mass)

    print '(A,F16.8,A,F16.8)', "Initial potential energy=", calc_vdw_pbc(pos=positions, cutoff=datablock%cutoff_set, boundary=datablock%box), " Initial kinetic energy=", calc_KE(velocities)

    call write_velocities(vel=velocities, unit_nr=vel_unit, step=0_i64)

    call cpu_time(time0)

    ! ~ Integrem les posicions dels atoms ~
    call velocity_verlet(vel=velocities, pos=positions, parambox=datablock, &
                         log_unit=log_unit, therm_ptr=thermostat_ptr)

    call cpu_time(time1)

    call write_velocities(velocities, vel_unit, stp)

    print '(A,F12.8,A,F12.8)', "Execution time for heating:", time1 - time0, " time/iteration: ", (time1 - time0) / datablock%n_steps

    close(log_unit)
    close(vel_unit)

    deallocate(positions)
    deallocate(velocities)

end program main