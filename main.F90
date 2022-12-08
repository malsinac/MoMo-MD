program main
    use, intrinsic :: iso_fortran_env, only: DP => REAL64, I64 => INT64
    use            :: therm_m,         only: calc_KE, calc_vdw_pbc, pbc, calc_vdw_force, compute_com_momenta
    use            :: init_cond_m,     only: init_positions_sc, bimodal_dist_velocities
    use            :: math_utils_m,    only: r8_normal_ab
    use            :: thermostats_m,   only: andersen_thermostat
    use            :: integrtors_m,    only: velocity_verlet, euler_integrator
    use            :: writers_m,       only: write_frame, write_system_information, write_velocities
    use            :: interface_m,     only: databloc_params_t
    implicit none

    ! Derived types
    type(databloc_params_t) :: datablock
    
    ! Scalar values 
    real(kind=DP)       :: time0, time1
    integer(kind=I64)   :: log_unit, vel_unit, rdf_unit
    character(len=2048) :: vel_name, log_name

    ! Arrays
    real(kind=DP), allocatable, dimension(:, :) :: positions, velocities

    
    ! ~ Definim parametres ~
    
    ! Particle related variables
    datablock%lj_epsilon = 0.998_dp
    datablock%lj_sigma = 3.4_dp
    datablock%mass = 1.0_dp
    ! Simulation related variables
    datablock%timestep = 0.1_dp
    datablock%cutoff_set = 2.5_dp
    datablock%n_particles = 216_I64
    datablock%density = 0.7_DP
    datablock%n_steps = 100000_I64
    datablock%n_steps_prod = 500000_i64
    ! I/O variables
    datablock%write_file = 4_I64
    datablock%write_stats = 0_i64
    datablock%sim_name = "vverlet_0-1"
    ! Thermostat variables
    datablock%ref_temp = 100_dp
    datablock%andersen_nu = 0.2_dp
    ! Simulation-dependent variables
    datablock%box = (real(datablock%n_particles, kind=DP)/datablock%density) ** (1.0_DP / 3.0_DP)


    ! ~ Realitzem l'equilibrat del sistema ~

    log_name = trim(datablock%sim_name) // ".log"
    vel_name = trim(datablock%sim_name) // ".vel"

    ! Files
    open(newunit=log_unit, file=trim(log_name), access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=vel_unit, file=trim(vel_name), access='sequential', action='write',&
    status='replace', form='formatted')

    ! ~ Inicialitzem les posicions i velocitats~
    allocate(positions(datablock%n_particles, 3), velocities(datablock%n_particles, 3))
    
    call init_positions_sc(rho=datablock%density, pos=positions)
    call bimodal_dist_velocities(vel=velocities, temp=datablock%ref_temp, mass=datablock%mass)

    print '(A,F16.8,A,F16.8)', "Initial potential energy=", calc_vdw_pbc(pos=positions, cutoff=datablock%cutoff_set, boundary=datablock%box), " Initial kinetic energy=", calc_KE(velocities)

    call write_velocities(vel=velocities, unit_nr=vel_unit, step=0_i64)

    call cpu_time(time0)

    call velocity_verlet(vel=velocities, pos=positions, parambox=datablock, &
                         log_unit=log_unit)

    call cpu_time(time1)

    call write_velocities(velocities, vel_unit, datablock%n_steps)

    print '(A,F12.8,A,F12.8)', "Execution time for heating:", time1 - time0, " time/iteration: ", (time1 - time0) / datablock%n_steps
    close(log_unit)
    close(vel_unit)

    ! ~ Realitzem la producció del sistema ~
    !datablock%write_stats = 10_i64

    !log_name = trim(datablock%sim_name) // "prod.log"
    !vel_name = trim(datablock%sim_name) // "prod.vel"
    !rdf_name = trim(datablock%sim_name) // "rdf.dat"

    ! Files
    !open(newunit=log_unit, file=trim(log_name), access='sequential', action='write',&
    !status='replace', form='formatted')
    !open(newunit=vel_unit, file=trim(vel_name), access='sequential', action='write',&
    !status='replace', form='formatted')
    !open(newunit=rdf_unit, file=trim(rdf_name), access='sequential', action='write',&
    !status='replace', form='formatted')
    
    !call cpu_time(time0)
    
    !call velocity_verlet(vel=velocities, pos=positions, parambox=datablock, &
    !                     log_unit=log_unit)
    !call cpu_time(time1)

    !print '(A,F12.8,A,F12.8)', "Execution time for production:", time1 - time0, " time/iteration: ", (time1 - time0) / datablock%n_steps_prod


    !close(log_unit)
    !close(vel_unit)

    deallocate(positions)
    deallocate(velocities)

end program main