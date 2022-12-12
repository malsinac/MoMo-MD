program main
    use, intrinsic :: iso_fortran_env, only: DP => REAL64, I64 => INT64
    use            :: therm_m,         only: calc_KE, calc_vdw_pbc
    use            :: init_cond_m,     only: init_positions_sc, bimodal_dist_velocities
    use            :: integrtors_m,    only: velocity_verlet, euler_integrator
    use            :: writers_m,       only: write_velocities
    use            :: interface_m,     only: databloc_params_t
    implicit none

    ! Derived types
    type(databloc_params_t) :: datablock
    
    ! Scalar values 
    real(kind=DP)       :: time0, time1
    integer(kind=I64)   :: log_unit, vel_unit, rdf_unit, msd_unit, frame_unit
    character(len=2048) :: vel_name, log_name, rdf_name, msd_name, frame_name

    ! Arrays
    real(kind=DP), allocatable, dimension(:, :) :: positions, velocities, init_position

    
    ! ~ Definim parametres ~
    
    ! Particle related variables
    datablock%lj_epsilon = 0.998_dp  ! [kJ/mol]
    datablock%lj_sigma = 3.4_dp      ! [A]
    datablock%mass = 40.0_dp         ! [g/mol]
    ! Simulation related variables
    datablock%timestep = 0.001_dp
    datablock%n_particles = 343_I64
    datablock%density = 0.1_DP
    datablock%n_steps = 100000_I64
    ! I/O variables
    datablock%write_file = 100_I64
    datablock%write_stats = 1000_i64
    datablock%write_frame = 0_i64
    datablock%sim_name = "dens_0-01"
    ! Thermostat variables
    datablock%ref_temp = 100_dp
    datablock%andersen_nu = 0.2_dp
    ! Simulation-dependent variables
    datablock%box = (real(datablock%n_particles, kind=DP)/datablock%density) ** (1.0_DP / 3.0_DP)
    datablock%cutoff_set = 0.5_dp * datablock%box
    ! Analysis dependent variables
    datablock%gdr_num_bins = 100_i64
    datablock%gdr_max_dist = datablock%box

    ! ~ Realitzem l'equilibrat del sistema ~

    log_name = trim(datablock%sim_name) // "_heating.log"
    vel_name = trim(datablock%sim_name) // "_heating.vel"

    ! Files
    open(newunit=log_unit, file=trim(log_name), access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=vel_unit, file=trim(vel_name), access='sequential', action='write',&
    status='replace', form='formatted')

    ! ~ Inicialitzem les posicions i velocitats~
    allocate(positions(datablock%n_particles, 3), velocities(datablock%n_particles, 3))
    allocate(init_position(datablock%n_particles, 3))
    
    call init_positions_sc(rho=datablock%density, pos=positions)
    call bimodal_dist_velocities(vel=velocities, temp=datablock%ref_temp)

    print '(A,F16.8,A,F16.8)', "Initial potential energy=", calc_vdw_pbc(pos=positions, cutoff=datablock%cutoff_set, boundary=datablock%box), &
                               " Initial kinetic energy=", calc_KE(velocities)

    call write_velocities(vel=velocities, unit_nr=vel_unit, step=0_i64)

    call cpu_time(time0)

    call velocity_verlet(vel=velocities, pos=positions, parambox=datablock, &
                         log_unit=log_unit, thermost=.TRUE.)

    call cpu_time(time1)

    call write_velocities(velocities, vel_unit, datablock%n_steps)

    print '(A,F12.8,A,F12.8)', "Execution time for heating:", time1 - time0, " time/iteration: ", (time1 - time0) / datablock%n_steps
    close(log_unit)
    close(vel_unit)

    ! ----------------------------------------------------------------------------------------------------------------------------------------
    
    ! ~ Realitzem la producció del sistema ~
    
    datablock%ref_temp = 1.2_dp
    datablock%n_steps = 500000_i64
    datablock%write_frame = 10000_i64
    init_position = positions

    ! Files
    log_name = trim(datablock%sim_name) // "_prod.log"
    vel_name = trim(datablock%sim_name) // "_prod.vel"
    frame_name = trim(datablock%sim_name) // "_prod.trr"
    rdf_name = trim(datablock%sim_name) // "_rdf.dat"
    msd_name = trim(datablock%sim_name) // "_msd.dat"

    open(newunit=log_unit, file=trim(log_name), access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=vel_unit, file=trim(vel_name), access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=rdf_unit, file=trim(rdf_name), access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=msd_unit, file=trim(msd_name), access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=frame_unit, file=trim(frame_name), access='sequential', action='write',&
    status='replace', form='formatted')
    
    call cpu_time(time0)
    
    call velocity_verlet(vel=velocities, pos=positions, parambox=datablock, &
                         log_unit=log_unit, thermost=.TRUE., rdf_unit=rdf_unit, &
                         msd_unit=msd_unit, init_pos=init_position, frame_unit=frame_unit)
    
    call cpu_time(time1)

    print '(A,F12.8,A,F12.8)', "Execution time for production: ", time1 - time0, " time/iteration: ", (time1 - time0) / datablock%n_steps

    close(log_unit)
    close(vel_unit)
    close(rdf_unit)
    close(msd_unit)
    close(frame_unit)

    deallocate(positions)
    deallocate(velocities)

end program main