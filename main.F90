program main
    use, intrinsic :: iso_fortran_env, only: DP => REAL64, I64 => INT64
    use            :: therm_m,         only: calc_KE, calc_vdw_pbc, pbc, calc_vdw_force, compute_com_momenta
    use            :: init_cond_m,     only: init_positions_sc, bimodal_dist_velocities
    use            :: math_utils_m,    only: r8_normal_ab
    use            :: thermostats_m,   only: andersen_thermostat
    use            :: integrtors_m,    only: velocity_verlet
    use            :: writers_m,       only: write_frame, write_system_information, write_velocities
    implicit none

    ! Scalar values 
    real(kind=DP)       :: box, density, time0, time1, mass, timestep, cutoff_set
    integer(kind=I64)   :: n_particles, n_steps, stp, write_file, log_unit, vel_unit

    ! Arrays
    real(kind=DP), allocatable, dimension(:, :) :: positions, velocities

    ! Files
    open(newunit=log_unit, file="part1.log", access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=vel_unit, file="part1_vel.log", access='sequential', action='write',&
    status='replace', form='formatted')

    ! ~ Definim parametres ~
    mass = 1.0_dp
    timestep = 0.0001_dp
    cutoff_set = 2.5_dp
    write_file = 4_I64
    n_particles = 125_I64
    density = 0.7_DP
    n_steps = 100000_I64
    box = (real(n_particles, kind=DP)/density) ** (1.0_DP / 3.0_DP)

    ! ~ Inicialitzem les posicions i velocitats~
    allocate(positions(n_particles, 3), velocities(n_particles, 3))
    call init_positions_sc(rho=density, pos=positions)
    call bimodal_dist_velocities(vel=velocities, temp=100.0_dp, mass=mass)

    print '(A,F16.8,A,F16.8)', "Initial potential energy=", calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box), " Initial kinetic energy=", calc_KE(velocities)

    call write_velocities(vel=velocities, unit_nr=vel_unit, step=0_i64)

    call cpu_time(time0)

    ! ~ Integrem les posicions dels atoms ~
    do stp = 1, n_steps

        ! Aplico integrador de velocity verlet
        call velocity_verlet(pos=positions, vel=velocities, dt=timestep, boundary=box, mass=mass, cutoff=cutoff_set)

        ! Donem informacio del sistema
        if (mod(stp, write_file) == 0_i64) then
            call write_system_information(pos=positions, vel=velocities,&
            cutoff=cutoff_set, frame=stp, unit=log_unit, boundary=box, mass=mass, dens=density)
        end if
    end do

    call cpu_time(time1)

    call write_velocities(velocities, vel_unit, stp)

    print '(A,F12.8,A,F12.8)', "Execution time for heating:", time1 - time0, " time/iteration: ", (time1 - time0) / n_steps

    close(log_unit)
    close(vel_unit)

    deallocate(positions)
    deallocate(velocities)

end program main