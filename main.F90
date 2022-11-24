program main
    use, intrinsic :: iso_fortran_env, only: DP => REAL64, I64 => INT64
    use            :: therm_m,         only: calc_KE, calc_vdw_pbc, pbc, calc_vdw_force
    use            :: init_cond_m,     only: init_positions_sc
    use            :: math_utils_m,    only: r8_normal_ab
    use            :: thermostats_m,   only: andersen_thermostat
    use            :: integrtors_m,    only: velocity_verlet
    use            :: writers_m,       only: write_frame
    implicit none

    ! Scalar values 
    real(kind=DP)       :: box, ke_calc, pe_calc, density, time0, time1, instant_temperature, mass, timestep, cutoff_set
    integer(kind=I64)   :: n_particles, n_steps, stp, traj_unit, write_file, log_unit

    ! Arrays
    real(kind=DP), allocatable, dimension(:, :) :: positions, velocities


    ! Files
    open(newunit=traj_unit, file="heating.xyz", access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=log_unit, file="heating.log", access='sequential', action='write',&
    status='replace', form='formatted')

    ! ~ Definim parametres ~
    mass = 1.0_dp
    timestep = 0.0001_dp
    cutoff_set = 2.5_dp
    write_file = 8_I64
    n_particles = 125_I64
    density = 0.7_DP
    n_steps = 10000_I64
    box = (real(n_particles, kind=DP)/density) ** (1.0_DP / 3.0_DP)

    ! ~ Inicialitzem les posicions i velocitats~
    allocate(positions(n_particles, 3), velocities(n_particles, 3))
    call init_positions_sc(rho=density, pos=positions)

    call init_velocities(vel=velocities)

    print '(A,F16.8,A,F16.8)', "Initial potential energy=", calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box), " Initial kinetic energy=", calc_KE(velocities)

    call write_frame(unit_nr=traj_unit, pos=positions, stp_c=0_I64)

    call cpu_time(time0)

    ! ~ Integrem les posicions dels atoms ~
    do stp = 1, n_steps

        ! Aplico integrador de velocity verlet
        call velocity_verlet(pos=positions, vel=velocities, dt=timestep, boundary=box, mass=mass, cutoff=cutoff_set)

        ! Aplico el termostat
        call andersen_thermostat(vel=velocities, temp=100.0_DP, nu=0.2_DP)

        ! Donem informacio del sistema
        call write_frame(unit_nr=traj_unit, pos=positions, stp_c=stp)
        pe_calc = calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box)
        ke_calc = calc_KE(velocities)
        instant_temperature = (2.0_DP / (3.0_DP * n_particles)) * ke_calc
        write(unit=log_unit, fmt='(I6,F25.8,F25.8,F25.8,F25.8)') stp, ke_calc, pe_calc, ke_calc + pe_calc, instant_temperature
    end do

    print '(A,F16.8,A,F16.8)', "Post-heating potential energy=", calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box), " kinetic energy=", calc_KE(velocities)

    call cpu_time(time1)

    print '(A,F12.8,A,F12.8)', "Execution time for heating:", time1 - time0, " time/iteration: ", (time1 - time0) / n_steps

    call write_frame(unit_nr=traj_unit, pos=positions, stp_c = stp)

    close(traj_unit)
    close(log_unit)

    ! Production trajectory
    open(newunit=traj_unit, file="thermodynamics.dat", access='sequential', action='write',&
    status='replace', form='formatted')

    n_steps = 500000

    call cpu_time(time0)

    print '(A,F16.8,A,F16.8)', "Initial traj potential energy=", calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box), " kinetic energy=", calc_KE(velocities)

    ! ~ Integrem les posicions dels atoms ~
    do stp = 1, n_steps

        ! Aplico integrador de velocity verlet
        call velocity_verlet(pos=positions, vel=velocities, dt=timestep, boundary=box, mass=mass, cutoff=cutoff_set)

        ! Aplico el termostat
        call andersen_thermostat(vel=velocities, temp=1.5_DP, nu=0.2_DP)

        ! Donem informacio del sistema
        if (mod(stp, write_file) == 0) then
            pe_calc = calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box)
            ke_calc = calc_KE(velocities)
            instant_temperature = (2.0_DP / (3.0_DP * n_particles)) * ke_calc
            write(unit=traj_unit, fmt='(I6,F25.8,F25.8,F25.8,F25.8)') stp, ke_calc, pe_calc, ke_calc + pe_calc, instant_temperature
        end if
    end do

    print '(A,F16.8,A,F16.8)', "Final traj potential energy=", calc_vdw_pbc(pos=positions, cutoff=cutoff_set, boundary=box), " kinetic energy=", calc_KE(velocities)

    call cpu_time(time1)

    print '(A,F12.8,A,F12.8)', "Execution time for production:", time1 - time0, " time/iteration: ", (time1 - time0) / n_steps

    close(traj_unit)

    deallocate(positions)
    deallocate(velocities)

contains

    pure subroutine init_velocities(vel)
        implicit none
        real(kind=DP), dimension(:, :), intent(inout) :: vel

        vel = 0.0_DP
    end subroutine init_velocities

end program main