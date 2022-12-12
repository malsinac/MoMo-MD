module integrtors_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: therm_m,         only: pbc, calc_vdw_force
    use            :: writers_m,       only: write_system_information, write_rdf, write_msd
    use            :: interface_m,     only: databloc_params_t
    use            :: thermostats_m,   only: andersen_thermostat

contains

    subroutine velocity_verlet(pos, vel, parambox, log_unit, rdf_unit, msd_unit, &
        init_pos, thermost)
        implicit none
        ! In/Out variables
        real(kind=dp), intent(inout), dimension(:,:)        :: pos, vel
        type(databloc_params_t), intent(in)                 :: parambox
        integer(kind=i64), intent(in)                       :: log_unit
        integer(kind=i64), intent(in), optional             :: rdf_unit, msd_unit
        real(kind=dp), intent(in), dimension(:,:), optional :: init_pos
        logical, intent(in)                                 :: thermost
        ! Internal_variables
        integer(kind=I64)                                   :: idx_, stp
        real(kind=DP), dimension(:, :), allocatable         :: actual_force, new_force
        logical                                             :: present_rdf, present_msd

        allocate(actual_force(parambox%n_particles, 3))
        allocate(new_force(parambox%n_particles, 3))

        present_rdf = present(rdf_unit)
        present_msd = present(msd_unit)

        ! Posem tot abans a la caixa unitaria
        do idx_ = 1, parambox%n_particles
            call pbc(pos(idx_, :), parambox%box)
        end do

        call calc_vdw_force(pos=pos, cutoff=parambox%cutoff_set, forces=actual_force, boundary=parambox%box)
        do stp = 1, parambox%n_steps

            ! Calculem r(t)
            pos = pos + (vel*parambox%timestep) + ((actual_force * 0.5_dp) * (parambox%timestep ** 2))

            ! Apliquem pbcs
            do idx_ = 1, parambox%n_particles
                call pbc(pos(idx_, :), parambox%box)
            end do

            ! Calculem forces a r(t+dt)
            call calc_vdw_force(pos=pos, cutoff=parambox%cutoff_set, forces=new_force, boundary=parambox%box)
            vel = vel + (((actual_force + new_force) * 0.5_dp) * parambox%timestep)

            ! Apliquem el thermostat
            if (thermost) then
                call andersen_thermostat(vel, parambox)
            end if

            ! We write information of the system to the log unit
            if (mod(stp, parambox%write_file) == 0_i64) then
                call write_system_information(pos=pos, vel=vel, frame=stp, unit=log_unit, parambox=parambox)
            end if

            if (present_rdf) then
                if (mod(stp, parambox%write_stats) == 0) then
                    call write_rdf(stp=stp, parambox=parambox, pos=pos, &
                    write_unit=rdf_unit)
                end if
            end if

            if (present_msd) then
                if (mod(stp, parambox%write_stats) == 0) then
                    call write_msd(stp=stp, parambox=parambox, pos=pos, initpos=init_pos,&
                    write_unit=msd_unit)
                end if
            end if

            ! Swapping de matrius
            actual_force = new_force
        end do

        deallocate(actual_force)
        deallocate(new_force)
    end subroutine velocity_verlet

    subroutine euler_integrator(pos, vel, parambox, log_unit)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(inout), dimension(:,:) :: pos, vel
        type(databloc_params_t), intent(in)          :: parambox
        integer(kind=i64), intent(in)                :: log_unit
        ! Function variables
        integer(kind=I64)                            :: idx_, stp
        real(kind=DP), dimension(:, :), allocatable  :: actual_force

        allocate(actual_force(parambox%n_particles, 3))


        ! Posem tot abans a la caixa unitaria
        do idx_ = 1, parambox%n_particles
            call pbc(pos(idx_, :), parambox%box)
        end do

        actual_force = 0.0_dp

        do stp = 1, parambox%n_steps
            
            ! Compute forces
            call calc_vdw_force(pos=pos, cutoff=parambox%cutoff_set, forces=actual_force, boundary=parambox%box)
        
            ! Update positions
            pos = pos + (vel * parambox%timestep) + (0.5_DP * actual_force * (parambox%timestep**2))
            
            ! Apliquem pbcs
            do idx_ = 1, parambox%n_particles
                call pbc(pos(idx_, :), parambox%box)
            end do
            
            vel = vel + (actual_force * parambox%timestep)
        
            ! Escribim al output
            if (mod(stp, parambox%write_file) == 0_i64) then
                call write_system_information(pos=pos, vel=vel, frame=stp, unit=log_unit, parambox=parambox)
            end if

        end do

        deallocate(actual_force)

    end subroutine euler_integrator

    subroutine compute_velocities(new_pos, prev_pos, dt, vel)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(in), dimension(:, :)   :: new_pos, prev_pos
        real(kind=DP), intent(inout), dimension(:,:) :: vel
        real(kind=DP), intent(in)                    :: dt

        vel = (new_pos - prev_pos) / (2.0_DP * dt)

    end subroutine compute_velocities

end module integrtors_m