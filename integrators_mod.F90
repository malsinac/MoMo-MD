#define MASS 1.0_DP
#define TIMESTEP 0.0001_DP
#define CUTOFF_SET 2.5_DP

module integrtors_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: therm_m, only: pbc, calc_vdw_force

contains

    subroutine velocity_verlet(pos, vel, dt, boundary)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(inout), dimension(:,:) :: pos, vel
        real(kind=DP), intent(in)                    :: dt, boundary
        ! Internal_variables
        integer(kind=I64)                             :: n_p, idx_
        real(kind=DP), dimension(:, :), allocatable   :: actual_force, new_force

        n_p = size(pos, dim=1, kind=I64)

        allocate(actual_force(n_p, 3))
        allocate(new_force(n_p, 3))

        ! Calculem forces a r(t)
        call calc_vdw_force(pos=pos, cutoff=CUTOFF_SET, forces=actual_force, boundary=boundary)

        pos = pos + (vel*dt) + ((actual_force/(2.0_DP * MASS)) * (dt ** 2))

        ! Apliquem pbcs
        do idx_ = 1, n_p
            call pbc(pos(idx_, :), boundary)
        end do

        ! Calculem forces a r(t+dt)
        call calc_vdw_force(pos=pos, cutoff=CUTOFF_SET, forces=new_force, boundary=boundary)
        vel = vel + (((actual_force + new_force) / (2.0_DP * MASS)) * dt)

        deallocate(actual_force)
        deallocate(new_force)
    end subroutine velocity_verlet

    subroutine euler_integrator(pos, new_pos, vel, new_vel, boundary, dt, forc)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(in), dimension(:,:)    :: pos, vel, forc
        real(kind=DP), intent(out), dimension(:,:)   :: new_pos, new_vel
        real(kind=DP), intent(in)                    :: boundary, dt
        ! Function variables
        integer(kind=I64)                           :: n_p, i

        n_p = size(pos, dim=1, kind=I64)
        
        ! Update positions
        new_pos = pos + (vel * dt) + (0.5_DP * forc * dt * dt)
        do i = 1, n_p
            call pbc(new_pos(i, :), boundary)
        end do

        new_vel = vel + ((forc/MASS) * dt)

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