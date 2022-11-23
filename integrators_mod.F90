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

end module integrtors_m