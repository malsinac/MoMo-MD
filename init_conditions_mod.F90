module init_cond_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: bimodal_dist_velocities, init_positions_sc

contains

    pure subroutine bimodal_dist_velocities(vel, temp)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: vel
        real(kind=dp), intent(in)                    :: temp
        ! Internal variables
        real(kind=dp)                                :: sqrt_temp, total_mass
        integer(kind=i64)                            :: i_pos, j_dim, n_p, mid_idx, k_counter
        real(kind=dp), dimension(3)                  :: v_cm

        vel = 0.0_dp
        n_p = size(vel, dim=2, kind=i64)
        sqrt_temp = sqrt(temp)
        mid_idx = ceiling(real(n_p, kind=dp) * 3.0_dp / 2.0_dp, kind=i64)
        k_counter = 1_i64
        v_cm = 0.0_dp
        total_mass = 0.0_dp

        ! Assignem velocitats a les particules
        do i_pos = 1, n_p
            do j_dim = 1, 3
                if (k_counter <= mid_idx) then 
                    vel(j_dim, i_pos) = sqrt_temp
                else if (k_counter > mid_idx) then
                    vel(j_dim, i_pos) = -sqrt_temp
                end if
                k_counter = k_counter + 1_i64
            end do
        end do

        ! Calculem la velocitat del centre de masses
        do i_pos = 1, n_p
            v_cm(1) = v_cm(1) + vel(1, i_pos)
            v_cm(2) = v_cm(2) + vel(2, i_pos)
            v_cm(3) = v_cm(3) + vel(3, i_pos)
            total_mass = total_mass + 1.0_dp
        end do

        v_cm = v_cm / total_mass

        ! Restem la velocitat de centre de masses a totes les velocitats perque aquesta sigui 0
        do i_pos = 1, n_p
            vel(1, i_pos) = vel(1, i_pos) - v_cm(1)
            vel(2, i_pos) = vel(2, i_pos) - v_cm(2)
            vel(3, i_pos) = vel(3, i_pos) - v_cm(3)
        end do

    end subroutine bimodal_dist_velocities

    subroutine init_positions_sc(rho, pos, gen_log_unt)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(inout), dimension(:, :) :: pos
        real(kind=DP), intent(in)                     :: rho
        integer(kind=i64), optional                   :: gen_log_unt
        ! Internal variables
        integer(kind=I64)                             :: M, i, j, k, p, N
        real(kind=DP)                                 :: L, a

        N = size(pos, dim=2, kind=I64)
        M = int(N ** (1.0_DP / 3.0_DP), kind=I64) + 1
        L = (real(N, kind=DP)/rho) ** (1.0_DP/3.0_DP)
        a = L/M

        if (present(gen_log_unt)) then
            write(unit=gen_log_unt, fmt='(A,I3,A,I3,A,F12.8,A,F12.8)') "sc lattice, parameters: N=", N, " M=", M, " L=", L, " a=", a
        end if

        do i = 0, M - 1
            do j = 0, M - 1
                do k = 0, M - 1
                    p = (k+1) + (j*M) + (i * (M**2))
                    pos(:, p) = [real(i, kind=DP)*a, real(j, kind=DP)*a, real(k, kind=DP)*a]
                end do
            end do
        end do

    end subroutine    

end module init_cond_m