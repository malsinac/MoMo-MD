module init_cond_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: bimodal_dist_velocities

contains

    pure subroutine bimodal_dist_velocities(vel, temp)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: vel
        real(kind=dp), intent(in)                    :: temp
        ! Internal variables
        real(kind=dp)                                :: sqrt_temp, u_num
        integer(kind=i64)                            :: i_pos, j_dim, n_p, mid_idx, k_counter
        real(kind=dp), dimension(3)                  :: v_cm

        vel = 0.0_dp
        n_p = size(vel, dim=1)
        sqrt_temp = sqrt(temp)
        mid_idx = ceiling(real(n_p, kind=dp) * 3.0_dp / 2.0_dp, kind=i64)
        k_counter = 0_i64
        v_cm = 0.0_dp

        ! Assignem velocitats a les particules
        do i_pos = 1, n_p
            do j_dim = 1, 3
                if (k_counter <= mid_idx) then 
                    vel(i_pos, j_dim) = sqrt_temp
                else if (k_counter > mid_idx) then
                    vel(i_pos, j_dim) = -sqrt_temp
                end if
                k_counter = k_counter + 1_i64
            end do
        end do

        ! Calculem la velocitat del centre de masses
        do i_pos = 1, n_p
            v_cm(1) = v_cm(1) + vel(i_pos, 1)
            v_cm(2) = v_cm(2) + vel(i_pos, 2)
            v_cm(3) = v_cm(3) + vel(i_pos, 3)
        end do

        ! Restem la velocitat de centre de masses a totes les velocitats perque aquesta sigui 0
        do i_pos = 1, n_p
            vel(i_pos, 1) = vel(i_pos, 1) - v_cm(1)
            vel(i_pos, 2) = vel(i_pos, 2) - v_cm(2)
            vel(i_pos, 3) = vel(i_pos, 3) - v_cm(3)
        end do

    end subroutine bimodal_dist_velocities

end module init_cond_m