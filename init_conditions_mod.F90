module init_cond_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: bimodal_dist_velocities

contains

    subroutine bimodal_dist_velocities(vel, temp)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: vel
        real(kind=dp), intent(in)                    :: temp
        ! Internal variables
        real(kind=dp)                                :: sqrt_temp, u_num
        integer(kind=i64)                            :: i_pos, j_dim, n_p, k_val

        n_p = size(vel, dim=1)
        sqrt_temp = sqrt(temp)
        
        do i_pos = 1, n_p
            do j_dim = 1, 3
                call random_number(u_num)
                k_val = -1 + floor(3.0_dp*u_num, kind=i64)
                vel(i_pos,j_dim) = k_val * sqrt_temp
            end do
        end do
    end subroutine bimodal_dist_velocities

end module init_cond_m