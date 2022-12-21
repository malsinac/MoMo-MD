module analysis_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: therm_m,         only: pbc
    use            :: interface_m,     only: databloc_params_t
    implicit none

    public :: g_r, calc_msd

contains

    subroutine g_r(gr_mat, pos, parambox)
        ! Notes
        ! gr_mat(1,:) -> valors de distancia
        ! gr_mat(2,:) -> numero de partícules a aquesta distancia
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: gr_mat
        real(kind=dp), dimension(:,:), intent(in)    :: pos
        type(databloc_params_t), intent(in)          :: parambox
        ! Internal variables
        integer(kind=i64)                            :: n_bins, i_ax, index_mat, j_ax, n_p
        real(kind=dp)                                :: dr, dist, dv, ndg
        real(kind=dp), parameter                     :: PI = 4.0_dp * atan(1.0_dp)
        real(kind=dp), dimension(3)                  :: rij

        n_bins = size(gr_mat, dim=2, kind=i64)
        n_p = parambox%n_particles
        dr = parambox%gdr_max_dist / real(n_bins, kind=dp)

        gr_mat(1,:) = [(real(i_ax, dp)*dr, i_ax=1, n_bins)]
        gr_mat(2,:) = 0.0_dp

        ! Calculem n(r)
        do j_ax = 1, n_p - 1
            do i_ax = j_ax + 1, n_p
                ! Calculem rij
                rij = pos(:, j_ax) - pos(:, i_ax)
                call pbc(x=rij, l_=parambox%box)
                dist = norm2(rij)
                
                ! Apliquem el cutoff de maxima distancia
                if (dist < parambox%cutoff_set) then
                    index_mat = int(dist/dr, kind=i64) + 1_i64
                    gr_mat(2, index_mat) = gr_mat(2, index_mat) + 2.0_dp
                end if
            end do
        end do

        ! Calculem g(r) en unitats reals
        do i_ax = 1, n_bins
            associate(r => gr_mat(1, i_ax), gdr => gr_mat(2, i_ax))
                dv = (((real(i_ax, kind=dp) + 1.0_dp) ** 3) - (real(i_ax, kind=dp) ** 3)) * (dr ** 3)
                ndg = (4.0_dp / 3.0_dp) * pi * dv * parambox%density
                gdr = gdr / (parambox%n_particles * ndg)
            end associate
        end do
        gr_mat(1,:) = gr_mat(1,:) * parambox%lj_sigma

    end subroutine g_r

    subroutine calc_msd(msd_vec, pos, init_pos, box, lj_sigma)
        implicit none
        ! In/Out variables
        real(kind=dp), intent(out)                  :: msd_vec
        real(kind=dp), dimension(:,:), intent(in)   :: pos, init_pos
        real(kind=dp), intent(in)                   :: box, lj_sigma
        ! Internal variables
        integer(kind=i64)                           :: j_part, n_p
        real(kind=dp)                               :: dd
        real(kind=dp), dimension(3)                 :: rij

        n_p = size(pos, dim=2, kind=i64)
        msd_vec = 0.0_dp

        do j_part = 1, n_p
            rij = pos(:, j_part) - init_pos(:, j_part)
            call pbc(rij, box)
            dd = (rij(1) ** 2) + (rij(2) ** 2) + (rij(3) ** 2)
            msd_vec = msd_vec + dd
        end do

        msd_vec = (msd_vec * lj_sigma) / real(n_p, kind=dp)

    end subroutine calc_msd

end module analysis_m

