module analysis_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: g_r, calc_msd

contains

    subroutine g_r(gr_mat, dens, pos, max_dist)
        ! Notes
        ! gr_mat(2,:) on la segona dimensió es el numero de bins que faré
        ! gr_mat(1,:) -> valors de distancia
        ! gr_mat(2,:) -> numero de partícules a aquesta distancia
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: gr_mat
        real(kind=dp), dimension(:,:)                :: pos
        real(kind=dp), intent(in)                    :: dens, max_dist
        ! Internal variables
        integer(kind=i64)                            :: n_bins, i_ax, index_mat, j_ax, n_p
        real(kind=dp)                                :: dr, dist
        real(kind=dp), parameter                     :: PI = acos(-1.0_dp)

        n_bins = size(gr_mat, dim=2, kind=i64)
        n_p = size(pos, dim=1, kind=i64)
        dr = max_dist / real(n_bins, kind=dp)

        gr_mat(1,:) = [(real(i_ax, dp)*dr, i_ax=1, n_bins)]
        gr_mat(2,:) = 0.0_dp

        ! Calculem n(r)
        do j_ax = 1, n_p
            do i_ax = 2, n_p
                dist = norm2(pos(j_ax,:) - pos(i_ax,:))
                index_mat = ceiling(dist/dr)
                gr_mat(2, index_mat) = gr_mat(2, index_mat) + 1.0_dp
            end do
        end do

        ! Calculem g(r)
        do i_ax = 1, n_bins
            gr_mat(2, i_ax) = gr_mat(2, i_ax) / (dr * (gr_mat(1,i_ax)**2))
        end do

        gr_mat(2,:) = gr_mat(2,:) * (1.0_dp / (dens*4.0_dp*PI))

    end subroutine g_r

    subroutine calc_msd(msd_vec, pos, init_pos)
        ! Necessito guardar les posicions per a cada particula per cada frame
        implicit none
        ! In/Out variables
        real(kind=dp), intent(out)                  :: msd_vec
        real(kind=dp), dimension(:,:), intent(in)   :: pos, init_pos
        ! Internal variables
        integer(kind=i64)                           :: num_frames, i_frame, j_part, n_p
        real(kind=dp)                               :: dd

        n_p = size(pos, dim=2)
        msd_vec = 0.0_dp

        do j_part = 1, n_p
            dd = norm2(pos(j_part, :) - init_pos) ** 2
            msd_vec = msd_vec + dd
        end do
        msd_vec = msd_vec / real(n_p, kind=dp)

    end subroutine calc_msd

end module analysis_m

