module analysis_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

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
        integer(kind=i64)                            :: n_bins, i_ax, index_mat
        real(kind=dp)                                :: dr, dist
        real(kind=dp), parameter                     :: PI = acos(-1.0_dp)

        n_bins = size(gr_mat, dim=2, kind=i64)
        dr = max_dist / real(n_bins, kind=dp)

        gr_mat(1,:) = [(real(i_ax, dp)*dr, i_ax=0, n_bins-1)]
        gr_mat(2,:) = 0.0_dp

        ! Calculem n(r)
        do i_ax = 2, size(pos, dim=1)
            dist = norm2(pos(1,:) - pos(i_ax,:))
            index_mat = ceiling(dist/dr)
            gr_mat(2, index_mat) = gr_mat(2, index_mat) + 1.0_dp
        end do

        ! Calculem g(r)
        do i_ax = 1, n_bins
            gr_mat(2, i_ax) = gr_mat(2, i_ax) / (dr * (gr_mat(1,i_ax)**2))
        end do

        gr_mat(2,:) = gr_mat(2,:) * (1.0_dp / (dens*4.0_dp*PI))

    end subroutine g_r

end module analysis_m

