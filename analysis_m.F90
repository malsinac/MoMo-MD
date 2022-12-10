module analysis_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: therm_m,         only: pbc
    use            :: interface_m,     only: databloc_params_t
    implicit none

    public :: g_r, calc_msd

contains

    subroutine g_r(gr_mat, dens, pos, parambox)
        ! Notes
        ! gr_mat(1,:) -> valors de distancia
        ! gr_mat(2,:) -> numero de partícules a aquesta distancia
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(inout) :: gr_mat
        real(kind=dp), dimension(:,:), intent(in)    :: pos
        real(kind=dp), value                         :: dens
        type(databloc_params_t), intent(in)          :: parambox
        ! Internal variables
        integer(kind=i64)                            :: n_bins, i_ax, index_mat, j_ax, n_p
        real(kind=dp)                                :: dr, dist, div_bin
        real(kind=dp), parameter                     :: PI = acos(-1.0_dp)
        real(kind=dp), dimension(3)                  :: rij

        n_bins = size(gr_mat, dim=2, kind=i64)
        n_p = size(pos, dim=1, kind=i64)
        dr = parambox%gdr_max_dist / real(n_bins, kind=dp)

        gr_mat(1,:) = [(real(i_ax, dp)*dr, i_ax=1, n_bins)]
        gr_mat(2,:) = 0.0_dp

        ! Calculem n(r)
        do j_ax = 1, n_p - 1
            do i_ax = j_ax + 1, n_p
                ! Calculem rij
                rij = pos(j_ax,:) - pos(i_ax,:)
                call pbc(rij, parambox%box)
                dist = norm2(rij)
                
                ! Apliquem el cutoff de maxima distancia
                if (dist > parambox%gdr_max_dist) cycle
                index_mat = int(dist/dr, kind=i64) + 1_i64
                gr_mat(2, index_mat) = gr_mat(2, index_mat) + 2.0_dp
            end do
        end do

        ! Calculem g(r) en unitats reals
        dens = dens * (parambox%lj_sigma ** 3) ! \rho' -> \rho
        dr = dr * parambox%lj_sigma            ! dr' -> dr
        do i_ax = 1, n_bins
            associate(r => gr_mat(1, i_ax), gdr => gr_mat(2, i_ax))
                r = r * parambox%lj_sigma              ! r' -> r
                div_bin = dens * 4.0_dp * pi * r * r * dr
                gdr = gdr / div_bin
            end associate
        end do
        gr_mat(1,:) = gr_mat(1,:) * parambox%lj_sigma

    end subroutine g_r

    subroutine calc_msd(msd_vec, pos, init_pos)
        ! Necessito guardar les posicions per a cada particula per cada frame
        implicit none
        ! In/Out variables
        real(kind=dp), intent(out)                  :: msd_vec
        real(kind=dp), dimension(:,:), intent(in)   :: pos, init_pos
        ! Internal variables
        integer(kind=i64)                           :: j_part, n_p
        real(kind=dp)                               :: dd

        n_p = size(pos, dim=2)
        msd_vec = 0.0_dp

        do j_part = 1, n_p
            dd = norm2(pos(j_part, :) - init_pos(j_part, :)) ** 2
            msd_vec = msd_vec + dd
        end do

    end subroutine calc_msd

end module analysis_m

