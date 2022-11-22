module therm_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

contains

    function calc_pressure(dens, lenth, forces, positions, temp, cutoff) result(press)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in) :: forces, positions
        real(kind=dp), intent(in)                 :: dens, lenth, temp, cutoff
        real(kind=dp)                             :: press
        ! Internal variables
        integer(kind=i64)                         :: i_part, j_part, n_p
        real(kind=dp), dimension(3)               :: rij, fij
        real(kind=dp)                             :: dij

        press = 0.0_dp
        n_p = size(positions, dim=1, kind=i64)

        do i_part = 1, n_p
            do j_part = i_part+1, n_p

                ! Calculem rij
                rij(1) = positions(i_part, 1) - positions(j_part, 1)
                rij(2) = positions(i_part, 2) - positions(j_part, 2)
                rij(3) = positions(i_part, 3) - positions(j_part, 3)

                ! TODO afegir les pbcs un cop hagis fet el modul per pbcs

                ! Calculem fij
                dij = norm2(rij)

                if (cutoff > dij) cycle

                fij(1) = (48.0_dp / dij**14.0_dp - 24.0_dp /dij**8)*rij(1)
                fij(2) = (48.0_dp / dij**14.0_dp - 24.0_dp /dij**8)*rij(2)
                fij(3) = (48.0_dp / dij**14.0_dp - 24.0_dp /dij**8)*rij(3)

                ! Upgredagem el valor de la pressió
                press = press + (dot_product(rij, fij))

            end do
        end do

        press = (press*dens*temp) + (1.0_dp / (3.0_dp * (lenth**3)))
    end function calc_pressure

end module therm_m