module therm_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: calc_pressure, calc_KE, calc_vdw_pbc, pbc, calc_vdw_force, compute_com_momenta

contains

    function calc_pressure(lenth, positions, temp, cutoff) result(press)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in) :: positions
        real(kind=dp), intent(in)                 :: lenth, temp, cutoff
        real(kind=dp)                             :: press
        ! Internal variables
        integer(kind=i64)                         :: i_part, j_part, n_p
        real(kind=dp), dimension(3)               :: rij, fij
        real(kind=dp)                             :: dij, cutoff2, virial, vol

        press = 0.0_dp
        virial = 0.0_dp
        cutoff2 = cutoff * cutoff
        n_p = size(positions, dim=2, kind=i64)
        vol = lenth ** 3
        rij = 0.0_dp
        fij = 0.0_dp

        do i_part = 1, n_p - 1
            do j_part = i_part + 1, n_p

                ! Calculem rij
                rij(1) = positions(1, i_part) - positions(1, j_part)
                rij(2) = positions(2, i_part) - positions(2, j_part)
                rij(3) = positions(3, i_part) - positions(3, j_part)

                call pbc(x=rij, l_=lenth)
                
                ! Calculem fij
                dij = (rij(1)**2) + (rij(2)**2) + (rij(3)**2)

                if (dij < cutoff2) then

                    fij(1) = (48.0_dp / dij**7 - 24.0_dp / dij**4) * rij(1)
                    fij(2) = (48.0_dp / dij**7 - 24.0_dp / dij**4) * rij(2)
                    fij(3) = (48.0_dp / dij**7 - 24.0_dp / dij**4) * rij(3)

                    ! Upgredagem el valor de la pressio
                    virial = virial + dot_product(rij, fij)
                end if

            end do
        end do
        press = (real(n_p, kind=dp)*temp / vol) + ((1.0_dp / (3.0_dp * vol)) * virial)
    end function calc_pressure

    pure function calc_KE(vel) result(ke)
        implicit none
        ! In/Out variables
        real(kind=DP)                             :: ke
        real(kind=DP), dimension(:,:), intent(in) :: vel
        ! Internal variables
        integer(kind=I64)                         :: n_p, i

        ! Variable initialization
        n_p = size(vel, dim=2, kind=I64)
        ke = 0.0_DP
        
        do i = 1, n_p
            ke = ke + sum(vel(:, i)**2)
        end do

        ! Final multiplications
        ke = ke * 0.5_DP
    end function calc_KE

    function calc_vdw_pbc(pos, cutoff, boundary) result(vdw_calc)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(in), dimension(:, :) :: pos
        real(kind=DP)                              :: vdw_calc
        real(kind=DP), intent(in)                  :: cutoff, boundary
        ! Internal variables
        real(kind=DP)                              :: dist, ecalc, cutoff2
        integer(kind=I64)                          :: i, j, num_particles
        real(kind=DP), dimension(3)                :: rij

        ! Initializing parameters
        vdw_calc = 0.0_DP
        num_particles = size(pos, dim=2, kind=I64)

        cutoff2 = cutoff ** 2
        rij = 0.0_dp

        do j = 1, num_particles - 1
            do i = j + 1, num_particles
                rij(1) = pos(1, j) - pos(1, i)
                rij(2) = pos(2, j) - pos(2, i)
                rij(3) = pos(3, j) - pos(3, i)
                
                call pbc(x=rij, l_=boundary)

                dist = (rij(1)**2) + (rij(2)**2) + (rij(3)**2)

                if (dist < cutoff2) then
                    
                    dist = dist ** 3
                    
                    ecalc = (4.0_DP * ((1.0_DP / (dist**2)) - (1.0_DP/(dist)))) - &
                            (4.0_DP * ((1.0_DP / (cutoff2**6)) - (1.0_DP/(cutoff2**3))))
            
                    vdw_calc = vdw_calc + ecalc
                
                end if
            end do
        end do
    end function calc_vdw_pbc

    subroutine calc_vdw_force(pos, cutoff, forces, boundary)
        implicit none
        !In/Out variables
        real(kind=DP), intent(in), dimension(:, :)    :: pos
        real(kind=DP), intent(in)                     :: cutoff, boundary
        real(kind=DP), intent(out), dimension(:, :)   :: forces
        ! Internal variables
        integer(kind=I64)           :: n_part, i, j
        real(kind=DP)               :: dist, cutoff2
        real(kind=DP), dimension(3) :: rij

        n_part = size(pos, dim=2, kind=I64)

        forces = 0.0_DP
        rij = 0.0_dp
        cutoff2 = cutoff * cutoff

        do i = 1, n_part - 1
            do j = i + 1, n_part
                
                rij(1) = pos(1, i) - pos(1, j)
                rij(2) = pos(2, i) - pos(2, j)
                rij(3) = pos(3, i) - pos(3, j)
                
                call pbc(x=rij, l_=boundary)
                
                dist = (rij(1)**2) + (rij(2)**2) + (rij(3)**2)
                
                if (dist < cutoff) then
                    
                    ! Calculem la forc entre particula i, j

                    forces(1, i) = forces(1, i) + (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(1)
                    forces(2, i) = forces(2, i) + (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(2)
                    forces(3, i) = forces(3, i) + (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(3)

                    forces(1, j) = forces(1, j) - (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(1)
                    forces(2, j) = forces(2, j) - (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(2)
                    forces(3, j) = forces(3, j) - (48.0_DP / dist**7 - 24.0_dp / dist**4) * rij(3)
                end if
            end do
        end do
    end subroutine calc_vdw_force

    pure subroutine pbc(x, l_)
        ! In/Out variables
        real(kind=DP), intent(inout), dimension(3) :: x
        real(kind=DP), intent(in)                  :: l_
        ! Internal variables
        integer(kind=I64) :: i_

        do i_ = 1, 3
            if (x(i_) > l_ / 2.0_DP) then
                x(i_) = x(i_) - l_
            end if
            if (x(i_) < - l_ / 2.0_DP) then
                x(i_) = x(i_) + l_
            end if
        end do
    end subroutine pbc

    pure subroutine compute_com_momenta(vel, com_momenta)
        implicit none
        ! In/Out variables
        real(kind=DP), dimension(:, :), intent(in) :: vel
        real(kind=DP), dimension(3), intent(out)   :: com_momenta
        ! Internal variables
        integer(kind=I64)                          :: i_aux, n_p
        real(kind=dp)                              :: total_mass

        com_momenta = 0.0_DP
        total_mass = 0.0_dp
        n_p = size(vel, dim=2, kind=I64)

        do i_aux = 1, n_p
            com_momenta(1) = com_momenta(1) + vel(1, i_aux)
            com_momenta(2) = com_momenta(2) + vel(2, i_aux)
            com_momenta(3) = com_momenta(3) + vel(3, i_aux)
            total_mass = total_mass + 1.0_dp
        end do

        com_momenta = com_momenta / total_mass

    end subroutine compute_com_momenta

end module therm_m