#define MASS 1.0_DP
#define TIMESTEP 0.0001_DP
#define CUTOFF_SET 2.5_DP

module therm_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: calc_pressure, calc_KE, calc_vdw_pbc, pbc, calc_vdw_force

contains

    function calc_pressure(dens, lenth, positions, temp, cutoff) result(press)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in) :: positions
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

                call pbc(rij, lenth)
                
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

    pure function calc_KE(vel) result(ke)
        implicit none
        ! In/Out variables
        real(kind=DP) :: ke
        real(kind=DP), dimension(:,:), intent(in) :: vel
        ! Internal variables
        integer(kind=I64) :: n_p, i

        ! Variable initialization
        n_p = size(vel, dim=1, kind=I64)
        ke = 0.0_DP
        
        do i = 1, n_p
            ke = ke + sum(vel(i, :)**2)
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
        num_particles = size(pos, dim=1, kind=I64)

        cutoff2 = cutoff ** 2.0_DP
        rij = 0

        do j = 1, num_particles
            do i = j+1, num_particles
                rij(1) = pos(i, 1) - pos(j, 1)
                rij(2) = pos(i, 2) - pos(j, 2)
                rij(3) = pos(i, 3) - pos(j, 3)
                
                call pbc(x=rij, l_=boundary)

                dist = norm2(rij) ** 2

                if (dist < cutoff2) then
                    
                    dist = dist ** 3
                    
                    ecalc = (4.0_DP * ((1.0_DP / (dist**2)) - (1.0_DP/(dist)))) - &
                        (4.0_DP * ((1.0_DP / (cutoff2**6))-(1.0_DP/(cutoff2**3))))
            
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
        real(kind=DP), intent(inout), dimension(:, :) :: forces
        ! Internal variables
        integer(kind=I64) :: n_part, i, j
        real(kind=DP)     :: dist
        real(kind=DP), dimension(3) :: rij

        n_part = size(pos, dim=1, kind=I64)

        if (size(forces, dim=1, kind=I64) /= n_part) stop 1

        forces = 0.0_DP
        rij = 0

        do i = 1, n_part
            do j = i+1, n_part
                
                rij(1) = pos(i, 1) - pos(j, 1)
                rij(2) = pos(i, 2) - pos(j, 2)
                rij(3) = pos(i, 3) - pos(j, 3)
                
                call pbc(x=rij, l_=boundary)
                
                dist = norm2(rij)
                if (dist < cutoff) then
                    ! Calculem la força entre particula i i j

                    forces(i, 1) = forces(i, 1) + (48.0_DP / dist**14 - 24 / dist**8) * rij(1)
                    forces(i, 2) = forces(i, 2) + (48.0_DP / dist**14 - 24 / dist**8) * rij(2)
                    forces(i, 3) = forces(i, 3) + (48.0_DP / dist**14 - 24 / dist**8) * rij(3)

                    forces(j, 1) = forces(j, 1) - (48.0_DP / dist**14 - 24 / dist**8) * rij(1)
                    forces(j, 2) = forces(j, 2) - (48.0_DP / dist**14 - 24 / dist**8) * rij(2)
                    forces(j, 3) = forces(j, 3) - (48.0_DP / dist**14 - 24 / dist**8) * rij(3)
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
        
        !do i_ = 1, 3
        !    if (x( i_ ) > l_ ) then
        !         x( i_ ) = x( i_ ) - l_
        !    else if (x( i_ ) < 0.0_DP) then
        !        x( i_ ) = x( i_ ) + l_
        !    end if
        !end do

        do i_ = 1, 3
            if (x(i_) > l_ / 2.0_DP) then
                x(i_) = x(i_) - l_
            end if
            if (x(i_) < - l_ / 2.0_DP) then
                x(i_) = x(i_) + l_
            end if
        end do
    end subroutine pbc

end module therm_m