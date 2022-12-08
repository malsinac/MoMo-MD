module writers_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: therm_m, only: calc_vdw_pbc, calc_KE, compute_com_momenta, calc_pressure
    implicit none

    public :: write_frame, write_system_information, write_velocities

contains

    subroutine write_frame(unit_nr, pos, stp_c)
        implicit none
        integer(kind=I64), intent(in)             :: unit_nr, stp_c
        real(kind=DP), dimension(:,:), intent(in) :: pos
        integer(kind=I64) :: n_p, i_aux

        n_p = size(array=pos, dim=1, kind=I64)
        write(unit=unit_nr, fmt='(I3)') n_p
        write(unit=unit_nr, fmt='(I6)') stp_c

        do i_aux=1, n_p
            write(unit=unit_nr, fmt='(A,F20.8,F20.8,F20.8)') "Ar ", pos(i_aux, 1), pos(i_aux, 2), &
        pos(i_aux, 3)
        end do
        write(unit=unit_nr, fmt='(A)') ""
    end subroutine write_frame

    subroutine write_system_information(pos, vel, cutoff, frame, unit, boundary, mass, dens)
        implicit none
        ! In/Out variables
        real(kind=dp), intent(in), dimension(:,:) :: pos, vel
        real(kind=dp), intent(in)                 :: cutoff, boundary, mass, dens
        integer(kind=i64), intent(in)             :: frame, unit
        ! Internal variables
        real(kind=dp)                 :: pe_calc, ke_calc, temper, com_vel_mod, calc_press, time
        real(kind=dp), dimension(3)   :: com_vector
        integer(kind=i64)             :: n_p

        n_p = size(pos, dim=1, kind=i64)

        pe_calc = calc_vdw_pbc(pos=pos, cutoff=cutoff, boundary=boundary)
        ke_calc = calc_KE(vel)
        temper = (2.0_DP / (3.0_DP * n_p)) * ke_calc
        call compute_com_momenta(vel=vel, com_momenta=com_vector, mass=mass)
        com_vel_mod = norm2(com_vector)
        calc_press = calc_pressure(dens=dens, lenth=boundary, positions=pos, temp=temper, cutoff=cutoff)
        ! time = real(frame, kind=dp) * 

        write(unit=unit, fmt='(A,I6,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4)') "t: ", frame, " KE=", ke_calc, " PE=", pe_calc, &
        " H=",ke_calc+pe_calc, " T=", temper, " P=", calc_press, " COM mod=", com_vel_mod

    end subroutine write_system_information

    subroutine write_velocities(vel, unit_nr, step)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in) :: vel
        integer(kind=i64), intent(in)             :: unit_nr, step
        ! Internal variables
        integer(kind=I64) :: n_p, i_aux

        n_p = size(vel, dim=1, kind=i64)

        write(unit=unit_nr, fmt='(I6)') step
        do i_aux=1, n_p
            write(unit=unit_nr, fmt='(F14.8,F14.8,F14.8,F14.8)') vel(i_aux,1), vel(i_aux,2), vel(i_aux,3), norm2(vel(i_aux,:))
        end do

    end subroutine write_velocities

    !subroutine write_rdf()

    !end subroutine write_rdf



end module writers_m