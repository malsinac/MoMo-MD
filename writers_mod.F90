module writers_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    implicit none

    public :: write_frame, write_system_information

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
    end subroutine write_frame

    subroutine write_system_information(kinetic_energy, potential_energy, frame, unit)
        implicit none
        real(kind=dp), intent(in)     :: kinetic_energy, potential_energy
        integer(kind=i64), intent(in) :: frame, unit

        write(unit=unit, fmt=*) "frame: ", frame, " KE=", kinetic_energy, " PE=", potential_energy, " H=",kinetic_energy+potential_energy

    end subroutine write_system_information

    subroutine write_velocities(vel, unit_nr)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in) :: vel
        integer(kind=i64), intent(in)             :: unit_nr
        ! Internal variables
        integer(kind=I64) :: n_p, i_aux

        do i_aux=1, n_p
            write(unit=unit_nr, fmt='(F12.8,F12.8,F12.8,F12.8)') vel(i_aux,1), vel(i_aux,2), vel(i_aux,3), norm2(vel(i_aux,:))
        end do

    end subroutine write_velocities



end module writers_m