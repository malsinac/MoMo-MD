module writers_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use            :: therm_m,         only: calc_vdw_pbc, calc_KE, compute_com_momenta, calc_pressure
    use            :: interface_m,     only: databloc_params_t
    use            :: analysis_m,      only: g_r, calc_msd
    implicit none

    public :: write_frame_file, write_system_information, write_velocities

contains

    subroutine write_frame_file(unit_nr, pos, stp_c)
        implicit none
        integer(kind=I64), intent(in)             :: unit_nr, stp_c
        real(kind=DP), dimension(:,:), intent(in) :: pos
        integer(kind=I64)                         :: n_p, i_aux

        n_p = size(array=pos, dim=2, kind=I64)
        write(unit=unit_nr, fmt='(I3)') n_p
        write(unit=unit_nr, fmt='(I6)') stp_c

        do i_aux=1, n_p
            write(unit=unit_nr, fmt='(A,F20.8,F20.8,F20.8)') "Ar ", pos(1, i_aux), pos(2, i_aux), pos(3, i_aux)
        end do
    end subroutine write_frame_file

    subroutine write_system_information(pos, vel, frame, unit, parambox)
        implicit none
        ! In/Out variables
        real(kind=dp), intent(in), dimension(:,:) :: pos, vel
        integer(kind=i64), intent(in)             :: frame, unit
        type(databloc_params_t), intent(in)       :: parambox
        ! Internal variables
        real(kind=dp)                             :: pe_calc, ke_calc, temper, com_vel_mod, calc_press, time
        real(kind=dp), dimension(3)               :: com_vector
        real(kind=dp), parameter                  :: kb = 8.314462618e-3_dp  ! [kJ / mol K]
        real(kind=dp), parameter                  :: na = 6.02214076e23_dp

        ! Energies
        pe_calc = calc_vdw_pbc(pos=pos, cutoff=parambox%cutoff_set, boundary=parambox%box)  ! u'
        ke_calc = calc_KE(vel)
        
        ! Temperatures
        temper = (2.0_DP / (3.0_DP * real(parambox%n_particles, kind=dp))) * ke_calc  ! T'
        call compute_com_momenta(vel=vel, com_momenta=com_vector)
        com_vel_mod = norm2(com_vector)
        
        ! Pressio
        calc_press = calc_pressure(lenth=parambox%box, positions=pos, temp=temper, cutoff=parambox%cutoff_set) ! P'
        
        ! Temps
        time = real(frame, kind=dp) * parambox%timestep

        ! We do the transformation to real units. Those are
        !    P = [Pa], E = [kJ/mol], dens = [g/cm??], T = [K], t = [ps]
        temper = temper * (parambox%lj_epsilon / kb)  ! [K]
        calc_press = (calc_press * parambox%lj_epsilon * 1.0e33_dp) / (na * (parambox%lj_sigma ** 3)) ! [Pa]
        pe_calc = pe_calc * parambox%lj_epsilon  ! [kJ/mol]
        ke_calc = ke_calc * parambox%lj_epsilon  ! [kJ/mol]
        time = (time / sqrt(parambox%lj_epsilon / ((parambox%mass * 1.0e-6_dp) * ((parambox%lj_sigma * 1.0e-10_dp)**2)))) * 1.0e12_dp  ! [ps]


        write(unit=unit, fmt='(A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4,A,ES18.8e4)') "t: ", time, " KE=", ke_calc, " PE=", pe_calc, &
        " H=", ke_calc+pe_calc, " T=", temper, " P=", calc_press, " COM mod=", com_vel_mod

    end subroutine write_system_information

    subroutine write_velocities(vel, unit_nr, step)
        implicit none
        ! In/Out variables
        real(kind=dp), dimension(:,:), intent(in) :: vel
        integer(kind=i64), intent(in)             :: unit_nr, step
        ! Internal variables
        integer(kind=I64) :: n_p, i_aux

        n_p = size(vel, dim=2, kind=i64)

        write(unit=unit_nr, fmt='(I6)') step
        do i_aux=1, n_p
            write(unit=unit_nr, fmt='(F14.8,F14.8,F14.8,F14.8)') vel(1, i_aux), vel(2, i_aux), vel(3, i_aux), norm2(vel(:, i_aux))
        end do

    end subroutine write_velocities

    subroutine write_rdf(stp, parambox, pos, write_unit)
        implicit none
        ! In/Out variables
        integer(kind=i64), intent(in)              :: stp
        type(databloc_params_t), intent(in)        :: parambox
        real(kind=dp), intent(in), dimension(:,:)  :: pos
        integer(kind=i64), intent(in)              :: write_unit
        ! Internal variables
        real(kind=dp), allocatable, dimension(:,:) :: gdr_mat
        integer(kind=i64)                          :: i_aux

        !if (.NOT. allocated(gdr_mat)) then
            allocate(gdr_mat(2, parambox%gdr_num_bins))
        !end if

        call g_r(gr_mat=gdr_mat, pos=pos, parambox=parambox)

        if (stp == parambox%write_stats) then
            ! Write the distance binning
            do i_aux = 1, parambox%gdr_num_bins
                write(unit=write_unit, fmt='(ES18.8e4)', advance='no') gdr_mat(1, i_aux)
            end do
            write(unit=write_unit, fmt='(A)') ""
        end if

        do i_aux = 1, parambox%gdr_num_bins
            write(unit=write_unit, fmt='(ES18.8e4)', advance='no') gdr_mat(2, i_aux)
        end do
        write(unit=write_unit, fmt='(A)') ""


        deallocate(gdr_mat)
    end subroutine write_rdf

    subroutine write_msd(stp, parambox, pos, initpos, write_unit)
        implicit none
        ! In/Out variables
        integer(kind=i64), intent(in)             :: stp, write_unit
        type(databloc_params_t), intent(in)       :: parambox
        real(kind=dp), dimension(:,:), intent(in) :: pos, initpos
        ! Internal variables
        real(kind=dp)                             :: msd_val, time

        call calc_msd(msd_vec=msd_val, pos=pos, init_pos=initpos, box=parambox%box, lj_sigma=parambox%lj_sigma)

        time = real(stp, kind=dp) * parambox%timestep
        time = (time / sqrt(parambox%lj_epsilon / ((parambox%mass * 1.0e-6_dp) * (parambox%lj_sigma * 1.0e-10_dp)**2))) * 1.0e12_dp  ! [ps]

        ! Output is ps vs A


        write(unit=write_unit, fmt='(ES18.8e4,ES18.8e4)') time, msd_val

    end subroutine write_msd



end module writers_m