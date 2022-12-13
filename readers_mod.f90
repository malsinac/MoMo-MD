module readers_m 
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64, error_unit
    use            :: interface_m    , only: databloc_params_t
    implicit none

contains

    subroutine read_nml(param_file, datablock)
        implicit none
        ! In/Out variables
        type(databloc_params_t), intent(out) :: datablock
        character(len=2048), intent(in)      :: param_file
        ! Internal variables
        integer(kind=I64)                :: unit_nr, iost
        character(len=1024)              :: msg
        ! Namelist variables
        real(kind=dp)                    :: lj_epsilon, lj_sigma, mass, timestep, density, andersen_nu
        integer(kind=i64)                :: n_particles, n_steps, write_file, write_stats, gdr_num_bins, write_frame
        character(len=2048)              :: sim_name

        namelist /PARAMS/ lj_epsilon, lj_sigma, mass, &
        timestep, n_particles, density, n_steps, &
        write_file, write_stats, write_frame, sim_name, &
        andersen_nu, &
        gdr_num_bins

        ! Checking if file is present
        inquire(file=trim(param_file), iostat=iost)
        if (iost /= 0_I64) then
            write(unit=error_unit, fmt='(A)') "FATAL ERROR: namelist file not found"
            stop 1
        end if

        ! Reading namelist
        open(newunit=unit_nr, file=trim(param_file), access='sequential', &
        action='read', iostat=iost, iomsg=msg, form='formatted')

        read(iostat=iost, unit=unit_nr, nml=PARAMS, iomsg=msg)
        if (iost /= 0_I64) then
            write(unit=error_unit, fmt='(A,A)') "Error when loading parameters on the namelist file: ", trim(msg)
            stop 1
        end if
        close(unit_nr)

        ! Assigning to the datablock and defining parameters
        datablock%mass = mass
        datablock%lj_epsilon = lj_epsilon
        datablock%lj_sigma = lj_sigma
        datablock%timestep = timestep
        datablock%n_particles = n_particles
        datablock%density = density
        datablock%write_file = write_file
        datablock%write_stats = write_stats
        datablock%write_frame = write_frame
        datablock%sim_name = sim_name
        datablock%andersen_nu = andersen_nu
        datablock%gdr_num_bins = gdr_num_bins
        datablock%box = (real(datablock%n_particles, kind=DP)/datablock%density) ** (1.0_DP / 3.0_DP)
        datablock%cutoff_set = 0.5_dp * datablock%box
        datablock%gdr_max_dist = datablock%box

    end subroutine read_nml

end module readers_m