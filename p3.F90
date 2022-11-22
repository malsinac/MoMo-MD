#define MASS 1.0_DP
#define TIMESTEP 0.0001_DP
#define CUTOFF_SET 2.5_DP

program main
    use, intrinsic :: iso_fortran_env, only: DP => REAL64, I64 => INT64
    implicit none

    ! Scalar values 
    real(kind=DP)       :: box, ke_calc, pe_calc, density, time0, time1, instant_temperature
    integer(kind=I64)   :: n_particles, n_steps, stp, traj_unit, write_file, log_unit

    ! Arrays
    real(kind=DP), allocatable, dimension(:, :) :: positions, velocities


    ! Files
    open(newunit=traj_unit, file="heating.xyz", access='sequential', action='write',&
    status='replace', form='formatted')
    open(newunit=log_unit, file="heating.log", access='sequential', action='write',&
    status='replace', form='formatted')

    ! ~ Definim parametres ~
    write_file = 8_I64
    n_particles = 125_I64
    density = 0.7_DP
    n_steps = 10000_I64
    box = (real(n_particles, kind=DP)/density) ** (1.0_DP / 3.0_DP)

    ! ~ Inicialitzem les posicions i velocitats~
    allocate(positions(n_particles, 3), velocities(n_particles, 3))
    call init_positions_sc(rho=density, pos=positions)

    call init_velocities(vel=velocities)

    print '(A,F16.8,A,F16.8)', "Initial potential energy=", calc_vdw_pbc(pos=positions, cutoff=CUTOFF_SET, boundary=box), " Initial kinetic energy=", calc_KE(velocities)

    call write_frame(unit_nr=traj_unit, pos=positions, stp_c=0_I64)

    call cpu_time(time0)

    ! ~ Integrem les posicions dels atoms ~
    do stp = 1, n_steps

        ! Aplico integrador de velocity verlet
        call velocity_verlet(pos=positions, vel=velocities, dt=TIMESTEP, boundary=box)

        ! Aplico el termostat
        call andersen_thermostat(vel=velocities, temp=100.0_DP, nu=0.2_DP)

        ! Donem informacio del sistema
        call write_frame(unit_nr=traj_unit, pos=positions, stp_c=stp)
        pe_calc = calc_vdw_pbc(pos=positions, cutoff=CUTOFF_SET, boundary=box)
        ke_calc = calc_KE(velocities)
        instant_temperature = (2.0_DP / (3.0_DP * n_particles)) * ke_calc
        write(unit=log_unit, fmt='(I6,F25.8,F25.8,F25.8,F25.8)') stp, ke_calc, pe_calc, ke_calc + pe_calc, instant_temperature
    end do

    print '(A,F16.8,A,F16.8)', "Post-heating potential energy=", calc_vdw_pbc(pos=positions, cutoff=CUTOFF_SET, boundary=box), " kinetic energy=", calc_KE(velocities)

    call cpu_time(time1)

    print '(A,F12.8,A,F12.8)', "Execution time for heating:", time1 - time0, " time/iteration: ", (time1 - time0) / n_steps

    call write_frame(unit_nr=traj_unit, pos=positions, stp_c = stp)

    close(traj_unit)
    close(log_unit)

    ! Production trajectory
    open(newunit=traj_unit, file="thermodynamics.dat", access='sequential', action='write',&
    status='replace', form='formatted')

    n_steps = 500000

    call cpu_time(time0)

    print '(A,F16.8,A,F16.8)', "Initial traj potential energy=", calc_vdw_pbc(pos=positions, cutoff=CUTOFF_SET, boundary=box), " kinetic energy=", calc_KE(velocities)

    ! ~ Integrem les posicions dels atoms ~
    do stp = 1, n_steps

        ! Aplico integrador de velocity verlet
        call velocity_verlet(pos=positions, vel=velocities, dt=TIMESTEP, boundary=box)

        ! Aplico el termostat
        call andersen_thermostat(vel=velocities, temp=1.5_DP, nu=0.2_DP)

        ! Donem informacio del sistema
        if (mod(stp, write_file) == 0) then
            pe_calc = calc_vdw_pbc(pos=positions, cutoff=CUTOFF_SET, boundary=box)
            ke_calc = calc_KE(velocities)
            instant_temperature = (2.0_DP / (3.0_DP * n_particles)) * ke_calc
            write(unit=traj_unit, fmt='(I6,F25.8,F25.8,F25.8,F25.8)') stp, ke_calc, pe_calc, ke_calc + pe_calc, instant_temperature
        end if
    end do

    print '(A,F16.8,A,F16.8)', "Final traj potential energy=", calc_vdw_pbc(pos=positions, cutoff=CUTOFF_SET, boundary=box), " kinetic energy=", calc_KE(velocities)

    call cpu_time(time1)

    print '(A,F12.8,A,F12.8)', "Execution time for production:", time1 - time0, " time/iteration: ", (time1 - time0) / n_steps

    close(traj_unit)

    deallocate(positions)
    deallocate(velocities)

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
            write(unit=unit_nr, fmt='(A,F20.8,F20.8,F20.8)') "C ", pos(i_aux, 1), pos(i_aux, 2), &
        pos(i_aux, 3)
        end do

    end subroutine write_frame

    subroutine velocity_verlet(pos, vel, dt, boundary)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(inout), dimension(:,:) :: pos, vel
        real(kind=DP), intent(in)                    :: dt, boundary
        ! Internal_variables
        integer(kind=I64)                             :: n_p, idx_
        real(kind=DP), dimension(:, :), allocatable   :: actual_force, new_force

        n_p = size(pos, dim=1, kind=I64)

        allocate(actual_force(n_p, 3))
        allocate(new_force(n_p, 3))

        ! Calculem forces a r(t)
        call calc_vdw_force(pos=pos, cutoff=CUTOFF_SET, forces=actual_force, boundary=boundary)

        pos = pos + (vel*dt) + ((actual_force/(2.0_DP * MASS)) * (dt ** 2))

        ! Apliquem pbcs
        do idx_ = 1, n_p
            call pbc(pos(idx_, :), boundary)
        end do

        ! Calculem forces a r(t+dt)
        call calc_vdw_force(pos=pos, cutoff=CUTOFF_SET, forces=new_force, boundary=boundary)
        vel = vel + (((actual_force + new_force) / (2.0_DP * MASS)) * dt)

        deallocate(actual_force)
        deallocate(new_force)
    end subroutine velocity_verlet

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

    subroutine init_positions_sc(rho, pos)
        implicit none
        ! In/Out variables
        real(kind=DP), intent(inout), dimension(:, :) :: pos
        real(kind=DP), intent(in)                     :: rho
        ! Internal variables
        integer(kind=I64)                             :: M, i, j, k, p, N
        real(kind=DP)                                 :: L, a

        N = size(pos, dim=1, kind=I64)
        M = int(N ** (1.0_DP / 3.0_DP), kind=I64) + 1
        L = (real(N, kind=DP)/rho) ** (1.0_DP/3.0_DP)
        a = L/M

        print '(A)', ""
        print '(A,I3,A,I3,A,F12.8,A,F12.8)', "sc lattice, parameters: N=", N, " M=", M, " L=", L, " a=", a

        p = 1
        do i = 0, M - 1
            do j = 0, M - 1
                do k = 0, M - 1
                    ! p = (k+1) + (j*M) + (i * (M**2))
                    pos(p, :) = [real(i, kind=DP)*a + a*0.5_DP, real(j, kind=DP)*a + a*0.5_DP, real(k, kind=DP)*a + a*0.5_DP]
                    p = p + 1
                end do
            end do
        end do

    end subroutine

    subroutine andersen_thermostat(vel, temp, nu)
        implicit none
        ! In/Out variables
        real(kind=DP), dimension(:,:), intent(out) :: vel
        real(kind=DP), intent(in)                    :: temp, nu
        ! Internal variables
        integer(kind=I64) :: j_part, i_coord
        real(kind=DP)     :: rand_numb, sig

        sig = sqrt(temp)

        do j_part = 1, size(vel, dim=1)
            call random_number( harvest = rand_numb )
            if (rand_numb < nu) then
                do i_coord = 1, 3
                    vel(j_part, i_coord) = r8_normal_ab(a=0.0_DP, b=sig)
                end do 
            end if
        end do 
    end subroutine andersen_thermostat

    pure subroutine init_velocities(vel)
        implicit none
        real(kind=DP), dimension(:, :), intent(inout) :: vel

        vel = 0.0_DP
    end subroutine init_velocities

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

    function r8_normal_ab ( a, b )

        !*****************************************************************************80
        !
        !! R8_NORMAL_AB returns a scaled pseudonormal R8.
        !
        !  Discussion:
        !
        !    The normal probability distribution function (PDF) is sampled,
        !    with mean A and standard deviation B.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    06 August 2013
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real ( kind = rk ) A, the mean of the PDF.
        !
        !    Input, real ( kind = rk ) B, the standard deviation of the PDF.
        !
        !    Output, real ( kind = rk ) R8_NORMAL_AB, a sample of the normal PDF.
        !
          implicit none
        
          integer, parameter :: rk = kind ( 1.0D+00 )
        
          real ( kind = rk ) a
          real ( kind = rk ) b
          real ( kind = rk ) r1
          real ( kind = rk ) r2
          real ( kind = rk ) r8_normal_ab
          real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
          real ( kind = rk ) x
        
          call random_number ( harvest = r1 )
          call random_number ( harvest = r2 )
          x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
        
          r8_normal_ab = a + b * x
        
          return
        end function r8_normal_ab

end program main