module thermostats_m
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use :: math_utils_m, only: r8_normal_ab
    implicit none

    public :: andersen_thermostat

contains

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

end module thermostats_m