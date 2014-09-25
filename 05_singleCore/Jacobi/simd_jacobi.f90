PROGRAM jacobi

  IMPLICIT NONE

  INTEGER, PARAMETER :: ndummy = 1

  DOUBLE PRECISION, ALLOCATABLE :: phi(:,:,:,:)
  DOUBLE PRECISION :: east, west, north, south, top, bottom

  INTEGER :: imax, jmax, kmax
  INTEGER :: niter

  DOUBLE PRECISION :: diff_max

  ! local variables
  INTEGER :: i, j, k
  INTEGER :: i_sol
  INTEGER :: ios, plength
  CHARACTER(len=255) :: cmdline ! command line argument
  CHARACTER(len=255) :: message

  NAMELIST /input/ imax, jmax, kmax, east, west, north, south, top, bottom, niter

  ! read the command line parameters
  CALL get_COMMAND(cmdline)
  CALL get_command_ARGUMENT(0, LENGTH=plength)
  ! remove the first command line argument which is the program name
  ! create a "namelist" string
  cmdline = "&INPUT "//cmdline(plength+2:)// " /"
  WRITE(*,*) cmdline
  ! read the variables from the namelist string
  READ(cmdline, NML=input, iostat=ios, iomsg=message)
  IF (ios /= 0) THEN
     WRITE(*,*) 'ERROR: ', ios
     WRITE(*,*) message
     STOP
  END IF

  WRITE(*, NML=input)

  ! allocate and initialize
  ALLOCATE(phi(1-ndummy:imax+ndummy, 1-ndummy:jmax+ndummy, 1-ndummy:kmax+ndummy, 2))

  phi(:,:,:,:) = 0.0D0

  CALL set_boundaries(imax, jmax, kmax, ndummy, east, west, north, south, top, bottom, phi(:,:,:,1))
  CALL set_boundaries(imax, jmax, kmax, ndummy, east, west, north, south, top, bottom, phi(:,:,:,2))

  ! now run the solver for niter iterations

  CALL solve_phi(imax, jmax, kmax, ndummy, niter, phi(:,:,:,:), i_sol, diff_max)


CONTAINS

  SUBROUTINE solve_phi(im, jm, km, nd, niter, p, i_sol, diff_max)
    INTEGER, INTENT(in) :: im, jm, km, nd, niter
    DOUBLE PRECISION, INTENT(inout) :: p(1-nd:im+nd, 1-nd:jm+nd, 1-nd:km+nd, 2)
    INTEGER, INTENT(out) :: i_sol
    DOUBLE PRECISION, INTENT(out) :: diff_max

    double precision, parameter :: flops = 6.0D0
    double precision :: problem_size
    double precision :: walltime, tot_walltime, max_walltime, min_walltime, avg_walltime
    double precision :: max_mflops, avg_mflops
    DOUBLE PRECISION :: oos = 1.d0/6.d0
    integer :: clock_start, clock_rate, clock_max, clock_end
    INTEGER :: i_old, i_new, i_tmp
    INTEGER :: nt

    i_old = 1
    i_new = 2
    
    walltime = 0.0D0
    tot_walltime = 0.0D0
    max_walltime = 0.0D0
    min_walltime = 1.0D10
    avg_walltime = 0.0D0
    problem_size = dble(im)*dble(jm)*dble(km)
    DO nt = 1, niter
       walltime = 0.0
       diff_max = 0.0
       call system_clock(clock_start, clock_rate, clock_max)
       DO k = 1, km
          DO j = 1, jm
             !DIR$ SIMD
             DO i = 1, im
                ! 6 flops, 6 loads one store
                p(i, j, k, i_new) = oos * ( &
                     & p(i-1, j, k, i_old) + p(i+1, j, k, i_old) + &
                     & p(i, j-1, k, i_old) + p(i, j+1, k, i_old) + &
                     & p(i, j, k-1, i_old) + p(i, j, k+1, i_old) )
                diff_max = MAX(ABS(p(i,j,k,i_new) - p(i,j,k,i_old)),diff_max)
             END DO
          END DO
       END DO
       call system_clock(clock_end, clock_rate, clock_max)
       walltime = walltime + real(clock_end - clock_start) / real(clock_rate)
       tot_walltime = tot_walltime + walltime
       max_walltime = max(walltime, max_walltime)
       min_walltime = min(walltime, min_walltime)
       WRITE(*,*) "Iteration: ", nt, "max difference: ", diff_max
       i_tmp = i_old
       i_old = i_new
       i_new = i_tmp
    END DO
    write(*,*) "Total wall clock time of solver (excluding I/O to stdout): ", tot_walltime
    write(*,*) "Average time: ", tot_walltime / dble(niter)
    write(*,*) "Min time: ", min_walltime
    write(*,*) "Max time: ", max_walltime
    write(*,*) "Wall clock time / point: ", tot_walltime / dble(niter) / problem_size

    avg_mflops = dble(niter)*problem_size*flops/(tot_walltime*1.0D06)
    max_mflops = problem_size*flops/(min_walltime*1.0D06)
    write(*,*) "Average MFLOP/s: ", avg_mflops, "Maximum MFLOP/s: ", max_mflops
    
    i_sol = i_old
       
  END SUBROUTINE solve_phi

  SUBROUTINE set_boundaries(im, jm, km, nd, e, w, n, s, t, b, p)
    INTEGER, INTENT(in) :: im, jm, km, nd
    DOUBLE PRECISION, INTENT(in) :: e, w, n, s, t, b
    DOUBLE PRECISION, INTENT(inout) :: p(1-nd:im+nd, 1-nd:jm+nd, 1-nd:km+nd)

    ! east
    p(im+1:im+nd, :, :) = e
    ! west
    p(1-nd:1, :, :) = w
    ! north
    p(:, jm+1:jm+nd, :) = n
    ! south
    p(:, 1-nd:1, :) = s
    ! top
    p(:, :, km+1:km+nd) = t
    ! bottom
    p(:, :, 1-nd:1) = b
  END SUBROUTINE set_boundaries


END PROGRAM jacobi
