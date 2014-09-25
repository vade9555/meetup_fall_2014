PROGRAM jacobi

  IMPLICIT NONE

  INTEGER, PARAMETER :: ndummy = 1

  DOUBLE PRECISION, ALLOCATABLE :: phi(:,:,:,:)
  DOUBLE PRECISION :: east, west, north, south, top, bottom

  INTEGER :: imax, jmax, kmax
  INTEGER :: bi, bj, bk         ! blocking factors
  INTEGER :: niter

  DOUBLE PRECISION :: diff_max

  ! local variables
  INTEGER :: i_sol
  INTEGER :: ios, plength
  CHARACTER(len=255) :: cmdline ! command line argument
  CHARACTER(len=255) :: message

  NAMELIST /input/ imax, jmax, kmax, bi, bj, bk, east, west, north, south, top, bottom, niter

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

  CALL solve_phi(imax, jmax, kmax, ndummy, bi, bj, bk, niter, phi(:,:,:,:), i_sol, diff_max)


CONTAINS


  SUBROUTINE solve_phi(im, jm, km, nd, bi, bj, bk, niter, p, i_sol, diff_max)
    INTEGER, INTENT(in) :: im, jm, km, nd, bi, bj, bk, niter
    DOUBLE PRECISION, INTENT(inout) :: p(1-nd:im+nd, 1-nd:jm+nd, 1-nd:km+nd, 2)
    INTEGER, INTENT(out) :: i_sol
    DOUBLE PRECISION, INTENT(out) :: diff_max

    DOUBLE PRECISION, PARAMETER :: flops = 6.0D0
    DOUBLE PRECISION :: problem_size
    DOUBLE PRECISION :: walltime, tot_walltime, max_walltime, min_walltime, avg_walltime
    DOUBLE PRECISION :: max_mflops, avg_mflops
    DOUBLE PRECISION :: oos = 1.d0/6.d0
    INTEGER :: clock_start, clock_rate, clock_max, clock_end
    INTEGER :: i_old, i_new, i_tmp
    INTEGER :: nt
    INTEGER :: i, j, k
    INTEGER :: is, js, ks, ie, je, ke
    i_old = 1
    i_new = 2
    
    walltime = 0.0D0
    tot_walltime = 0.0D0
    max_walltime = 0.0D0
    min_walltime = 1.0D10
    avg_walltime = 0.0D0
    problem_size = DBLE(im)*DBLE(jm)*DBLE(km)
    DO nt = 1, niter
       walltime = 0.0
       diff_max = 0.0
       CALL system_CLOCK(clock_start, clock_rate, clock_max)
       DO ks = 1, km, bk
          ke = MIN(km, ks+bk-1)
          DO js = 1, jm, bj
             je = MIN(jm, js+bj-1)
             DO is = 1, im, bi
                ie = MIN(im, is+bi-1)
                ! 6 flops, 6 loads one store
                !loop on a block
                DO k=ks, ke
                   DO j=js, je
                      !DIR$ SIMD
                      DO i=is, ie
                         p(i, j, k, i_new) = oos * ( &
                              & p(i-1, j, k, i_old) + p(i+1, j, k, i_old) + &
                              & p(i, j-1, k, i_old) + p(i, j+1, k, i_old) + &
                              & p(i, j, k-1, i_old) + p(i, j, k+1, i_old) )
                         diff_max = MAX(ABS(p(i,j,k,i_new) - p(i,j,k, i_old)), diff_max)
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
       CALL system_CLOCK(clock_end, clock_rate, clock_max)
       walltime = walltime + REAL(clock_end - clock_start) / REAL(clock_rate)
       tot_walltime = tot_walltime + walltime
       max_walltime = MAX(walltime, max_walltime)
       min_walltime = MIN(walltime, min_walltime)
       WRITE(*,*) "Iteration: ", nt, "max difference: ", diff_max
       i_tmp = i_old
       i_old = i_new
       i_new = i_tmp
    END DO
    WRITE(*,*) "Total wall clock time of solver (excluding I/O to stdout): ", tot_walltime
    WRITE(*,*) "Average time: ", tot_walltime / DBLE(niter)
    WRITE(*,*) "Min time: ", min_walltime
    WRITE(*,*) "Max time: ", max_walltime
    WRITE(*,*) "Wall clock time / point: ", tot_walltime / DBLE(niter) / problem_size

    avg_mflops = DBLE(niter)*problem_size*flops/(tot_walltime*1.0D06)
    max_mflops = problem_size*flops/(min_walltime*1.0D06)
    WRITE(*,*) "Average MFLOP/s: ", avg_mflops, "Maximum MFLOP/s: ", max_mflops
    
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
