PROGRAM testspline
  USE fastspline

  IMPLICIT NONE

  INTEGER, PARAMETER :: nn = 12, nx = 2*12 + 2 - 1

  INTEGER :: n
  REAL*8,  DIMENSION(nn) :: y 
  REAL*8 :: lam, T
  INTEGER :: JJ, r

  REAL*8 :: score
  REAL*8, DIMENSION(nx) :: x

! Input
  y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
  y = [-1, 2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2]
  n = nn
  r = 2
  T = 0.1
  JJ = 6
  lam = 5.8


  print *, 'Call cspline'
  CALL cspline(x, score, y, n, r, T, JJ, lam)

! Output
  WRITE(*,*)x

END PROGRAM testspline

