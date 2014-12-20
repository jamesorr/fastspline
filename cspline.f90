MODULE mcspline
CONTAINS

SUBROUTINE cspline(x, score, y, n, r, T, JJ, lam)

  IMPLICIT NONE 

! Input arguments
  INTEGER, INTENT(in) :: n
  REAL*8,  INTENT(in), DIMENSION(N) :: y 
!f2py optional , depend(sal) :: n=len(y)
  !> lambda parameter
  REAL*8,  INTENT(in) :: lam
  REAL*8,  INTENT(in) :: T 
  INTEGER, INTENT(in) :: JJ
  INTEGER, INTENT(in) :: r  
! Output arguments
  REAL*8, INTENT(out) :: score
! REAL*8, INTENT(out), ALLOCATABLE :: x(:)
! INTEGER, PARAMETER :: nx = r*n+r-1
! REAL*8, INTENT(out) :: x(nx)
  REAL*8, INTENT(out) :: x(r*n+r-1)
! REAL*8, INTENT(out) :: x(25)


! Local variables
  REAL*8, DIMENSION(:), ALLOCATABLE :: c, e, f, w, z
  REAL*8 :: a0, a1, d

! REAL (kind=8)    :: a2, x1, x2, alpha, beta
  COMPLEX (kind=8) :: a2, x1, x2, alpha, beta
  COMPLEX (kind=8) :: atmp

  REAL*8 :: fac1, fac2
  REAL*8 :: elim, flim, glim, hlim, qlim
  REAL*8 :: g1, g2, h
  REAL*8 :: j1, j2, j3, j4
  REAL*8 :: lamc, lamr, mu, q, sq
  REAL*8 :: Tcu
  REAL*8 :: tmp1, tmp2, tmp3, tmp4
  REAL*8 :: tr, tr1, tr2, tr3
  REAL*8 :: v
  INTEGER :: ir
  INTEGER :: nc   
  INTEGER :: r2, rn, rsq
  INTEGER :: i, j, k
  INTEGER :: NN

! print *, 'Initialized variables'

! n = SIZE(y)
  nc = CEILING(REAL(n)/REAL(2.))
  Tcu = T**3
  rn = r*n

! IF (ALLOCATED(x)) DEALLOCATE(x)
! ALLOCATE(x(rn+r-1))
  x = 0.0d0 

! print *, '2) Allocate c, w, z'

! IF (ALLOCATED(c)) DEALLOCATE(c)
  ALLOCATE(c(n))
  c = 0.0d0 

! IF (ALLOCATED(w)) DEALLOCATE(w)
  ALLOCATE(w(n-2))
  w = 0.0d0 

! IF (ALLOCATED(z)) DEALLOCATE(z)
  ALLOCATE(z(n))
  z = 0.0d0 

  DO j=1,n-2 
     w(j) = y(j)-2*y(j+1)+y(j+2)
  END DO

! print *, '3) determine lamr, a0, a1, a2, x1, x2'
  lamr = lam*Tcu
  a0 = 6.d0 + lamr*2.d0/3.d0
  a1 = 4.d0 - lamr/6.d0
  atmp = a1**2 - 4.d0*(a0-2.d0)
  a2 = SQRT(atmp)
  x1 =  (a1 + a2) / 2.d0
  x2 = -(a2 - a1) / 2.d0
! print *,' a1**2 - 4.d0*(a0-2.d0) = ', atmp
! print *, 'a0, a1, a2 =', a0, a1, a2
! print *, 'x1, x2 =', x1, x2

! print *, '4) determine alpha, beta'
  IF (lamr > 24.d0) THEN 
     alpha = 0.5d0*(x1 + SQRT(x1**2-4.d0))
     beta  = 0.5d0*(x2 + SQRT(x2**2-4.d0))
  ELSE IF (lamr < 24.d0) THEN 
     alpha = 0.5d0*(x1 - SQRT(x1**2-4.d0))
     beta  = 0.5d0*(x2 - SQRT(x2**2-4))
  ELSE
     alpha = 0.5d0*(x1 - SQRT(x1**2-4.d0))
     beta  = 0.5d0*(x2 + SQRT(x2**2-4.d0))
  END IF
! print *, 'alpha, beta = ', alpha, beta

! IF (JJ > LOG10(alpha*beta) - (nc-1)*2*LOG10(ABS(alpha))) THEN 
  IF (JJ > LOG10(DBLE(alpha*beta)) - (nc-1)*2*LOG10(ABS(alpha))) THEN 
     !print *, '5) JJ if block: GREATER THAN'

     !Use untruncated algorithm since NN > nc-2 
     !Factor coefficient matrix, solve triangular systems, find trace 

!    IF (ALLOCATED(e)) DEALLOCATE(e)
     ALLOCATE(e(n-2))
     e = 0.0d0 

!    IF (ALLOCATED(f)) DEALLOCATE(f)
     ALLOCATE(f(n-2))
     f = 0.0d0 

     d = a0
     f(1) = 1.d0 / d
     c(2) = f(1)*w(1)
     mu = a1
     e(1) = mu*f(1)
     d = a0-mu*e(1)
     f(2) = 1.d0 / d
     c(3) = f(2)*(w(2)+mu*c(2))
     mu = a1 - e(1)
     e(2) = mu*f(2)

     DO j=3,n-2 
        d = a0-mu*e(j-1)-f(j-2)
        f(j) = 1.d0 / d
        c(j+1) = f(j)*(w(j)+mu*c(j)-c(j-1))
        mu = a1 - e(j-1)
        e(j) = mu * f(j)
     END DO
     c(n-2) = c(n-2)+e(n-3)*c(n-1)

     DO j=n-4,1,-1 
        c(j+1) = c(j+1) + e(j)*c(j+2) - f(j)*c(j+3)
     END DO

     g2 = f(n-2)
     tr1 = g2
     h = e(n-3)*g2
     tr2 = h
     g1 = f(n-3) + e(n-3)*h
     tr1 = tr1 + g1
     tr3=0

     DO k=n-4,n-nc,-1 
        q = e(k)*h - f(k)*g2
        tr3 = tr3 + q
        h = e(k)*g1 - f(k)*h
        tr2 = tr2 + h
        g2 = g1
        g1 = f(k)*(1-q) + e(k)*h
        tr1 = tr1 + g1
     END DO

     q = e(n-nc-1)*h - f(n-nc-1)*g2
     tr3 = tr3 + q
     h = e(n-nc-1)*g1 - f(n-nc-1)*h
     tr2 = tr2 + h

     tr1 =  6.d0*(2.d0*tr1 - DBLE(MOD(n,2))*g1)
     tr2 = -8.d0*(2.d0*tr2 - (1 + DBLE(MOD(n,2)) )*h)
     tr3 =  2.d0*(2.d0*tr3 - DBLE(MOD(n,2))*q)
     tr = (tr1+tr2+tr3)/DBLE(n)

  ELSE 
     !print *, '5) JJ if block: LESS than'
     !Use truncated algorithm since NN < nc-1
     !Factor coefficient matrix, solve triangular systems, find trace 
     flim = alpha*beta
     elim = alpha + beta
     glim = flim*(1+flim) / ((1 - flim)*((1.d0+flim)**2 - elim**2))
     hlim = elim*glim/(1.d0+flim)
     qlim = elim*hlim - flim*glim
     NN = CEILING((LOG10(flim) - JJ) / (2*LOG10(ABS(alpha))))

     !print *, '6) Allocate e, f'
     IF (ALLOCATED(e)) DEALLOCATE(e)
     ALLOCATE(e(NN))
     e = 0.0d0 

     IF (ALLOCATED(f)) DEALLOCATE(f)
     ALLOCATE(f(NN))
     f = 0.0d0 

     !print *, '7) produce tr%, g%, etc.'
     d = a0
     f(1) = 1.d0 / d
     c(2) = f(1)*w(1)
     mu = a1
     e(1) = mu*f(1)
     d = a0 - mu*e(1)
     f(2) = 1.d0 / d
     c(3) = f(2)*(w(2) + mu*c(2))
     mu = a1 - e(1)
     e(2) = mu*f(2)
     g2 = flim
     tr1 = g2
     h = elim*g2
     tr2 = h
     g1 = flim+elim*h
     tr1 = tr1+g1
     tr3 = 0.d0

     !print *, '8) Loop to produce tr1 etc'
     DO j=3,NN 
        d = a0 - mu*e(j-1) - f(j-2)
        f(j) = 1.d0 / d
        c(j+1) = f(j)*(w(j) + mu*c(j) - c(j-1))
        mu = a1 - e(j-1)
        e(j) = mu*f(j)
        q = elim*h - flim*g2
        tr3 = tr3 + q
        h = elim*g1 - flim*h
        tr2 = tr2+h
        g2 = g1
        g1 = flim*(1-q) + elim*h
        tr1 = tr1 + g1
     END DO

     !print *, '9) Compute tr1, tr2, tr3, tr, mu'
     tr1 = tr1+(nc-NN-1)*glim
     tr2 = tr2+(nc-NN)*hlim
     tr3 = tr3+(nc-NN)*qlim
     tr1 =  6.d0*(2.d0*tr1 - DBLE(MOD(n,2))*glim)
     tr2 = -8.d0*(2.d0*tr2 - (1.d0 + DBLE(MOD(n,2)) )*hlim)
     tr3 =  2.d0*(2.d0*tr3 - DBLE(MOD(n,2))*qlim)
     tr = (tr1 + tr2 + tr3)/DBLE(n)
     mu = a1 - elim

     !print *, '10) Compute c(j+1): first'
     !print *, 'c =', c
     !print *, 'NN, n =', NN, n
     !print *, 'flim, mu =', flim, mu
     DO j=NN+1,n-2 
        c(j+1) = flim*(w(j) + mu*c(j)-c(j-1))
     END DO
     print *, c
     print *, '10) Compute c(n-2)'
     c(n-2) = c(n-2) + elim*c(n-1)

     print *, '11) Compute c(j)'
     DO j=n-3,NN+2,-1 
        c(j) = c(j) + elim*c(j+1) - flim*c(j+2)
     END DO

     print *, '12) Compute c(j+1) : second'
     DO j=NN,1,-1 
        c(j+1) = c(j+1) + e(j)*c(j+2) - f(j)*c(j+3)
     END DO

  END IF

  print *, '13) Compute score'
  !Compute GCV score 
  z(1) = c(2)
  z(2) = c(3)-2*c(2)
  DO j=3,n-2 
     z(j) = c(j-1)-2*c(j)+c(j+1)
  END DO

  z(n-1) = c(n-2)-2*c(n-1)
  z(n) = c(n-1)

  !m2f: sq=(z*z')/n
! sq = (z*TRANSPOSE(z))/n
! sq = MatMul(z,TRANSPOSE(z)) / DBLE(n)
! Compute the dot product of the vector z
  sq = sum (z * z)
  score = sq / tr**2

  print *, 'Funny x(r:rn:r) operation'
  !Compute estimates 
  !x(r:r:rn)=y-z
  x(r:rn:r) = y-z
  IF (r < 8) THEN 
     fac1 = x(2*r) - x(r)    - lamr*c(2)/6.d0
     fac2 = x(rn)  - x(rn-r) + lamr*c(n-1)/6.d0
     DO j=1,r-1 
        j1 = DBLE(j)/DBLE(r)
        j2 = 1.d0 - j1
        v = lamr*j1*j2 / 6.d0
        j3 = v*(1.d0 + j1)
        j4 = v*(2.d0 - j1)
        DO i=1,n-1 
           ir = i*r
           x(ir+j) = j2*x(ir) + j1*x(ir+r) - j3*c(i+1) - j4*c(i)
        END DO
        x(j)    = x(r)  - j2*fac1
        x(rn+j) = x(rn) + j1*fac2
     END DO
  ELSE
     lamc = lamr/(6*r**3)
     r2 = 2*r
     rsq = r**2
     DO i=1,n-1 
        ir = i*r
        tmp1 = x(ir)   / r
        tmp2 = x(ir+r) / r
        tmp3 = lamc*c(i+1)
        tmp4 = lamc*c(i)
        DO j=1,r-1 
           x(ir+j) = DBLE(r-j)*tmp1 + DBLE(j)*tmp2 - DBLE(j)*DBLE(rsq-j*j)*tmp3 &
                   - DBLE(j)*DBLE(r-j)*DBLE(r2-j)*tmp4
        END DO
     END DO
     tmp1 = x(r) - x(r+1)
     tmp2 = x(rn) - x(rn-1)
     DO j = 1,r-1 
        x(j)    = x(r) + (r-j)*tmp1
        x(rn+j) = x(rn) + j*tmp2
     END DO
  END IF
END SUBROUTINE cspline
END MODULE mcspline
