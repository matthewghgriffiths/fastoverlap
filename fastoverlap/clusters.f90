!    FASTOVERLAP
!
!    FORTRAN Module for calculating Fast SO(3) Fourier transforms (SOFTs)
!    Copyright (C) 2017  Matthew Griffiths
!    
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!    
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    
!    You should have received a copy of the GNU General Public License along
!    with this program; if not, write to the Free Software Foundation, Inc.,
!    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


!    Includes code from https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
!
!    Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.

MODULE CLUSTERFASTOVERLAP

DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0

CONTAINS

SUBROUTINE CALCLSCALE(L, RJ, SIGMA, LSCALE)

IMPLICIT NONE
INTEGER, INTENT(IN) :: L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA
DOUBLE PRECISION, INTENT(OUT) :: LSCALE(0:L)

INTEGER J
DOUBLE PRECISION H, TMP, RJSIGMA, FACTOR, G1, G2, G3, GA, GB, GC, MAXH
LOGICAL APPROXH

MAXH = HUGE(H)
APPROXH = .FALSE.

RJSIGMA = RJ/SIGMA/SQRT(2.D0)
TMP = (2*PI*SIGMA**2)**1.5D0

DO J=0,L

CALL GAMMA((3.D0+J)*0.5D0, G1)
CALL GAMMA(1.5D0 + J, G2)

IF(.NOT.APPROXH) THEN
    CALL HYP1F1(0.5D0*J, 1.5D0+J, -RJSIGMA**2, H)
    H = H/G2
    IF ((H.NE.H).OR.(H.GE.MAXH)) THEN
        APPROXH = .TRUE.
        CALL HYP1F1(0.5D0*(J-1), 0.5D0+J, -RJSIGMA**2, H)
        CALL GAMMA((2.D0+J)*0.5D0, GA)
        CALL GAMMA(0.5D0 + J, GB)       
        CALL GAMMA(0.5D0 * (J-1), GC)
        H = H/GB
        FACTOR = H/((EXP(-RJSIGMA**2))*(-RJSIGMA)**(-(2+J))/GC + (RJSIGMA**(-J+1)/GA))
!WRITE(*,*) J, H, (EXP(-RJSIGMA**2))*(-RJSIGMA)**(-(2+J))/GC + (RJSIGMA**(-J+1)/GA)
!WRITE(*,*) GA, GB, GC, factor
    ENDIF
ENDIF

IF(APPROXH) THEN
    CALL GAMMA(0.5D0 * L, G3)
    H = FACTOR * ((EXP(-RJSIGMA**2))*(-RJSIGMA)**(-3-J)/G3 + (RJSIGMA**(-J))/G1)
ENDIF

!WRITE(*,*) J, G1, G2, H, (EXP(-RJSIGMA**2))*(-RJSIGMA)**(3+J)/G3 + (RJSIGMA**(-J))/G1

LSCALE(J) = H*G2
!LSCALE(J) = TMP * (RJSIGMA**J*H) * G1

ENDDO

CALL HYP1F1(0.5D0*L, 1.5D0+L, -RJSIGMA**2, H)

IF(H.NE.H) THEN
    
ENDIF

END SUBROUTINE CALCLSCALE

SUBROUTINE HARMSCALE(N, L, RJ, SIGMA, R0, DSCALE)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: DSCALE

DOUBLE PRECISION RJSIGMA, SIGMAR0, LSCALE, H, G1, G2, G3, TMP

RJSIGMA = RJ/SIGMA/SQRT(2.D0)
SIGMAR0 = (SIGMA/R0)**2
TMP = (2*PI*SIGMA**2)**1.5D0

CALL GAMMA((3.D0+L)*0.5, G1)
CALL GAMMA(1.5D0 + L, G2)

IF( EXPONENT(RJSIGMA)*2*L .LT. 20) THEN
    CALL HYP1F1(0.5D0*L, 1.5D0+L, -RJSIGMA**2, H)
    LSCALE = TMP * ((RJSIGMA)**L * H) * (G1 / G2)
ELSE
    CALL GAMMA(0.5D0 * L, G3)
    H = G1 / G3 * EXP(-RJSIGMA**2) * (-RJSIGMA)**(-(3+L)) + RJSIGMA**(-l)   
    LSCALE = TMP * (RJSIGMA)**L * H
ENDIF

DSCALE = EXP(-2*N*SIGMAR0) * LSCALE

END SUBROUTINE HARMSCALE

SUBROUTINE SOLVEBANDED(N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
INTEGER            IPIV( * )
DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )

EXTERNAL DGBSV

CALL DGBSV(N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

END SUBROUTINE SOLVEBANDED

SUBROUTINE BANDEDTODENSE(N, KL, KU, A, AB)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, KL, KU
DOUBLE PRECISION, INTENT(IN) :: AB(0:2*KL+KU, 0:N-1)
DOUBLE PRECISION, INTENT(OUT) :: A(0:N-1,0:N-1)

INTEGER I,J,K

A=0.D0
DO K=0,KU
    DO J=K,N-1
        ! AB(KL+KU+1+i-j,j) = A(i,j)
        A(J-K, J) = AB(KL+KU-K, J)
    ENDDO
ENDDO

DO K=1,KL
    DO J=0,N-K-1
        A(J+K, J) = AB(KL+KU+K, J)
    ENDDO
ENDDO

END SUBROUTINE BANDEDTODENSE

SUBROUTINE DENSETOBANDED(N, KL, KU, A, AB)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, KL, KU
DOUBLE PRECISION, INTENT(OUT) :: AB(0:2*KL+KU, 0:N-1)
DOUBLE PRECISION, INTENT(IN) :: A(0:N-1,0:N-1)

INTEGER I,J,K

AB = 0.D0
DO K=0,KU
    DO J=K,N-1
!WRITE(*,*) J-K, J, KL+KU-K+1, J, A(J-K, J)
        ! AB(KL+KU+1+i-j,j) = A(i,j)
        AB(KL+KU-K, J) = A(J-K, J)
    ENDDO
ENDDO

DO K=1,KL
    DO J=0,N-K-1
!WRITE(*,*) J+K, J, KL+KU+K+1, J, A(J+K, J)
        ! AB(KL+KU+1+i-j,j) = A(i,j)
        AB(KL+KU+K, J) = A(J+K, J)
    ENDDO
ENDDO

END SUBROUTINE DENSETOBANDED

SUBROUTINE RECURRSIVEARRAY(N, L, RJ, SIGMA, R0, AB)

! AB(KL+KU+i-j,j) = A(i,j) (0-indexing)
! Calculates matrix associated with following recurrence:
! 0 = \eta_{n,l} d_{n-1,l+2} + \mu_{n,l} d_{n,l} + \nu_{n,l} d_{n,l+1} + 
!      \xi_{n,l} d_{n,l+2} + \tau_{n,l} d_{n+1,l}

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: AB(0:4*(2*N+L-4),0:(N-2)*(L+N-2)-1)

DOUBLE PRECISION R0SIGMA,   SQRTI
INTEGER I,J,K,IND, KL, KU, NL, LDAB

KU = 0
KL = 2*(2*N+L-4)
NL = (N-2)*(L+N-2)
LDAB = 2*KL+KU+1
R0SIGMA = SIGMA**2/RJ/R0


write(*,*) SHAPE(AB)

! Set \tau
DO I=3,N
    DO J=0,L+2*N-2*I
        ! Put on 0 diagonal
        IND = (I-3)*(L+2*N-I-1) + J 
        AB(KL+KU, IND) = 1.D0 !SQRT(I+J+1.5D0) / (2*J+3)
    ENDDO
ENDDO

! Set \xi
DO I=4,N
    DO J=0,L+2*N-2*I
        IND = (I-4)*(L+2*N-I-1) + J + I - 2
        SQRTI = SQRT(1.D0*I)
        AB(KL+KU+L+2*N-2*I+1, IND) =  SQRT(I+J+1.5D0)/SQRTI
        AB(KL+KU+L+2*N-2*I+2, IND-1) = (2*J+3)*R0SIGMA/SQRTI
        AB(KL+KU+L+2*N-2*I+3, IND-2) = -SQRT(I+J+0.5D0)/SQRTI
    ENDDO
ENDDO

! Set \eta
DO I=5,N
    DO J=0,L+2*N-2*I
        IND = (I-5)*(L+2*N-I-1) + J + 2*I - 8
        AB(KL+KU+2*L+4*N-4*I+6, IND) = -SQRT((I-1.D0)/(I))
    ENDDO
ENDDO

END SUBROUTINE RECURRSIVEARRAY

SUBROUTINE RECURRSIVEARRAY3(N, L, RJ, SIGMA, R0, AB)

! AB(KL+KU+i-j,j) = A(i,j) (0-indexing)
! Calculates matrix associated with following recurrence:
! 0 = \eta_{n,l} d_{n-1,l+2} + \mu_{n,l} d_{n,l} + \nu_{n,l} d_{n,l+1} + 
!      \xi_{n,l} d_{n,l+2} + \tau_{n,l} d_{n+1,l}

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: AB(0:8*N-4, 0: N*N - N)

DOUBLE PRECISION R0SIGMA, SQRTI
INTEGER I,J,K,IND, KL, KU, NL, LDAB

KU = 0
KL = 4*N-2
NL = N*N - N
LDAB = 2*KL+KU+1

write(*,*) SHAPE(AB)

R0SIGMA = SQRT(2.D0) * R0 * RJ / (R0**2+SIGMA**2)

! Setting value for n,l = (0,0)
AB(KL+KU, 0) = 1.D0
! Setting values for n,l = 0,1<=l<=2n-5
IND = 1
DO J=1,2*N-5
    AB(KL+KU, IND) = 1.D0
    AB(KL+KU+1, IND-1) = - R0SIGMA / SQRT(1.D0 + 2.D0*J)
    IND = IND + 1
ENDDO

! Setting value for n,l = (1,0)
AB(KL+KU, IND) = 1.D0
! Setting value for n,l = (1,0)
AB(KL+KU+2*N-4, IND-(2*N-4)) = -(3.D0*R0**4 - 2*R0**2 * RJ**2 - 3.D0*SIGMA**2) / SQRT(6.D0) / (R0**2+SIGMA**2)**2
! Setting values for n,l = 1,1<=l<=2n-5
IND = IND + 1
DO J=1,2*N-5
    AB(KL+KU, IND) = 1.D0
    AB(KL+KU+1, IND-1) = - R0SIGMA* ((3.D0 + 2.D0*J)*R0**4 - 2*R0**2 * RJ**2 - (3.D0 + 2.D0*J)*SIGMA**2)/&
            ((1.D0 + 2.D0*J)*R0**4 - 2*R0**2 * RJ**2 - (1.D0 + 2.D0*J)*SIGMA**2) / SQRT(3.D0 + 2.D0*J)
    IND = IND + 1
ENDDO

R0SIGMA = SIGMA**2/RJ/R0
I=2
SQRTI = SQRT(1.D0*I)
DO J=0,2*N-2*I-3
    ! Put on 0 diagonal, n,l -> i,j
    AB(KL+KU, IND) = 1.D0! I + 0.1*J
    ! n,l -> i-1,j+2
    AB(KL+KU + 2*N-2*I-2, IND - (2*N-2*I-2)) = SQRT(I+J+1.5D0)/SQRTI
    ! n,l -> i-1,j+1
    AB(KL+KU + 2*N-2*I-1, IND - (2*N-2*I-1)) = (2*J+3)*R0SIGMA/SQRTI
    ! n,l -> i-1,j
    AB(KL+KU + 2*N-2*I  , IND - (2*N-2*I  )) = -SQRT(I+J+0.5D0)/SQRTI
    ! n,l -> i-2,j+2
    AB(KL+KU + 4*N-4*I-2, IND - (4*N-4*I-2)) = -SQRT(I-1.D0)/SQRTI
    IND = IND + 1
ENDDO

DO I=3,N-2
    SQRTI = SQRT(1.D0*I)
    DO J=0,2*N-2*I-3
        ! Put on 0 diagonal, n,l -> i,j
        AB(KL+KU, IND) = 1.D0
        ! n,l -> i-1,j+2
        AB(KL+KU + 2*N-2*I-2, IND - (2*N-2*I-2)) = SQRT(I+J+1.5D0)/SQRTI
        ! n,l -> i-1,j+1
        AB(KL+KU + 2*N-2*I-1, IND - (2*N-2*I-1)) = (2*J+3)*R0SIGMA/SQRTI
        ! n,l -> i-1,j
        AB(KL+KU + 2*N-2*I  , IND - (2*N-2*I  )) = -SQRT(I+J+0.5D0)/SQRTI
        ! n,l -> i-2,j+2
        AB(KL+KU + 4*N-4*I  , IND - (4*N-4*I  )) = -SQRT(I-1.D0)/SQRTI
        IND = IND + 1
    ENDDO
ENDDO

! (n,l) = (n-1,0)
I = N-1
SQRTI = SQRT(1.D0*I)
J=0
AB(KL+KU, IND) = 1.D0 
! n,l -> i-1,j+1
AB(KL+KU+1, IND-1) = (2*J+3)*R0SIGMA/SQRTI
! n,l -> i-1,j
AB(KL+KU+2, IND-2) = -SQRT(I+J+0.5D0)/SQRTI
! n,l -> i-2,j+2
AB(KL+KU+4, IND-4) = -SQRT(I-1.D0)/SQRTI
IND = IND + 1

J=1
! (n,l) = (n-1,1)
AB(KL+KU, IND) = 1.D0
! n,l -> i-1,j
AB(KL+KU+2, IND-2) = -SQRT(I+J+0.5D0)/SQRTI
! n,l -> i-2,j+2
AB(KL+KU+5, IND-5) = -SQRT(I-1.D0)/SQRTI
IND = IND + 1

I = N
SQRTI = SQRT(1.D0*I)
J=0
! (n,l) = (n,0)
AB(KL+KU, IND) = 1.D0
AB(KL+KU+1, IND-1) = (2*J+3)*R0SIGMA/SQRTI
AB(KL+KU+2, IND-2) = -SQRT(I+J+0.5D0)/SQRTI

END SUBROUTINE RECURRSIVEARRAY3

SUBROUTINE RECURRSIVEARRAY4(N, L, RJ, SIGMA, R0, AB)

! AB(KL+KU+i-j,j) = A(i,j) (0-indexing)
! Calculates matrix associated with following recurrence:
! 0 = \eta_{n,l} d_{n-1,l+2} + \mu_{n,l} d_{n,l} + \nu_{n,l} d_{n,l+1} + 
!      \xi_{n,l} d_{n,l+2} + \tau_{n,l} d_{n+1,l}

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: AB(0:8*N-4, 0: N*N - N)

DOUBLE PRECISION R0SIGMA, SQRTI, EXPSR
INTEGER I,J,K,IND, KL, KU, NL, LDAB

KU = 0
KL = 4*N-2
NL = N*N - N
LDAB = 2*KL+KU+1

EXPSR = EXP(2.D0 * SIGMA**2 / R0**2)
R0SIGMA = R0 * RJ / (R0**2+SIGMA**2)

! Setting value for n,l = (0,0)
AB(KL+KU, 0) = 1.D0
! Setting values for n,l = 0,1<=l<=2n-5
IND = 1
DO J=1,2*N-5
    AB(KL+KU, IND) = 1.D0
    AB(KL+KU+1, IND-1) = - 1.D0
    IND = IND + 1
ENDDO

! Setting value for n,l = (1,0)
AB(KL+KU, IND) = 1.D0
! Setting value for n,l = (1,0)
AB(KL+KU+2*N-4, IND-(2*N-4)) = -0.5D0 * (3.D0*R0**4 - 2*R0**2 * RJ**2 - 3.D0*SIGMA**2) * EXPSR / (R0**2+SIGMA**2)**2
! Setting values for n,l = 1,1<=l<=2n-5
IND = IND + 1
DO J=1,2*N-5
    AB(KL+KU, IND) = 1.D0
    AB(KL+KU+1, IND-1) = - ((3.D0 + 2.D0*J)*R0**4 - 2*R0**2 * RJ**2 - (3.D0 + 2.D0*J)*SIGMA**2)/&
        ((1.D0 + 2.D0*J)*R0**4 - 2*R0**2 * RJ**2 - (1.D0 + 2.D0*J)*SIGMA**2) !*R0SIGMA**2/(1.5D0 + J)
    IND = IND + 1
ENDDO

R0SIGMA = R0*RJ/(R0**2+SIGMA**2)
I=2
SQRTI = SQRT(1.D0*I)
DO J=0,2*N-2*I-3
    ! Put on 0 diagonal, n,l -> i,j
    AB(KL+KU, IND) = 1.D0! I + 0.1*J
    ! n,l -> i-1,j+2
    AB(KL+KU + 2*N-2*I-2, IND - (2*N-2*I-2)) = R0SIGMA**2 * EXPSR / I
    ! n,l -> i-1,j+1
    AB(KL+KU + 2*N-2*I-1, IND - (2*N-2*I-1)) = SIGMA**2 * (2.D0*J+3.D0) /(R0**2+SIGMA**2) /I * EXPSR
    ! n,l -> i-1,j
    AB(KL+KU + 2*N-2*I  , IND - (2*N-2*I  )) = - (I+J+0.5D0)/I * EXPSR
    ! n,l -> i-2,j+2
    AB(KL+KU + 4*N-4*I-2, IND - (4*N-4*I-2)) = - R0SIGMA**2 / I * EXPSR**2
    IND = IND + 1
ENDDO

DO I=3,N-2
    SQRTI = SQRT(1.D0*I)
    DO J=0,2*N-2*I-3
        ! Put on 0 diagonal, n,l -> i,j
        AB(KL+KU, IND) = 1.D0! I + 0.1*J
        ! n,l -> i-1,j+2
        AB(KL+KU + 2*N-2*I-2, IND - (2*N-2*I-2)) = R0SIGMA**2 * EXPSR / I
        ! n,l -> i-1,j+1
        AB(KL+KU + 2*N-2*I-1, IND - (2*N-2*I-1)) = SIGMA**2 * (2.D0*J+3.D0) /(R0**2+SIGMA**2) /I * EXPSR
        ! n,l -> i-1,j
        AB(KL+KU + 2*N-2*I  , IND - (2*N-2*I  )) = - (I+J+0.5D0)/I * EXPSR
        ! n,l -> i-2,j+2
        AB(KL+KU + 4*N-4*I, IND - (4*N-4*I)) = - R0SIGMA**2 / I * EXPSR**2
        IND = IND + 1
    ENDDO
ENDDO

! (n,l) = (n-1,0)
I = N-1
SQRTI = SQRT(1.D0*I)
J=0
AB(KL+KU, IND) = 1.D0 
! n,l -> i-1,j+1
AB(KL+KU+1, IND-1) = 2*SIGMA**2 * SQRT((2.D0*J+3.D0)/I) /(R0**2+SIGMA**2)
! n,l -> i-1,j
AB(KL+KU+2, IND-2) = -SQRT(I+J+0.5D0)/SQRTI
! n,l -> i-2,j+2
AB(KL+KU+4, IND-4) = -SQRT((I-1.D0)/(I*(3.D0+2.D0*L)*(5+2.D0*L))) * 2 *R0SIGMA**2
IND = IND + 1

J=1
! (n,l) = (n-1,1)
AB(KL+KU, IND) = 1.D0
! n,l -> i-1,j
AB(KL+KU+2, IND-2) = -SQRT(I+J+0.5D0)/SQRTI
! n,l -> i-2,j+2
AB(KL+KU+5, IND-5) = -SQRT((I-1.D0)/(I*(3.D0+2.D0*L)*(5+2.D0*L))) * 2 *R0SIGMA**2
IND = IND + 1

I = N
SQRTI = SQRT(1.D0*I)
J=0
! (n,l) = (n,0)
AB(KL+KU, IND) = 1.D0
AB(KL+KU+1, IND-1) = -SQRT(I+J+0.5D0)/SQRTI
AB(KL+KU+2, IND-2) = -SQRT(I+J+0.5D0)/SQRTI

END SUBROUTINE RECURRSIVEARRAY4

SUBROUTINE RECURRSIVEARRAY2(N, L, RJ, SIGMA, R0, AB)

! AB(KL+KU+i-j,j) = A(i,j) (0-indexing)
! Calculates matrix associated with following recurrence:
! 0 = \eta_{n,l} d_{n-1,l+2} + \mu_{n,l} d_{n,l} + \nu_{n,l} d_{n,l+1} + 
!      \xi_{n,l} d_{n,l+2} + \tau_{n,l} d_{n+1,l}

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: AB(0:3*L+6*N-19,0:(N-2)*(L+N-2)-1)

DOUBLE PRECISION R0SIGMA
INTEGER I,J,K,IND, KL, KU, NL, LDAB

KU = L+1 + 2*N - 6
KL = L-1 + 2*N - 6
NL = (N-2)*(L+N-2)
LDAB = 2*KL+KU+1


! Set \nu_{n,l}
R0SIGMA = SIGMA**2/RJ/R0
DO I=3,N
    DO J=0,L+2*N-2*I-1
        IND = (I-3)*(L+2*N-I-1) + J + 1
        ! Put on +1 diagonal
        AB(KL+KU-1, IND) = I+0.1*J! (2*J+3) * R0SIGMA / SQRT(I+J+1.5D0) !I+0.1*J! 
    ENDDO
ENDDO

! Set \xi
DO I=3,N
    DO J=0,L+2*N-2*I-2
        IND = (I-3)*(L+2*N-I-1)  + J + 2
        ! Put on +2 diagonal
        AB(KL+KU-2, IND) = I+0.1*J! SQRT(I+J+2.5D0) / SQRT(I+J+1.5D0)
    ENDDO
ENDDO

write(*,*) SHAPE(AB)
! Set \mu
DO I=3,N
    DO J=0,L+2*N-2*I
        ! Put on 0 diagonal
        IND = (I-3)*(L+2*N-I-1) + J 
        AB(KL+KU, IND) = I+0.1*J! - 1.D0 !SQRT(I+J+1.5D0) / (2*J+3)
    ENDDO
ENDDO

! Set \tau
DO I=3,N-1
    DO J=0,L+2*N-2*I-2
        IND = (I-2)*(L+2*N-I-1) + J + 2 -I
        AB(KL+KU-L-1-2*N+2*I, IND) = I+0.1*J ! SQRT(I+1.D0) / SQRT(I+J+1.5D0)
    ENDDO
ENDDO

! Set \eta
DO I=4,N
    DO J=0,L+2*N-2*I
        IND = (I-4)*(L+2*N-I-1) + J + I - 2
        write(*,*) I, J, IND, -L-1-2*N+2*I
        AB(KL+KU+L+2*N-2*I+1, IND) = I+0.1*J! -SQRT(I*1.D0) / SQRT(I+J+1.5D0)
    ENDDO
ENDDO

END SUBROUTINE RECURRSIVEARRAY2

SUBROUTINE BANDHARMONICNL(N, L, RJ, SIGMA, R0, RET)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: RET(0:N,0:L)

DOUBLE PRECISION R0SIGMA, AB(0:3*L,0:(N-2)*(L-1)-1), B(0:N*L-1)
INTEGER I,J,K,IND, KL, KU, NL, LDAB, IPIV(N*(L+1)), INFO

EXTERNAL DGBSV

KU = L
KL = L
NL = N*(L+1)
LDAB = 2*KL+KU+1


! Calculate values of integral when n=0 and
! Populate B with appropriate values
B = 0.D0
R0SIGMA = 1.D0/(R0**2+SIGMA**2)
RET(0,0) = SQRT(2.D0*SQRT(PI)*(R0*R0SIGMA)**3) * SIGMA**3 * EXP(-0.5D0*RJ**2*R0SIGMA)*4*PI
R0SIGMA = SQRT(2.D0) * R0 * RJ * R0SIGMA
DO J=1,L
    RET(0,J) = R0SIGMA / SQRT(1.D0+2.D0*J) * RET(0,J-1)
    B(J-1) = - RET(0,J)
ENDDO

! Calculate values of integral when l=0 and
! Populate B with appropriate values
RET(1,0) = SQRT(2.D0/3.D0) * RET(0,0)
B(0) = B(0) - SQRT(2.5D0)*RET(1,0)
DO I=2,N
    RET(I,0) = SQRT(2.D0*I/(1+2.D0*I)) * RET(I-1,0)
    B((I-1)*L) = - SQRT(I+1.5D0) * RET(I,0)
ENDDO


CALL RECURRSIVEARRAY(N, L, RJ, SIGMA, R0, AB)
CALL DGBSV(NL, KL, KU, 1, AB, LDAB, IPIV, B, NL, INFO )

!write(*,*) info, shape(ret), shape(b)
!
!DO I=1,N
!    DO J=1,L
!        IND = (I-1)*L + J - 1
!
!IF(IND.GT.N*L-1) THEN
!        write(*,*) I, J, IND
!ENDIF
!        RET(I,J) = B(IND)
!    ENDDO
!ENDDO

END SUBROUTINE BANDHARMONICNL

SUBROUTINE HARMONIC0L(N, RJ, SIGMA, R0, RET)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: RET(0:N)

DOUBLE PRECISION R0SIGMA
INTEGER I,J,K

R0SIGMA = 1.D0/(R0**2+SIGMA**2)
RET(0) = SQRT(2.D0*SQRT(PI)*(R0*R0SIGMA)**3) * SIGMA**3 * EXP(-0.5D0*RJ**2*R0SIGMA)*4*PI

R0SIGMA = SQRT(2.D0) * R0 * RJ * R0SIGMA
DO I=1,N
    RET(I) = R0SIGMA / SQRT(1.D0+2.D0*I) * RET(I-1)
ENDDO

END SUBROUTINE HARMONIC0L

SUBROUTINE HARMONICNL(N,L,RJ,SIGMA,R0,RET)

!
! Calculates the value of the overlap integral up to N and L
!
! 4\pi \int_0^{\infty} g_{nl}(r)\exp{\left(-\frac{r^2+{r^p_j}^2}{2\sigma^2}\right)} 
! i_l \left( \frac{r r^p_{j}}{\sigma^2} \right) r^2\; \mathrm{d}r
!
! N is the maximum quantum number of the Harmonic basis to calculate up to
! L is the maximum angular moment number to calculate
! SIGMA is the width of the Gaussian Kernels
! R0 is the length scale of the Harmonic Basis
! RET is the matrix of calculate values of the overlap integral
!

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
DOUBLE PRECISION, INTENT(IN) :: RJ, SIGMA, R0
DOUBLE PRECISION, INTENT(OUT) :: RET(0:N,0:L)

DOUBLE PRECISION R0SIGMA, RET2, SQRTI
INTEGER I,J,K

! Initiate Recurrence
R0SIGMA = 1.D0/(R0**2+SIGMA**2)
RET(0,0) = SQRT(2.D0*SQRT(PI)*(R0*R0SIGMA)**3) * SIGMA**3 * EXP(-0.5D0*RJ**2*R0SIGMA)*4*PI
R0SIGMA = SQRT(2.D0) * R0 * RJ * R0SIGMA
DO J=1,L
    RET(0,J) = R0SIGMA / SQRT(1.D0+2.D0*J) * RET(0,J-1)
ENDDO

R0SIGMA = SIGMA**2/RJ/R0
! When I=1 don't calculate RET(I-2,J)
I = 1
SQRTI = 1.D0 !SQRT(REAL(I,8))
DO J=0,L-2
!write(*,*) J, SQRT(I+J+0.5D0), (2.D0*J+3.D0)*SIGMA**2/RJ/R0, SQRT(I+J+1.5D0)
!write(*,*) J, RET(I-1,J), RET(I-1,J+1), RET(I-1,J+2)
!write(*,*) J, SQRT(I+J+0.5D0)*RET(I-1,J), (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1), SQRT(I+J+1.5D0) * RET(I-1,J+2)
RET(I,J) = (SQRT(I+J+0.5D0)*RET(I-1,J) - (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1) -&
    SQRT(I+J+1.5D0) * RET(I-1,J+2))/SQRTI
!write(*,*) J, RET(I,J)
ENDDO
! Assuming that integral for J>L = 0
!RET(I,L-1) = (SQRT(I+J+0.5D0)*RET(I-1,J) - (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1))/SQRTI
!RET(I,L) = (SQRT(I+J+0.5D0)*RET(I-1,J))/SQRTI


DO I=2,N
SQRTI = SQRT(REAL(I,8))
DO J=0,L-2*I
RET(I,J) = (SQRT(I+J+0.5D0)*RET(I-1,J) - (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1) -&
    SQRT(I+J+1.5D0) * RET(I-1,J+2) + SQRT(I-1.D0) * RET(I-2,J+2))/SQRTI
ENDDO
! Assuming that integral for J>L = 0
!RET(I,L-1) = (SQRT(I+J+0.5D0)*RET(I-1,J) - (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1))/SQRTI
!RET(I,L) = (SQRT(I+J+0.5D0)*RET(I-1,J))/SQRTI
ENDDO

END SUBROUTINE HARMONICNL

SUBROUTINE RYLM(COORD, R, YLM, L)

! Calculates the Spherical Harmonics associated with coordinate COORD
! up to L, returns R, the distance COORD is from origin
! Calculates value of Legendre Polynomial Recursively

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: COORD(3)
INTEGER, INTENT(IN) :: L
DOUBLE PRECISION, INTENT(OUT) :: R
COMPLEX*16, INTENT(OUT) :: YLM(0:L,0:2*L)

INTEGER J, M, INDM1, INDM0, INDM2
DOUBLE PRECISION THETA, PHI, Z, FACTORIALS(0:2*L)
COMPLEX*16 EXPIM(0:2*L)
DOUBLE PRECISION, EXTERNAL :: RLEGENDREL0, RLEGENDREM0, RLEGENDREM1


R = (COORD(1)**2+COORD(2)**2+COORD(3)**2)**0.5
PHI = ATAN2(COORD(2), COORD(1))
Z = COORD(3)/R

!Calculating Associate Legendre Function
YLM = CMPLX(0.D0,0.D0, 8)
YLM(0,0) = (4*PI)**-0.5

! Initialising Recurrence for Associated Legendre Polynomials
DO J=0, L-1
    YLM(J+1,J+1) = RLEGENDREL0(J, Z) * YLM(J,J) !* ((2.D0*J+3,D0)/(2.D0*J+1.D0)/(2.D0*J+1.D0)/(2.D0*J+2.D0))
ENDDO

! Recurrence for Associated Legendre Polynomials
DO J=1,L
    DO M=J,-J+1,-1
        INDM1 = MODULO(M+1, 2*L+1)
        INDM2 = MODULO(M-1, 2*L+1)
        INDM0 = MODULO(M, 2*L+1)
        YLM(J,INDM2) = RLEGENDREM0(M,J,Z) * YLM(J,INDM0) + RLEGENDREM1(M,J,Z) * YLM(J,INDM1)
    ENDDO
ENDDO

! Calculating exp(imPHI) component
DO M=-L,L
    INDM0 = MODULO(M, 2*L+1)
    EXPIM(INDM0) = EXP(CMPLX(0.D0, M*PHI, 8))
ENDDO

! Calculating Factorials
FACTORIALS(0) = 1.D0
DO M=1,2*L
    FACTORIALS(M) = M*FACTORIALS(M-1)
ENDDO

! Calculate Spherical Harmonics
DO J=1,L
    DO M=-J,J
        INDM0 = MODULO(M, 2*L+1)
        ! Could probably calculate the prefactor through another recurrence relation...
        YLM(J,INDM0) = EXPIM(INDM0)*YLM(J,INDM0) * ((2.D0*J+1.D0)*FACTORIALS(J-M)/FACTORIALS(J+M))**0.5
    ENDDO
ENDDO

END SUBROUTINE RYLM

SUBROUTINE HARMONICCOEFFS(COORDS, NATOMS, CNLM, N, L, HWIDTH, KWIDTH)

!
! For a set of Gaussian Kernels of width KWIDTH at COORDS, 
! this will calculate the coefficients of the isotropic quantum harmonic basis
! cnlm with length scale HWIDTH up to N and L.
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, N, L
DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), HWIDTH, KWIDTH
COMPLEX*16, INTENT(OUT) :: CNLM(0:N,0:L,0:2*L)

COMPLEX*16 :: YLM(0:L,0:2*L)
DOUBLE PRECISION HARMCOEFFS(0:2*N+L,0:N,0:L), DNL(0:N,0:L+2*N), RJ
INTEGER I,J,K,SI,M,INDM, S

CNLM = CMPLX(0.D0,0.D0,8)

DO K=1,NATOMS
    CALL RYLM(COORDS(3*K-2:3*K), RJ, YLM, L)
    CALL HARMONICNL(N,L+2*N,RJ,KWIDTH,HWIDTH,DNL)
    DO J=0,L
        DO M=-J,J
            INDM = MODULO(M,2*L+1)
            DO I=0,N
            CNLM(I,J,INDM) = CNLM(I,J,INDM) + DNL(I,J) * CONJG(YLM(J,INDM))
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE HARMONICCOEFFS

SUBROUTINE DOTHARMONICCOEFFS(C1NLM, C2NLM, N, L, ILMM)

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, L
COMPLEX*16, INTENT(IN) :: C1NLM(0:N,0:L,0:2*L), C2NLM(0:N,0:L,0:2*L)
COMPLEX*16, INTENT(OUT) :: ILMM(0:L,0:2*L,0:2*L)

INTEGER I, J, M1, M2, INDM1, INDM2

ILMM = CMPLX(0.D0,0.D0,8)

DO J=0,L
    DO M1=-J,J
        INDM1 = MODULO(M1, 2*L+1)
        DO M2=-J,J
            INDM2 = MODULO(M2, 2*L+1)
            DO I=0,N
                ILMM(J,INDM1,INDM2) = ILMM(J,INDM1,INDM2) + CONJG(C1NLM(I,J,INDM1))*C2NLM(I,J,INDM2)
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE DOTHARMONICCOEFFS

SUBROUTINE FOURIERCOEFFS(COORDSB, COORDSA, NATOMS, L, KWIDTH, ILMM, YLMA, YLMB)
!
! Calculates S03 Coefficients of the overlap integral of two structures
! does this calculation by direct calculation of the overlap between every pair
! of atoms, slower than the Harmonic basis, but slightly more accurate.
!

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), KWIDTH
INTEGER, INTENT(IN) :: NATOMS, L
COMPLEX*16, INTENT(OUT) :: ILMM(0:L,0:2*L,0:2*L)

COMPLEX*16, INTENT(OUT) :: YLMA(0:L,0:2*L,NATOMS), YLMB(0:L,0:2*L,NATOMS)
DOUBLE PRECISION RA(NATOMS), RB(NATOMS), IL(0:L), R1R2, EXPRA(NATOMS), EXPRB(NATOMS), FACT, TMP

INTEGER IA,IB,I,J,K,M1,M2,INDM1,INDM2

! Precalculate some values
DO I=1,NATOMS
CALL RYLM(COORDSA(3*I-2:3*I), RA(I), YLMA(:,:,I), L)
CALL RYLM(COORDSB(3*I-2:3*I), RB(I), YLMB(:,:,I), L)
EXPRA(I) = EXP(-0.25D0 * RA(I)**2 / KWIDTH**2)
EXPRB(I) = EXP(-0.25D0 * RB(I)**2 / KWIDTH**2)
ENDDO


FACT = 4.D0 * PI**2.5 * KWIDTH**3

ILMM = CMPLX(0.D0,0.D0,8)
DO IA=1,NATOMS
    DO IB=1,NATOMS
        R1R2 = 0.5D0 * RA(IA)*RB(IB)/KWIDTH**2
        CALL SPHI(L, R1R2, K, IL)
        TMP = FACT*EXPRA(IA)*EXPRB(IB)
        DO J=0,L
            DO M2=-L,L
                INDM2 = MODULO(M2, 2*L+1)
                DO M1=-L,L
                    INDM1 = MODULO(M1, 2*L+1)
                    ILMM(J,INDM1,INDM2) = ILMM(J,INDM1,INDM2) + IL(J)*YLMA(J,INDM1,IA)*CONJG(YLMA(J,INDM2,IB))*TMP
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE FOURIERCOEFFS

END MODULE CLUSTERFASTOVERLAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION RLEGENDREL0(L, Z)

! Calcualates recurrence factor M1 for associated legendre polynomials@
! P^{L+1}_{L+1} (Z) = L0*P^L_L (Z)

IMPLICIT NONE
INTEGER, INTENT(IN) :: L
DOUBLE PRECISION, INTENT(IN) :: Z

RLEGENDREL0 = - (2.D0*L+1) * (1-Z**2)**0.5 

END FUNCTION RLEGENDREL0


DOUBLE PRECISION FUNCTION RLEGENDREM0(M, L, Z)
! Calcualates recurrence factor M1 for associated legendre polynomials@
! P^{M-1}_L (Z) = M0*P^M_L (Z) + M1*P^{M+1}_L (Z)

IMPLICIT NONE
INTEGER, INTENT(IN) :: M, L
DOUBLE PRECISION, INTENT(IN) :: Z

RLEGENDREM0 = - 2.D0 * M * Z / (1.D0-Z**2)**0.5 / (L+M) / (L-M+1.D0)

END FUNCTION RLEGENDREM0

DOUBLE PRECISION FUNCTION RLEGENDREM1(M, L, Z)
! Calcualates recurrence factor M1 for associated legendre polynomials@
! P^{M-1}_L (Z) = M0*P^M_L (Z) + M1*P^{M+1}_L (Z)

IMPLICIT NONE
INTEGER, INTENT(IN) :: M, L
DOUBLE PRECISION, INTENT(IN) :: Z

RLEGENDREM1 = - 1.D0 / (L+M) / (L-M+1.D0)

END FUNCTION RLEGENDREM1

function envj ( n, x )

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Discussion:
!
!    ENVJ estimates -log(Jn(x)) from the estimate
!    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modifications suggested by Vincent Lafage, 11 January 2016.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bessel function.
!
!    Input, real ( kind = 8 ) X, the absolute value of the argument.
!
!    Output, real ( kind = 8 ) ENVJ, the value.
!
  implicit none

  real ( kind = 8 ) envj
  real ( kind = 8 ) logten
  integer ( kind = 4 ) n
  real ( kind = 8 ) n_r8
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x
!
!  Original code
!
  if ( .true. ) then

    envj = 0.5D+00 * log10 ( 6.28D+00 * n ) &
      - n * log10 ( 1.36D+00 * x / n )
!
!  Modification suggested by Vincent Lafage.
!
  else

    n_r8 = real ( n, kind = 8 )
    logten = log ( 10.0D+00 )
    envj = r8_gamma_log ( n_r8 + 1.0D+00 ) / logten - n_r8 * log10 ( x )

  end if

  return
end



function msta1 ( x, mp )

!*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward  
!    recurrence such that the magnitude of    
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) MP, the negative logarithm of the 
!    desired magnitude.
!
!    Output, integer ( kind = 4 ) MSTA1, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) x

  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
    nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )                  
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta1 = nn

  return
end function msta1

function msta2 ( x, n, mp )

!*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!    Jianming Jin supplied a modification to this code on 12 January 2016.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Jn(x).
!
!    Input, integer ( kind = 4 ) N, the order of Jn(x).
!
!    Input, integer ( kind = 4 ) MP, the number of significant digits.
!
!    Output, integer ( kind = 4 ) MSTA2, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) ejn
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) hmp
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) obj
  real ( kind = 8 ) x

  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1D+00 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1D+00 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
    nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta2 = nn + 10

  return
end function msta2

subroutine sphi ( n, x, nm, si)

!*****************************************************************************80
!
!! SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    18 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of In(X).
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) SI(0:N), DI(0:N), the values and derivatives
!    of the function of orders 0 through N.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: n

  real ( kind = 8 ) cs
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ), intent(out) :: nm
  real ( kind = 8 ), intent(out) :: si(0:n)
  real ( kind = 8 ) si0
  real ( kind = 8 ), intent(in) :: x

  nm = n

  if ( abs ( x ) < 1.0D-100 ) then
    do k = 0, n
      si(k) = 0.0D+00
    end do
    si(0) = 1.0D+00
    return
  end if

  si(0) = sinh ( x ) / x
  si(1) = -( sinh ( x ) / x - cosh ( x ) ) / x
  si0 = si(0)

  if ( 2 <= n ) then

    m = msta1 ( x, 200 )
    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if
    f0 = 0.0D+00
    f1 = 1.0D+00-100
    do k = m, 0, -1
      f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x + f0
      if ( k <= nm ) then
        si(k) = f
      end if
      f0 = f1
      f1 = f
    end do
    cs = si0 / f
    do k = 0, nm
      si(k) = cs * si(k)
    end do

  end if

  return
end subroutine sphi


subroutine HYP1F1 ( ain, bin, xin, hg )

!*****************************************************************************80
!
!! CHGM computes the confluent hypergeometric function M(a,b,x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    27 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, parameters.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) HG, the value of M(a,b,x).
!
  implicit none

  real ( kind = 8 ), intent(in) :: ain
  real ( kind = 8 ), intent(in) :: bin
  real ( kind = 8 ), intent(in) :: xin
  real ( kind = 8 ), intent(out) :: hg

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) x

  real ( kind = 8 ) a0
  real ( kind = 8 ) a1
  real ( kind = 8 ) aa

  real ( kind = 8 ) hg1
  real ( kind = 8 ) hg2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nl
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) rg
  real ( kind = 8 ) sum1
  real ( kind = 8 ) sum2
  real ( kind = 8 ) ta
  real ( kind = 8 ) tb
  real ( kind = 8 ) tba
  real ( kind = 8 ) x0
  real ( kind = 8 ) xg
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1

  a=ain
  b=bin
  x=xin
  pi = 3.141592653589793D+00
  a0 = a
  a1 = a
  x0 = x
  hg = 0.0D+00

  y1 = hg

  if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
    hg = 1.0D+300
  else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
    hg = 1.0D+00
  else if ( a == -1.0D+00 ) then
    hg = 1.0D+00 - x / b
  else if ( a == b ) then
    hg = exp ( x )
  else if ( a - b == 1.0D+00 ) then
    hg = ( 1.0D+00 + x / b ) * exp ( x )
  else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
    hg = ( exp ( x ) - 1.0D+00 ) / x
  else if ( a == int ( a ) .and. a < 0.0D+00 ) then
    m = int ( - a )
    r = 1.0D+00
    hg = 1.0D+00
    do k = 1, m
      r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
      hg = hg + r
    end do
  end if

  if ( hg /= 0.0D+00 ) then
    return
  end if

  if ( x < 0.0D+00 ) then
    a = b - a
    a0 = a
    x = abs ( x )
  end if

  if ( a < 2.0D+00 ) then
    nl = 0
  end if

  if ( 2.0D+00 <= a ) then
    nl = 1
    la = int ( a )
    a = a - la - 1.0D+00
  end if

  do n = 0, nl

    if ( 2.0D+00 <= a0 ) then
      a = a + 1.0D+00
    end if

    if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then

      hg = 1.0D+00
      rg = 1.0D+00
      do j = 1, 500
        rg = rg * ( a + j - 1.0D+00 ) &
          / ( j * ( b + j - 1.0D+00 ) ) * x
        hg = hg + rg
        if ( abs ( rg / hg ) < 1.0D-15 ) then
          exit
        end if
      end do

    else

      call gamma ( a, ta )
      call gamma ( b, tb )
      xg = b - a
      call gamma ( xg, tba )
      sum1 = 1.0D+00
      sum2 = 1.0D+00
      r1 = 1.0D+00
      r2 = 1.0D+00
      do i = 1, 8
        r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
        r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
        sum1 = sum1 + r1
        sum2 = sum2 + r2
      end do
      hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
      hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
      hg = hg1 + hg2

    end if

    if ( n == 0 ) then
      y0 = hg
    else if ( n == 1 ) then
      y1 = hg
    end if

  end do

  if ( 2.0D+00 <= a0 ) then
    do i = 1, la - 1
      hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
      y0 = y1
      y1 = hg
      a = a + 1.0D+00
    end do
  end if

  if ( x0 < 0.0D+00 ) then
    hg = hg * exp ( x0 )
  end if

  a = a1
  x = x0

  return
end

subroutine gamma ( x, ga )

!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
  implicit none

  real ( kind = 8 ), intent(in) :: x
  real ( kind = 8 ), intent(out) :: ga

  real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
    1.0D+00, &
    0.5772156649015329D+00, &
   -0.6558780715202538D+00, &
   -0.420026350340952D-01, &
    0.1665386113822915D+00, &
   -0.421977345555443D-01, &
   -0.96219715278770D-02, &
    0.72189432466630D-02, &
   -0.11651675918591D-02, &
   -0.2152416741149D-03, &
    0.1280502823882D-03, & 
   -0.201348547807D-04, &
   -0.12504934821D-05, &
    0.11330272320D-05, &
   -0.2056338417D-06, & 
    0.61160950D-08, &
    0.50020075D-08, &
   -0.11812746D-08, &
    0.1043427D-09, & 
    0.77823D-11, &
   -0.36968D-11, &
    0.51D-12, &
   -0.206D-13, &
   -0.54D-14, &
    0.14D-14, &
    0.1D-15 /)

  real ( kind = 8 ) gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) z

  if ( x == aint ( x ) ) then

    if ( 0.0D+00 < x ) then
      ga = 1.0D+00
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0D+00 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0D+00
      do k = 1, m
        r = r * ( z - real ( k, kind = 8 ) )
      end do
      z = z - real ( m, kind = 8 )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0D+00 / ( gr * z )

    if ( 1.0D+00 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0D+00 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return
end