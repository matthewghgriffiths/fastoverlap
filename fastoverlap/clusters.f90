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

! Including functions used to calculate orthogonal polynomials
INCLUDE "polynomials.f90"

! Module for performing discrete SO(3) transforms, depends on fftw.
INCLUDE "DSOFT.f90"

INCLUDE "utils.f90"

MODULE CLUSTERFASTOVERLAP

DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0

CONTAINS

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

SUBROUTINE RYML(COORD, R, YML, L)

! Calculates the Spherical Harmonics associated with coordinate COORD
! up to L, returns R, the distance COORD is from origin
! Calculates value of Legendre Polynomial Recursively

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: COORD(3)
INTEGER, INTENT(IN) :: L
DOUBLE PRECISION, INTENT(OUT) :: R
COMPLEX*16, INTENT(OUT) :: YML(-L:L,0:L)

INTEGER J, M, INDM1, INDM0, INDM2
DOUBLE PRECISION THETA, PHI, Z, FACTORIALS(0:2*L)
COMPLEX*16 EXPIM(-L:L)
DOUBLE PRECISION, EXTERNAL :: RLEGENDREL0, RLEGENDREM0, RLEGENDREM1


R = (COORD(1)**2+COORD(2)**2+COORD(3)**2)**0.5
PHI = ATAN2(COORD(2), COORD(1))
Z = COORD(3)/R

!Calculating Associate Legendre Function
YML = CMPLX(0.D0,0.D0, 8)
YML(0,0) = (4*PI)**-0.5

! Initialising Recurrence for Associated Legendre Polynomials
DO J=0, L-1
    YML(J+1,J+1) = RLEGENDREL0(J, Z) * YML(J,J) !* ((2.D0*J+3,D0)/(2.D0*J+1.D0)/(2.D0*J+1.D0)/(2.D0*J+2.D0))
ENDDO

! Recurrence for Associated Legendre Polynomials
DO J=1,L
    DO M=J,-J+1,-1
        YML(M-1, J) = RLEGENDREM0(M,J,Z) * YML(M, J) + RLEGENDREM1(M,J,Z) * YML(M+1,J)
    ENDDO
ENDDO

! Calculating exp(imPHI) component
DO M=-L,L
    EXPIM(M) = EXP(CMPLX(0.D0, M*PHI, 8))
ENDDO

! Calculating Factorials
FACTORIALS(0) = 1.D0
DO J=1,2*L
    FACTORIALS(J) = J*FACTORIALS(J-1)
ENDDO

! Calculate Spherical Harmonics
DO J=1,L
    DO M=-J,J
        INDM0 = MODULO(M, 2*L+1)
        ! Could probably calculate the prefactor through another recurrence relation...
        YML(M,J) = EXPIM(M)*YML(M,J) * ((2.D0*J+1.D0)*FACTORIALS(J-M)/FACTORIALS(J+M))**0.5
    ENDDO
ENDDO

END SUBROUTINE RYML

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

SUBROUTINE HARMONICCOEFFSNML(COORDS, NATOMS, CNML, N, L, HWIDTH, KWIDTH)

!
! For a set of Gaussian Kernels of width KWIDTH at COORDS, 
! this will calculate the coefficients of the isotropic quantum harmonic basis
! cnlm with length scale HWIDTH up to N and L.
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, N, L
DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), HWIDTH, KWIDTH
COMPLEX*16, INTENT(OUT) :: CNML(0:N,-L:L,0:L)

COMPLEX*16 :: YML(-L:L,0:L)
DOUBLE PRECISION HARMCOEFFS(0:2*N+L,0:N,0:L), DNL(0:N,0:L+2*N), RJ
INTEGER I,J,K,SI,M,INDM, S

CNML = CMPLX(0.D0,0.D0,8)

DO K=1,NATOMS
    CALL RYML(COORDS(3*K-2:3*K), RJ, YML, L)
    CALL HARMONICNL(N,L+2*N,RJ,KWIDTH,HWIDTH,DNL)
    DO J=0,L
        DO M=-J,J
            INDM = MODULO(M,2*L+1)
            DO I=0,N
                CNML(I,M,J) = CNML(I,M,J) + DNL(I,J) * CONJG(YML(M,J))
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE HARMONICCOEFFSNML

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

SUBROUTINE DOTHARMONICCOEFFSNML(C1NML, C2NML, N, L, IMML)

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, L
COMPLEX*16, INTENT(IN) :: C1NML(0:N,-L:L,0:L), C2NML(0:N,-L:L,0:L)
COMPLEX*16, INTENT(OUT) :: IMML(-L:L,-L:L,0:L)

INTEGER I, J, M1, M2, INDM1, INDM2

IMML = CMPLX(0.D0,0.D0,8)

DO J=0,L
    DO M2=-J,J
        DO M1=-J,J
            DO I=0,N
                IMML(M1,M2,J) = IMML(M1,M2,J) + CONJG(C1NML(I,M1,J))*C2NML(I,M2,J)
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE DOTHARMONICCOEFFSNML

SUBROUTINE FOURIERCOEFFS(COORDSB, COORDSA, NATOMS, L, KWIDTH, ILMM)
!
! Calculates S03 Coefficients of the overlap integral of two structures
! does this calculation by direct calculation of the overlap between every pair
! of atoms, slower than the Harmonic basis, but slightly more accurate.
!

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), KWIDTH
INTEGER, INTENT(IN) :: NATOMS, L
COMPLEX*16, INTENT(OUT) :: ILMM(0:L,0:2*L,0:2*L)

COMPLEX*16 YLMA(0:L,0:2*L,NATOMS), YLMB(0:L,0:2*L,NATOMS)
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

SUBROUTINE FOURIERCOEFFSMML(COORDSB, COORDSA, NATOMS, L, KWIDTH, IMML)
!
! Calculates S03 Coefficients of the overlap integral of two structures
! does this calculation by direct calculation of the overlap between every pair
! of atoms, slower than the Harmonic basis, but slightly more accurate.
!

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), KWIDTH
INTEGER, INTENT(IN) :: NATOMS, L
COMPLEX*16, INTENT(OUT) :: IMML(-L:L,-L:L,0:L)

COMPLEX*16 YMLA(-L:L,0:L,NATOMS), YMLB(-L:L,0:L,NATOMS)
DOUBLE PRECISION RA(NATOMS), RB(NATOMS), IL(0:L), R1R2, EXPRA(NATOMS), EXPRB(NATOMS), FACT, TMP

INTEGER IA,IB,I,J,K,M1,M2,INDM1,INDM2

! Precalculate some values
DO I=1,NATOMS
    CALL RYML(COORDSA(3*I-2:3*I), RA(I), YMLA(:,:,I), L)
    CALL RYML(COORDSB(3*I-2:3*I), RB(I), YMLB(:,:,I), L)
    EXPRA(I) = EXP(-0.25D0 * RA(I)**2 / KWIDTH**2)
    EXPRB(I) = EXP(-0.25D0 * RB(I)**2 / KWIDTH**2)
ENDDO


FACT = 4.D0 * PI**2.5 * KWIDTH**3

IMML = CMPLX(0.D0,0.D0,8)
DO IA=1,NATOMS
    DO IB=1,NATOMS
        R1R2 = 0.5D0 * RA(IA)*RB(IB)/KWIDTH**2
        CALL SPHI(L, R1R2, K, IL)
        TMP = FACT*EXPRA(IA)*EXPRB(IB)
        DO J=0,L
            DO M2=-L,L
                DO M1=-L,L
                    IMML(M1,M2,J) = IMML(M1,M2,J) + IL(J)*YMLA(M1,J,IA)*CONJG(YMLB(M2,J,IB))*TMP
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE FOURIERCOEFFSMML

SUBROUTINE FOURIERCOEFFSMEM(COORDSB, COORDSA, NATOMS, L, KWIDTH, IMML, IL, YMLB, YMLA)
!
! Calculates S03 Coefficients of the overlap integral of two structures
! does this calculation by direct calculation of the overlap between every pair
! of atoms, slower than the Harmonic basis, but slightly more accurate.
!

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), KWIDTH
INTEGER, INTENT(IN) :: NATOMS, L
COMPLEX*16, INTENT(OUT) :: IMML(-L:L,-L:L,0:L)

COMPLEX*16, INTENT(OUT) ::  YMLA(-L:L,0:L,NATOMS), YMLB(-L:L,0:L,NATOMS)
DOUBLE PRECISION RA(NATOMS), RB(NATOMS), R1R2, EXPRA(NATOMS), EXPRB(NATOMS), FACT, TMP

DOUBLE PRECISION, INTENT(OUT) ::  IL(0:L, NATOMS, NATOMS)

INTEGER IA,IB,I,J,K,M1,M2,INDM1,INDM2

YMLA = CMPLX(0.D0,0.D0,8)
YMLB = CMPLX(0.D0,0.D0,8)
! Precalculate some values
DO I=1,NATOMS
    CALL RYML(COORDSA(3*I-2:3*I), RA(I), YMLA(:,:,I), L)
    CALL RYML(COORDSB(3*I-2:3*I), RB(I), YMLB(:,:,I), L)
    EXPRA(I) = EXP(-0.25D0 * RA(I)**2 / KWIDTH**2)
    EXPRB(I) = EXP(-0.25D0 * RB(I)**2 / KWIDTH**2)
ENDDO

DO IA=1,NATOMS
    DO IB=1,NATOMS
        R1R2 = 0.5D0 * RA(IA)*RB(IB)/KWIDTH**2
        CALL SPHI(L, R1R2, K, IL(:,IA,IB))
    ENDDO
ENDDO


FACT = 4.D0 * PI**2.5 * KWIDTH**3

IMML = CMPLX(0.D0,0.D0,8)
DO IA=1,NATOMS
    DO IB=1,NATOMS
!        R1R2 = 0.5D0 * RA(IA)*RB(IB)/KWIDTH**2
!        CALL SPHI(L, R1R2, K, IL)
        TMP = FACT*EXPRA(IA)*EXPRB(IB)
        DO J=0,L
            DO M2=-L,L
                DO M1=-L,L
                    IMML(M1,M2,J) = IMML(M1,M2,J) + IL(J,IA,IB)*YMLA(M1,J,IA)*CONJG(YMLB(M2,J,IB))*TMP
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE FOURIERCOEFFSMEM

SUBROUTINE CALCOVERLAP(IMML, OVERLAP, L)
USE DSOFT, ONLY : ISOFT

IMPLICIT NONE

INTEGER, INTENT(IN) :: L
COMPLEX*16, INTENT(IN) :: IMML(-L:L,-L:L,0:L)
COMPLEX*16, INTENT(OUT) :: OVERLAP(2*L+2,2*L+2,2*L+2)

END SUBROUTINE CALCOVERLAP

END MODULE CLUSTERFASTOVERLAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

