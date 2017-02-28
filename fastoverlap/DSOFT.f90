!    SOFT
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

!    This code is a FORTRAN reimplementation of the SOFT C++ library from
!    http://www.cs.dartmouth.edu/~geelong/soft/ under the GNU GPL licence

!    Citation:
!    Kostelec, P. J., & Rockmore, D. N. (2008). FFTs on the rotation group. 
!    Journal of Fourier Analysis and Applications, 14(2), 145–179. 
!    http://doi.org/10.1007/s00041-008-9013-5

!    Dependencies:
!        1. FFTW

MODULE DSOFT

INTEGER*8, SAVE :: BW
DOUBLE PRECISION, SAVE, ALLOCATABLE :: WEIGHTS(:), WIGNERD(:,:,:,:)

CONTAINS

SUBROUTINE SETBANDWIDTH(BANDWIDTH)

IMPLICIT NONE
INTEGER*8, INTENT(IN) :: BANDWIDTH

!WRITE(*,*) BW, BANDWIDTH
! Check if bandwidth has already been calculated
IF (BW.NE.BANDWIDTH) THEN
    IF (ALLOCATED(WEIGHTS)) THEN
        DEALLOCATE(WEIGHTS)
    ENDIF
    IF (ALLOCATED(WIGNERD)) THEN
        DEALLOCATE(WIGNERD)
    ENDIF
    ALLOCATE(WEIGHTS(2*BANDWIDTH))
    ALLOCATE(WIGNERD(2*BANDWIDTH,BANDWIDTH,2*BANDWIDTH-1,2*BANDWIDTH-1))
    CALL MAKEWEIGHTS(BANDWIDTH)
    CALL CALCWIGNERD(BANDWIDTH)
ENDIF

BW = BANDWIDTH

END SUBROUTINE SETBANDWIDTH

SUBROUTINE MAKEWEIGHTS(BANDWIDTH)

IMPLICIT NONE
INTEGER*8, INTENT(IN) :: BANDWIDTH

DOUBLE PRECISION FUDGE, SINJ
INTEGER*8 J, K

FUDGE = 3.141592653589793D0 / 4 / BANDWIDTH

DO J=1, BANDWIDTH*2
    WEIGHTS(J) = 0
    SINJ = 2.D0 * SIN((2*J-1)*FUDGE) / BANDWIDTH  
    DO K=1,BANDWIDTH  
    WEIGHTS(J) = WEIGHTS(J) + SINJ * SIN((2*J-1)*(2*K-1)*FUDGE) / (2*K - 1) 
    ENDDO
ENDDO

END SUBROUTINE MAKEWEIGHTS

SUBROUTINE RECURRTERMS(J,M1,M2,A,B,C)

! The Wigner little d elements are calculated with a recurrence relation
! This subroutine calculates the appropriate coefficients of the recurrent 
! relation. For more information see:
!
!    Kostelec, P. J., & Rockmore, D. N. (2008). FFTs on the rotation group. 
!    Journal of Fourier Analysis and Applications, 14(2), 145–179. 
!    http://doi.org/10.1007/s00041-008-9013-5

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: J,M1,M2
DOUBLE PRECISION, INTENT(OUT) :: A, B, C

DOUBLE PRECISION T1,T2,T3,T4,T5,DJ,DM1,DM2

DJ = REAL(J,8)
DM1 = REAL(M1,8)
DM2 = REAL(M2,8)

T1 = ((2.D0*DJ +3.D0)/(2.D0*DJ + 1.D0))**0.5D0
T3 = (DJ+1.D0)*(2.D0*DJ+1.D0)
T5 = (((DJ+1.D0)**2-DM1**2)*((DJ+1.D0)**2-DM2**2))**(-0.5)
B = T1*T3*T5

IF (J.EQ.0) THEN
    A=0.D0
    C=0.D0
ELSE
    T2 = ( (2.D0*DJ +3.D0)/(2.D0*DJ-1.D0) )**0.5D0 * (DJ+1.D0)/DJ
    T4 = ( (DJ**2-DM1**2)*(DJ**2-DM2**2) )**(0.5)
    A = T2*T4*T5
    C = M1*M2 / (DJ*(DJ+1.D0))
ENDIF

!WRITE(*,*) J, T1

END SUBROUTINE RECURRTERMS

SUBROUTINE CALCWIGNERD(BANDWIDTH)

!
! Calculates normalised Wigner little-d matrix coefficients for euler angles
! $\beta_k = \frac{\pi(2k+1)}{4 B}, 0\leq k < 2B$ and B = the bandwidth
! stores result in WIGNERD(k, l, m1, m2) in the SOFT module.
!
! Follows method described in:
!
!    Kostelec, P. J., & Rockmore, D. N. (2008). FFTs on the rotation group. 
!    Journal of Fourier Analysis and Applications, 14(2), 145–179. 
!    http://doi.org/10.1007/s00041-008-9013-5

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: BANDWIDTH

DOUBLE PRECISION COSB(2*BANDWIDTH+1), COSB2(2*BANDWIDTH+1), SINB2(2*BANDWIDTH+1)
DOUBLE PRECISION SINCOSB2(2*BANDWIDTH+1), SINDIVCOSB2(2*BANDWIDTH+1)
DOUBLE PRECISION, DIMENSION(0:3*BANDWIDTH-1) :: FACTORIALS
DOUBLE PRECISION FACTOR, FUDGE, BETA, A, B, C, JM1(2*BANDWIDTH+1), T1,T2,T3,T4
INTEGER*8 I,J,M1,M2,IND1,IND2,MAXM

FUDGE = 3.141592653589793D0 / 4 / BANDWIDTH
DO I=1,2*BANDWIDTH
    BETA = FUDGE * (2*I-1)
    COSB(I) = COS(BETA)
    COSB2(I) = COS(BETA/2)
    SINB2(I) = SIN(BETA/2)
    SINCOSB2 = SINB2*COSB2
    SINDIVCOSB2 = SINB2/COSB2
ENDDO 

FACTORIALS(0) = 1.D0
DO I=1, 3*BANDWIDTH-1
    FACTORIALS(I) = I*FACTORIALS(I-1)
ENDDO

! Initialise recurrence
WIGNERD(:,:,:,:) = 0.D0
DO M1=-BANDWIDTH-1,BANDWIDTH-1
    IND1 = MODULO(M1, 2*BANDWIDTH-1) + 1
    DO J=ABS(M1), BANDWIDTH-1
        FACTOR = ((2.D0*J+1.D0)*FACTORIALS(2*J)/FACTORIALS(J+M1)/FACTORIALS(J-M1)/2.D0)**0.5
        IND2 = MODULO(-J, 2*BANDWIDTH-1) + 1
        DO I=1,2*BANDWIDTH
        WIGNERD(I,J+1,J+1,IND1)  = FACTOR * COSB2(I)**(J+M1) * (-SINB2(I))**(J-M1) 
        WIGNERD(I,J+1,IND2,IND1) = FACTOR * COSB2(I)**(J-M1) * (SINB2(I))**(J+M1)
        WIGNERD(I,J+1,IND1,J+1)  = FACTOR * COSB2(I)**(J+M1) * (SINB2(I))**(J-M1) 
        WIGNERD(I,J+1,IND1,IND2) = FACTOR * COSB2(I)**(J-M1) * (-SINB2(I))**(J+M1)
        ENDDO
    ENDDO
ENDDO

! Perform recurrence to calculate Wigner Matrix elements
DO M2=-BANDWIDTH-2,BANDWIDTH-2
    IND2 = MODULO(M2, 2*BANDWIDTH-1) + 1
    DO M1=-BANDWIDTH-2,BANDWIDTH-2
        IND1 = MODULO(M1, 2*BANDWIDTH-1) + 1
        MAXM = MAX(ABS(M1),ABS(M2))
        DO J=MAXM, BANDWIDTH-2
            CALL RECURRTERMS(J,M1,M2,A,B,C)
            DO I=1,2*BANDWIDTH
                WIGNERD(I,J+2,IND1,IND2) = B * (COSB(I) - C) * WIGNERD(I,J+1,IND1,IND2)
            ENDDO
            IF (J.GT.0) THEN
                DO I=1,2*BANDWIDTH
                    WIGNERD(I,J+2,IND1,IND2) = WIGNERD(I,J+2,IND1,IND2) - A*WIGNERD(I,J,IND1,IND2)         
                ENDDO
            ENDIF
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE CALCWIGNERD


SUBROUTINE SOFT(INPUT, OUTPUT, BANDWIDTH)

! Performs discrete SO3 Fourier Analysis for a real input array for a function
! defined on SO(3) returns a complex array of the Fourier Coefficients.

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: BANDWIDTH
DOUBLE PRECISION, INTENT(IN) :: INPUT(2*BANDWIDTH,2*BANDWIDTH,2*BANDWIDTH)
COMPLEX*16, INTENT(OUT) :: OUTPUT(BANDWIDTH, 2*BANDWIDTH-1, 2*BANDWIDTH-1)
!COMPLEX*16, INTENT(OUT) :: TEMP(2*BANDWIDTH, 2*BANDWIDTH, 2*BANDWIDTH)

INCLUDE "fftw3.f90"
COMPLEX*16 IN1D(2*BANDWIDTH), OUT1D(2*BANDWIDTH), TEMP(2*BANDWIDTH, 2*BANDWIDTH, 2*BANDWIDTH)
INTEGER*8 PLAN, K1,K2,K3,M1,M2,I1,I2,IND1,IND2,J,MAXM


CALL SETBANDWIDTH(BANDWIDTH)

CALL DFFTW_PLAN_DFT_1D(PLAN, (2*BANDWIDTH), IN1D, OUT1D, FFTW_FORWARD, FFTW_ESTIMATE)

! Do FFT on axis 1
DO K1=1,2*BANDWIDTH
    DO K2=1,2*BANDWIDTH
        DO K3=1,2*BANDWIDTH
            IN1D(K3) = CMPLX(INPUT(K3,K2,K1),0.D0, 16)
        ENDDO
        CALL DFFTW_EXECUTE_(PLAN, IN1D, OUT1D)
        DO K3=1,2*BANDWIDTH
            TEMP(K3,K2,K1) = OUT1D(K3)
        ENDDO
    ENDDO
ENDDO

! Do FFT on axis 3
DO K1=1,2*BANDWIDTH
    DO K2=1,2*BANDWIDTH
        DO K3=1,2*BANDWIDTH
            IN1D(K3) = TEMP(K2,K1,K3)
        ENDDO
        CALL DFFTW_EXECUTE_(PLAN, IN1D, OUT1D)
        DO K3=1,2*BANDWIDTH
            TEMP(K2,K1,K3) = OUT1D(K3)/(2*BANDWIDTH)**2
        ENDDO
    ENDDO
ENDDO

! Perform Discrete Wigner Transform
OUTPUT = 0
DO M2=-BANDWIDTH-1,BANDWIDTH-1
    I2 = MODULO(M2, 2*BANDWIDTH) + 1
    IND2 = MODULO(M2, 2*BANDWIDTH-1) + 1
    DO M1=-BANDWIDTH-1,BANDWIDTH-1
        I1 = MODULO(M1, 2*BANDWIDTH) + 1
        IND1 = MODULO(M1, 2*BANDWIDTH-1) + 1
        MAXM = MAX(ABS(M1),ABS(M2))
        DO J=MAXM, BANDWIDTH-1
            DO K1=1,2*BANDWIDTH
                OUTPUT(J+1,IND1,IND2) = OUTPUT(J+1,IND1,IND2) + WIGNERD(K1,J+1,IND1,IND2)*WEIGHTS(K1)*TEMP(I1,K1,I2)
            ENDDO
        ENDDO
    ENDDO
ENDDO

CALL DFFTW_DESTROY_PLAN_(PLAN)

END SUBROUTINE SOFT

SUBROUTINE ISOFT(INPUT, OUTPUT, BANDWIDTH)

! Performs SO3 Fourier Synthesis for a complex input array of Fourier Coefficients
! Generates a complex output array.

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: BANDWIDTH
COMPLEX*16, INTENT(IN) :: INPUT(BANDWIDTH, 2*BANDWIDTH-1, 2*BANDWIDTH-1)
COMPLEX*16, INTENT(OUT) :: OUTPUT(2*BANDWIDTH,2*BANDWIDTH,2*BANDWIDTH)
!COMPLEX*16, INTENT(OUT) :: TEMP(2*BANDWIDTH, 2*BANDWIDTH, 2*BANDWIDTH)

INCLUDE "fftw3.f90"
COMPLEX*16 IN1D(2*BANDWIDTH), OUT1D(2*BANDWIDTH), TEMP(2*BANDWIDTH, 2*BANDWIDTH, 2*BANDWIDTH)
INTEGER*8 PLAN, K1,K2,K3,M1,M2,I1,I2,IND1,IND2,J,MAXM


CALL SETBANDWIDTH(BANDWIDTH)

CALL DFFTW_PLAN_DFT_1D(PLAN, (2*BANDWIDTH), IN1D, OUT1D, FFTW_BACKWARD, FFTW_ESTIMATE)

! Discrete inverse Wigner Transform
TEMP = 0
DO M2=-BANDWIDTH-1,BANDWIDTH-1
    I2 = MODULO(M2, 2*BANDWIDTH) + 1
    IND2 = MODULO(M2, 2*BANDWIDTH-1) + 1
    DO M1=-BANDWIDTH-1,BANDWIDTH-1
        I1 = MODULO(M1, 2*BANDWIDTH) + 1
        IND1 = MODULO(M1, 2*BANDWIDTH-1) + 1
        MAXM = MAX(ABS(M1),ABS(M2))
        DO K1=1,2*BANDWIDTH
            DO J=MAXM, BANDWIDTH-1
                TEMP(I1,K1,I2) = TEMP(I1,K1,I2) + WIGNERD(K1,J+1,IND1,IND2)*INPUT(J+1,IND1,IND2)
            ENDDO
        ENDDO
    ENDDO
ENDDO

! Inverse Fourier Transform on axis 3
DO K1=1,2*BANDWIDTH
    DO K2=1,2*BANDWIDTH
        DO K3=1,2*BANDWIDTH
            IN1D(K3) = TEMP(K2,K1,K3)
        ENDDO
        CALL DFFTW_EXECUTE_(PLAN, IN1D, OUT1D)
        DO K3=1,2*BANDWIDTH
            TEMP(K2,K1,K3) = OUT1D(K3)
        ENDDO
    ENDDO
ENDDO

! Inverse Fourier Transform on axis 1
DO K1=1,2*BANDWIDTH
    DO K2=1,2*BANDWIDTH
        DO K3=1,2*BANDWIDTH
            IN1D(K3) = TEMP(K3,K2,K1)
        ENDDO
        CALL DFFTW_EXECUTE_(PLAN, IN1D, OUT1D)
        DO K3=1,2*BANDWIDTH
            OUTPUT(K3,K2,K1) = OUT1D(K3)!/(2*BANDWIDTH)**2
        ENDDO
    ENDDO
ENDDO

CALL DFFTW_DESTROY_PLAN_(PLAN)

ENDSUBROUTINE ISOFT

! TODO Implement version of these algorithms that take advantage of the symmetries
! imposed by a real input array.

END MODULE DSOFT

















