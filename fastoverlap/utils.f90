!    FASTOVERLAP
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


MODULE FASTOVERLAPUTILS

!***********************************************************************
! This module contains some subroutines that are useful for FASTOVERLAP 
! alignment for both periodic and isolated structures
!***********************************************************************
! Subroutines:
!     Peakfinding subroutines:
!         SETINDEXES
!         SETFSPACE
!         GAUSSIAN
!         FCN
!         FIT
!         FINDPEAK
!         FINDPEAKS
!     FFT subroutines
!         FFT3D
!         IFFT3D
!         FFT1D
!         IFFT1D
!***********************************************************************

IMPLICIT NONE

! Variables and arrays needed for peakfinding
INTEGER, PARAMETER :: N=11, DEFAULTWIDTH=2
DOUBLE PRECISION, PARAMETER :: DEFAULTTOL=1.D-6
INTEGER, SAVE :: FSIZE, FSHAPE(3)
DOUBLE PRECISION, SAVE, ALLOCATABLE :: FSPACE(:,:,:), FVEC(:), FJAC(:,:)

! An array of the fastest length arrays on which to perform FFTs
INTEGER*8, SAVE :: FASTLEN(200) = (/1, 2, 3, 4, 5, 6, 8, 8, 9, 10, 12, 12, 15, &
    15, 15, 16, 18, 18, 20, 20, 24, 24, 24, 24, 25, 27, 27, 30, 30, 30, 32, &
    32, 36, 36, 36, 36, 40, 40, 40, 40, 45, 45, 45, 45, 45, 48, 48, 48, 50, &
    50, 54, 54, 54, 54, 60, 60, 60, 60, 60, 60, 64, 64, 64, 64, 72, 72, 72, &
    72, 72, 72, 72, 72, 75, 75, 75, 80, 80, 80, 80, 80, 81, 90, 90, 90, 90, &
    90, 90, 90, 90, 90, 96, 96, 96, 96, 96, 96, 100, 100, 100, 100, 108, 108, &
    108, 108, 108, 108, 108, 108, 120, 120, 120, 120, 120, 120, 120, 120, 120,&
    120, 120, 120, 125, 125, 125, 125, 125, 128, 128, 128, 135, 135, 135, 135,&
    135, 135, 135, 144, 144, 144, 144, 144, 144, 144, 144, 144, 150, 150, 150,&
    150, 150, 150, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 162, 162,&
    180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180,&
    180, 180, 180, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192,&
    200, 200, 200, 200, 200, 200, 200, 200/)

CONTAINS

SUBROUTINE SETINDEXES(NEWSHAPE)

! Helper routine to allocate memory to appropriate arrays needed to perform
! Levenberg-Marquardt non-linear least-squares curve fitting to find peaks

IMPLICIT NONE

INTEGER, INTENT(IN) :: NEWSHAPE(3)

IF (.NOT.ALL(FSHAPE.EQ.NEWSHAPE)) THEN
    FSHAPE = NEWSHAPE    
    IF(ALLOCATED(FSPACE)) THEN
        DEALLOCATE(FSPACE)
    ENDIF
    IF(ALLOCATED(FVEC)) THEN
        DEALLOCATE(FVEC)
    ENDIF
    IF(ALLOCATED(FJAC)) THEN
        DEALLOCATE(FJAC)
    ENDIF
    
    ALLOCATE( FSPACE( FSHAPE(1),FSHAPE(2),FSHAPE(3) ) )
    FSIZE = SIZE(FSPACE)
    
    ALLOCATE(FVEC(FSIZE))
    ALLOCATE(FJAC(N,FSIZE))
ENDIF

END SUBROUTINE SETINDEXES

!***********************************************************************

SUBROUTINE SETFSPACE(NEWFSPACE)

IMPLICIT NONE

!INTEGER, INTENT(IN) :: NX,NY,NZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:,:) :: NEWFSPACE
!INTEGER NSHAPE(3)

!NSHAPE=(/NX,NY,NZ/)
CALL SETINDEXES(SHAPE(NEWFSPACE))

FSPACE = NEWFSPACE

END SUBROUTINE SETFSPACE

!***********************************************************************

SUBROUTINE GAUSSIAN(X,NX,NY,NZ,FOUT)

! Routine to calculate the values of a 3-D gaussian
! FOUT(IX, IY, IZ) = A * Exp(-(I-I0)^T SIGMA (I-I0))
! I = (/IX, IY, IZ/)
!specified by the parameter vector X:
! X = (\A, mean, SIGMA(1,1), SIGMA(2,2), SIGMA(3,3), SIGMA(1,2),SIGMA(2,3),SIGMA(1,3), I0(1), I0(2), I0(3) \)

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: NX, NY, NZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X
DOUBLE PRECISION, INTENT(OUT) :: FOUT(NX,NY,NZ)

INTEGER*8 IX,IY,IZ,J
DOUBLE PRECISION SIGMA(3,3), A, MEAN, Y, EXPY, FY, IND0(3), DY(3)

A = X(1)
MEAN = X(2)
SIGMA(1,1) = X(3)
SIGMA(2,2) = X(4)
SIGMA(3,3) = X(5)
SIGMA(1,2) = X(6)
SIGMA(2,1) = 0.D0!X(6)
SIGMA(2,3) = X(7)
SIGMA(3,2) = 0.D0!X(7)
SIGMA(1,3) = X(8)
SIGMA(3,1) = 0.D0!X(8)
!IND0 = X(9:11)

DO IZ=1,NZ
    DO IY=1,NY
        DO IX=1,NX
            IND0 = (/IX,IY,IZ/) - X(9:11)
            DO J=1,3
                DY(J) = SUM(SIGMA(J,:)*IND0)
            ENDDO
            Y = SUM(IND0*DY)
            EXPY = EXP(-Y)
            FOUT(IX,IY,IZ) =  (A*EXPY + MEAN)
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE GAUSSIAN

!***********************************************************************

SUBROUTINE FCN(M,N,X,FVEC,FJAC,LDFJAC,IFLAG)

! 
! subroutine passed to lmder1 to perform least squares regression, minimizing
! SUM((FOUT - FSPACE)**2)
! where  FOUT(IX, IY, IZ) = A * Exp(-(I-I0)^T SIGMA (I-I0))
! I = (/IX, IY, IZ/)
!specified by the parameter vector X:
! X = (\A, mean, SIGMA(1,1), SIGMA(2,2), SIGMA(3,3), SIGMA(1,2),SIGMA(2,3),SIGMA(1,3), I0(1), I0(2), I0(3) \)
! M = SIZE(FSPACE) is the number of observations
! LDFJAC = N specifies the dimension of the jacobian matrix
! N = 11 is the number of parameters to optimise
! If IFLAG=1 then calculates FVEC, the vector of square difference of each observation
! If IFLAG=2 then calculates FVEC and FJAC, the jacobian maxtrix of FVEC

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: LDFJAC, N, M, IFLAG
DOUBLE PRECISION, INTENT(OUT) :: FJAC(LDFJAC, N), FVEC(M)
DOUBLE PRECISION, INTENT(INOUT) :: X(N)

DOUBLE PRECISION SIGMA(3,3), A, MEAN, Y, EXPY, FY, DIFF, DY(3), IND0(3)
INTEGER*8 :: I,J,K,IND(3)!,IX,IY,IZ!,S(2)=(/3,1/)

! if IFLAG =/= 1/2 then do nothing...
IF(IFLAG.EQ.1 .OR. IFLAG.EQ.2) THEN
A = X(1)
MEAN = X(2)
SIGMA(1,1) = X(3)
SIGMA(2,2) = X(4)
SIGMA(3,3) = X(5)
SIGMA(1,2) = X(6)
SIGMA(2,1) = 0.D0!X(6)
SIGMA(2,3) = X(7)
SIGMA(3,2) = 0.D0!X(7)
SIGMA(1,3) = X(8)
SIGMA(3,1) = 0.D0!X(8)
!IND0 = X(9:11)

DO I=1,M
    !Some pointer arithmetic to get the 3D index location
    !I miss 0-indexing
    IND(1) = (I-1)/FSHAPE(2)/FSHAPE(3) + 1
    IND(2) = MOD((I-1)/FSHAPE(3), FSHAPE(2)) + 1
    IND(3) = MOD(I-1, FSHAPE(3)) + 1
    IND0 = IND - X(9:11)
    !Y = 0.D0
    DO J=1,3
        DY(J) = SUM(SIGMA(J,:)*IND0)
    ENDDO
    Y = SUM(IND0*DY)
    EXPY = EXP(-Y)
    FY = (A*EXPY + MEAN)
    DIFF = (FY - FSPACE(IND(1),IND(2),IND(3)))
    FVEC(I) = DIFF**2
    IF(IFLAG.EQ.2) THEN
        ! Calculating Jacobian
        FJAC(I,1) = 2 * EXPY * DIFF
        FJAC(I,2) = 2 * DIFF
        FJAC(I,3) = -(IND0(1)*IND0(1))*A*EXPY * DIFF * 2
        FJAC(I,4) = -(IND0(2)*IND0(2))*A*EXPY * DIFF * 2
        FJAC(I,5) = -(IND0(3)*IND0(3))*A*EXPY * DIFF * 2
        FJAC(I,6) = -(IND0(1)*IND0(2))*A*EXPY * DIFF * 2
        FJAC(I,7) = -(IND0(2)*IND0(3))*A*EXPY * DIFF * 2
        FJAC(I,8) = -(IND0(1)*IND0(3))*A*EXPY * DIFF * 2
        FJAC(I,9:11) = 4 * DY * A * EXPY * DIFF
    ENDIF
ENDDO
ENDIF

END SUBROUTINE FCN

!***********************************************************************

SUBROUTINE FIT(X, NEWFSPACE, NX, NY, NZ, INFO, TOL)

! This fits a 3 dimensional gaussian of the form
! A exp (- (I-I0)T Sigma (I-I0) ) + mean
! Where I is the 3-D vector of the indexes
! To the 3 dimensional array specified by FSPACE
! This uses the Levenberg-Marquardt method. 
! Usage:
! CALL FIT(X0, FSPACE, INFO, TOL(optional))
! X0 = (\A, mean, SIGMA(1,1), SIGMA(2,2), SIGMA(3,3), SIGMA(1,2),SIGMA(2,3),SIGMA(1,3), I0(1), I0(2), I0(3) \)
!INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares is
!       possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: NX,NY,NZ
DOUBLE PRECISION, INTENT(IN) :: NEWFSPACE(NX,NY,NZ)
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: TOL
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: X
INTEGER*8, INTENT(OUT) :: INFO

DOUBLE PRECISION USETOL

IF (PRESENT(TOL)) THEN
    USETOL = TOL
ELSE
    USETOL = DEFAULTTOL
ENDIF

CALL SETFSPACE(NEWFSPACE)
!Perform Levenberg-Marquardt non-linear least square regression
CALL LMDER1 (FCN, FSIZE, 11, X, FVEC, FJAC, FSIZE, USETOL, INFO)

END SUBROUTINE FIT

!***********************************************************************

SUBROUTINE FINDPEAK (A, WIDTH, X, INFO, TOL, AMAX)

! Finds maximum value of 3D array A Selects the indexes within WIDTH
! Fits Gaussian to these indexes, then outputs the fit as X

! ASSUMES PERIODIC BOUNDARY CONDITIONS

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:,:) :: A
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: TOL
INTEGER*8, INTENT(IN) :: WIDTH
DOUBLE PRECISION, INTENT(OUT) :: X(11)
INTEGER*8, INTENT(OUT) :: INFO, AMAX(3)

DOUBLE PRECISION FSPACE(WIDTH*2+1,WIDTH*2+1,WIDTH*2+1)
DOUBLE PRECISION MAXA, MEANA
INTEGER*8 ASHAPE(3),I1,I2,I3,IND(3) !AMAX(3)

AMAX = MAXLOC(A)
MEANA = SUM(A)/SIZE(A)
MAXA = MAXVAL(A) - MEANA
! initialise guess for parameter array
X = (/MAXA,MEANA,1.D0,1.D0,1.D0,0.D0,0.D0,0.D0,WIDTH+1.D0,WIDTH+1.D0,WIDTH+1.D0/)
ASHAPE = SHAPE(A)

! selecting subarray to fit peak to
DO I3=1,2*WIDTH+1
    DO I2=1,2*WIDTH+1
        DO I1=1,2*WIDTH+1
            ! Ensures periodic boundary conditions
            IND = MODULO(AMAX+(/I1,I2,I3/)-2-WIDTH,ASHAPE) + 1
            FSPACE(I1,I2,I3) = A(IND(1),IND(2),IND(3))
        ENDDO
    ENDDO
ENDDO

IF(PRESENT(TOL)) THEN
    CALL FIT(X, FSPACE, WIDTH*2+1, WIDTH*2+1, WIDTH*2+1,INFO, TOL)
ELSE
    CALL FIT(X, FSPACE, WIDTH*2+1, WIDTH*2+1, WIDTH*2+1, INFO)
ENDIF

END SUBROUTINE FINDPEAK

!***********************************************************************

SUBROUTINE FINDPEAKS(FSPACE, PEAKS, AMPLITUDES, NPEAKS)

! This finds up to npeaks of a 3D periodic array
! The locations are returned in peaks as fractional index coordinates
! Amplitude gives the relative amplitude of each of the peaks
! NPEAKS gives the actual number of peaks found

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:,:) :: FSPACE
INTEGER*8, INTENT(INOUT) :: NPEAKS
!INTEGER, INTENT(IN), OPTIONAL :: WIDTH
DOUBLE PRECISION, INTENT(OUT) :: PEAKS(NPEAKS,3), AMPLITUDES(NPEAKS)

INTEGER*8 WIDTH, NFOUND, FSHAPE(3), INFO, N, FMAX(3)
DOUBLE PRECISION T, X(11), PEAK(3)
DOUBLE PRECISION, ALLOCATABLE :: FSPACECOPY(:,:,:), GAUSARRAY(:,:,:)

WIDTH = DEFAULTWIDTH
FSHAPE = SHAPE(FSPACE)
ALLOCATE(FSPACECOPY(FSHAPE(1),FSHAPE(2),FSHAPE(3)))
ALLOCATE(GAUSARRAY(FSHAPE(1),FSHAPE(2),FSHAPE(3)))
FSPACECOPY = FSPACE

NFOUND = 0
DO WHILE(NFOUND.EQ.0)
    DO N=1,NPEAKS
        CALL FINDPEAK(FSPACECOPY, WIDTH, X, INFO, DEFAULTTOL, FMAX)

        IF(INFO.EQ.4.OR.INFO.EQ.5) THEN
            EXIT
        ELSE
            ! Find the location of the peak and subtract this peak from the
            ! copy of the data
            NFOUND = NFOUND + 1
            PEAK = (X(9:11) - WIDTH - 1 + FMAX)
            PEAKS(N,:) = PEAK
            AMPLITUDES(N) = X(1)
            X(9:11) = PEAK
            CALL GAUSSIAN(X,FSHAPE(1),FSHAPE(2),FSHAPE(3),GAUSARRAY)
            FSPACECOPY = FSPACECOPY - GAUSARRAY
        ENDIF
    ENDDO
    ! If we've failed to find any peaks, increase the size of the box and start again
    IF (NFOUND.EQ.0) THEN
        WIDTH = WIDTH + 1
    ENDIF
ENDDO

NPEAKS = NFOUND

END SUBROUTINE FINDPEAKS

!***********************************************************************
! FFT subroutines
!***********************************************************************
    
SUBROUTINE FFT3D(NX, NY, NZ, IN, OUT)
! calculates forward FFT in 3D
IMPLICIT NONE

INTEGER*8, INTENT(IN) :: NX, NY, NZ
COMPLEX*16, INTENT(IN) :: IN(NX, NY, NZ)
COMPLEX*16, INTENT(OUT) :: OUT(NX, NY, NZ)

INCLUDE "fftw3.f90"
INTEGER*8 PLAN_FORWARD

CALL DFFTW_PLAN_DFT_3D_(PLAN_FORWARD, NX, NY, NZ, IN, OUT, FFTW_FORWARD, FFTW_ESTIMATE )
CALL DFFTW_EXECUTE_(PLAN_FORWARD)
!CALL DFFTW_DESTROY_PLAN(PLAN_FORWARD)

END SUBROUTINE FFT3D

!***********************************************************************

SUBROUTINE IFFT3D(NX, NY, NZ, IN, OUT)

! calculates UNNORMALISED inverse fourier transform so,
! IN == IFFT3D(NX,NY,NZ, FFT3D(NX,NY,NZ, IN))/(NX*NY*NZ)

IMPLICIT NONE

INTEGER*8, INTENT(IN) :: NX, NY, NZ
COMPLEX*16, INTENT(IN) :: IN(NX, NY, NZ)
COMPLEX*16, INTENT(OUT) :: OUT(NX, NY, NZ)

INCLUDE "fftw3.f90"
INTEGER*8 PLAN_BACKWARD

CALL DFFTW_PLAN_DFT_3D_(PLAN_BACKWARD,NX,NY,NZ,IN,OUT,FFTW_BACKWARD,FFTW_ESTIMATE)
CALL DFFTW_EXECUTE_(PLAN_BACKWARD)
CALL DFFTW_DESTROY_PLAN_(PLAN_BACKWARD)

END SUBROUTINE IFFT3D

SUBROUTINE FFT1D(N, IN, OUT)
! calculates forward FFT in 1D

IMPLICIT NONE

INTEGER*4, INTENT(IN) :: N
COMPLEX*16, INTENT(IN) :: IN(N)
COMPLEX*16, INTENT(OUT) :: OUT(N)

INCLUDE "fftw3.f90"
INTEGER*8 PLAN_FORWARD

CALL DFFTW_PLAN_DFT_1D_(PLAN_FORWARD, N, IN, OUT, FFTW_FORWARD, FFTW_ESTIMATE )
CALL DFFTW_EXECUTE_(PLAN_FORWARD)
CALL DFFTW_DESTROY_PLAN_(PLAN_FORWARD)

END SUBROUTINE FFT1D

!***********************************************************************

SUBROUTINE IFFT1D(N, IN, OUT)

! calculates UNNORMALISED inverse fourier transform so,
! IN == IFFT1D(N, FFT1D(N, IN))/N

IMPLICIT NONE

INTEGER*4, INTENT(IN) :: N
COMPLEX*16, INTENT(IN) :: IN(N)
COMPLEX*16, INTENT(OUT) :: OUT(N)

INCLUDE "fftw3.f90"
INTEGER*8 PLAN_BACKWARD

CALL DFFTW_PLAN_DFT_1D_(PLAN_BACKWARD, N, IN, OUT, FFTW_BACKWARD, FFTW_ESTIMATE )
CALL DFFTW_EXECUTE_(PLAN_BACKWARD)
CALL DFFTW_DESTROY_PLAN_(PLAN_BACKWARD)

END SUBROUTINE IFFT1D

END MODULE FASTOVERLAPUTILS