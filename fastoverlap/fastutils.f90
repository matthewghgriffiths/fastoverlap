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


!    Fortran 90/95 modules:
!      fastoverlaputils --- fshape,fspace,fvec,defaulttol,fsize,n,fastlen,defaultwidth,fjac,setindexes(),setfspace(),gaussian(),fcn(),fit(),findpeak(),findpeaks(),fft3d(),ifft3d(),fft1d(),ifft1d().
!    Functions:
!      rlegendrel0 = rlegendrel0(l,z)
!      rlegendrem0 = rlegendrem0(m,l,z)
!      rlegendrem1 = rlegendrem1(m,l,z)
!      envj = envj(n,x)
!      msta1 = msta1(x,mp)
!      msta2 = msta2(x,n,mp)
!      nm,si = sphi(n,x)
!      hg = hyp1f1(ain,bin,xin)
!      ga = gamma(x)
!      fvec,fjac,info,nfev,njev,qtf = lmder(fcn,m,x,ldfjac,ftol,xtol,gtol,maxfev,diag,mode,factor,nprint,ipvt,n=len(x),fcn_extra_args=())
!      fvec,fjac,info = lmder1(fcn,m,x,ldfjac,tol,n=len(x),fcn_extra_args=())
!      enorm = enorm(x,n=len(x))
!      enorm2 = enorm2(x,n=len(x))
!      lmpar(r,ipvt,diag,qtb,delta,par,x,sdiag,n=shape(r,1),ldr=shape(r,0))
!      qrsolv(r,ipvt,diag,qtb,x,sdiag,n=shape(r,1),ldr=shape(r,0))
!      qrfac(m,a,pivot,ipvt,rdiag,acnorm,n=shape(a,1),lda=shape(a,0),lipvt=len(ipvt))
!      xmed = median(x,n=len(x))

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
INTEGER, SAVE :: FASTLEN(200) = (/1, 2, 3, 4, 5, 6, 8, 8, 9, 10, 12, 12, 15, &
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

INTEGER, INTENT(IN) :: NX, NY, NZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: X
DOUBLE PRECISION, INTENT(OUT) :: FOUT(NX,NY,NZ)

INTEGER IX,IY,IZ,J
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

INTEGER, INTENT(IN) :: LDFJAC, N, M, IFLAG
DOUBLE PRECISION, INTENT(OUT) :: FJAC(LDFJAC, N), FVEC(M)
DOUBLE PRECISION, INTENT(INOUT) :: X(N)

DOUBLE PRECISION SIGMA(3,3), A, MEAN, Y, EXPY, FY, DIFF, DY(3), IND0(3)
INTEGER :: I,J,K,IND(3)!,IX,IY,IZ!,S(2)=(/3,1/)

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

INTEGER, INTENT(IN) :: NX,NY,NZ
DOUBLE PRECISION, INTENT(IN) :: NEWFSPACE(NX,NY,NZ)
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: TOL
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: X
INTEGER, INTENT(OUT) :: INFO

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
INTEGER, INTENT(IN) :: WIDTH
DOUBLE PRECISION, INTENT(OUT) :: X(11)
INTEGER, INTENT(OUT) :: INFO, AMAX(3)

DOUBLE PRECISION FSPACE(WIDTH*2+1,WIDTH*2+1,WIDTH*2+1)
DOUBLE PRECISION MAXA, MEANA
INTEGER ASHAPE(3),I1,I2,I3,IND(3) !AMAX(3)

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
INTEGER, INTENT(INOUT) :: NPEAKS
!INTEGER, INTENT(IN), OPTIONAL :: WIDTH
DOUBLE PRECISION, INTENT(OUT) :: PEAKS(NPEAKS,3), AMPLITUDES(NPEAKS)

INTEGER WIDTH, NFOUND, FSHAPE(3), INFO, N, FMAX(3)
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

INTEGER, INTENT(IN) :: NX, NY, NZ
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

INTEGER, INTENT(IN) :: NX, NY, NZ
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

!***********************************************************************

! Some helper functions for calculating various orthogonal polynomials

!***********************************************************************


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



!    CODE REPRODUCED FROM MINPACK UNDER THE GNU LPGL LICENCE:

!    REFERENCES:

!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
    
!    Jorge More, Danny Sorenson, Burton Garbow, Kenneth Hillstrom,
!    The MINPACK Project,
!    in Sources and Development of Mathematical Software,
!    edited by Wayne Cowell,
!    Prentice-Hall, 1984,
!    ISBN: 0-13-823501-5,
!    LC: QA76.95.S68.



subroutine lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
  diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

!*****************************************************************************80
!
!! LMDER minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, real ( kind = 8 ) GTOL.  Termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN with IFLAG = 1 is at least MAXFEV by the end of an iteration.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of
!       squares are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN with
!    IFLAG = 1.
!
!    Output, integer ( kind = 4 ) NJEV, the number of calls to FCN with
!    IFLAG = 2.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P
!    such that JAC*P = Q*R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column
!    IPVT(J) of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), contains the first N elements of Q'*FVEC.
!
  implicit none

  integer ( kind = 4 ), INTENT(IN) :: ldfjac
  integer ( kind = 4 ), INTENT(IN) ::  m
  integer ( kind = 4 ), INTENT(IN) ::  n

  real ( kind = 8 ) actred
  real ( kind = 8 ) delta
  real ( kind = 8 ), INTENT(INOUT) :: diag(n)
  real ( kind = 8 ) dirder
  real ( kind = 8 ) enorm
  real ( kind = 8 ) epsmch
  real ( kind = 8 ), INTENT(IN) :: factor
  external  fcn
  real ( kind = 8 ), INTENT(OUT) :: fjac(ldfjac,n)
  real ( kind = 8 ) fnorm
  real ( kind = 8 ) fnorm1
  real ( kind = 8 ), INTENT(IN) :: ftol
  real ( kind = 8 ), INTENT(OUT) :: fvec(m)
  real ( kind = 8 ) gnorm
  real ( kind = 8 ), INTENT(IN) :: gtol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), INTENT(OUT) :: info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ), INTENT(IN) :: maxfev
  integer ( kind = 4 ), INTENT(IN) :: mode
  integer ( kind = 4 ), INTENT(OUT) :: nfev
  integer ( kind = 4 ), INTENT(OUT) :: njev
  integer ( kind = 4 ), INTENT(IN) :: nprint
  real ( kind = 8 ) par
  logical pivot
  real ( kind = 8 ) pnorm
  real ( kind = 8 ) prered
  real ( kind = 8 ), INTENT(OUT) :: qtf(n)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) wa3(n)
  real ( kind = 8 ) wa4(m)
  real ( kind = 8 ) xnorm
  real ( kind = 8 ), INTENT(INOUT) ::  x(n)
  real ( kind = 8 ), INTENT(IN) :: xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    go to 300
  end if

  if ( m < n ) then
    go to 300
  end if

  if ( ldfjac < m &
    .or. ftol < 0.0D+00 .or. xtol < 0.0D+00 .or. gtol < 0.0D+00 &
     .or. maxfev <= 0 .or. factor <= 0.0D+00 ) then
    go to 300
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
  nfev = 1
  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
30   continue
!
!  Calculate the jacobian matrix.
!
    iflag = 2
    call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )

    njev = njev + 1

    if ( iflag < 0 ) then
      go to 300
    end if
!
!  If requested, call FCN to enable printing of iterates.
!
    if ( 0 < nprint ) then
      iflag = 0
      if ( mod ( iter - 1, nprint ) == 0 ) then
        call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
      end if
      if ( iflag < 0 ) then
        go to 300
      end if
    end if
!
!  Compute the QR factorization of the jacobian.
!
    pivot = .true.
    call qrfac ( m, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2 )
!
!  On the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then
        diag(1:n) = wa2(1:n)
        do j = 1, n
          if ( wa2(j) == 0.0D+00 ) then
            diag(j) = 1.0D+00
          end if
        end do
      end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n) * x(1:n)

      xnorm = enorm ( n, wa3 )
      delta = factor * xnorm
      if ( delta == 0.0D+00 ) then
        delta = factor
      end if

    end if
!
!  Form Q'*FVEC and store the first N components in QTF.
!
    wa4(1:m) = fvec(1:m)

    do j = 1, n

      if ( fjac(j,j) /= 0.0D+00 ) then
        sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
        temp = - sum2 / fjac(j,j)
        wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
      end if

      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)

    end do
!
!  Compute the norm of the scaled gradient.
!
    gnorm = 0.0D+00

    if ( fnorm /= 0.0D+00 ) then

      do j = 1, n
        l = ipvt(j)
        if ( wa2(l) /= 0.0D+00 ) then
          sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
          gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
        end if
      end do

    end if
!
!  Test for convergence of the gradient norm.
!
    if ( gnorm <= gtol ) then
      info = 4
      go to 300
    end if
!
!  Rescale if necessary.
!
    if ( mode /= 2 ) then
      do j = 1, n
        diag(j) = max ( diag(j), wa2(j) )
      end do
    end if
!
!  Beginning of the inner loop.
!
200    continue
!
!  Determine the Levenberg-Marquardt parameter.
!
    call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction p and x + p. calculate the norm of p.
!
    wa1(1:n) = - wa1(1:n)
    wa2(1:n) = x(1:n) + wa1(1:n)
    wa3(1:n) = diag(1:n) * wa1(1:n)

    pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
    if ( iter == 1 ) then
      delta = min ( delta, pnorm )
    end if
!
!  Evaluate the function at x + p and calculate its norm.
!
    iflag = 1
    call fcn ( m, n, wa2, wa4, fjac, ldfjac, iflag )

    nfev = nfev + 1

    if ( iflag < 0 ) then
      go to 300
    end if

    fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
    actred = -1.0D+00
    if ( 0.1D+00 * fnorm1 < fnorm ) then
      actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
    end if
!
!  Compute the scaled predicted reduction and
!  the scaled directional derivative.
!
    do j = 1, n
      wa3(j) = 0.0D+00
      l = ipvt(j)
      temp = wa1(l)
      wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
    end do

    temp1 = enorm ( n, wa3 ) / fnorm
    temp2 = ( sqrt ( par ) * pnorm ) / fnorm
    prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
    dirder = - ( temp1 ** 2 + temp2 ** 2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
    if ( prered /= 0.0D+00 ) then
      ratio = actred / prered
    else
      ratio = 0.0D+00
    end if
!
!  Update the step bound.
!
    if ( ratio <= 0.25D+00 ) then

      if ( 0.0D+00 <= actred ) then
        temp = 0.5D+00
      end if

      if ( actred < 0.0D+00 ) then
        temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
      end if

      if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
        temp = 0.1D+00
      end if

      delta = temp * min ( delta, pnorm / 0.1D+00 )
      par = par / temp

    else

      if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
        delta = 2.0D+00 * pnorm
        par = 0.5D+00 * par
      end if

    end if
!
!  Successful iteration.
!
!  Update X, FVEC, and their norms.
!
    if ( 0.0001D+00 <= ratio ) then
      x(1:n) = wa2(1:n)
      wa2(1:n) = diag(1:n) * x(1:n)
      fvec(1:m) = wa4(1:m)
      xnorm = enorm ( n, wa2 )
      fnorm = fnorm1
      iter = iter + 1
    end if
!
!  Tests for convergence.
!
    if ( abs ( actred) <= ftol .and. &
      prered <= ftol .and. &
      0.5D+00 * ratio <= 1.0D+00 ) then
      info = 1
    end if

    if ( delta <= xtol * xnorm ) then
      info = 2
    end if

    if ( abs ( actred) <= ftol .and. prered <= ftol &
      .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) then
      info = 3
    end if

    if ( info /= 0 ) then
      go to 300
    end if
!
!  Tests for termination and stringent tolerances.
!
    if ( nfev >= maxfev ) then
      info = 5
    end if

    if ( abs ( actred ) <= epsmch .and. prered <= epsmch &
      .and. 0.5D+00 * ratio <= 1.0D+00 ) then
      info = 6
    end if

    if ( delta <= epsmch * xnorm ) then
      info = 7
    end if

    if ( gnorm <= epsmch ) then
      info = 8
    end if

    if ( info /= 0 ) then
      go to 300
    end if
!
!  End of the inner loop. repeat if iteration unsuccessful.
!
    if ( ratio < 0.0001D+00 ) then
      go to 200
    end if
!
!  End of the outer loop.
!
    go to 30

  300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
  end if

  return
end

subroutine lmder1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! LMDER1 minimizes M functions in N variables by Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDER1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDER.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, is the leading dimension of FJAC,
!    which must be no less than M.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
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
!
  implicit none

  integer ( kind = 4 ), INTENT(IN) ::  ldfjac
  integer ( kind = 4 ), INTENT(IN) ::  m
  integer ( kind = 4 ), INTENT(IN) ::  n

  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) factor
  external fcn
  real ( kind = 8 ), INTENT(OUT) ::  fjac(ldfjac,n)
  real ( kind = 8 ) ftol
  real ( kind = 8 ), INTENT(OUT) ::  fvec(m)
  real ( kind = 8 ) gtol
  integer ( kind = 4 ), INTENT(OUT) ::  info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) qtf(n)
  real ( kind = 8 ), INTENT(IN) ::  tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( ldfjac < m ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  factor = 100.0D+00
  maxfev = 100 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  mode = 1
  nprint = 0

  !hack to get f2py to work, not needed for pure FORTRAN compilation
  iflag = 0
  call fcn ( m, n, x, fvec, fjac, ldfjac, iflag )

  call lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end

function enorm ( n, x )

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) enorm

  enorm = sqrt ( sum ( x(1:n) ** 2 ))

  return
end

function enorm2 ( n, x )

!*****************************************************************************80
!
!! ENORM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This routine was named ENORM.  It has been renamed "ENORM2",
!    and a simplified routine has been substituted.
!
!    The Euclidean norm is computed by accumulating the sum of
!    squares in three different sums.  The sums of squares for the
!    small and large components are scaled so that no overflows
!    occur.  Non-destructive underflows are permitted.  Underflows
!    and overflows do not occur in the computation of the unscaled
!    sum of squares for the intermediate components.
!
!    The definitions of small, intermediate and large components
!    depend on two constants, RDWARF and RGIANT.  The main
!    restrictions on these constants are that RDWARF^2 not
!    underflow and RGIANT^2 not overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM2, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) agiant
  real ( kind = 8 ) enorm2
  integer ( kind = 4 ) i
  real ( kind = 8 ) rdwarf
  real ( kind = 8 ) rgiant
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) s3
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xabs
  real ( kind = 8 ) x1max
  real ( kind = 8 ) x3max

  rdwarf = sqrt ( tiny ( rdwarf ) )
  rgiant = sqrt ( huge ( rgiant ) )

  s1 = 0.0D+00
  s2 = 0.0D+00
  s3 = 0.0D+00
  x1max = 0.0D+00
  x3max = 0.0D+00
  agiant = rgiant / real ( n, kind = 8 )

  do i = 1, n

    xabs = abs ( x(i) )

    if ( xabs <= rdwarf ) then

      if ( x3max < xabs ) then
        s3 = 1.0D+00 + s3 * ( x3max / xabs ) ** 2
        x3max = xabs
      else if ( xabs /= 0.0D+00 ) then
        s3 = s3 + ( xabs / x3max ) ** 2
      end if

    else if ( agiant <= xabs ) then

      if ( x1max < xabs ) then
        s1 = 1.0D+00 + s1 * ( x1max / xabs ) ** 2
        x1max = xabs
      else
        s1 = s1 + ( xabs / x1max ) ** 2
      end if

    else

      s2 = s2 + xabs ** 2

    end if

  end do
!
!  Calculation of norm.
!
  if ( s1 /= 0.0D+00 ) then

    enorm2 = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )

  else if ( s2 /= 0.0D+00 ) then

    if ( x3max <= s2 ) then
      enorm2 = sqrt ( s2 * ( 1.0D+00 + ( x3max / s2 ) * ( x3max * s3 ) ) )
    else
      enorm2 = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
    end if

  else

    enorm2 = x3max * sqrt ( s3 )

  end if

  return
end

subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )

!*****************************************************************************80
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  If, however, the limit of 10 iterations
!    is reached, then the output PAR will contain the best
!    value obtained so far.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the 
!    identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real ( kind = 8 ) DELTA, an upper bound on the euclidean norm
!    of D*X.  DELTA should be positive.
!
!    Input/output, real ( kind = 8 ) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real ( kind = 8 ) delta
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) dwarf
  real ( kind = 8 ) dxnorm
  real ( kind = 8 ) enorm
  real ( kind = 8 ) gnorm
  real ( kind = 8 ) fp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nsing
  real ( kind = 8 ) par
  real ( kind = 8 ) parc
  real ( kind = 8 ) parl
  real ( kind = 8 ) paru
  real ( kind = 8 ) qnorm
  real ( kind = 8 ) qtb(n)
  real ( kind = 8 ) r(ldr,n)
  real ( kind = 8 ) sdiag(n)
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) x(n)
!
!  DWARF is the smallest positive magnitude.
!
  dwarf = tiny ( dwarf )
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  do j = 1, n
    wa1(j) = qtb(j)
    if ( r(j,j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if
    if ( nsing < n ) then
      wa1(j) = 0.0D+00
    end if
  end do

  do k = 1, nsing
    j = nsing - k + 1
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, n
    l = ipvt(j)
    x(l) = wa1(j)
  end do
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm ( n, wa2 )
  fp = dxnorm - delta

  if ( fp <= 0.1D+00 * delta ) then
    if ( iter == 0 ) then
      par = 0.0D+00
    end if
    return
  end if
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
  parl = 0.0D+00

  if ( n <= nsing ) then

    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    end do

    temp = enorm ( n, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp

  end if
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
  do j = 1, n
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  end do

  gnorm = enorm ( n, wa1 )
  paru = gnorm / delta

  if ( paru == 0.0D+00 ) then
    paru = dwarf / min ( delta, 0.1D+00 )
  end if
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
  par = max ( par, parl )
  par = min ( par, paru )
  if ( par == 0.0D+00 ) then
    par = gnorm / dxnorm
  end if
!
!  Beginning of an iteration.
!
  do
 
    iter = iter + 1
!
!  Evaluate the function at the current value of PAR.
!
    if ( par == 0.0D+00 ) then
      par = max ( dwarf, 0.001D+00 * paru )
    end if

    wa1(1:n) = sqrt ( par ) * diag(1:n)

    call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

    wa2(1:n) = diag(1:n) * x(1:n)
    dxnorm = enorm ( n, wa2 )
    temp = fp
    fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
    if ( abs ( fp ) <= 0.1D+00 * delta ) then
      exit
    end if
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
    if ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) then
      exit
    else if ( iter == 10 ) then
      exit
    end if
!
!  Compute the Newton correction.
!
    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      wa1(j) = wa1(j) / sdiag(j)
      temp = wa1(j)
      wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
    end do

    temp = enorm ( n, wa1 )
    parc = ( ( fp / delta ) / temp ) / temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
    if ( 0.0D+00 < fp ) then
      parl = max ( parl, par )
    else if ( fp < 0.0D+00 ) then
      paru = min ( paru, par )
    end if
!
!  Compute an improved estimate for PAR.
!
    par = max ( parl, par + parc )
!
!  End of an iteration.
!
  end do
!
!  Termination.
!
  if ( iter == 0 ) then
    par = 0.0D+00
  end if

  return
end

subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )

!*****************************************************************************80
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P such 
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) cotan
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nsing
  real ( kind = 8 ) qtb(n)
  real ( kind = 8 ) qtbpj
  real ( kind = 8 ) r(ldr,n)
  real ( kind = 8 ) s
  real ( kind = 8 ) sdiag(n)
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa(n)
  real ( kind = 8 ) x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  do j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0D+00 ) then

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        if ( sdiag(k) /= 0.0D+00 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * t ** 2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k+1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n

    if ( sdiag(j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if

    if ( nsing < n ) then
      wa(j) = 0.0D+00
    end if

  end do

  do j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  end do
!
!  Permute the components of Z back to components of X.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  end do

  return
end
subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This subroutine uses Householder transformations with column
!    pivoting (optional) to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  The Householder transformation for
!    column K, K = 1,2,...,min(M,N), is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  The form of
!    this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK subroutine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q (the non-trivial elements of
!    the U vectors described above).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer ( kind = 4 ) IPVT(LIPVT), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity 
!    matrix.  If PIVOT is false, IPVT is not referenced.
!
!    Input, integer ( kind = 4 ) LIPVT, the dimension of IPVT, which should 
!    be N if pivoting is used.
!
!    Output, real ( kind = 8 ) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real ( kind = 8 ) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) lipvt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) acnorm(n)
  real ( kind = 8 ) ajnorm
  real ( kind = 8 ) enorm
  real ( kind = 8 ) epsmch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_temp
  integer ( kind = 4 ) ipvt(lipvt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) minmn
  logical pivot
  real ( kind = 8 ) r8_temp(m)
  real ( kind = 8 ) rdiag(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa(n)

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1:m,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if ( rdiag(kmax) < rdiag(k) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        r8_temp(1:m) = a(1:m,j)
        a(1:m,j)     = a(1:m,kmax)
        a(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, a(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( a(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j + 1, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0D+00 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp ** 2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) ) ** 2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = - ajnorm

  end do

  return
end

! Code from Alan Miller released into the public domain
! http://jblevins.org/mirror/amiller/ 2017-2-24


!PROGRAM t_median
!
!IMPLICIT NONE
!
!INTEGER           :: n
!DOUBLE PRECISION, ALLOCATABLE :: x(:)
!DOUBLE PRECISION              :: xmed
!
!INTERFACE
!  SUBROUTINE median(x, n, xmed)
!    IMPLICIT NONE
!    INTEGER, INTENT(IN)                :: n
!    DOUBLE PRECISION, INTENT(IN OUT), DIMENSION(:) :: x
!    DOUBLE PRECISION, INTENT(OUT)                  :: xmed
!  END SUBROUTINE median
!END INTERFACE
!
!DO
!  WRITE(*, *)'Enter n: '
!  READ(*, *) n
!  ALLOCATE( x(n) )
!  CALL RANDOM_NUMBER(x)
!  CALL median(x, n, xmed)
!  WRITE(*, 900) x(1), x(n), xmed
!  900 FORMAT(' First & last = ', 2f10.4, '    Median = ', f10.4/)
!  DEALLOCATE( x )
!END DO
!
!STOP
!END PROGRAM t_median



SUBROUTINE median(x, n, xmed)

! Find the median of X(1), ... , X(N), using as much of the quicksort
! algorithm as is needed to isolate it.
! N.B. On exit, the array X is partially ordered.

!     Latest revision - 26 November 1996
IMPLICIT NONE

INTEGER, INTENT(IN)                :: n
DOUBLE PRECISION, INTENT(IN OUT) :: x(n)
DOUBLE PRECISION, INTENT(OUT)                  :: xmed

! Local variables

DOUBLE PRECISION    :: temp, xhi, xlo, xmax, xmin
LOGICAL :: odd
INTEGER :: hi, lo, nby2, nby2p1, mid, i, j, k

nby2 = n / 2
nby2p1 = nby2 + 1
odd = .true.


!     HI & LO are position limits encompassing the median.

IF (n == 2 * nby2) odd = .false.
lo = 1
hi = n
IF (n < 3) THEN
  IF (n < 1) THEN
    xmed = 0.0
    RETURN
  END IF
  xmed = x(1)
  IF (n == 1) RETURN
  xmed = 0.5*(xmed + x(2))
  RETURN
END IF

!     Find median of 1st, middle & last values.

10 mid = (lo + hi)/2
xmed = x(mid)
xlo = x(lo)
xhi = x(hi)
IF (xhi < xlo) THEN          ! Swap xhi & xlo
  temp = xhi
  xhi = xlo
  xlo = temp
END IF
IF (xmed > xhi) THEN
  xmed = xhi
ELSE IF (xmed < xlo) THEN
  xmed = xlo
END IF

! The basic quicksort algorithm to move all values <= the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

i = lo
j = hi
50 DO
  IF (x(i) >= xmed) EXIT
  i = i + 1
END DO
DO
  IF (x(j) <= xmed) EXIT
  j = j - 1
END DO
IF (i < j) THEN
  temp = x(i)
  x(i) = x(j)
  x(j) = temp
  i = i + 1
  j = j - 1

!     Decide which half the median is in.

  IF (i <= j) GO TO 50
END IF

IF (.NOT. odd) THEN
  IF (j == nby2 .AND. i == nby2p1) GO TO 130
  IF (j < nby2) lo = i
  IF (i > nby2p1) hi = j
  IF (i /= j) GO TO 100
  IF (i == nby2) lo = nby2
  IF (j == nby2p1) hi = nby2p1
ELSE
  IF (j < nby2p1) lo = i
  IF (i > nby2p1) hi = j
  IF (i /= j) GO TO 100

! Test whether median has been isolated.

  IF (i == nby2p1) RETURN
END IF
100 IF (lo < hi - 1) GO TO 10

IF (.NOT. odd) THEN
  xmed = 0.5*(x(nby2) + x(nby2p1))
  RETURN
END IF
temp = x(lo)
IF (temp > x(hi)) THEN
  x(lo) = x(hi)
  x(hi) = temp
END IF
xmed = x(nby2p1)
RETURN

! Special case, N even, J = N/2 & I = J + 1, so the median is
! between the two halves of the series.   Find max. of the first
! half & min. of the second half, then average.

130 xmax = x(1)
DO k = lo, j
  xmax = MAX(xmax, x(k))
END DO
xmin = x(n)
DO k = i, hi
  xmin = MIN(xmin, x(k))
END DO
xmed = 0.5*(xmin + xmax)

RETURN
END SUBROUTINE median
