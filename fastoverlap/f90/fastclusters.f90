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

!***********************************************************************
! CLUSTERFASTOVERLAP MODULE
!***********************************************************************

! Subroutines:
!
!    ALIGN(COORDSB, COORDSA, NATOMS, DEBUG, L, KWIDTH, DISTANCE, DIST2, RMATBEST, NROTATIONS)
!        MAIN ALIGNMENT ALGORITHM ROUTINE
!        KWIDTH is the Gaussian Kernel width, this should probably be set to ~1/3 interatomic separation.
!        Performs alignment using SO(3) Coefficients calculated directly. 
!        Needs PERMGROUP, NPERMSIZE, NPERMGROUP, BESTPERM to be set and properly allocated
!   
!    ALIGNHARM(COORDSB, COORDSA, NATOMS, DEBUG, N, L, HWIDTH, KWIDTH, DISTANCE, DIST2, RMATBEST, NROTATIONS)
!        Performs alignment using SO(3) Coefficients calculated using Quantum Harmonic Oscillator Basis 
!        KWIDTH is the Gaussian Kernel width,  HWIDTH is the Quantum Harmonic Oscillator Basis length scale
!        These need to be carefully chosen along with N and L to ensure calculation is stable and accurate.
!        Needs PERMGROUP, NPERMSIZE, NPERMGROUP, BESTPERM to be set and properly allocated
! 
!    ALIGNCOEFFS(COORDSB,COORDSA,NATOMS,IMML,L,DEBUG,DISTANCE,DIST2,RMATBEST,NROTATIONS,ANGLES)
!        Primary alignment routine, called by ALIGN1
!        Needs PERMGROUP, NPERMSIZE, NPERMGROUP, BESTPERM to be set and properly allocated
!
!    HARMONIC0L(L, RJ, SIGMA, R0, RET)
!        Calculates the Harmonic integral when n=0
!
!    HARMONICNL(N,L,RJ,SIGMA,R0,RET)
!        Calculates Harmonic integral up to N,L
!        Note calculation unstable, so SIGMA must be > 10 RJ to get good results
!    
!    RYML(COORD, R, YML, L)
!        Calculates |COORD| and the Spherical Harmonic associated with COORD up to l
!    
!    HARMONICCOEFFS(COORDS, NATOMS, CNML, N, L, HWIDTH, KWIDTH)
!        Projects structure into Quantum Harmonic Oscillator Basis with scale HWIDTH and
!        Gaussian kernel width KWIDTH up to order N and angular moment degree L
!    
!    DOTHARMONICCOEFFS(C1NML, C2NML, N, L, IMML)
!        Calculates the SO(3) Fourier Coefficients of the overlap integral of two 
!        structures with coefficient arrays C1NML and C2NML
!    
!    FOURIERCOEFFS(COORDSB, COORDSA, NATOMS, L, KWIDTH, IMML, YMLB, YMLA)
!        Calculates the SO(3) Fourier Coefficients of the overlap integral of two 
!        structures directly by calculating the coefficients of the NATOMS**2
!        Gaussian overlap functions.
!    
!    CALCOVERLAP(IMML, OVERLAP, L, ILMM)
!        Calculates the overlap integral array from SO(3) Fourier Coefficients IMML
!        Also returns ILMM, the transposed and rolled version of IMML used by DSOFT
!    
!    FINDROTATIONS(OVERLAP, L, ANGLES, AMPLITUDES, NROTATIONS, DEBUG)
!        Finds the maximum overlap Euler angles of an overlap integral array
!    
!    EULERM(A,B,G,ROTM)
!        Calculates rotation matrix, ROTM, corresponding to  Euler angles, a,b,g
!    
!    EULERINVM(A,B,G,ROTM)
!        Calculates transpose/inverse of rotation matrix corresponding to Euler angles, a,b,g
!    
!    SETCLUSTER()
!        Used to set keywords if they're not set already
!    
!    CHECKKEYWORDS()
!        Sanity checks for the keywords

!***********************************************************************

! EXTERNAL SUBROUTINES
!    MINPERMDIST (minpermdist.f90) depends on (bulkmindist.f90,minperm.f90,newmindist.f90,orient.f90)
!    XDNRMP (legendre.f90)
!        Needed to calculate Legendre polynomials

!***********************************************************************

! EXTERNAL MODULES
!    COMMONS (commons.f90)
!    FASTOVERLAPUTILS (fastutils.f90) depends on (minperm.f90)
!        Helper Module Needed for Peak Fitting and FFT routines
!    DSOFT (DSOFT.f90) 
!        Module for performing discrete SO(3) transforms, depends on fftw.

!***********************************************************************

INCLUDE "commons.f90"
INCLUDE "fastutils.f90"

! Module for performing discrete SO(3) transforms, depends on fftw.
INCLUDE "DSOFT.f90"

MODULE CLUSTERFASTOVERLAP

USE COMMONS, ONLY : PERMGROUP, NPERMSIZE, NPERMGROUP, NATOMS, BESTPERM, MYUNIT
USE FASTOVERLAPUTILS, ONLY : DUMMYA, DUMMYB, XBESTA, XBESTASAVE

LOGICAL, SAVE :: PERMINVOPTSAVE, NOINVERSIONSAVE

DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0

CONTAINS

SUBROUTINE ALIGN(COORDSB, COORDSA, NATOMS, DEBUG, L, KWIDTH, DISTANCE, DIST2, RMATBEST, NROTATIONS)

!  COORDSA becomes the optimal alignment of the optimal permutation(-inversion)
!  isomer. DISTANCE is the residual square distance for the best alignment with 
!  respect to permutation(-inversion)s as well as orientation and centre of mass.
!  COORDSA and COORDSB are both centred on the ORIGIN

!  KWIDTH is the width of the Gaussian kernels that are centered on each of the
!  atomic coordinates, whose overlap integral is maximised to find the optimal
!  rotations

!  RMATBEST gives the optimal rotation matrix

!  L is the maximum angular momentum degree up to which the SO(3) coefficients 
!  are calculated number of coefficients that will be calculated = 1/3 (L+1)(2L+1)(2L+3)

!  Number of Calculations for SO(3) calculations ~ O(1/3 (L+1)(2L+1)(2L+3) * NATOMS**2)

USE COMMONS, ONLY: BESTPERM, PERMOPT, PERMINVOPT, NOINVERSION, CHRMMT, AMBERT, AMBER12T
IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, L
INTEGER, INTENT(IN) :: NROTATIONS
LOGICAL, INTENT(IN) :: DEBUG
DOUBLE PRECISION, INTENT(IN) :: KWIDTH ! Gaussian Kernel width
DOUBLE PRECISION, INTENT(INOUT) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, DIST2, RMATBEST(3,3)

COMPLEX*16 PIMML(-L:L,-L:L,0:L)
COMPLEX*16 IMML(-L:L,-L:L,0:L), YMLA(-L:L,0:L,NATOMS), YMLB(-L:L,0:L,NATOMS)

DOUBLE PRECISION SAVEA(3*NATOMS),SAVEB(3*NATOMS),COMA(3),COMB(3)
DOUBLE PRECISION ANGLES(NROTATIONS,3), DISTSAVE, RMATSAVE(3,3), WORSTRAD, DIST2SAVE
INTEGER J,J1,J2,M1,M2,IND2,NROT,NDUMMY,INVERT,PATOMS
INTEGER SAVEPERM(NATOMS), KEEPPERM(NATOMS)

! Checking keywords are set properly
CALL CHECKKEYWORDS()

! Setting keywords for fastoverlap use of minpermdist, will be reset when exiting program
PERMINVOPTSAVE = PERMINVOPT
NOINVERSIONSAVE = NOINVERSION
PERMINVOPT = .FALSE.
NOINVERSION = .TRUE.

! Centering COORDSA and COORDSB on the origin
COMA = 0.D0
COMB = 0.D0
DO J=1,NATOMS
    COMA = COMA + COORDSA(3*J-2:3*J)
    COMB = COMB + COORDSB(3*J-2:3*J)
ENDDO
COMA = COMA/NATOMS
COMB = COMB/NATOMS
DO J=1,NATOMS
    COORDSA(3*J-2:3*J) = COORDSA(3*J-2:3*J) - COMA
    COORDSB(3*J-2:3*J) = COORDSB(3*J-2:3*J) - COMB
ENDDO


! Calculating overlap integral separately for each permutation group
IMML = CMPLX(0.D0,0.D0,8)
NDUMMY=1
DO J1=1,NPERMGROUP
    PATOMS=INT(NPERMSIZE(J1),4)
    DO J2=1,PATOMS
        IND2 = PERMGROUP(NDUMMY+J2-1)
        SAVEA(3*J2-2:3*J2)=COORDSA(3*IND2-2:3*IND2)
        SAVEB(3*J2-2:3*J2)=COORDSB(3*IND2-2:3*IND2)
    ENDDO
    CALL FOURIERCOEFFS(SAVEB,SAVEA,PATOMS,L,KWIDTH,PIMML,YMLB,YMLA)
    DO J=0,L
        DO M2=-J,J
            DO M1=-J,J
            IMML(M1,M2,J) = IMML(M1,M2,J) + PIMML(M1,M2,J)
            ENDDO
        ENDDO
    ENDDO
    NDUMMY=NDUMMY+NPERMSIZE(J1)
ENDDO

SAVEA(1:3*NATOMS) = COORDSA(1:3*NATOMS)
SAVEB(1:3*NATOMS) = COORDSB(1:3*NATOMS)

NROT = NROTATIONS
CALL ALIGNCOEFFS(SAVEB,SAVEA,NATOMS,IMML,L,DEBUG,DISTSAVE,DIST2SAVE,RMATSAVE,NROT,ANGLES)

IF (PERMINVOPTSAVE.AND.(.NOT.(CHRMMT.OR.AMBERT.OR.AMBER12T))) THEN 
    IF (DEBUG) WRITE(MYUNIT,'(A)') 'fastoverlap> inverting geometry for comparison with target'
    ! Saving non inverted configuration
    XBESTASAVE(1:3*NATOMS) = SAVEA(1:3*NATOMS)

    ! Calculating overlap integral for inverted configuration
    NDUMMY=1
    DO J1=1,NPERMGROUP
        PATOMS=INT(NPERMSIZE(J1),4)
        DO J2=1,PATOMS
            IND2 = PERMGROUP(NDUMMY+J2-1)
            SAVEA(3*J2-2:3*J2)=-COORDSA(3*IND2-2:3*IND2)
            SAVEB(3*J2-2:3*J2)=COORDSB(3*IND2-2:3*IND2)
        ENDDO
        CALL FOURIERCOEFFS(SAVEB,SAVEA,PATOMS,L,KWIDTH,PIMML,YMLB,YMLA)
        DO J=0,L
            DO M2=-J,J
                DO M1=-J,J
                    IMML(M1,M2,J) = IMML(M1,M2,J) + PIMML(M1,M2,J)
                ENDDO
            ENDDO
        ENDDO
        NDUMMY=NDUMMY+NPERMSIZE(J1)
    ENDDO
    SAVEA(1:3*NATOMS) = -COORDSA(1:3*NATOMS)
    SAVEB(1:3*NATOMS) = COORDSB(1:3*NATOMS)

    NROT = NROTATIONS
    CALL ALIGNCOEFFS(SAVEB,SAVEA,NATOMS,IMML,L,DEBUG,DISTANCE,DIST2,RMATBEST,NROT,ANGLES)
    IF (DISTANCE.LT.DISTSAVE) THEN
        IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') &
    &   'fastoverlap> inversion found better alignment, distance=', distance
        COORDSA(1:3*NATOMS) = SAVEA(1:3*NATOMS)
    ELSE
        COORDSA(1:3*NATOMS) = XBESTASAVE(1:3*NATOMS)
        DISTANCE = DISTSAVE
        DIST2 = DIST2SAVE
        RMATBEST = RMATSAVE
        IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') &
    &   'fastoverlap> better alignment with no-inversion, distance=', distance
    ENDIF
ELSE
    IF (DEBUG) WRITE(MYUNIT,'(A)') 'fastoverlap> not inverting geometry for comparison with target'
    COORDSA(1:3*NATOMS) = SAVEA(1:3*NATOMS)
    DISTANCE = DISTSAVE
    DIST2 = DIST2SAVE
    RMATBEST = RMATSAVE
ENDIF

IF (DEBUG) THEN
    WRITE(MYUNIT,'(A,G20.10)') 'fastoverlap> overall best distance=', distance
    WRITE(MYUNIT,'(A)') 'fastoverlap> overall best rotation matrix:'
    WRITE(MYUNIT, '(3F20.10)') RMATBEST(1:3,1:3)
ENDIF

PERMINVOPT = PERMINVOPTSAVE
NOINVERSION = NOINVERSIONSAVE

END SUBROUTINE ALIGN

SUBROUTINE ALIGNHARM(COORDSB, COORDSA, NATOMS, DEBUG, N, L, HWIDTH, KWIDTH, DISTANCE, DIST2, RMATBEST, NROTATIONS)
!  COORDSA becomes the optimal alignment of the optimal permutation(-inversion)
!  isomer. DISTANCE is the residual square distance for the best alignment with 
!  respect to permutation(-inversion)s as well as orientation and centre of mass.
!  COORDSA and COORDSB are both centred on the ORIGIN

!  RMATBEST gives the optimal rotation matrix

!  KWIDTH is the width of the Gaussian kernels that are centered on each of the
!  atomic coordinates, whose overlap integral is maximised to find the optimal
!  rotations
!  L is the maximum angular momentum degree up to which the SO(3) coefficients 
!  are calculated number of coefficients that will be calculated = 1/3 (L+1)(2L+1)(2L+3)

!  HWIDTH is the lengthscale of the Quantum Harmonic Oscillator Basis
!  N is the maximum order of the Quantum Harmonic Oscillator basis

!  Number of Calculations for SO(3) calculations ~ O(1/3 (L+1)(2L+1)(2L+3) * NATOMS**2)

USE COMMONS, ONLY: BESTPERM, PERMOPT, PERMINVOPT, NOINVERSION, CHRMMT, AMBERT, AMBER12T
IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, N, L
INTEGER, INTENT(IN) :: NROTATIONS
LOGICAL, INTENT(IN) :: DEBUG
DOUBLE PRECISION, INTENT(IN) :: HWIDTH, KWIDTH
DOUBLE PRECISION, INTENT(INOUT) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, DIST2, RMATBEST(3,3)

COMPLEX*16 PIMML(-L:L,-L:L,0:L)
COMPLEX*16 IMML(-L:L,-L:L,0:L), YMLA(-L:L,0:L,NATOMS), YMLB(-L:L,0:L,NATOMS)
COMPLEX*16 COEFFSA(0:N,-L:L,0:L,NPERMGROUP), COEFFSB(0:N,-L:L,0:L,NPERMGROUP)

DOUBLE PRECISION SAVEA(3*NATOMS),SAVEB(3*NATOMS)
DOUBLE PRECISION ANGLES(NROTATIONS,3), DISTSAVE, RMATSAVE(3,3), WORSTRAD, DIST2SAVE
INTEGER J,J1,J2,M1,M2,IND2,NROT,NDUMMY,INVERT,PATOMS
INTEGER SAVEPERM(NATOMS), KEEPPERM(NATOMS)


! Checking keywords are set properly
CALL CHECKKEYWORDS()

! Setting keywords for fastoverlap use of minpermdist, will be reset when exiting program
PERMINVOPTSAVE = PERMINVOPT
NOINVERSIONSAVE = NOINVERSION
PERMINVOPT = .FALSE.
NOINVERSION = .TRUE.

! Calculating overlap integral separately for each permutation group
IMML = CMPLX(0.D0,0.D0,8)
NDUMMY=1
DO J1=1,NPERMGROUP
    PATOMS=INT(NPERMSIZE(J1),4)
    DO J2=1,PATOMS
        IND2 = PERMGROUP(NDUMMY+J2-1)
        SAVEA(3*J2-2:3*J2)=COORDSA(3*IND2-2:3*IND2)
        SAVEB(3*J2-2:3*J2)=COORDSB(3*IND2-2:3*IND2)
    ENDDO
    CALL HARMONICCOEFFS(SAVEA, PATOMS, COEFFSA(:,:,:,J1), N, L, HWIDTH, KWIDTH)
    CALL HARMONICCOEFFS(SAVEB, PATOMS, COEFFSB(:,:,:,J1), N, L, HWIDTH, KWIDTH)
    CALL DOTHARMONICCOEFFS(COEFFSB(:,:,:,J1), COEFFSA(:,:,:,J1), N, L, PIMML)
    DO J=0,L
        DO M2=-J,J
            DO M1=-J,J
            IMML(M1,M2,J) = IMML(M1,M2,J) + PIMML(M1,M2,J)
            ENDDO
        ENDDO
    ENDDO
    NDUMMY=NDUMMY+NPERMSIZE(J1)
ENDDO

NROT = NROTATIONS
CALL ALIGNCOEFFS(SAVEB,SAVEA,NATOMS,IMML,L,DEBUG,DISTSAVE,DIST2SAVE,RMATSAVE,NROT,ANGLES)

IF (PERMINVOPTSAVE.AND.(.NOT.(CHRMMT.OR.AMBERT.OR.AMBER12T))) THEN 
    IF (DEBUG) WRITE(MYUNIT,'(A)') 'fastoverlap> inverting geometry for comparison with target'
    ! Saving non inverted configuration
    XBESTASAVE(1:3*NATOMS) = SAVEA(1:3*NATOMS)
    KEEPPERM(1:NATOMS) = BESTPERM(1:NATOMS)
    SAVEA = -COORDSA(1:3*NATOMS)
    NROT = NROTATIONS

    ! Recalculating Fourier Coefficients for inverted COORDSA
    IMML = CMPLX(0.D0,0.D0,8)
    NDUMMY=1
    DO J1=1,NPERMGROUP
        DO J=0,L
            COEFFSA(:,:,J,J1) = COEFFSA(:,:,J,J1) * (-1)**(J)
        ENDDO
        CALL DOTHARMONICCOEFFS(COEFFSB(:,:,:,J1), COEFFSA(:,:,:,J1), N, L, PIMML)
        DO J=0,L
            DO M2=-J,J
                DO M1=-J,J
                IMML(M1,M2,J) = IMML(M1,M2,J) + PIMML(M1,M2,J)
                ENDDO
            ENDDO
        ENDDO
        NDUMMY=NDUMMY+NPERMSIZE(J1)
    ENDDO
    CALL ALIGNCOEFFS(SAVEB,SAVEA,NATOMS,IMML,L,DEBUG,DISTANCE,DIST2,RMATBEST,NROT,ANGLES)
    
    IF (DISTANCE.LT.DISTSAVE) THEN
        IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') &
    &   'fastoverlap> inversion found better alignment, distance=', distance
        COORDSA(1:3*NATOMS) = SAVEA(1:3*NATOMS)
        RMATBEST = RMATSAVE
    ELSE
        COORDSA(1:3*NATOMS) = XBESTASAVE(1:3*NATOMS)
        DISTANCE = DISTSAVE
        DIST2 = DIST2SAVE
        RMATBEST = RMATSAVE
    ENDIF
ELSE
    IF (DEBUG) WRITE(MYUNIT,'(A)') 'fastoverlap> not inverting geometry for comparison with target'
    COORDSA(1:3*NATOMS) = SAVEA(1:3*NATOMS)
    DISTANCE = DISTSAVE
    DIST2 = DIST2SAVE
    RMATBEST = RMATSAVE
ENDIF

IF (DEBUG) THEN
    WRITE(MYUNIT,'(A,G20.10)') 'fastoverlap> overall best distance=', distance
    WRITE(MYUNIT,'(A)') 'fastoverlap> overall best rotation matrix:'
    WRITE(MYUNIT, '(3F20.10)') RMATBEST(1:3,1:3)
ENDIF

PERMINVOPT = PERMINVOPTSAVE
NOINVERSION = NOINVERSIONSAVE

END SUBROUTINE ALIGNHARM

SUBROUTINE ALIGNCOEFFS(COORDSB,COORDSA,NATOMS,IMML,L,DEBUG,DISTANCE,DIST2,RMATBEST,NROTATIONS,ANGLES)
! Aligns two structures, specified by COORDSA and COORDSB, aligns COORDSA so it most
! closely matches COORDSB. 
! Assumes that COORDSA and COORDSB are both centered on their Centers of Mass
! Uses precalculated Fourier Coefficients, IMML
! Uses minpermdist to refine alignment

! Low-level routine, better to use ALIGN or ALIGNHARM

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, L
INTEGER, INTENT(INOUT) :: NROTATIONS
LOGICAL, INTENT(IN) :: DEBUG
DOUBLE PRECISION, INTENT(INOUT) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)
DOUBLE PRECISION, INTENT(OUT) :: ANGLES(NROTATIONS,3)
DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, DIST2, RMATBEST(3,3)
COMPLEX*16, INTENT(IN) :: IMML(-L:L,-L:L,0:L)

COMPLEX*16 ILMM(0:L,0:2*L,0:2*L)
DOUBLE PRECISION OVERLAP(2*L+2,2*L+2,2*L+2)
DOUBLE PRECISION AMPLITUDES(NROTATIONS), BESTDIST, RMATSAVE(3,3), RMAT(3,3), WORSTRAD
INTEGER J, J1


CALL CALCOVERLAP(IMML, OVERLAP, L, ILMM)
CALL FINDROTATIONS(OVERLAP, L, ANGLES, AMPLITUDES, NROTATIONS, DEBUG)
IF (DEBUG) WRITE(MYUNIT,'(A,I3,A)') 'fastoverlap> found ', NROTATIONS, ' candidate rotations'


BESTDIST = HUGE(BESTDIST)
DUMMYB(:) = COORDSB(:3*NATOMS)

DO J=1,NROTATIONS

    CALL EULERM(ANGLES(J,1),ANGLES(J,2),ANGLES(J,3),RMATSAVE)
    DO J1=1,NATOMS
        DUMMYA(J1*3-2:J1*3) = MATMUL(RMATSAVE, COORDSA(J1*3-2:J1*3))
    ENDDO

    IF (DEBUG) THEN
        WRITE(MYUNIT,'(A,I3,A)') 'fastoverlap> testing rotation', J, ' with Euler angles:'
        WRITE(MYUNIT, '(3F20.10)') ANGLES(J,:)
        WRITE(MYUNIT,'(A)') 'fastoverlap> testing rotation matrix:'
        WRITE(MYUNIT, '(3F20.10)') RMATSAVE(1:3,1:3)
    ENDIF

    CALL MINPERMDIST(DUMMYB,DUMMYA,NATOMS,DEBUG,0.D0,0.D0,0.D0,.FALSE.,.FALSE.,DISTANCE,DIST2,.FALSE.,RMAT)
    IF (DISTANCE.LT.BESTDIST) THEN
        BESTDIST = DISTANCE
        XBESTA(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
        RMATBEST = MATMUL(RMAT,RMATSAVE)

        IF (DEBUG) THEN
            WRITE(MYUNIT,'(A,G20.10)') 'fastoverlap> new best alignment distance=', BESTDIST
            WRITE(MYUNIT,'(A)') 'fastoverlap> new best rotation matrix:'
            WRITE(MYUNIT, '(3F20.10)') RMATBEST(1:3,1:3)
        END IF

    ELSE IF (DEBUG) THEN
        WRITE(MYUNIT,'(A,G20.10)') 'fastoverlap> best aligment distance found=', BESTDIST
        WRITE(MYUNIT,'(A)') 'fastoverlap> best rotation matrix found:'
        WRITE(MYUNIT, '(3F20.10)') RMATBEST(1:3,1:3)
    ENDIF
ENDDO


! Returning Best Coordinates
COORDSA(1:3*NATOMS) = XBESTA(1:3*NATOMS)

DISTANCE = BESTDIST
DIST2 = BESTDIST**2

END SUBROUTINE ALIGNCOEFFS

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
SQRTI = 1.D0
DO J=0,L-2
    RET(I,J) = (SQRT(I+J+0.5D0)*RET(I-1,J) - (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1) -&
        SQRT(I+J+1.5D0) * RET(I-1,J+2))/SQRTI
ENDDO

DO I=2,N
    SQRTI = SQRT(REAL(I,8))
    DO J=0,L-2*I
    RET(I,J) = (SQRT(I+J+0.5D0)*RET(I-1,J) - (2.D0*J+3.D0)*SIGMA**2/RJ/R0 * RET(I-1,J+1) -&
        SQRT(I+J+1.5D0) * RET(I-1,J+2) + SQRT(I-1.D0) * RET(I-2,J+2))/SQRTI
    ENDDO
ENDDO

END SUBROUTINE HARMONICNL

SUBROUTINE RYML2(COORD, R, YML, L)

! Calculates the Spherical Harmonics associated with coordinate COORD
! up to L, returns R, the distance COORD is from origin
! Calculates value of Legendre Polynomial Recursively

! UNSTABLE WHEN Z CLOSE TO 0 OR 1

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: COORD(3)
INTEGER, INTENT(IN) :: L
DOUBLE PRECISION, INTENT(OUT) :: R
COMPLEX*16, INTENT(OUT) :: YML(-L:L,0:L)

INTEGER J, M, INDM1, INDM0, INDM2
DOUBLE PRECISION THETA, PHI, Z, FACTORIALS(0:2*L), SQRTZ, SQRTMJ
COMPLEX*16 EXPIM(-L:L)

R = (COORD(1)**2+COORD(2)**2+COORD(3)**2)**0.5
PHI = ATAN2(COORD(2), COORD(1))
Z = COORD(3)/R
SQRTZ = SQRT(1.D0-Z**2)

!Calculating Associate Legendre Function
YML = CMPLX(0.D0,0.D0, 8)
YML(0,0) = (4*PI)**-0.5

! Initialising Recurrence for Associated Legendre Polynomials
! Calculating normalised Legendre Polynomials for better numerical stability
! Pnorm^m_l = \sqrt{(l-m)!/(l+m)!} P^m_l
DO J=0, L-1
    YML(J+1,J+1) = - SQRT((2.D0*J+1.D0)/(2.D0*J+2.D0)) * SQRTZ* YML(J,J)
    ! Calculating first recurrence term
    YML(J, J+1) = -SQRT(2.D0*(J+1))*Z/SQRTZ * YML(J+1, J+1)
ENDDO

! Recurrence for normalised Associated Legendre Polynomials
DO J=1,L
    DO M=J-1,-J+1,-1
        SQRTMJ = SQRT((J+M)*(J-M+1.D0))
        YML(M-1, J) = -2*M*Z/SQRTMJ/SQRTZ * YML(M, J) - SQRT((J-M)*(J+M+1.D0))/SQRTMJ * YML(M+1,J)
    ENDDO
ENDDO

! Calculating exp(imPHI) component
DO M=-L,L
    EXPIM(M) = EXP(CMPLX(0.D0, M*PHI, 8))
ENDDO

! Calculate Spherical Harmonics
DO J=1,L
    DO M=-J,J
        INDM0 = MODULO(M, 2*L+1)
        YML(M,J) = EXPIM(M)*YML(M,J) * SQRT((2.D0*J+1.D0))
    ENDDO
ENDDO

END SUBROUTINE RYML2

SUBROUTINE RYML(COORD, R, YML, L)

! Calculates the Spherical Harmonics associated with coordinate COORD
! up to L, returns R, the distance COORD is from origin
! Calculates value of Legendre Polynomial Recursively

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: COORD(3)
INTEGER, INTENT(IN) :: L
DOUBLE PRECISION, INTENT(OUT) :: R
COMPLEX*16, INTENT(OUT) :: YML(-L:L,0:L)

INTEGER J, M, INDM1, INDM0, INDM2, ISIG
DOUBLE PRECISION THETA, PHI, Z, FACTORIALS(0:2*L), SQRTZ, SQRTMJ, PLM(0:L), IPN(0:L), FACT
COMPLEX*16 EXPIM(-L:L)

R = (COORD(1)**2+COORD(2)**2+COORD(3)**2)**0.5
PHI = ATAN2(COORD(2), COORD(1))
Z = COORD(3)/R

!Calculating Associate Legendre Function
YML = CMPLX(0.D0,0.D0, 8)
YML(0,0) = (4*PI)**-0.5

FACT = (2*PI)**-0.5

DO J=0, L
    ! Calculate Normalised Legendre Polynomial
    CALL XDNRMP(J,0,J,Z,1,PLM(0:J),IPN(0:J),ISIG)
    YML(0:J,J) = PLM(0:J) * FACT
    DO M=1,J
        YML(-M,J) = YML(M,J)
        YML(M,J) = YML(-M,J) * (-1)**M
    ENDDO
ENDDO

! Calculating exp(imPHI) component
DO M=-L,L
    EXPIM(M) = EXP(CMPLX(0.D0, M*PHI, 8))
ENDDO

! Calculate Spherical Harmonics
DO J=1,L
    DO M=-J,J
        INDM0 = MODULO(M, 2*L+1)
        YML(M,J) = EXPIM(M)*YML(M,J) !* SQRT((2.D0*J+1.D0))
    ENDDO
ENDDO

END SUBROUTINE RYML

SUBROUTINE HARMONICCOEFFS(COORDS, NATOMS, CNML, N, L, HWIDTH, KWIDTH)

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

END SUBROUTINE HARMONICCOEFFS

SUBROUTINE HARMONICCOEFFSPERM(COORDS, NATOMS, CNML, N, L, HWIDTH, KWIDTH, NPERMGROUP)

!
! For a set of Gaussian Kernels of width KWIDTH at COORDS, 
! this will calculate the coefficients of the isotropic quantum harmonic basis
! cnlm with length scale HWIDTH up to N and L.
! Returns coefficients of the different permutations groups
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, N, L, NPERMGROUP
DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS), HWIDTH, KWIDTH
COMPLEX*16, INTENT(OUT) :: CNML(0:N,-L:L,0:L,1:NPERMGROUP)

DOUBLE PRECISION DUMMY(3*NATOMS)
INTEGER J1, J2, IND2, NDUMMY, PATOMS

! Calculating overlap integral separately for each permutation group
NDUMMY=1
DO J1=1,NPERMGROUP
    PATOMS=NPERMSIZE(J1)
    DO J2=1,PATOMS
        IND2 = PERMGROUP(NDUMMY+J2-1)
        DUMMY(3*J2-2:3*J2)=COORDS(3*IND2-2:3*IND2)
    ENDDO
    CALL HARMONICCOEFFS(DUMMY, PATOMS, CNML(:,:,:,J1), N, L, HWIDTH, KWIDTH)
    NDUMMY=NDUMMY+PATOMS
ENDDO

END SUBROUTINE HARMONICCOEFFSPERM

SUBROUTINE HARMONICCOEFFSMULTI(COORDSLIST,NATOMS,NLIST,CNMLLIST,N,L,HWIDTH,KWIDTH,NPERMGROUP)

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, NLIST, N, L, NPERMGROUP
DOUBLE PRECISION, INTENT(IN) :: COORDSLIST(3*NATOMS, NLIST), HWIDTH, KWIDTH
COMPLEX*16, INTENT(OUT) :: CNMLLIST(0:N,-L:L,0:L,1:NPERMGROUP, NLIST)

INTEGER I

!write(*,*) NATOMS, NLIST, N, L, NPERMGROUP
!WRITE(*,*) SHAPE(CNMLLIST), SHAPE(COORDSLIST)

DO I=1,NLIST
    CALL HARMONICCOEFFSPERM(COORDSLIST(:,I),NATOMS,CNMLLIST(:,:,:,:,I),N,L,HWIDTH,KWIDTH,NPERMGROUP)
ENDDO

END SUBROUTINE HARMONICCOEFFSMULTI

SUBROUTINE DOTHARMONICCOEFFS(C1NML, C2NML, N, L, IMML)

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

END SUBROUTINE DOTHARMONICCOEFFS

SUBROUTINE DOTHARMONICCOEFFSPERM(C1NML, C2NML, N, L, IMML, NPERMGROUP)

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, L, NPERMGROUP
COMPLEX*16, INTENT(IN) :: C1NML(0:N,-L:L,0:L,NPERMGROUP), C2NML(0:N,-L:L,0:L,NPERMGROUP)
COMPLEX*16, INTENT(OUT) :: IMML(-L:L,-L:L,0:L)

INTEGER I, J, M1, M2, K, INDM1, INDM2

IMML = CMPLX(0.D0,0.D0,8)

DO K=1,NPERMGROUP
    DO J=0,L
        DO M2=-J,J
            DO M1=-J,J
                DO I=0,N
                    IMML(M1,M2,J) = IMML(M1,M2,J) + CONJG(C1NML(I,M1,J,K))*C2NML(I,M2,J,K)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE DOTHARMONICCOEFFSPERM

SUBROUTINE CALCSIMILARITY(C1NML, C2NML, N, L, NPERMGROUP, NORM, MAXOVER)

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, L, NPERMGROUP
COMPLEX*16, INTENT(IN) :: C1NML(0:N,-L:L,0:L,NPERMGROUP), C2NML(0:N,-L:L,0:L,NPERMGROUP)
DOUBLE PRECISION, INTENT(OUT) :: NORM, MAXOVER

COMPLEX*16 IMML(-L:L,-L:L,0:L), ILMM(0:L,0:2*L,0:2*L)
DOUBLE PRECISION OVERLAP(2*L+2,2*L+2,2*L+2)

INTEGER J,M1,M2

CALL DOTHARMONICCOEFFSPERM(C1NML, C2NML, N, L, IMML, NPERMGROUP)

! Calculated average overlap
DO J=0,L
    DO M2=-J,J
        DO M1=-J,J
            NORM = NORM + REAL(IMML(M1,M2,J),8)**2 + AIMAG(IMML(M1,M2,J))**2
        ENDDO
    ENDDO
ENDDO

! Calculate max overlap
CALL CALCOVERLAP(IMML, OVERLAP, L, ILMM)
MAXOVER = MAXVAL(OVERLAP)

END SUBROUTINE CALCSIMILARITY

SUBROUTINE CALCSIMILARITIES(C1NMLLIST,N1LIST,C2NMLLIST,N2LIST,N,L,NPERMGROUP,NORMS,MAXOVERS,SYM)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N1LIST, N2LIST, N, L, NPERMGROUP
COMPLEX*16, INTENT(IN) :: C1NMLLIST(0:N,-L:L,0:L,NPERMGROUP,N1LIST), &
    & C2NMLLIST(0:N,-L:L,0:L,NPERMGROUP,N2LIST)
LOGICAL, INTENT(IN) :: SYM
DOUBLE PRECISION, INTENT(OUT) :: NORMS(N1LIST,N2LIST), MAXOVERS(N1LIST,N2LIST)

INTEGER I1, I2

IF (SYM) THEN
    ! if C1NMLLIST == C2NMLLIST then only need to calculate half the values
    DO I1=1,N1LIST
        DO I2=I1,N1LIST
            CALL CALCSIMILARITY(C1NMLLIST(:,:,:,:,I1), C2NMLLIST(:,:,:,:,I2), N, L, NPERMGROUP, &
                & NORMS(I1,I2), MAXOVERS(I1,I2))
            NORMS(I2,I1) = NORMS(I1,I2)
            MAXOVERS(I2,I1) = MAXOVERS(I1,I2)
        ENDDO
    ENDDO
ELSE
    ! Calculate all values
    DO I1=1,N1LIST
        DO I2=1,N1LIST
            CALL CALCSIMILARITY(C1NMLLIST(:,:,:,:,I1), C2NMLLIST(:,:,:,:,I2), N, L, NPERMGROUP, &
                & NORMS(I1,I2), MAXOVERS(I1,I2))
        ENDDO
    ENDDO
ENDIF

END SUBROUTINE CALCSIMILARITIES

SUBROUTINE CALCOVERLAPMATRICES(COORDSLIST,NATOMS,NLIST,N,L,HWIDTH,KWIDTH,NORMS,MAXOVERS)

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, NLIST, N, L
DOUBLE PRECISION, INTENT(IN) :: COORDSLIST(3*NATOMS, NLIST), HWIDTH, KWIDTH
DOUBLE PRECISION, INTENT(OUT) :: NORMS(NLIST,NLIST), MAXOVERS(NLIST,NLIST)

COMPLEX*16 CNMLLIST(0:N,-L:L,0:L,1:NPERMGROUP, NLIST)

CALL HARMONICCOEFFSMULTI(COORDSLIST,NATOMS,NLIST,CNMLLIST,N,L,HWIDTH,KWIDTH,NPERMGROUP)
CALL CALCSIMILARITIES(CNMLLIST,NLIST,CNMLLIST,NLIST,N,L,NPERMGROUP,NORMS,MAXOVERS,.TRUE.)

END SUBROUTINE CALCOVERLAPMATRICES

SUBROUTINE FOURIERCOEFFS(COORDSB, COORDSA, NATOMS, L, KWIDTH, IMML, YMLB, YMLA)
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
DOUBLE PRECISION RA(NATOMS), RB(NATOMS), IL(0:L), R1R2, EXPRA(NATOMS), EXPRB(NATOMS), FACT, TMP

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

FACT = 4.D0 * PI**2.5 * KWIDTH**3

IMML = CMPLX(0.D0,0.D0,8)
DO IA=1,NATOMS
    DO IB=1,NATOMS
        ! Don't calculate cross terms for points separated by 4 kwidths to speed up calculation
        IF (ABS(RA(IA)-RB(IB)).LT.(4*KWIDTH)) THEN
            R1R2 = 0.5D0 * RA(IA)*RB(IB)/KWIDTH**2
            CALL SPHI(L, R1R2, K, IL)
            TMP = FACT*EXPRA(IA)*EXPRB(IB)!*SQRT(PI/2/R1R2)
            DO J=0,L
                DO M2=-L,L
                    DO M1=-L,L
                        IMML(M1,M2,J) = IMML(M1,M2,J) + IL(J)*YMLB(M1,J,IB)*CONJG(YMLA(M2,J,IA))*TMP
                    ENDDO
                ENDDO
            ENDDO
        END IF
    ENDDO
ENDDO

END SUBROUTINE FOURIERCOEFFS

SUBROUTINE CALCOVERLAP(IMML, OVERLAP, L, ILMM)
! Converts an array of SO(3) Fourier Coefficients to a discrete
! overlap array using a fast discrete SO(3) Fourier Transform (DSOFT)

USE DSOFT, ONLY : ISOFT

IMPLICIT NONE
INTEGER, INTENT(IN) :: L
COMPLEX*16, INTENT(IN) :: IMML(-L:L,-L:L,0:L)
DOUBLE PRECISION, INTENT(OUT) :: OVERLAP(2*L+2,2*L+2,2*L+2)

COMPLEX*16, INTENT(OUT) :: ILMM(0:L,0:2*L,0:2*L)
COMPLEX*16 FROT(2*L+2,2*L+2,2*L+2)
INTEGER I,J,M1,M2, NJ
INTEGER*8 BW

! Convert array into format usable by DSOFT:
BW = INT(L+1,8)
NJ = 2*L + 1

ILMM = CMPLX(0.D0, 0.D0, 8)
DO J=0,L
    ILMM(J,0,0) = IMML(0,0,J)
    DO M2=1,J
        ILMM(J,0,M2) = IMML(0,M2,J)
        ILMM(J,0,NJ-M2) = IMML(0,-M2,J)
        ILMM(J,M2,0) = IMML(M2,0,J)
        ILMM(J,NJ-M2,0) = IMML(-M2,0,J)
        DO M1=1,J
            ILMM(J,M1,M2) = IMML(M1,M2,J)
            ILMM(J,NJ-M1,M2) = IMML(-M1,M2,J)
            ILMM(J,M1,NJ-M2) = IMML(M1,-M2,J)
            ILMM(J,NJ-M1,NJ-M2) = IMML(-M1,-M2,J)
        ENDDO
    ENDDO
ENDDO

! Perform inverse discrete SO(3) Fourier Transform (DSOFT)
CALL ISOFT(ILMM, FROT, BW)
! Output is complex so must be converted back to real
OVERLAP = REAL(FROT, 8)

END SUBROUTINE CALCOVERLAP

SUBROUTINE FINDROTATIONS(OVERLAP, L, ANGLES, AMPLITUDES, NROTATIONS, DEBUG)
! Fits a set of Gaussians to the overlap integral and calculates the Euler angles these correspond to

USE FASTOVERLAPUTILS, ONLY: FINDPEAKS

IMPLICIT NONE

INTEGER, INTENT(IN) :: L
INTEGER, INTENT(INOUT) :: NROTATIONS
LOGICAL, INTENT(IN) :: DEBUG
DOUBLE PRECISION, INTENT(IN) :: OVERLAP(2*L+2,2*L+2,2*L+2)
DOUBLE PRECISION, INTENT(OUT) :: ANGLES(NROTATIONS,3), AMPLITUDES(NROTATIONS)

DOUBLE PRECISION CONVERT
INTEGER J

ANGLES=0.D0

CALL FINDPEAKS(OVERLAP, ANGLES, AMPLITUDES, NROTATIONS, DEBUG)

! Convert index locations to Euler Angles
CONVERT = PI / (2*L+2)
ANGLES(:NROTATIONS,1) = (ANGLES(:NROTATIONS,1)-1.0D0) * 2 * CONVERT
ANGLES(:NROTATIONS,2) = (ANGLES(:NROTATIONS,2)-0.5D0) * CONVERT
ANGLES(:NROTATIONS,3) = (ANGLES(:NROTATIONS,3)-1.0D0) * 2 * CONVERT

END SUBROUTINE FINDROTATIONS

SUBROUTINE EULERM(A,B,G,ROTM)
! Calculates rotation matrix of the Euler angles A,B,G
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: A,B,G
DOUBLE PRECISION, INTENT(OUT) :: ROTM(3,3)

DOUBLE PRECISION  COSA, SINA, COSB, SINB, COSG, SING

COSA = COS(A)
SINA = SIN(A)
COSB = COS(B)
SINB = SIN(B)
COSG = COS(G)
SING = SIN(G)

  ROTM (1,1) =   COSG * COSB * COSA  -  SING * SINA
  ROTM (1,2) = + SING * COSB * COSA  +  COSG * SINA
  ROTM (1,3) =          SINB * COSA
  ROTM (2,1) = - COSG * COSB * SINA  -  SING * COSA
  ROTM (2,2) = - SING * COSB * SINA  +  COSG * COSA
  ROTM (2,3) = -        SINB * SINA
  ROTM (3,1) = - COSG * SINB
  ROTM (3,2) = - SING * SINB
  ROTM (3,3) =          COSB

END SUBROUTINE EULERM

SUBROUTINE EULERINVM(A,B,G,ROTM)
! Calculates inverse (transposed) rotation matrix of the Euler angles A,B,G
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: A,B,G
DOUBLE PRECISION, INTENT(OUT) :: ROTM(3,3)

DOUBLE PRECISION  COSA, SINA, COSB, SINB, COSG, SING

COSA = COS(A)
SINA = SIN(A)
COSB = COS(B)
SINB = SIN(B)
COSG = COS(G)
SING = SIN(G)

  ROTM (1,1) =   COSG * COSB * COSA  -  SING * SINA
  ROTM (2,1) =   SING * COSB * COSA  +  COSG * SINA
  ROTM (3,1) =          SINB * COSA
  ROTM (1,2) = - COSG * COSB * SINA  -  SING * COSA
  ROTM (2,2) = - SING * COSB * SINA  +  COSG * COSA
  ROTM (3,2) = -        SINB * SINA
  ROTM (1,3) = - COSG * SINB
  ROTM (2,3) = - SING * SINB
  ROTM (3,3) =          COSB

END SUBROUTINE EULERINVM

SUBROUTINE SETCLUSTER()

USE COMMONS, ONLY : MYUNIT,NFREEZE,GEOMDIFFTOL,ORBITTOL,FREEZE,PULLT,TWOD,  &
    &   EFIELDT,AMBERT,QCIAMBERT,AMBER12T,CHRMMT,STOCKT,CSMT,PERMDIST,      &
    &   LOCALPERMDIST,LPERMDIST,OHCELLT,QCIPERMCHECK,PERMOPT,PERMINVOPT,    &
    &   NOINVERSION,GTHOMSONT,MKTRAPT,MULLERBROWNT,RIGID,OHCELLT

IMPLICIT NONE

MYUNIT = 6
NFREEZE = 0
GEOMDIFFTOL = 0.5D0
ORBITTOL = 1.0D-3

FREEZE = .FALSE.
PULLT = .FALSE.
TWOD = .FALSE.
EFIELDT = .FALSE.
AMBERT = .FALSE.
QCIAMBERT = .FALSE.
AMBER12T = .FALSE.
CHRMMT = .FALSE.
STOCKT = .FALSE.
CSMT = .FALSE.
PERMDIST = .TRUE.
LOCALPERMDIST = .FALSE.
LPERMDIST = .FALSE.
QCIPERMCHECK = .FALSE.
PERMOPT = .TRUE.
PERMINVOPT = .TRUE.
NOINVERSION = .FALSE.
GTHOMSONT = .FALSE.
MKTRAPT = .FALSE.
MULLERBROWNT = .FALSE.
RIGID = .FALSE.
OHCELLT = .FALSE.

END SUBROUTINE SETCLUSTER

SUBROUTINE CHECKKEYWORDS()

USE COMMONS, ONLY : MYUNIT,NFREEZE,GEOMDIFFTOL,ORBITTOL,FREEZE,PULLT,TWOD,  &
    &   EFIELDT,AMBERT,QCIAMBERT,AMBER12T,CHRMMT,STOCKT,CSMT,PERMDIST,      &
    &   LOCALPERMDIST,LPERMDIST,OHCELLT,QCIPERMCHECK,PERMOPT,PERMINVOPT,    &
    &   NOINVERSION,GTHOMSONT,MKTRAPT,MULLERBROWNT,RIGID, OHCELLT

IMPLICIT NONE

IF(OHCELLT) THEN
    WRITE(*,'(A)') 'ERROR - cluster fastoverlap not compatible with OHCELL keyword'
    STOP
ENDIF

IF(STOCKT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with STOCK keyword'
    STOP
ENDIF

IF(CSMT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with CSM keyword'
    STOP
ENDIF

IF(PULLT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with PULL keyword'
    STOP
ENDIF

IF(EFIELDT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with EFIELD keyword'
    STOP
ENDIF

IF(RIGID) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with RIGID keyword'
    STOP
ENDIF

IF(QCIPERMCHECK) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with QCIPERMCHECK keyword'
    STOP
ENDIF

IF(QCIAMBERT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with QCIAMBER keyword'
    STOP
ENDIF

IF(GTHOMSONT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with GTHOMSON keyword'
    STOP
ENDIF

IF(MKTRAPT) THEN
    WRITE(*,'(A)') 'ERROR - fastoverlap not compatible with MKTRAP keyword'
    STOP
ENDIF

END SUBROUTINE CHECKKEYWORDS

END MODULE CLUSTERFASTOVERLAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INCLUDE "bulkmindist.f90"
INCLUDE "minpermdist.f90"
INCLUDE "minperm.f90"
INCLUDE "newmindist.f90"
INCLUDE "orient.f90"
INCLUDE "legendre.f90"
