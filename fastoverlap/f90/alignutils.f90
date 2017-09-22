
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

! Subroutines:

!    ITERATIVEALIGN(COORDSB,COORDSA,NCOORDS,NDEBUG,NBOXLX,NBOXLY,NBOXLZ,NBULKT, &
!     & DISTANCE,DIST2,RMATBEST,DISPBEST,PERMBEST)
!        Main alignment algorithm
!        SAFE TO CALL AS LONG AS NPERMGROUP, NPERMSIZE and PERMGROUP exist
!        iteratively permutes then moves coordsa to best match coordsb
!        returns the rotation matrix RMATBEST or displacement vector DISPBEST  that
!        best maps coordsa onto coordsb along with the permutation PERMBEST
!        along with the distance DIST2 and the distance squared DISTANCE

!    MINIMISESEPARATION(COORDSB,COORDSA,NCOORDS,DISTANCE,RMATBEST,DISPBEST)
!        Moves coordsa to best match coordsb

!    FINDROTATION(COORDSB,COORDSA,NCOORDS,DIST,RMAT)
!        rotates coordsa around the origin to match coordsb

!    FINDDISPLACEMENT(COORDSB,COORDSA,NCOORDS,DIST,DISP)
!        minimizes the average displacement between points
!        (whilst applying periodic BC)

!    FINDBESTPERMUTATION(COORDSB,COORDSA,NCOORDS,NEWPERM,DISTANCE,DIST2)
!        finds the best permutational alignment between coordsa and coordsb

!    PERMPAIRDISTS(COORDSB,COORDSA,NCOORDS,MAXNEI,NDISTS,NIDX,NPERMGROUP)
!        calculates the value of the distance matrix between coordsa and coordsb
!        only up to the PMAXNEI nearest neighbour distances are stored

!    FINDBESTPERM(NDISTS,NIDX,NCOORDS,MAXNEI,PERM,DIST,NPERMGROUP,INFO)
!        solves the permutation problem given the results of PERMPAIRDISTS

!    PAIRDISTS(n, p, q, sx, sy, sz, pbc, cc, kk, maxnei)
!        calculates the pairwise distance matrix for a homoatomix pair of structures
!        p and q

!    REALLOCATEARRAYS()
!        this allocates the arrays needed by the algorithm

!    SETPERM(NEWNATOMS, NEWPERMGROUP, NEWNPERMSIZE)
!        this allocates the permutation arrays in Commons, not needed in GMIN or OPTIM

!    JOVOSAP(N,SZ,CC,KK,FIRST,X,Y,U,V,H)
!        this code finds the minimal permutation alignment between two structures,
!        abandon all hope all ye who enter this code

! functions:
!    PAIRDIST(C1, C2)
!        calculates the distance between points C1 and C2
!        includes periodic boundary conditions


! INCLUDE "commons.f90"

MODULE ALIGNUTILS

USE COMMONS, ONLY : PERMGROUP, NPERMSIZE, NPERMGROUP, BESTPERM, MYUNIT, &
 & NSETS, SETS, PERMINVOPT, NOINVERSION, BOXLX, BOXLY, BOXLZ, OHCELLT, TWOD!, PERMDIST, PERMOPT

IMPLICIT NONE

INTEGER, SAVE :: NATOMS, NLAP, NPERM, PATOMS, NTRIES, INFO
INTEGER, SAVE :: PMAXNEI = 60
DOUBLE PRECISION, PARAMETER :: PSCALE = 1.D6 ! Scale for linear assignment problem
INTEGER, PARAMETER :: MAXIMUMTRIES=20 ! Maximum number of iterations

! Arrays of distances and nearest neighbour distances
DOUBLE PRECISION, SAVE, ALLOCATABLE :: DUMMYDISTS(:,:), DUMMYNEARDISTS(:)

INTEGER, SAVE, ALLOCATABLE :: DUMMYIDX(:,:)
INTEGER, SAVE, ALLOCATABLE :: INVPERMGROUP(:)

DOUBLE PRECISION, SAVE :: ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), TMAT(3,3)
DOUBLE PRECISION, SAVE :: CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ, RMATCUMUL(3,3), RMATNEW(3,3)
DOUBLE PRECISION, SAVE :: NEWDISTANCE, NEWDIST2, PDIST2

!DOUBLE PRECISION, SAVE :: BOXLY, BOXLY, BOXLZ
DOUBLE PRECISION, SAVE :: BOXVEC(3), DISPCUMUL(3), DISPNEW(3)

! Used when solving assignment problem
DOUBLE PRECISION, SAVE, ALLOCATABLE :: PDUMMYA(:), PDUMMYB(:), DUMMYA(:), &
    & DUMMYB(:), DUMMY(:)
INTEGER, SAVE, ALLOCATABLE :: NEWPERM(:), LPERM(:), ALLPERM(:), SAVEPERM(:)

LOGICAL, SAVE :: DEBUG = .TRUE., SAVECOORDS = .TRUE., BULKT

! For saving alignments
INTEGER, SAVE :: NSTORED, NSAVE=20
DOUBLE PRECISION, SAVE :: DTOL=1E-3
DOUBLE PRECISION, SAVE, ALLOCATABLE ::  BESTDISTS(:), BESTCOORDS(:,:)
DOUBLE PRECISION, SAVE, ALLOCATABLE ::  BESTRMATS(:,:,:), BESTDISPS(:,:)

CONTAINS

SUBROUTINE ITERATIVEALIGN(COORDSB,COORDSA,NCOORDS,NDEBUG,NBOXLX,NBOXLY,NBOXLZ,NBULKT, &
 & DISTANCE,DIST2,RMATBEST,DISPBEST,PERMBEST)

INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION, INTENT(IN) :: NBOXLX, NBOXLY, NBOXLZ
LOGICAL, INTENT(IN) :: NDEBUG, NBULKT

DOUBLE PRECISION, INTENT(INOUT) :: COORDSA(3*NCOORDS), COORDSB(3*NCOORDS)
DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, DIST2, RMATBEST(3,3), DISPBEST(3)
INTEGER, INTENT(OUT) :: PERMBEST(NCOORDS)

INTEGER J1, J2, J3

! Setting module variables
DEBUG = NDEBUG
BULKT = NBULKT
NATOMS = NCOORDS
BOXLX = NBOXLX; BOXLY = NBOXLY; BOXLZ = NBOXLZ
BOXVEC = (/BOXLX,BOXLY,BOXLZ/)

CALL REALLOCATEARRAYS()

IF (BULKT) THEN
    DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)
    DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)

    DISPBEST(1:3) = 0.D0
ELSE
    ! Calculating centres of mass of coordinates
    ! Superimposing centre of mass of COORDSA with COORDSB
    ! Sets centres of mass of both structures to origin
    CMAX=0.0D0; CMAY=0.0D0; CMAZ=0.0D0
    DO J1=1,NATOMS
        CMAX=CMAX+COORDSA(3*(J1-1)+1)
        CMAY=CMAY+COORDSA(3*(J1-1)+2)
        CMAZ=CMAZ+COORDSA(3*(J1-1)+3)
    ENDDO
    CMAX=CMAX/NATOMS; CMAY=CMAY/NATOMS; CMAZ=CMAZ/NATOMS

    CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
    DO J1=1,NATOMS
        CMBX=CMBX+COORDSB(3*(J1-1)+1)
        CMBY=CMBY+COORDSB(3*(J1-1)+2)
        CMBZ=CMBZ+COORDSB(3*(J1-1)+3)
    ENDDO
    CMBX=CMBX/NATOMS; CMBY=CMBY/NATOMS; CMBZ=CMBZ/NATOMS

    DO J1=1,NATOMS
        DUMMYA(3*(J1-1)+1) = COORDSA(3*(J1-1)+1) - CMAX
        DUMMYA(3*(J1-1)+2) = COORDSA(3*(J1-1)+2) - CMAY
        DUMMYA(3*(J1-1)+3) = COORDSA(3*(J1-1)+3) - CMAZ

        DUMMYB(3*(J1-1)+1) = COORDSB(3*(J1-1)+1) - CMBX
        DUMMYB(3*(J1-1)+2) = COORDSB(3*(J1-1)+2) - CMBY
        DUMMYB(3*(J1-1)+3) = COORDSB(3*(J1-1)+3) - CMBZ
    ENDDO

    RMATBEST(1:3,1:3) = 0.0D0
    RMATBEST(1,1) = 1.0D0; RMATBEST(2,2) = 1.0D0; RMATBEST(3,3) = 1.0D0
END IF

DO J1=1,NATOMS
!    BESTPERM(J1)  = J1
    PERMBEST(J1) = J1
!    SAVEPERM(J1) = J1
ENDDO

NTRIES = 0
NPERM = NCOORDS
DO WHILE(NPERM.GT.0)

    IF (DEBUG) WRITE(MYUNIT,'(A,I2)') 'alignutils> beginning iteration ', NTRIES+1

    ! Saving unpermuted coordinates
    DUMMY(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
    SAVEPERM(1:NATOMS) = PERMBEST(1:NATOMS)

    CALL FINDBESTPERMUTATION(DUMMYB,DUMMYA,NATOMS,NEWPERM,NEWDISTANCE,PDIST2)

    ! Applying permutation
    NPERM = 0
    DO J1=1,NATOMS
        DUMMYA(3*(J1-1)+1)=DUMMY(3*(NEWPERM(J1)-1)+1)
        DUMMYA(3*(J1-1)+2)=DUMMY(3*(NEWPERM(J1)-1)+2)
        DUMMYA(3*(J1-1)+3)=DUMMY(3*(NEWPERM(J1)-1)+3)
        PERMBEST(J1) = SAVEPERM(NEWPERM(J1))
        IF (J1.NE.NEWPERM(J1)) THEN
            NPERM=NPERM+1
        ENDIF
    ENDDO

    IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,G20.10)') &
    & 'alignutils> distance after permuting ',NPERM,' pairs of atoms=', PDIST2

    CALL MINIMISESEPARATION(DUMMYB,DUMMYA,NATOMS,NEWDIST2,RMATNEW,DISPNEW)

    IF (DEBUG.AND.BULKT) THEN
        WRITE(MYUNIT,'(A,G20.10)') &
        & 'alignutils> distance after minimising displacement', NEWDIST2
    ELSE IF (DEBUG) THEN
        WRITE(MYUNIT,'(A,G20.10)') &
        & 'alignutils> distance after minimising rotation', NEWDIST2
    ENDIF

    ! Updating coordinates
    IF (BULKT) THEN
        DISPBEST = DISPBEST + DISPNEW
        DO J1=1,NATOMS
            DUMMYA(3*(J1-1)+1) = COORDSA(3*(PERMBEST(J1)-1)+1) + DISPBEST(1)
            DUMMYA(3*(J1-1)+2) = COORDSA(3*(PERMBEST(J1)-1)+2) + DISPBEST(2)
            DUMMYA(3*(J1-1)+3) = COORDSA(3*(PERMBEST(J1)-1)+3) + DISPBEST(3)
        ENDDO
    ELSE
        RMATBEST = MATMUL(RMATNEW, RMATBEST)
        DO J1=1,NATOMS
            DUMMYA(3*(J1-1)+1) = COORDSA(3*(PERMBEST(J1)-1)+1) - CMAX
            DUMMYA(3*(J1-1)+2) = COORDSA(3*(PERMBEST(J1)-1)+2) - CMAY
            DUMMYA(3*(J1-1)+3) = COORDSA(3*(PERMBEST(J1)-1)+3) - CMAZ

            DUMMYA(3*J1-2:3*J1) = MATMUL(RMATBEST, DUMMYA(3*J1-2:3*J1))
        ENDDO
    ENDIF

    NTRIES = NTRIES + 1

    IF (((NEWDIST2-PDIST2)/NEWDIST2).GT.(SQRT(1.D0*NCOORDS)/PSCALE)) THEN
        IF (DEBUG) WRITE(MYUNIT, '(A)') 'alignutils> WARNING - distance increased with nonzero permutations'
        EXIT
    ENDIF
    IF (NTRIES.GT.MAXIMUMTRIES) THEN
        IF (DEBUG) WRITE(MYUNIT, '(A)') 'alignutils> WARNING - number of tries exceeded'
        EXIT
    ENDIF
ENDDO

! Assigning solution to COORDSA
IF (BULKT) THEN
    COORDSA(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
ELSE
    COORDSA(1:3*NATOMS-2:3) = DUMMYA(1:3*NATOMS-2:3) + CMBX
    COORDSA(2:3*NATOMS-1:3) = DUMMYA(2:3*NATOMS-1:3) + CMBY
    COORDSA(3:3*NATOMS  :3) = DUMMYA(3:3*NATOMS  :3) + CMBZ
ENDIF

DISTANCE = NEWDIST2**2
DIST2 = NEWDIST2

IF (SAVECOORDS) CALL ADDCOORDS(COORDSA, NATOMS, BULKT, DIST2, RMATBEST, DISPBEST)

IF (DEBUG) THEN
    WRITE(MYUNIT, '(A,G20.10,A,I2,A)') 'alignutils> best distance found=', NEWDIST2, ' after ', NTRIES, ' iterations'
    IF (BULKT) THEN
        WRITE(MYUNIT, '(A)') 'alignutils> best displacement found:'
        WRITE(MYUNIT, '(3G20.10)') DISPBEST(1:3)
    ELSE
        WRITE(MYUNIT, '(A)') 'alignutils> best rotation found:'
        WRITE(MYUNIT, '(3G20.10)') RMATBEST(1:3,1:3)
    ENDIF
ENDIF

END SUBROUTINE ITERATIVEALIGN

SUBROUTINE MINIMISESEPARATION(COORDSB,COORDSA,NCOORDS,DISTANCE,RMATBEST,DISPBEST)

IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NCOORDS), COORDSB(3*NCOORDS)

DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, RMATBEST(3,3), DISPBEST(3)

IF (BULKT) THEN
    CALL FINDDISPLACEMENT(COORDSB,COORDSA,NCOORDS,DISTANCE,DISPBEST)
ELSE
    CALL FINDROTATION(COORDSB,COORDSA,NCOORDS,DISTANCE,RMATBEST)
ENDIF

END SUBROUTINE MINIMISESEPARATION

SUBROUTINE FINDROTATION(COORDSB,COORDSA,NCOORDS,DIST,RMAT)
! Finds the rotation that minimises the Euclidean distance between
! COORDSA onto COORDSB around the origin

IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NCOORDS), COORDSB(3*NCOORDS)

DOUBLE PRECISION, INTENT(OUT) :: RMAT(3,3), DIST

INTEGER, PARAMETER :: LWORK=12
INTEGER J1, JMIN, INFO
DOUBLE PRECISION QMAT(4,4), MINV, DIAG(4), TEMPA(LWORK), XM, YM, ZM, XP, YP, ZP
DOUBLE PRECISION Q1, Q2, Q3, Q4

!  The formula below is not invariant to overall translation because XP, YP, ZP
!  involve a sum of coordinates! We need to have COORDSA and COORDSB coordinate
!  centres both at the origin!!

QMAT(1:4,1:4)=0.0D0
DO J1=1,NCOORDS
      XM=COORDSB(3*(J1-1)+1)-COORDSA(3*(J1-1)+1)
      YM=COORDSB(3*(J1-1)+2)-COORDSA(3*(J1-1)+2)
      ZM=COORDSB(3*(J1-1)+3)-COORDSA(3*(J1-1)+3)
      XP=COORDSB(3*(J1-1)+1)+COORDSA(3*(J1-1)+1)
      YP=COORDSB(3*(J1-1)+2)+COORDSA(3*(J1-1)+2)
      ZP=COORDSB(3*(J1-1)+3)+COORDSA(3*(J1-1)+3)
      QMAT(1,1)=QMAT(1,1)+XM**2+YM**2+ZM**2
      QMAT(1,2)=QMAT(1,2)+YP*ZM-YM*ZP
      QMAT(1,3)=QMAT(1,3)+XM*ZP-XP*ZM
      QMAT(1,4)=QMAT(1,4)+XP*YM-XM*YP
      QMAT(2,2)=QMAT(2,2)+YP**2+ZP**2+XM**2
      QMAT(2,3)=QMAT(2,3)+XM*YM-XP*YP
      QMAT(2,4)=QMAT(2,4)+XM*ZM-XP*ZP
      QMAT(3,3)=QMAT(3,3)+XP**2+ZP**2+YM**2
      QMAT(3,4)=QMAT(3,4)+YM*ZM-YP*ZP
      QMAT(4,4)=QMAT(4,4)+XP**2+YP**2+ZM**2
ENDDO

QMAT(2,1)=QMAT(1,2); QMAT(3,1)=QMAT(1,3); QMAT(3,2)=QMAT(2,3)
QMAT(4,1)=QMAT(1,4); QMAT(4,2)=QMAT(2,4); QMAT(4,3)=QMAT(3,4)

CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,LWORK,INFO)
IF (INFO.NE.0) WRITE(MYUNIT,'(A,I6,A)') 'alignutils> FINDROTATION WARNING - INFO=',INFO,' in DSYEV'

MINV=1.0D100
DO J1=1,4
    IF (DIAG(J1).LT.MINV) THEN
    JMIN=J1
    MINV=DIAG(J1)
    ENDIF
ENDDO
IF (MINV.LT.0.0D0) THEN
    IF (ABS(MINV).LT.1.0D-6) THEN
        MINV=0.0D0
    ELSE
        WRITE(MYUNIT,'(A,G20.10,A)') 'alignutils> FINDROTATION WARNING MINV is ',MINV,' change to absolute value'
        MINV=-MINV
    ENDIF
ENDIF
DIST=SQRT(MINV)

!IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,I6)') 'alignutils> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)

RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2

END SUBROUTINE FINDROTATION

SUBROUTINE FINDDISPLACEMENT(COORDSB,COORDSA,NCOORDS,DIST,DISP)

IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NCOORDS), COORDSB(3*NCOORDS)

DOUBLE PRECISION, INTENT(OUT) :: DISP(3), DIST

INTEGER J1
DOUBLE PRECISION XM, YM, ZM

! Calculate average displacement
DO J1=1,NCOORDS
    XM = COORDSB(3*J1-2) - COORDSA(3*J1-2)
    YM = COORDSB(3*J1-1) - COORDSA(3*J1-1)
    DISP(1) = DISP(1) + XM - BOXLX*NINT(XM/BOXLX)
    DISP(2) = DISP(2) + YM - BOXLY*NINT(YM/BOXLY)
ENDDO

IF (TWOD) THEN
    DISP(3) = 0.D0
ELSE
    DO J1=1,NCOORDS
        ZM = COORDSB(3*J1  ) - COORDSA(3*J1  )
        DISP(3) = DISP(3) + ZM - BOXLZ*NINT(ZM/BOXLZ)
    ENDDO
END IF

DISP = DISP/NCOORDS

! Calculate new distance
DIST = 0.D0
DO J1=1,NCOORDS
    DIST = DIST + PAIRDIST(COORDSB(3*J1-2:3*J1),COORDSA(3*J1-2:3*J1)+DISP)
ENDDO
DIST = SQRT(DIST)

END SUBROUTINE FINDDISPLACEMENT

SUBROUTINE FINDBESTPERMUTATION(COORDSB,COORDSA,NCOORDS,NEWPERM,DISTANCE,DIST2)

IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NCOORDS), COORDSB(3*NCOORDS)

INTEGER, INTENT(OUT) :: NEWPERM(NCOORDS)
DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, DIST2

CALL PERMPAIRDISTS(COORDSB,COORDSA,NCOORDS,PMAXNEI,DUMMYDISTS,DUMMYIDX,NPERMGROUP)
CALL FINDBESTPERM(DUMMYDISTS,DUMMYIDX,NCOORDS,PMAXNEI,NEWPERM,DISTANCE,NPERMGROUP,INFO)

DIST2 = SQRT(DISTANCE)

IF ((INFO.GT.0).AND.DEBUG) WRITE(MYUNIT, "(A,I3)") &
 & "alignutils> WARNING LAP algorithm failed to align npoints= ", INFO

END SUBROUTINE FINDBESTPERMUTATION

SUBROUTINE PERMPAIRDISTS(COORDSB,COORDSA,NCOORDS,MAXNEI,NDISTS,NIDX,NPERMGROUP)

! Calculates the maxtrix of closest distances between COORDSB and COORDSA
! Only stores up to MAXNEI nearest neighbours
! NIDX returns the indexes of the nearest neighbour distances, contained in NDISTS
! Uses module variables BOXLX, BOXLY, BOXLZ, BULKT when calculating periodic distances

IMPLICIT NONE

INTEGER, INTENT(IN) :: NCOORDS, NPERMGROUP, MAXNEI
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NCOORDS), COORDSB(3*NCOORDS)

INTEGER, INTENT(OUT) :: NIDX(MAXNEI*NCOORDS,NPERMGROUP)
DOUBLE PRECISION, INTENT(OUT) :: NDISTS(MAXNEI*NCOORDS,NPERMGROUP)

INTEGER NDUMMY,J1,J2,NPERM

NATOMS = NCOORDS
CALL REALLOCATEARRAYS()

NDUMMY = 0

NIDX   = -1
NDISTS = HUGE(1.D0)

DO J1=1,NPERMGROUP
    NPERM=NPERMSIZE(J1)
    DO J2=1,NPERM
        PDUMMYA(3*(J2-1)+1)=COORDSA(3*((PERMGROUP(NDUMMY+J2))-1)+1)
        PDUMMYA(3*(J2-1)+2)=COORDSA(3*((PERMGROUP(NDUMMY+J2))-1)+2)
        PDUMMYA(3*(J2-1)+3)=COORDSA(3*((PERMGROUP(NDUMMY+J2))-1)+3)
        PDUMMYB(3*(J2-1)+1)=COORDSB(3*((PERMGROUP(NDUMMY+J2))-1)+1)
        PDUMMYB(3*(J2-1)+2)=COORDSB(3*((PERMGROUP(NDUMMY+J2))-1)+2)
        PDUMMYB(3*(J2-1)+3)=COORDSB(3*((PERMGROUP(NDUMMY+J2))-1)+3)
    ENDDO
    CALL PAIRDISTS(NPERM,PDUMMYB(1:3*NPERM),PDUMMYA(1:3*NPERM),BOXLX,BOXLY, &
 & BOXLZ,BULKT,NDISTS(1:MAXNEI*NPERM,J1),NIDX(1:MAXNEI*NPERM,J1),MAXNEI)
    NDUMMY = NDUMMY + NPERM
ENDDO

END SUBROUTINE PERMPAIRDISTS

SUBROUTINE FINDBESTPERM(NDISTS,NIDX,NCOORDS,MAXNEI,PERM,DIST,NPERMGROUP,INFO)

! Solves assignment problem using the shortest augmenting path algorithm:
! Jonker, R., & Volgenant, A. (1987).
! A shortest augmenting path algorithm for dense and sparse linear assignment problems.
! Computing, 38(4), 325–340. http://doi.org/10.1007/BF02278710

! This calculates the exact distance as well!

! Code copied from GMIN/source/minperm.f90

IMPLICIT NONE

INTEGER, INTENT(IN) :: NCOORDS,NPERMGROUP,MAXNEI,NIDX(MAXNEI*NCOORDS,NPERMGROUP)
DOUBLE PRECISION, INTENT(IN) :: NDISTS(MAXNEI*NCOORDS,NPERMGROUP)

DOUBLE PRECISION, INTENT(OUT) :: DIST
INTEGER, INTENT(OUT) :: PERM(NCOORDS), INFO

! COULD SET THESE AS MODULE VARIABLES
INTEGER*8 :: KK(NCOORDS*MAXNEI), CC(NCOORDS*MAXNEI)
INTEGER*8 :: FIRST(NCOORDS+1), X(NCOORDS), Y(NCOORDS)
INTEGER*8 :: U(NCOORDS), V(NCOORDS), N8, SZ8, H
!INTEGER(KIND=INT64) :: KK(NATOMS*MAXNEI), CC(NATOMS*MAXNEI)
!INTEGER(KIND=INT64) :: FIRST(NATOMS+1), X(NATOMS), Y(NATOMS)
!INTEGER(KIND=INT64) :: U(NATOMS), V(NATOMS), N8, SZ8, H
INTEGER N,M,I,J,K,K1,I1,J1,J2,NDUMMY

DOUBLE PRECISION D2

NATOMS = NCOORDS
CALL REALLOCATEARRAYS()

D2=0.D0
DIST=0.D0
INFO=0

NDUMMY=0

DO J1=1,NPERMGROUP

    N = NPERMSIZE(J1)
    M = MAXNEI
    IF(N.LE.MAXNEI) M=N
    SZ8 = M*N
    N8 = N

    DO I=0,N
        FIRST(I+1) = I*M +1
    ENDDO
    KK = -1
    CC = HUGE(1)
    DO J=1,N
        K = FIRST(J)-1
        DO I=1,M
            KK(I+K) = NIDX(I+K,J1)
            CC(I+K) = INT(NDISTS(I+K,J1)*PSCALE, 8)
        ENDDO
    ENDDO

    ! Solving the assignment problem to deduce the correct permutation
    CALL JOVOSAP(N8, SZ8, CC(:M*N), KK(:M*N), FIRST(:N+1), Y(:N), X(:N), U(:N), V(:N), H)
    NLAP = NLAP + 1

    DO J=1,N
        IF (Y(J).GT.N) THEN
            Y(J)=N
            INFO = INFO + 1
        END IF
        IF (Y(J).LT.1) THEN
            Y(J)=1
            INFO = INFO + 1
        END IF
        PERM(PERMGROUP(NDUMMY+J)) = PERMGROUP(NDUMMY+Y(J))

        ! Calculating exact distance
        K = FIRST(J)-1
        J2 = MIN(Y(J),M)
        IF (Y(J).NE.NIDX(J2+K,J1)) THEN
            DO J2=1,M !If N>MAXNEI then we must search the list
                IF (Y(J).EQ.NIDX(J2+K,J1)) EXIT
            ENDDO
        END IF
        DIST = DIST + NDISTS(J2+K,J1)
    ENDDO

    ! untested!!
    IF (NSETS(J1).GT.0) THEN
        DO I=1,N
            DO K=1,NSETS(J1)
                PERM(SETS(PERMGROUP(NDUMMY+I),K))=SETS(PERM(PERMGROUP(NDUMMY+Y(I))),K)
            ENDDO
        ENDDO
    ENDIF

    NDUMMY = NDUMMY + NPERMSIZE(J1)
ENDDO

END SUBROUTINE FINDBESTPERM

SUBROUTINE PAIRDISTS(n, p, q, sx, sy, sz, pbc, cc, kk, maxnei)
      implicit none

!     Input
!       n  : System size
!       p,q: Coordinate vectors (n particles)
!       s  : Box lengths (or dummy if open B.C.)
!       pbc: Periodic boundary conditions?
      integer, intent(in) :: n, maxnei
      double precision, intent(in) :: p(3*n), q(3*n), sx, sy, sz
      logical, intent(in) :: pbc
      double precision s(3)

!     Output
!       perm: Permutation so that p(i) <--> q(perm(i))
!       dist: Minimum attainable distance
!     We have
      double precision, intent(out) :: cc(n*maxnei)
      integer, intent(out) :: kk(n*maxnei)
      double precision DUMMY

!     Parameters
!       scale : Precision
!       maxnei: Maximum number of closest neighbourspa
      double precision scale, d, h

      parameter (scale = 1.0d6   )
!      parameter (maxnei = 60     )

      integer*8 first(n+1)!, x(n), y(n)
!      integer*8 u(n), v(n)
      integer   m, i, j, k, l, l2, t, a
      integer*8 n8, sz8
      integer J1

      BOXVEC = (/sx,sy,sz/)
      s(1)=sx
      s(2)=sy
      s(3)=sz
      m = maxnei
      if(n .le. maxnei) m = n
      sz8 = m*n
      n8 = n

      do i=0,n
         first(i+1) = i*m + 1
      enddo

      if(m .eq. n) then
!     Compute the full matrix...
         do i=1,n
            k = first(i)-1
            do j=1,n
               cc(k+j) = PAIRDIST(p(3*i-2), q(3*j-2))
               kk(k+j) = j
!              write(*,*) i, j, '-->', cc(k+j)
            enddo
         enddo
      else
!     We need to store the distances of the maxnei closeest neighbors
!     of each particle. The following builds a heap to keep track of
!     the maxnei closest neighbours seen so far. It might be more
!     efficient to use quick-select instead... (This is definitely
!     true in the limit of infinite systems.)
        do i=1,n
           k = first(i)-1
           do j=1,m
              cc(k+j) = PAIRDIST(p(3*i-2), q(3*j-2))
              kk(k+j) = j
              l = j
10            if(l .le. 1) goto 11
              l2 = l/2
              if(cc(k+l2) .lt. cc(k+l)) then
                 h = cc(k+l2)
                 cc(k+l2) = cc(k+l)
                 cc(k+l) = h
                 t = kk(k+l2)
                 kk(k+l2) = kk(k+l)
                 kk(k+l) = t
                 l = l2
                 goto 10
              endif
11         enddo

           do j=m+1,n
              d = PAIRDIST(p(3*i-2), q(3*j-2))
              if(d .lt. cc(k+1)) then
                 cc(k+1) = d
                 kk(k+1) = j
                 l = 1
20               l2 = 2*l
                 if(l2+1 .gt. m) goto 21
                 if(cc(k+l2+1) .gt. cc(k+l2)) then
                    a = k+l2+1
                 else
                    a = k+l2
                 endif
                 if(cc(a) .gt. cc(k+l)) then
                    h = cc(a)
                    cc(a) = cc(k+l)
                    cc(k+l) = h
                    t = kk(a)
                    kk(a) = kk(k+l)
                    kk(k+l) = t
                    l = a-k
                    goto 20
                 endif
21               if (l2 .le. m) THEN ! split IF statements to avoid a segmentation fault
                    IF (cc(k+l2) .gt. cc(k+l)) then
                       h = cc(k+l2)
                       cc(k+l2) = cc(k+l)
                       cc(k+l) = h
                       t = kk(k+l2)
                       kk(k+l2) = kk(k+l)
                       kk(k+l) = t
                    ENDIF
                 endif
              endif
           enddo
        enddo
      ENDIF

END SUBROUTINE PAIRDISTS

FUNCTION PAIRDIST(C1, C2) RESULT(DIST)

! Calculates distance^2 between points C1 and C2
! Requires BULKT and BOXVEC variables to be set

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: C1(3), C2(3)
DOUBLE PRECISION T, DIST

INTEGER I

IF (BULKT) THEN
    DO I=1,3
        IF (BOXVEC(I).NE.0.0D0) THEN
            T = C1(i) - C2(i)
            T = T - BOXVEC(i)*anint(T/BOXVEC(I))
            DIST = DIST + T*T
        ENDIF
    ENDDO
ELSE
    DIST = (C1(1) - C2(1))**2+(C1(2) - C2(2))**2+(C1(3) - C2(3))**2
ENDIF

END FUNCTION PAIRDIST

SUBROUTINE REALLOCATEARRAYS()

IMPLICIT NONE

IF (SIZE(DUMMYDISTS).NE.(PMAXNEI*NATOMS*NPERMGROUP)) THEN
    IF (DEBUG) WRITE(MYUNIT,"(A)") 'alignutils> reallocating distance arrays'
    IF(ALLOCATED(DUMMYDISTS)) DEALLOCATE(DUMMYDISTS,DUMMYNEARDISTS, &
     & DUMMYIDX)
    ALLOCATE(DUMMYDISTS(PMAXNEI*NATOMS,NPERMGROUP),DUMMYNEARDISTS(NATOMS), &
     & DUMMYIDX(PMAXNEI*NATOMS,NPERMGROUP))
END IF

IF (SIZE(LPERM).NE.NATOMS) THEN
    IF (DEBUG) WRITE(MYUNIT,"(A)") 'alignutils> reallocating coordinate arrays'
    IF(ALLOCATED(PDUMMYA)) DEALLOCATE(PDUMMYA,PDUMMYB,DUMMYA,DUMMYB)
    IF(ALLOCATED(NEWPERM)) DEALLOCATE(NEWPERM,LPERM,ALLPERM,SAVEPERM,INVPERMGROUP,DUMMY)
    ALLOCATE(PDUMMYA(3*NATOMS),PDUMMYB(3*NATOMS),DUMMYA(3*NATOMS),DUMMY(3*NATOMS), &
     & DUMMYB(3*NATOMS),NEWPERM(NATOMS),LPERM(NATOMS),SAVEPERM(NATOMS), &
     & ALLPERM(NATOMS),INVPERMGROUP(NATOMS))
END IF

IF (SAVECOORDS.AND.(SIZE(BESTCOORDS).NE.(NSAVE*NATOMS*3))) THEN
    IF (DEBUG) WRITE(MYUNIT, "(A,I3,A)") "alignutils> reallocating arrays to save ", NSAVE, " coordinates"
    NSTORED = 0
    IF (ALLOCATED(BESTDISTS)) DEALLOCATE(BESTDISTS,BESTCOORDS,BESTRMATS,BESTDISPS)
    ALLOCATE(BESTDISTS(NSAVE),BESTCOORDS(3*NATOMS,NSAVE),BESTRMATS(3,3,NSAVE),BESTDISPS(3,NSAVE))
END IF

END SUBROUTINE REALLOCATEARRAYS

SUBROUTINE SETPERM(NEWNATOMS, NEWPERMGROUP, NEWNPERMSIZE)
! Not needed for GMIN/OPTIM/PATHSAMPLE
! (Re)allocates arrays that define allowed permuations
IMPLICIT NONE

INTEGER, INTENT(IN) :: NEWNATOMS, NEWPERMGROUP(:), NEWNPERMSIZE(:)

IF(.NOT.SIZE(PERMGROUP).EQ.SIZE(NEWPERMGROUP)) THEN
    IF(ALLOCATED(PERMGROUP)) THEN
        DEALLOCATE(PERMGROUP)
    ENDIF
    ALLOCATE(PERMGROUP(SIZE(NEWPERMGROUP)))
ENDIF

NPERMGROUP = SIZE(NEWNPERMSIZE)
IF(.NOT.SIZE(NPERMSIZE).EQ.SIZE(NEWNPERMSIZE)) THEN
    IF(ALLOCATED(NPERMSIZE)) THEN
        DEALLOCATE(NPERMSIZE)
    ENDIF
    ALLOCATE(NPERMSIZE(NPERMGROUP))
ENDIF

IF(.NOT.SIZE(BESTPERM).EQ.NEWNATOMS) THEN
    IF(ALLOCATED(BESTPERM)) THEN
        DEALLOCATE(BESTPERM)
    ENDIF
    ALLOCATE(BESTPERM(NEWNATOMS))
ENDIF

IF(.NOT.SIZE(NSETS).EQ.(3*NEWNATOMS)) THEN
    IF(ALLOCATED(NSETS)) THEN
        DEALLOCATE(NSETS)
    ENDIF
    ALLOCATE(NSETS(3*NEWNATOMS))
ENDIF

IF(.NOT.SIZE(SETS).EQ.(3*NEWNATOMS*70)) THEN
    IF(ALLOCATED(SETS)) THEN
        DEALLOCATE(SETS)
    ENDIF
    ALLOCATE(SETS(3*NEWNATOMS,70))
ENDIF

NATOMS = NEWNATOMS
PERMGROUP = NEWPERMGROUP
NPERMSIZE = NEWNPERMSIZE
NSETS = 0

CALL REALLOCATEARRAYS()

END SUBROUTINE SETPERM

SUBROUTINE OHOPS(X,Y,OPNUM,NLOCAL)
IMPLICIT NONE
INTEGER OPNUM, J2, J3, NLOCAL
DOUBLE PRECISION RMAT(3,3,48), X(3*NLOCAL), Y(3*NLOCAL)
DATA RMAT / &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  -1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & 0.0D0,  1.00000000000D0,  0.0D0,   &
 & 1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0,   &
 & 0.0D0,  -1.00000000000D0,  0.0D0,   &
 & -1.00000000000D0,  0.0D0,  0.0D0,   &
 & 0.0D0,  0.0D0,  1.00000000000D0 /

DO J2=1,NLOCAL
   J3=3*(J2-1)
   Y(J3+1)=RMAT(1,1,OPNUM)*X(J3+1)+RMAT(1,2,OPNUM)*X(J3+2)+RMAT(1,3,OPNUM)*X(J3+3)
   Y(J3+2)=RMAT(2,1,OPNUM)*X(J3+1)+RMAT(2,2,OPNUM)*X(J3+2)+RMAT(2,3,OPNUM)*X(J3+3)
   Y(J3+3)=RMAT(3,1,OPNUM)*X(J3+1)+RMAT(3,2,OPNUM)*X(J3+2)+RMAT(3,3,OPNUM)*X(J3+3)
ENDDO

END SUBROUTINE OHOPS

SUBROUTINE JOVOSAP(N,SZ,CC,KK,FIRST,X,Y,U,V,H)
      IMPLICIT NONE
      INTEGER*8, INTENT(IN)  :: N, SZ
      INTEGER*8, INTENT(IN)  :: CC(SZ),KK(SZ),FIRST(N+1)
      INTEGER*8, INTENT(OUT) :: X(N),Y(N),U(N),V(N),H
      INTEGER*8 CNT,L0,T,T0,TD,V0,VJ,DJ
      INTEGER*8 LAB(N),D(N),FREE(N),TODO(N)
      LOGICAL OK(N)
      INTEGER*8 J, I, J0, L, J1, MIN, K, I0
      INTEGER*8 BIGINT
      J1 = -1
      J0 = -1

!     I don't know how to make g77 read integer*8 constants/parameters.
!       PARAMETER (BIGINT = 10**12) does not work(!)
!     nor does
!       PARAMETER (BIGINT = 1000000000000)
!     but this seems to be ok:
      BIGINT = 10**9
      BIGINT = BIGINT * 1000

!
! THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
! ACCORDING
!
!   "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
!    Assignment Problems," Computing 38, 325-340, 1987
!
!   by
!
!   R. Jonker and A. Volgenant, University of Amsterdam.
!
!
! INPUT PARAMETERS :
! N = NUMBER OF ROWS AND COLUMNS
! C = WEIGHT MATRIX
!
! OUTPUT PARAMETERS
! X = COL ASSIGNED TO ROW
! Y = ROW ASSIGNED TO COL
! U = DUAL ROW VARIABLE
! V = DUAL COLUMN VARIABLE
! H = VALUE OF OPTIMAL SOLUTION
!
! INITIALIZATION

!     Next line added by tomaso@nada.kth.se, to enable detection
!     of solutions being equivalent to the initial guess

!
!  If Y(:) is initialised to zero then we see segmentation faults if
!  a Y element is unset, etc.
!

      Y(1:N) = 0
      X(1:N) = 0
      TODO(1:N)=0
      h = -1
      DO 10 J=1,N
         V(J)=BIGINT
   10 CONTINUE
      DO 20 I=1,N
         X(I)=0
         DO 15 T=FIRST(I),FIRST(I+1)-1
            J=KK(T)
            IF (CC(T).LT.V(J)) THEN
              V(J)=CC(T)
              Y(J)=I
            END IF
   15    CONTINUE
   20 CONTINUE
      DO 30 J=1,N
         J0=N-J+1
         I=Y(J0)
         IF (I.EQ.0) THEN
!           PRINT '(A,I6,A)','minperm> WARNING B - matching failed'
            RETURN
         ENDIF
         IF (X(I).NE.0) THEN
           X(I)=-ABS(X(I))
           Y(J0)=0
         ELSE
           X(I)=J0
         END IF
   30 CONTINUE
      L=0
      DO 40 I=1,N
         IF (X(I).EQ.0) THEN
           L=L+1
           FREE(L)=I
           GOTO 40
         END IF
         IF (X(I).LT.0) THEN
           X(I)=-X(I)
         ELSE
           J1=X(I)
           MIN=BIGINT
           DO 31 T=FIRST(I),FIRST(I+1)-1
              J=KK(T)
              IF (J.EQ.J1) GOTO 31
              IF (CC(T)-V(J).LT.MIN) MIN=CC(T)-V(J)
   31      CONTINUE
           V(J1)=V(J1)-MIN
         END IF
   40 CONTINUE
! IMPROVE THE INITIAL SOLUTION
      CNT=0
      IF (L.EQ.0) RETURN
   41 L0=L
      K=1
      L=0
   50 I=FREE(K)
      K=K+1
      V0=BIGINT
      VJ=BIGINT
      DO 42 T=FIRST(I),FIRST(I+1)-1
         J=KK(T)
         H=CC(T)-V(J)
         IF (H.LT.VJ) THEN
           IF (H.GE.V0) THEN
             VJ=H
             J1=J
           ELSE
             VJ=V0
             V0=H
             J1=J0
             J0=J
           END IF
         END IF
   42 CONTINUE
      I0=Y(J0)
      IF (V0.LT.VJ) THEN
        V(J0)=V(J0)-VJ+V0
      ELSE
         if (j1 .lt. 0) then
            write(*,*) "error j1 is being used uninitialized"
            stop
         endif
        IF (I0.EQ.0) GOTO 43
        J0=J1
        I0=Y(J1)
      END IF
      IF (I0.EQ.0) GOTO 43
      IF (V0.LT.VJ) THEN
        K=K-1
        FREE(K)=I0
      ELSE
        L=L+1
        FREE(L)=I0
      END IF
   43 X(I)=J0
      Y(J0)=I
      IF (K.LE.L0) GOTO 50
      CNT=CNT+1
      IF ((L.GT.0).AND.(CNT.LT.2)) GOTO 41
! AUGMENTATION PART
      L0=L
      DO 90 L=1,L0
         DO 51 J=1,N
            OK(J)=.FALSE.
            D(J)=BIGINT
   51    CONTINUE
         MIN=BIGINT
         I0=FREE(L)
         TD=N
         DO 52 T=FIRST(I0),FIRST(I0+1)-1
            J=KK(T)
            DJ=CC(T)-V(J)
            D(J)=DJ
            LAB(J)=I0
            IF (DJ.LE.MIN) THEN
              IF (DJ.LT.MIN) THEN
                MIN=DJ
                K=1
                TODO(1)=J
              ELSE
                K=K+1
                TODO(K)=J
              END IF
            END IF
   52    CONTINUE
         DO 53 H=1,K
            J=TODO(H)
            IF (J.EQ.0) THEN
!              PRINT '(A,I6,A)','minperm> WARNING C - matching failed'
               RETURN
            ENDIF
            IF (Y(J).EQ.0) GOTO 80
            OK(J)=.TRUE.
   53    CONTINUE
! REPEAT UNTIL A FREE ROW HAS BEEN FOUND
   60    IF (K.EQ.0) THEN
!           PRINT '(A,I6,A)','minperm> WARNING D - matching failed'
            RETURN
         ENDIF
         J0=TODO(K)
         K=K-1
         I=Y(J0)
         TODO(TD)=J0
         TD=TD-1
         T0=FIRST(I)
         T=T0-1
   61    T=T+1
         IF (KK(T).NE.J0) GOTO 61
         H=CC(T)-V(J0)-MIN
         DO 62 T=T0,FIRST(I+1)-1
            J=KK(T)
            IF (.NOT. OK(J)) THEN
              VJ=CC(T)-H-V(J)
              IF (VJ.LT.D(J)) THEN
                D(J)=VJ
                LAB(J)=I
                IF (VJ.EQ.MIN) THEN
                  IF (Y(J).EQ.0) GOTO 70
                  K=K+1
                  TODO(K)=J
                  OK(J)=.TRUE.
                END IF
              END IF
            END IF
   62    CONTINUE
         IF (K.NE.0) GOTO 60
         MIN=BIGINT-1
         DO 63 J=1,N
            IF (D(J).LE.MIN) THEN
              IF (.NOT. OK(J)) THEN
                IF (D(J).LT.MIN) THEN
                  MIN=D(J)
                  K=1
                  TODO(1)=J
                ELSE
                  K=K+1
                  TODO(K)=J
                END IF
              END IF
            END IF
   63    CONTINUE
         DO 64 J0=1,K
            J=TODO(J0)
            IF (Y(J).EQ.0) GOTO 70
            OK(J)=.TRUE.
   64    CONTINUE
         GOTO 60
   70    IF (MIN.EQ.0) GOTO 80
         DO 71 K=TD+1,N
            J0=TODO(K)
            V(J0)=V(J0)+D(J0)-MIN
   71    CONTINUE
   80    I=LAB(J)
         Y(J)=I
         K=J
         J=X(I)
         X(I)=K
         IF (I0.NE.I) GOTO 80
   90 CONTINUE
      H=0
      DO 100 I=1,N
         J=X(I)
         T=FIRST(I)-1
  101    T=T+1
         IF (T.GT.SZ) THEN
            PRINT '(A,I6,A)','alignutils> WARNING D - atom ',I,' not matched - maximum number of neighbours too small?'
            RETURN
         ENDIF
         IF (KK(T).NE.J) GOTO 101
         DJ=CC(T)
         U(I)=DJ-V(J)
         H=H+DJ
  100 CONTINUE

END SUBROUTINE JOVOSAP

SUBROUTINE ADDCOORDS(COORDS, NCOORDS, NBULKT, DIST, RMAT, DISP)

IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NCOORDS), DIST, RMAT(3,3), DISP(3)
LOGICAL, INTENT(IN) :: NBULKT

INTEGER J, STARTSHIFT
DOUBLE PRECISION DIFF

BULKT = NBULKT

NATOMS = NCOORDS
CALL REALLOCATEARRAYS()

IF (NSTORED.EQ.0) THEN
    STARTSHIFT = 1
ENDIF

DO STARTSHIFT=1,NSTORED
    IF (ABS(DIST-BESTDISTS(STARTSHIFT)).LT.DTOL) THEN
        ! Testing whether structure identical to one already stored
        DIFF = 0.D0
        DO J=1,NCOORDS
            DIFF = DIFF + PAIRDIST(BESTCOORDS(3*J-2:3*J,STARTSHIFT),COORDS(3*J-2:3*J))
        ENDDO
        IF (SQRT(DIFF).LT.DTOL) THEN
            IF (DEBUG) WRITE(MYUNIT, "(A,I3)") &
     & "alignutils> structure being added identical to structure=", STARTSHIFT
            RETURN
        END IF
    END IF
    IF (DIST.LT.BESTDISTS(STARTSHIFT)) EXIT
END DO

IF (STARTSHIFT.LE.(NSTORED+1).AND.(STARTSHIFT.LE.NSAVE)) THEN
    IF (DEBUG) WRITE(MYUNIT, "(A,I3,A,I3)") &
     & "alignutils> saving coords, added at=",STARTSHIFT, " total stored=", NSTORED
    CALL SHIFTCOORDS(STARTSHIFT)
    BESTDISTS(STARTSHIFT) = DIST
    BESTCOORDS(:,STARTSHIFT) = COORDS(:)
    BESTRMATS(:,:,STARTSHIFT) = RMAT(:,:)
    BESTDISPS(:,STARTSHIFT) = DISP(:)
ENDIF

END SUBROUTINE ADDCOORDS

SUBROUTINE PRINTDISTANCES()

IMPLICIT NONE
INTEGER J

WRITE(MYUNIT, "(A,I3,A)") "alignutils> found", NSTORED, " candidate alignments with distances:"
DO J=1,NSTORED
    WRITE(MYUNIT, "(G20.10)") BESTDISTS(J)
END DO

END SUBROUTINE PRINTDISTANCES

SUBROUTINE SHIFTCOORDS(STARTSHIFT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: STARTSHIFT

INTEGER J,MAXJ

MAXJ = MIN(NSTORED,NSAVE-1)
DO J=MAXJ,STARTSHIFT,-1
    BESTDISTS(J+1) = BESTDISTS(J)
    BESTCOORDS(:,J+1) = BESTCOORDS(:,J)
    BESTRMATS(:,:,J+1) = BESTRMATS(:,:,J)
    BESTDISPS(:,J+1) = BESTDISPS(:,J)
END DO

NSTORED = MIN(NSTORED+1,NSAVE)

END SUBROUTINE SHIFTCOORDS

END MODULE

!MODULE SAVEALIGNMENTS
!
!USE ALIGNUTILS, ONLY : BULKT, NSAVE, DEBUG, MYUNIT
!IMPLICIT NONE
!
!INTEGER, SAVE :: NATOMS=0,NSTORED=0
!
!DOUBLE PRECISION, SAVE :: DTOL=1E-5
!DOUBLE PRECISION, SAVE, ALLOCATABLE ::  BESTDISTS(:), BESTCOORDS(:,:)
!DOUBLE PRECISION, SAVE, ALLOCATABLE ::  BESTRMATS(:,:,:), BESTDISPS(:,:)
!
!CONTAINS
!
!SUBROUTINE REALLOCATEARRAYS(NCOORDS)
!
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: NCOORDS
!
!IF ((NCOORDS.NE.NATOMS).OR.(SIZE(BESTDISTS).NE.NSAVE)) THEN
!    IF (DEBUG) WRITE(MYUNIT, "(A,I3,A)") "savealignments> reallocating arrays to save ", NSAVE, " coordinates"
!    NATOMS = NCOORDS
!    NSTORED = 0
!    IF (ALLOCATED(BESTDISTS)) DEALLOCATE(BESTDISTS,BESTCOORDS,BESTRMATS,BESTDISPS)
!    ALLOCATE(BESTDISTS(NSAVE),BESTCOORDS(3*NCOORDS,NSAVE),BESTRMATS(3,3,NSAVE),BESTDISPS(3,NSAVE))
!END IF
!
!END SUBROUTINE REALLOCATEARRAYS
!
!SUBROUTINE ADDCOORDS(COORDS, NCOORDS, NBULKT, DIST, RMAT, DISP)
!
!USE ALIGNUTILS, ONLY : PAIRDIST
!
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: NCOORDS
!DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NCOORDS), DIST, RMAT(3,3), DISP(3)
!LOGICAL, INTENT(IN) :: NBULKT
!
!INTEGER J, STARTSHIFT
!DOUBLE PRECISION DIFF
!
!BULKT = NBULKT
!
!CALL REALLOCATEARRAYS(NCOORDS)
!
!DO STARTSHIFT=1,NSTORED
!    IF (ABS(DIST-BESTDISTS(STARTSHIFT)).LT.DTOL) THEN
!        ! Testing whether structure identical to one already stored
!        DIFF = 0.D0
!        DO J=1,NCOORDS
!            DIFF = DIFF + PAIRDIST(BESTCOORDS(3*J-2:3*J,STARTSHIFT),COORDS(3*J-2:3*J))
!        ENDDO
!        IF (SQRT(DIFF).LT.DTOL) RETURN
!    END IF
!    IF (DIST.LT.BESTDISTS(STARTSHIFT)) EXIT
!END DO
!
!IF (STARTSHIFT.LT.NSTORED) THEN
!    IF (DEBUG) WRITE(MYUNIT, "(A,I3,A,I3)") &
!     & "savealignments> saving coords, added at=",STARTSHIFT, " total stored=", NSTORED
!    CALL SHIFTCOORDS(STARTSHIFT)
!    BESTDISTS(STARTSHIFT) = DIST
!    BESTCOORDS(:,STARTSHIFT) = COORDS(1:NCOORDS)
!    BESTRMATS(:,:,STARTSHIFT) = RMAT(:,:)
!    BESTDISPS(:,STARTSHIFT) = DISP(:)
!ENDIF
!
!END SUBROUTINE ADDCOORDS
!
!SUBROUTINE SHIFTCOORDS(STARTSHIFT)
!
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: STARTSHIFT
!
!INTEGER J,MAXJ
!
!MAXJ = MIN(NSTORED,NSAVE-1)
!DO J=MAXJ,STARTSHIFT,-1
!    BESTDISTS(J+1) = BESTDISTS(J)
!    BESTCOORDS(:,J+1) = BESTCOORDS(:,J)
!    BESTRMATS(:,:,J+1) = BESTRMATS(:,:,J)
!    BESTDISPS(:,J+1) = BESTDISPS(:,J)
!END DO
!
!NSTORED = MIN(NSTORED+1,NSAVE)
!
!END SUBROUTINE SHIFTCOORDS
!
!END MODULE