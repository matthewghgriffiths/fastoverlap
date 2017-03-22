!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

MODULE COMMONS
    IMPLICIT NONE
    
    INTEGER, SAVE ::  NFREEZE, NPERMGROUP, MYUNIT, NTSITES, BESTINVERT, NATOMS
    INTEGER, SAVE, ALLOCATABLE :: NPERMSIZE(:), PERMGROUP(:), BESTPERM(:), NSETS(:), SETS(:,:)
    DOUBLE PRECISION, SAVE :: ORBITTOL, GEOMDIFFTOL, BOXLX, BOXLY, BOXLZ
    LOGICAL, SAVE ::  FREEZE, PULLT, TWOD, EFIELDT, AMBERT, QCIAMBERT, AMBER12T, CHRMMT, &
        & STOCKT, CSMT, PERMDIST, LOCALPERMDIST,  LPERMDIST, OHCELLT, QCIPERMCHECK, &
        & PERMOPT, PERMINVOPT, NOINVERSION, GTHOMSONT, MKTRAPT, MULLERBROWNT

CONTAINS

SUBROUTINE SETLOGIC()

IMPLICIT NONE

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
OHCELLT = .FALSE.
QCIPERMCHECK = .FALSE.
PERMOPT = .FALSE.
PERMINVOPT = .FALSE.
NOINVERSION = .FALSE.
GTHOMSONT = .FALSE.
MKTRAPT = .FALSE.
MULLERBROWNT = .FALSE.

END SUBROUTINE SETLOGIC

SUBROUTINE INITIALISE(NEWNATOMS, NEWPERMGROUP, NEWNPERMSIZE)

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

END SUBROUTINE INITIALISE

END MODULE COMMONS

MODULE GENRIGID

      USE COMMONS, ONLY : AMBERT, AMBER12T 

      INTEGER :: NRIGIDBODY, DEGFREEDOMS, MAXSITE, NRELAXRIGIDR, NRELAXRIGIDA
      INTEGER :: XNATOMS
      INTEGER, ALLOCATABLE :: NSITEPERBODY(:), REFVECTOR(:), RIGIDSINGLES(:)
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: RIGIDGROUPS
      DOUBLE PRECISION, ALLOCATABLE :: RIGIDCOORDS(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: SITESRIGIDBODY
      DOUBLE PRECISION, ALLOCATABLE :: GR_WEIGHTS(:) ! weights for com calculation, e.g. masses
      LOGICAL :: RIGIDINIT, ATOMRIGIDCOORDT, RELAXRIGIDT
      LOGICAL :: GENRIGIDT

      DOUBLE PRECISION, ALLOCATABLE :: IINVERSE(:,:,:)

      LOGICAL :: RIGIDOPTIMROTAT, FREEZERIGIDBODYT
      DOUBLE PRECISION :: OPTIMROTAVALUES(3)
      LOGICAL :: AACONVERGENCET
      INTEGER, ALLOCATABLE :: LRBNPERMGROUP(:), LRBNPERMSIZE(:,:), LRBPERMGROUP(:,:), LRBNSETS(:,:), LRBSETS(:,:,:)

!   vr274:  added lattice coordinates
!           if HAS_LATTICE_COORDS is true, the last two atoms are treated
!           as lattice coordintes and rigidcoords is in reduced lattice units
      LOGICAL HAS_LATTICE_COORDS

!   jdf43:  RIGIDISRIGID = logical array, size NATOMS, TRUE if atom is part of
!           RB
      LOGICAL, ALLOCATABLE :: RIGIDISRIGID(:)
      INTEGER, ALLOCATABLE :: RB_BY_ATOM(:) ! sn402: records which RB an atom belongs to
!   jdf43:  FROZENRIGIDBODY = logical array, size NATOMS, TRUE if RB is frozen
      LOGICAL, ALLOCATABLE :: FROZENRIGIDBODY(:)

END MODULE GENRIGID

SUBROUTINE CONVERT(X,Y,Z,A,B,C,OVEC,H1VEC,H2VEC)
!      USE PORFUNCS
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,A,B,C,OVEC(3),H1VEC(3),H2VEC(3),PI,OL,CL,TH,PH,PS,SINA,SINB,SINC,COSA,COSB,COSC,SP12,SP13,SP22,&
                      SP23,SP32,SP33,SY12,SY22,SY32,SZ13,SZ23,SZ33,ANGLE,HOLEN,RANGLE,OHZ,HYL,HZL
      DATA ANGLE,HOLEN/104.52D0,0.9572D0/

      PI=4.0D0*ATAN(1.0D0)
      RANGLE=PI*ANGLE/360.0D0
      OHZ=HOLEN*COS(RANGLE)
      HYL=HOLEN*SIN(RANGLE)
      HZL=16.0D0*OHZ/18.0D0
      OL=-OHZ+HZL
!C
!CCCCCC    A,B,C ARE EULER ANGLES
!C
      TH=A
      PH=B
      PS=C
      SINA=SIN(TH)
      SINB=SIN(PH)
      SINC=SIN(PS)
      COSA=COS(TH)
      COSB=COS(PH)
      COSC=COS(PS)
      SP12=-(SINC*COSB+COSA*SINB*COSC)
      SP13=SINA*SINB
      SP22=-SINC*SINB+COSA*COSB*COSC
      SP23=-SINA*COSB
      SP32=SINA*COSC
      SP33=COSA
      SY12=SP12*HYL
      SY22=SP22*HYL
      SY32=SP32*HYL
      SZ13=SP13*HZL
      SZ23=SP23*HZL
      SZ33=SP33*HZL
!C
!CCCCC HYDROGEN POSITIONS
!C
      H1VEC(1)=SY12+SZ13+X
      H1VEC(2)=SY22+SZ23+Y
      H1VEC(3)=SY32+SZ33+Z

      H2VEC(1)=-SY12+SZ13+X
      H2VEC(2)=-SY22+SZ23+Y
      H2VEC(3)=-SY32+SZ33+Z
!C
!CCCC OXYGEN POSITION
!C
      OVEC(1)=SP13*OL   +X
      OVEC(2)=SP23*OL   +Y
      OVEC(3)=SP33*OL   +Z

      RETURN
END SUBROUTINE CONVERT

SUBROUTINE CONVERT2(OVEC,H1VEC,H2VEC,X,Y,Z,A,B,C)
!      USE PORFUNCS
!C
!C  Convert H, H, O positions to Eulers - needs H, H, O
!C  positions and centres of mass. Here we allow for the possibility that the
!C  ideal rigid water geometry may be broken and take the best fit.
!C
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,A,B,C,OVEC(3),H1VEC(3),H2VEC(3),PI,OL,CL,TH,PH,PS,SINA,SINB,SINC,COSA,COSB,COSC,SP12,SP13,SP22,&
                      SP23,SP32,SP33,SY12,SY22,SY32,SZ13,SZ23,SZ33,ANGLE,HOLEN,RANGLE,OHZ,HYL,HZL,SP13T,SP12T,SP22T,PI2,TEMP,&
                      SP23T,SP33T,ASAVE,BSAVE,CSAVE,SY12T,SY22T,SY32T,SZ13T,SZ23T,SZ33T,DMIN,ABEST,BBEST,CBEST
      INTEGER JA,JB,JC
      DATA ANGLE,HOLEN/104.52D0,0.9572D0/

      PI=4.0D0*ATAN(1.0D0) 
      PI2=2.0D0*PI
      RANGLE=PI*ANGLE/360.0D0 
      OHZ=HOLEN*COS(RANGLE) 
      HYL=HOLEN*SIN(RANGLE) 
      HZL=16.0D0*OHZ/18.0D0 
      OL=-OHZ+HZL 
!C
      X=(H1VEC(1)+H2VEC(1)+16.0D0*OVEC(1))/18.0D0
      Y=(H1VEC(2)+H2VEC(2)+16.0D0*OVEC(2))/18.0D0
      Z=(H1VEC(3)+H2VEC(3)+16.0D0*OVEC(3))/18.0D0


      COSA=((OVEC(3)-Z)/OL)
      IF (ABS(COSA).GT.1.0D0) THEN
         COSA=ABS(COSA)/COSA
         A=DACOS(-COSA)
      ELSE
         A=DACOS((OVEC(3)-Z)/OL)
      ENDIF
      SINA=DSIN(A)
      IF (SINA.EQ.0.0D0) SINA=1.0D-10

      COSB=-(OVEC(2)-Y)/(OL*SINA)
      IF (ABS(COSB).GT.1.0D0) THEN
         COSB=ABS(COSB)/COSB
         B=DACOS(-COSB)
      ELSE
         B=DACOS(-(OVEC(2)-Y)/(OL*SINA))
      ENDIF
      SINB=DSIN(B)

      COSC=(H1VEC(3)-H2VEC(3))/(2.0D0*HYL*SINA)
      IF (ABS(COSC).GT.1.0D0) THEN
         COSC=ABS(COSC)/COSC
         C=DACOS(-COSC)
      ELSE
         C=DACOS((H1VEC(3)-H2VEC(3))/(2.0D0*HYL*SINA))
      ENDIF
      SINC=DSIN(C)

      SP13T=(OVEC(1)-X)/OL
      SP23T=(OVEC(2)-Y)/OL
      SP33T=(OVEC(3)-Z)/OL

      SY12T=(H1VEC(1)-H2VEC(1))/(2.0D0)
      SY22T=(H1VEC(2)-H2VEC(2))/(2.0D0)
      SY32T=(H1VEC(3)-H2VEC(3))/(2.0D0)

      SZ13T=(H1VEC(1)+H2VEC(1)-2.0D0*X)/(2.0D0)
      SZ23T=(H1VEC(2)+H2VEC(2)-2.0D0*Y)/(2.0D0)
      SZ33T=(H1VEC(3)+H2VEC(3)-2.0D0*Z)/(2.0D0)
 
      ASAVE=A
      BSAVE=B
      CSAVE=C

      DMIN=1.0D100
      DO JA=1,4
         IF (JA.EQ.1) A=ASAVE
         IF (JA.EQ.2) A=PI2-ASAVE
         IF (JA.EQ.3) A=PI-ASAVE
         IF (JA.EQ.4) A=PI+ASAVE
         SINA=SIN(A)
         COSA=COS(A)
         DO JB=1,4
            IF (JB.EQ.1) B=BSAVE
            IF (JB.EQ.2) B=PI2-BSAVE
            IF (JB.EQ.3) B=PI-BSAVE
            IF (JB.EQ.4) B=PI+BSAVE
            SINB=SIN(B)
            COSB=COS(B)
            DO JC=1,4
               IF (JC.EQ.1) C=CSAVE
               IF (JC.EQ.2) C=PI2-CSAVE
               IF (JC.EQ.3) C=PI-CSAVE
               IF (JC.EQ.4) C=PI+CSAVE
               SINC=SIN(C)
               COSC=COS(C)

               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP13=SINA*SINB
               SP22=-SINC*SINB+COSA*COSB*COSC
               SP23=-SINA*COSB
               SP32=SINA*COSC
               SP33=COSA
               SY12=SP12*HYL
               SY22=SP22*HYL
               SY32=SP32*HYL
               SZ13=SP13*HZL
               SZ23=SP23*HZL
               SZ33=SP33*HZL

               IF ( DABS(SP13T-SP13)+DABS(SP23T-SP23)+DABS(SP33T-SP33)+DABS(SY12T-SY12)+&
                    DABS(SY22T-SY22)+DABS(SY32T-SY32)+DABS(SZ13T-SZ13)+DABS(SZ23T-SZ23)+&
                    DABS(SZ33T-SZ33) .LT. DMIN) THEN
                   DMIN= DABS(SP13T-SP13)+DABS(SP23T-SP23)+DABS(SP33T-SP33)+DABS(SY12T-SY12)+&
                    DABS(SY22T-SY22)+DABS(SY32T-SY32)+DABS(SZ13T-SZ13)+DABS(SZ23T-SZ23)+DABS(SZ33T-SZ33) 
                  ABEST=A
                  BBEST=B
                  CBEST=C
                  IF (DMIN.LT.1.0D-10) GOTO 20
               ENDIF
            ENDDO
         ENDDO
      ENDDO

20    A=ABEST
      B=BBEST
      C=CBEST
!C     IF (DMIN.GT.0.1D0) WRITE(*,'(A,F15.5)') 'WARNING, deviation from rigid body geometry detected, best fit is ',DMIN

      RETURN
END SUBROUTINE CONVERT2

!MODULE PORFUNCS
!     implicit none
!     contains
!          subroutine getarg_subr(position,value) ! wraps getarg function so it can be use-associated
!               implicit none
! 
!               integer,intent(in) :: position
!               character(len=*),intent(out) :: value
! 
! 
!!               call getarg(position,value)
!          end subroutine getarg_subr
! 
!          subroutine iargc_subr(n) ! wraps iargc function so it can be use-associated
!               implicit none
!               integer,intent(out) :: n
! 
!               integer iargc
! 
!!               n = iargc()
!          end subroutine iargc_subr
! 
!          subroutine fork_subr(pid) ! returns zero in the child process, PID of child in parent process
!               implicit none
!               integer, intent(inout) :: pid
! 
!               integer fork
! 
!!               pid=fork()
!          end subroutine fork_subr
! 
!          subroutine system_subr(JobString,ExitStatus)
!               implicit none
! 
!               character(len=*),intent(in) :: JobString
!               integer,intent(out) :: ExitStatus
! 
!               integer system
! 
!!               ExitStatus=system(JobString)
!!               ExitStatus=ishft(ExitStatus,-8)
!          end subroutine system_subr
! 
!          subroutine wait_subr(pid,ExitStatus)
!               implicit none
! 
!               integer,intent(inout) :: pid,ExitStatus
! 
!               integer wait
! 
!!               pid=wait(ExitStatus)
!!               ExitStatus=ishft(ExitStatus,-8)
!          end subroutine wait_subr
! 
!          subroutine getpid_subr(pid)
!               implicit none
! 
!               integer,intent(out) :: pid
! 
!               integer getpid
! 
!!               pid=getpid()
!          end subroutine getpid_subr
!END MODULE PORFUNCS
