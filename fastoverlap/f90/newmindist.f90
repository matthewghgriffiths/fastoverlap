!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!   Finds the minimum distance between two geometries.
!   Geometry in RA should not change. RB is returned as the
!   closest geometry to RA if PRESERVET is .FALSE.
!
!   New analytic method based on quaterions from
!   Kearsley, Acta Cryst. A, 45, 208-210, 1989.
!
! jmc As long as zsym isn't 'W' (in which case mind does something special) mind
! doesn't care what atomic symbol we give it.
!
!     ----------------------------------------------------------------------------------------------
! jdf43>        Modified for general angle-axis 30/01/12
!     ----------------------------------------------------------------------------------------------

SUBROUTINE NEWMINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET,RIGIDBODY,DEBUG,RMAT)
USE COMMONS,ONLY : MYUNIT, MULLERBROWNT, BOXLX, BOXLY, BOXLZ, STOCKT, CSMT, MKTRAPT
IMPLICIT NONE
INTEGER J1, NATOMS, NSIZE, INFO, JINFO, JMIN
INTEGER,PARAMETER :: LWORK=12
DOUBLE PRECISION RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), XM, YM, ZM, XP, YP, ZP, OVEC(3), H1VEC(3), H2VEC(3), &
  &              DIAG(4), TEMPA(LWORK), RMAT(3,3), MINV, Q1, Q2, Q3, Q4, CMXA, CMYA, CMZA, CMXB, CMYB, CMZB, &
  &              NCMXB, NCMYB, NCMZB
DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
LOGICAL BULKT, TWOD, RIGIDBODY, PRESERVET, DEBUG
CHARACTER(LEN=5) ZUSE
INTEGER NCIT
DOUBLE PRECISION XSHIFT, YSHIFT, ZSHIFT, XSHIFTNEW, YSHIFTNEW, ZSHIFTNEW
DOUBLE PRECISION MYROTMAT(3,3), OMEGATOT(3,3)
COMMON /MINDOM/ MYROTMAT, OMEGATOT

! mg542 not implementing RIGIDBODY
!IF (RIGIDBODY) THEN
!   CALL RBMINDIST(RA,RB,NATOMS,DIST,RMAT,DEBUG)
!   RETURN
!ELSEIF (MKTRAPT) THEN
IF (MKTRAPT) THEN
   DIST=0.0D0
   DO J1=1,3*NATOMS
      DIST=DIST+(RA(J1)-RB(J1))**2
   ENDDO
   DIST=SQRT(DIST)
   RETURN
!
! Convert rigid body coordinates to Cartesians for rigid bodies.
!
ELSE IF (ZUSE(1:1).EQ.'W') THEN
   ALLOCATE(XA(3*3*(NATOMS/2)),XB(3*3*(NATOMS/2)))
   NSIZE=3*(NATOMS/2)
   DO J1=1,NATOMS/2
      CALL CONVERT(RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3), &
     &        RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
      XA(3*(J1-1)+1+1)=OVEC(1)
      XA(3*(J1-1)+1+2)=OVEC(2)
      XA(3*(J1-1)+1+3)=OVEC(3)
      XA(3*(J1-1)+2+1)=H1VEC(1)
      XA(3*(J1-1)+2+2)=H1VEC(2)
      XA(3*(J1-1)+2+3)=H1VEC(3)
      XA(3*(J1-1)+3+1)=H2VEC(1)
      XA(3*(J1-1)+3+2)=H2VEC(2)
      XA(3*(J1-1)+3+3)=H2VEC(3)
      CALL CONVERT(RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3), &
     &      RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
      XB(3*(J1-1)+1+1)=OVEC(1)
      XB(3*(J1-1)+1+2)=OVEC(2)
      XB(3*(J1-1)+1+3)=OVEC(3)
      XB(3*(J1-1)+2+1)=H1VEC(1)
      XB(3*(J1-1)+2+2)=H1VEC(2)
      XB(3*(J1-1)+2+3)=H1VEC(3)
      XB(3*(J1-1)+3+1)=H2VEC(1)
      XB(3*(J1-1)+3+2)=H2VEC(2)
      XB(3*(J1-1)+3+3)=H2VEC(3)
   ENDDO
!ELSEIF (RIGIDBODY) THEN
!   WRITE(MYUNIT,'(A)') 'newmindist> New quaterion procedure not yet coded for general angle-axis variables'
!   STOP
ELSEIF (STOCKT) THEN
   ALLOCATE(XA(3*(NATOMS/2)),XB(3*(NATOMS/2)))
   NSIZE=(NATOMS/2)
   XA(1:3*NSIZE)=RA(1:3*NSIZE)
   XB(1:3*NSIZE)=RB(1:3*NSIZE)
ELSEIF (MULLERBROWNT) THEN
   DIST=SQRT((RA(1)-RB(1))**2+(RA(2)-RB(2))**2)
   RETURN
ELSEIF (CSMT) THEN
   DIST=0.0D0
   DO J1=1,NATOMS
      DIST=DIST + (RA(3*(J1-1)+1)-RB(3*(J1-1)+1))**2 &
   &            + (RA(3*(J1-1)+2)-RB(3*(J1-1)+2))**2 &
   &            + (RA(3*(J1-1)+3)-RB(3*(J1-1)+3))**2
   ENDDO
   DIST=SQRT(DIST)
   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
   RETURN
ELSE
   ALLOCATE(XA(3*NATOMS),XB(3*NATOMS))
   NSIZE=NATOMS
   XA(1:3*NATOMS)=RA(1:3*NATOMS)
   XB(1:3*NATOMS)=RB(1:3*NATOMS)
ENDIF

!
! Move centre of coordinates of XA and XB to the origin.
!
!CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
!DO J1=1,NSIZE
!   CMXA=CMXA+XA(3*(J1-1)+1)
!   CMYA=CMYA+XA(3*(J1-1)+2)
!   CMZA=CMZA+XA(3*(J1-1)+3)
!ENDDO
!CMXA=CMXA/NSIZE; CMYA=CMYA/NSIZE; CMZA=CMZA/NSIZE
!DO J1=1,NSIZE
!   XA(3*(J1-1)+1)=XA(3*(J1-1)+1)-CMXA
!   XA(3*(J1-1)+2)=XA(3*(J1-1)+2)-CMYA
!   XA(3*(J1-1)+3)=XA(3*(J1-1)+3)-CMZA
!ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mg542!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mg542 superimposing COM means adding relative displacement meaningless
CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
CMXB=0.0D0; CMYB=0.0D0; CMZB=0.0D0
IF(BULKT) THEN
   XA(1:3*NATOMS)=RA(1:3*NATOMS)
   XB(1:3*NATOMS)=RB(1:3*NATOMS)
   DO J1=1,NSIZE
      XB(3*J1-2)=XB(3*J1-2) + BOXLX*NINT((XA(3*J1-2)-XB(3*J1-2))/BOXLX)
      XB(3*J1-1)=XB(3*J1-1) + BOXLY*NINT((XA(3*J1-1)-XB(3*J1-1))/BOXLY)
      XB(3*J1  )=XB(3*J1  ) + BOXLZ*NINT((XA(3*J1  )-XB(3*J1  ))/BOXLZ)
   ENDDO
ELSE
   !
   ! Move centre of coordinates of XA and XB to the origin.
   !
   CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
   DO J1=1,NSIZE
      CMXA=CMXA+XA(3*(J1-1)+1)
      CMYA=CMYA+XA(3*(J1-1)+2)
      CMZA=CMZA+XA(3*(J1-1)+3)
   ENDDO
   CMXA=CMXA/NSIZE; CMYA=CMYA/NSIZE; CMZA=CMZA/NSIZE
   DO J1=1,NSIZE
      XA(3*(J1-1)+1)=XA(3*(J1-1)+1)-CMXA
      XA(3*(J1-1)+2)=XA(3*(J1-1)+2)-CMYA
      XA(3*(J1-1)+3)=XA(3*(J1-1)+3)-CMZA
   ENDDO
      CMXB=0.0D0; CMYB=0.0D0; CMZB=0.0D0
      DO J1=1,NSIZE
         CMXB=CMXB+XB(3*(J1-1)+1)
         CMYB=CMYB+XB(3*(J1-1)+2)
         CMZB=CMZB+XB(3*(J1-1)+3)
      ENDDO
      CMXB=CMXB/NSIZE; CMYB=CMYB/NSIZE; CMZB=CMZB/NSIZE
      DO J1=1,NSIZE
         XB(3*(J1-1)+1)=XB(3*(J1-1)+1)-CMXB
         XB(3*(J1-1)+2)=XB(3*(J1-1)+2)-CMYB
         XB(3*(J1-1)+3)=XB(3*(J1-1)+3)-CMZB
      ENDDO
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mg542!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!  MG542 not implementing TWOD for time being.
!! Deal with TWO2 non-bulk here.
!!
!IF (TWOD.AND.(.NOT.BULKT)) THEN
!   CALL MINDIST(XA,XB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
!   RMAT(1:3,1:3)=OMEGATOT(1:3,1:3)
!!
!!  To align the structures for minpermdist we use the orientation in XA
!!  and translate the B coordinates to the centre of coordinates of RA.
!!
!   IF (.NOT.PRESERVET) THEN
!      DO J1=1,NATOMS
!         RB(3*(J1-1)+1)=XB(3*(J1-1)+1)-CMXB+CMXA+XSHIFT
!         RB(3*(J1-1)+2)=XB(3*(J1-1)+2)-CMYB+CMYA+YSHIFT
!         RB(3*(J1-1)+3)=XB(3*(J1-1)+3)-CMZB+CMZA+ZSHIFT
!      ENDDO
!   ENDIF
!   RETURN
!ENDIF


XSHIFT=0.0D0; YSHIFT=0.0D0; ZSHIFT=0.0D0
NCIT=0
IF (BULKT) THEN
! 1  NCIT=NCIT+1
!    IF (NCIT.GT.1000) THEN
!       PRINT '(A)','inertia> WARNING - iterative calculation of centre of mass shift did not converge'
!    ENDIF
!    XSHIFTNEW=0.0D0
!    YSHIFTNEW=0.0D0
!    ZSHIFTNEW=0.0D0
!    DO J1=1,NSIZE
!       XSHIFTNEW=XSHIFTNEW + XA(3*(J1-1)+1)-XB(3*(J1-1)+1) - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT)/BOXLX)
!       YSHIFTNEW=YSHIFTNEW + XA(3*(J1-1)+2)-XB(3*(J1-1)+2) - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT)/BOXLY)
!       ZSHIFTNEW=ZSHIFTNEW + XA(3*(J1-1)+3)-XB(3*(J1-1)+3) - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT)/BOXLZ)
!    ENDDO
!    XSHIFTNEW=XSHIFTNEW/NSIZE; YSHIFTNEW=YSHIFTNEW/NSIZE; ZSHIFTNEW=ZSHIFTNEW/NSIZE
!    IF ((ABS(XSHIFTNEW-XSHIFT).GT.1.0D-6).OR.(ABS(YSHIFTNEW-YSHIFT).GT.1.0D-6).OR.(ABS(ZSHIFTNEW-ZSHIFT).GT.1.0D-6)) THEN
!       IF (DEBUG) PRINT '(A,I6,6F15.7)',' newmindist> ',NCIT,XSHIFTNEW,YSHIFTNEW,ZSHIFTNEW,XSHIFT,YSHIFT,ZSHIFT
!       XSHIFT=0.05D0*XSHIFT+0.95D0*XSHIFTNEW
!       YSHIFT=0.05D0*YSHIFT+0.95D0*YSHIFTNEW
!       ZSHIFT=0.05D0*ZSHIFT+0.95D0*ZSHIFTNEW
!       GOTO 1
!    ENDIF
!    IF (DEBUG) PRINT '(A,I6,3F15.7)',' newmindist> coordinate shift converged. Cycles and values: ',NCIT,XSHIFT,YSHIFT,ZSHIFT
!    XSHIFT=XSHIFTNEW
!    YSHIFT=YSHIFTNEW
!    ZSHIFT=ZSHIFTNEW
!
! Actually, the iterative solution seems to be worse than simply putting the centre of mass
! at the origin.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mg542!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO J1=1,NSIZE
      XSHIFT = XSHIFT + XA(3*(J1-1)+1)-XB(3*(J1-1)+1)
      YSHIFT = YSHIFT + XA(3*(J1-1)+2)-XB(3*(J1-1)+2)
      ZSHIFT = ZSHIFT + XA(3*(J1-1)+3)-XB(3*(J1-1)+3)
   ENDDO
   XSHIFT = XSHIFT/NSIZE; YSHIFT = YSHIFT/NSIZE; ZSHIFT = ZSHIFT/NSIZE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!mg542!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DIST=0.0D0
   DO J1=1,NSIZE
      DIST=DIST + (XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT)/BOXLX))**2 &
   &            + (XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT)/BOXLY))**2 &
   &            + (XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT)/BOXLZ))**2
   ENDDO

   DIST=SQRT(DIST)
   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
ELSE
!
!  The formula below is not invariant to overall translation because XP, YP, ZP
!  involve a sum of coordinates! We need to have XA and XB coordinate centres both
!  at the origin!!
!
   QMAT(1:4,1:4)=0.0D0
!  PRINT *,'XA:'
!  PRINT '(6G20.10)',XA(1:3*NATOMS)
!  PRINT *,'XB:'
!  PRINT '(6G20.10)',XB(1:3*NATOMS)
   DO J1=1,NSIZE
      XM=XA(3*(J1-1)+1)-XB(3*(J1-1)+1)
      YM=XA(3*(J1-1)+2)-XB(3*(J1-1)+2)
      ZM=XA(3*(J1-1)+3)-XB(3*(J1-1)+3)
      XP=XA(3*(J1-1)+1)+XB(3*(J1-1)+1)
      YP=XA(3*(J1-1)+2)+XB(3*(J1-1)+2)
      ZP=XA(3*(J1-1)+3)+XB(3*(J1-1)+3)
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
   QMAT(2,1)=QMAT(1,2); QMAT(3,1)=QMAT(1,3); QMAT(3,2)=QMAT(2,3); QMAT(4,1)=QMAT(1,4); QMAT(4,2)=QMAT(2,4); QMAT(4,3)=QMAT(3,4)

   CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,LWORK,INFO)
   IF (INFO.NE.0) WRITE(MYUNIT,'(A,I6,A)') 'newmindist> WARNING - INFO=',INFO,' in DSYEV'

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
         WRITE(MYUNIT,'(A,G20.10,A)') 'newmindist> WARNING MINV is ',MINV,' change to absolute value'
         MINV=-MINV
      ENDIF
   ENDIF
   DIST=SQRT(MINV)

   IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,I6)') 'newmindist> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
   Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
!
! RMAT will contain the matrix that maps RB onto the best correspondence with RA
!
   RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
   RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
   RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
   RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
   RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
   RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
   RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
   RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
   RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2
ENDIF

!
!  Needs some thought for the angle/axis rigid body formulation.
!

IF (.NOT.PRESERVET) THEN
   IF (ZUSE(1:1).EQ.'W') THEN
      DO J1=1,NATOMS/2
         OVEC(1)=XB(1+(J1-1)*3+1)
         OVEC(2)=XB(2+(J1-1)*3+1)
         OVEC(3)=XB(3+(J1-1)*3+1)
         H1VEC(1)=XB(1+(J1-1)*3+2)
         H1VEC(2)=XB(2+(J1-1)*3+2)
         H1VEC(3)=XB(3+(J1-1)*3+2)
         H2VEC(1)=XB(1+(J1-1)*3+3)
         H2VEC(2)=XB(2+(J1-1)*3+3)
         H2VEC(3)=XB(3+(J1-1)*3+3)
         CALL CONVERT2(OVEC,H1VEC,H2VEC,RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3), &
  &                    RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3))
      ENDDO
   ELSEIF (RIGIDBODY) THEN
      WRITE(MYUNIT,'(A)') 'newmindist> back transformation not programmed yet for rigid bodies'
   ENDIF
!
!  Translate the RB coordinates to the centre of coordinates of RA.
!


      DO J1=1,NATOMS
         RB(3*(J1-1)+1)=RB(3*(J1-1)+1)-CMXB+CMXA+XSHIFT
         RB(3*(J1-1)+2)=RB(3*(J1-1)+2)-CMYB+CMYA+YSHIFT
         RB(3*(J1-1)+3)=RB(3*(J1-1)+3)-CMZB+CMZA+ZSHIFT
      ENDDO
      IF (STOCKT) THEN
         CALL NEWROTGEOMSTOCK(NATOMS,RB,RMAT,CMXA,CMYA,CMZA)
      ELSE
         CALL NEWROTGEOM(NSIZE,RB,RMAT,CMXA,CMYA,CMZA)
      ENDIF
ENDIF
! WRITE(MYUNIT,'(A,6G20.10)') 'CMA,CMB=',CMXA,CMYA,CMZA,CMXB,CMYB,CMZB

DEALLOCATE(XA,XB)

END SUBROUTINE NEWMINDIST

SUBROUTINE NEWROTGEOM(NATOMS,COORDS,ROTMAT,CX,CY,CZ)
IMPLICIT NONE
INTEGER I, J, K, NATOMS
DOUBLE PRECISION COORDS(*), R1, R0(3), ROTMAT(3,3), CX, CY, CZ

DO I=1,NATOMS
   R0(1)=COORDS(3*(I-1)+1)-CX
   R0(2)=COORDS(3*(I-1)+2)-CY
   R0(3)=COORDS(3*(I-1)+3)-CZ
   DO J=1,3
      R1=0.0D0
      DO K=1,3
         R1=R1+ROTMAT(J,K)*R0(K)
      ENDDO
      IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
      IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
      IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
   ENDDO
ENDDO

RETURN
END SUBROUTINE NEWROTGEOM

SUBROUTINE NEWROTGEOMSTOCK(NATOMS,COORDS,ROTMAT,CX,CY,CZ)
IMPLICIT NONE
INTEGER I, J, K, NATOMS, REALNATOMS, J3, J1, OFFSET
DOUBLE PRECISION COORDS(*), R1, R0(3), ROTMAT(3,3), CX, CY, CZ, X1, Y1, Z1, X2, Y2, Z2, CT1, ST1, P1, CP1, SP1, T1B, P1B, T1
DOUBLE PRECISION START(3), FINISH(3), DIFF, DIFFBEST
DOUBLE PRECISION, PARAMETER ::  PI=3.141592654D0

REALNATOMS=(NATOMS/2)
OFFSET = 3*REALNATOMS
!
! First rotate the dipoles.
!
DO J1=1,REALNATOMS
   J3=3*J1
   X1=COORDS(J3-2)-CX
   Y1=COORDS(J3-1)-CY
   Z1=COORDS(J3)-CZ
   T1=COORDS(OFFSET+J3-2)
   CT1=COS(T1)
   ST1=SIN(T1)
   P1=COORDS(OFFSET+J3-1)
   CP1=COS(P1)
   SP1=SIN(P1)
   X2=X1+ST1*CP1
   Y2=Y1+ST1*SP1
   Z2=Z1+CT1
   START(1)=ROTMAT(1,1)*X1+ROTMAT(1,2)*Y1+ROTMAT(1,3)*Z1
   START(2)=ROTMAT(2,1)*X1+ROTMAT(2,2)*Y1+ROTMAT(2,3)*Z1
   START(3)=ROTMAT(3,1)*X1+ROTMAT(3,2)*Y1+ROTMAT(3,3)*Z1
   FINISH(1)=ROTMAT(1,1)*X2+ROTMAT(1,2)*Y2+ROTMAT(1,3)*Z2
   FINISH(2)=ROTMAT(2,1)*X2+ROTMAT(2,2)*Y2+ROTMAT(2,3)*Z2
   FINISH(3)=ROTMAT(3,1)*X2+ROTMAT(3,2)*Y2+ROTMAT(3,3)*Z2
   T1=ACOS(FINISH(3)-START(3))
   DIFFBEST=1.0D100
   IF (SIN(T1).NE.0.0D0) THEN
      P1=ACOS((FINISH(1)-START(1))/SIN(T1))
   ELSE
      P1=1.0D0
   ENDIF
   DIFFBEST=(FINISH(1)-START(1)-SIN(T1)*COS(P1))**2+(FINISH(2)-START(2)-SIN(T1)*SIN(P1))**2+(FINISH(3)-START(3)-COS(T1))**2
   T1B=T1; P1B=P1
   DIFF=(FINISH(1)-START(1)-SIN(2*PI-T1)*COS(P1))**2+(FINISH(2)-START(2)-SIN(2*PI-T1)*SIN(P1))**2+ &
  &                (FINISH(3)-START(3)-COS(2*PI-T1))**2
   IF (DIFF.LT.DIFFBEST) THEN
      T1B=2*PI-T1; P1B=P1
      DIFFBEST=DIFF
   ENDIF
   DIFF=(FINISH(1)-START(1)-SIN(2*PI-T1)*COS(2*PI-P1))**2+(FINISH(2)-START(2)-SIN(2*PI-T1)*SIN(2*PI-P1))**2+ &
  &                (FINISH(3)-START(3)-COS(2*PI-T1))**2
   IF (DIFF.LT.DIFFBEST) THEN
      T1B=2*PI-T1; P1B=2*PI-P1
      DIFFBEST=DIFF
   ENDIF
   DIFF=(FINISH(1)-START(1)-SIN(T1)*COS(2*PI-P1))**2+(FINISH(2)-START(2)-SIN(T1)*SIN(2*PI-P1))**2+ &
  &                (FINISH(3)-START(3)-COS(T1))**2
   IF (DIFF.LT.DIFFBEST) THEN
      T1B=T1; P1B=2*PI-P1
      DIFFBEST=DIFF
   ENDIF
   IF (DIFFBEST.GT.1.0D-10) THEN
      PRINT '(A,G20.10)','newrotgeomstock> WARNING - angle rotation failed - DIFFBEST=',DIFFBEST
   ENDIF
!  PRINT '(A,G20.10)','newrotgeomstock> DIFFBEST=',DIFFBEST
!
!  Inverse cos gives us an angle between 0 and pi. However, 2*pi - angle gives
!  the same cos. There are therefore two possibilities for theta and two for phi,
!  and only one should regenerate the correct displacements. Find it!
!
   COORDS(OFFSET+J3-2)=T1B; COORDS(OFFSET+J3-1)=P1B
ENDDO

DO I=1,(NATOMS/2)
   R0(1)=COORDS(3*(I-1)+1)-CX
   R0(2)=COORDS(3*(I-1)+2)-CY
   R0(3)=COORDS(3*(I-1)+3)-CZ
   DO J=1,3
      R1=0.0D0
      DO K=1,3
         R1=R1+ROTMAT(J,K)*R0(K)
      ENDDO
      IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
      IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
      IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
   ENDDO
ENDDO

RETURN
END SUBROUTINE NEWROTGEOMSTOCK

