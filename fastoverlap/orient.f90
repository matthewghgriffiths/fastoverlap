
!C   Copyright (C) 1999-2006 David J. Wales
!C   This file is part of GMIN.
!C
!C   GMIN is free software; you can redistribute it and/or modify
!C   it under the terms of the GNU General Public License as published by
!C   the Free Software Foundation; either version 2 of the License, or
!C   (at your option) any later version.
!C
!C   GMIN is distributed in the hope that it will be useful,
!C   but WITHOUT ANY WARRANTY; without even the implied warranty of
!C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!C   GNU General Public License for more details.
!C
!C   You should have received a copy of the GNU General Public License
!C   along with this program; if not, write to the Free Software
!C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!C
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  Put permutational isomers into a standard orientation.
!C
      SUBROUTINE MYORIENT(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,DEBUG,ROT,ROTINV,STOCKT)
      USE COMMONS, ONLY : ORBITTOL, TWOD, PULLT, EFIELDT
      IMPLICIT NONE

!     mg542 modifying to work with f2py
      INTEGER, INTENT(IN) :: NATOMS, NCHOOSE1, NCHOOSE2
      INTEGER, INTENT(OUT) :: NORBIT1, NORBIT2
      DOUBLE PRECISION, INTENT(INOUT) :: Q1(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: T1(3*NATOMS), ROT(3,3), ROTINV(3,3)
      LOGICAL, INTENT(IN) :: DEBUG, STOCKT


      INTEGER J1, I, JMAX1, JMAX2, J2, J3
      DOUBLE PRECISION DIST(NATOMS), DMAX, RVEC(3), CMX, CMY, CMZ, T2(3*NATOMS), &
                       COST, SINT, RDOTN, TX, TY, TZ, IT(3,3), IV(3,3), PROJ, DMAX2, CUTOFF1
      DOUBLE PRECISION XVEC(3), YVEC(3), ZVEC(3), TEMP(3,3), TVEC(3*NATOMS), DUMMY
      DOUBLE PRECISION TXVEC(3), TYVEC(3), TZVEC(3)
      LOGICAL RBAAT
      DOUBLE PRECISION ENERGY, RMS, VNEW(3*NATOMS), DTEMP

!      INTEGER J1, I, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, J2, J3
!      DOUBLE PRECISION Q1(3*NATOMS), DIST(NATOMS), DMAX, RVEC(3), T1(3*NATOMS), CMX, CMY, CMZ, T2(3*NATOMS), &
!                      COST, SINT, RDOTN, TX, TY, TZ, IT(3,3), IV(3,3), PROJ, DMAX2, CUTOFF1
!      DOUBLE PRECISION ROT(3,3), XVEC(3), YVEC(3), ZVEC(3), ROTINV(3,3), TEMP(3,3), TVEC(3*NATOMS), DUMMY
!      DOUBLE PRECISION TXVEC(3), TYVEC(3), TZVEC(3)
!      LOGICAL DEBUG, STOCKT, RBAAT
!      DOUBLE PRECISION ENERGY, RMS, VNEW(3*NATOMS), DTEMP
!     mg542

!      IF (RBAAT) THEN

!         CALL RBORIENT(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,DEBUG,ROT,ROTINV)
!         RETURN

!      ENDIF
 
      CUTOFF1=ORBITTOL
      XVEC(1)=1.0D0; XVEC(2)=0.0D0; XVEC(3)=0.0D0
      YVEC(1)=0.0D0; YVEC(2)=1.0D0; YVEC(3)=0.0D0
      ZVEC(1)=0.0D0; ZVEC(2)=0.0D0; ZVEC(3)=1.0D0
!C
!C  Move centre of mass to the origin.
!C
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+Q1(3*(I-1)+1)
         CMY=CMY+Q1(3*(I-1)+2)
         CMZ=CMZ+Q1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
!C     PRINT*,'CMX,CMY,CMZ=',CMX,CMY,CMZ
      DO I=1,NATOMS
         Q1(3*(I-1)+1)=Q1(3*(I-1)+1)-CMX
         Q1(3*(I-1)+2)=Q1(3*(I-1)+2)-CMY
         Q1(3*(I-1)+3)=Q1(3*(I-1)+3)-CMZ
      ENDDO

      DMAX=-1.0D0
      NORBIT1=1
!
! We can only rotate around the z axis with a pulling potential!
! And for TWOD!
! GMIN does not have the GTHOMSONT logical, so this one is omitted here.
!
! hk286
      IF (PULLT.OR.EFIELDT.OR.TWOD) THEN
         T1(1:3*NATOMS)=Q1(1:3*NATOMS)
         GOTO 10
      ENDIF

      DO J1=1,NATOMS
         DIST(J1)=SQRT(Q1(3*(J1-1)+1)**2+Q1(3*(J1-1)+2)**2+Q1(3*(J1-1)+3)**2)
!        PRINT '(A,I8,4F20.10)','J1,DMAX,DIST(J1),ABS(DIST(J1)-DMAX),CUTOFF1=',J1,DMAX,DIST(J1),ABS(DIST(J1)-DMAX),CUTOFF1
         IF (ABS(DIST(J1)-DMAX).LT.CUTOFF1) THEN
            NORBIT1=NORBIT1+1
            IF (NORBIT1.EQ.NCHOOSE1) THEN
               JMAX1=J1
            ENDIF
!           IF (DIST(J1).GT.DMAX) THEN ! try to catch the case where 0.001 isn;t small enough!
!              DMAX=DIST(J1)
!              JMAX1=J1
!           ENDIF
         ELSE IF (DIST(J1).GT.DMAX) THEN
            DMAX=DIST(J1)
            NORBIT1=1
            JMAX1=J1
         ENDIF
!        PRINT '(A,I8,2G20.10,3I8)','J1,DIST,DMAX,JMAX1,NCHOOSE1,NORBIT1=',J1,DIST(J1),DMAX,JMAX1,NCHOOSE1,NORBIT1
      ENDDO
!C
!C  For tagged atoms the choice of the first atom matters if it belongs to an orbit of size > 1.
!C
!     WRITE(*,'(A,G20.10,A,G20.10)') 'atom ',JMAX1,' will be moved onto the z axis, distance=',DIST(JMAX1) 
      IF ((ABS(Q1(3*(JMAX1-1)+1)).LT.1.0D-8).AND.(ABS(Q1(3*(JMAX1-1)+2)).LT.1.0D-8)) THEN
         IF (Q1(3*(JMAX1-1)+3).GT.0.0D0) THEN
            T1(1:3*NATOMS)=Q1(1:3*NATOMS)
         ELSE  ! rotate about the x axis DO NOT INVERT!!
            DO J1=1,NATOMS
               T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)
               T1(3*(J1-1)+2)=-Q1(3*(J1-1)+2)
               T1(3*(J1-1)+3)=-Q1(3*(J1-1)+3)
            ENDDO
            YVEC(2)=-1.0D0
            ZVEC(3)=-1.0D0
         ENDIF
         GOTO 10
      ENDIF
!
! For sloppy cutoffs we cannot assume that DMAX is exactly the same for members of the same orbit!!
!
!     COST=Q1(3*(JMAX1-1)+3)/DMAX
!     SINT=SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)/DMAX
!
      COST=Q1(3*(JMAX1-1)+3)/DIST(JMAX1)
      SINT=SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)/DIST(JMAX1)
!C
!C  Rotate atom JMAX1 onto the z axis.
!C  Rotate all the atoms through ANGLE about RVEC. Use rotation formula
!C  from Goldstein p. 165.
!C
      RVEC(1)= Q1(3*(JMAX1-1)+2)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
      RVEC(2)=-Q1(3*(JMAX1-1)+1)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
      RVEC(3)=0.0D0
!C
      DO J1=1,NATOMS
!C        IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=Q1(3*(J1-1)+1)*RVEC(1)+Q1(3*(J1-1)+2)*RVEC(2)+Q1(3*(J1-1)+3)*RVEC(3)
            T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)-(Q1(3*(J1-1)+2)*RVEC(3)-Q1(3*(J1-1)+3)*RVEC(2))*SINT
            T1(3*(J1-1)+2)=Q1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)-(Q1(3*(J1-1)+3)*RVEC(1)-Q1(3*(J1-1)+1)*RVEC(3))*SINT
            T1(3*(J1-1)+3)=Q1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)-(Q1(3*(J1-1)+1)*RVEC(2)-Q1(3*(J1-1)+2)*RVEC(1))*SINT
!C        ENDIF
      ENDDO
!C
!C
!C  Track transformation of unit vectors.
!C
      RDOTN=RVEC(1)
      XVEC(1)=COST + RVEC(1)*RDOTN*(1-COST)
      XVEC(2)=       RVEC(2)*RDOTN*(1-COST)+RVEC(3)*SINT
      XVEC(3)=       RVEC(3)*RDOTN*(1-COST)-RVEC(2)*SINT

      RDOTN=RVEC(2)
      YVEC(1)=       RVEC(1)*RDOTN*(1-COST)-RVEC(3)*SINT
      YVEC(2)=COST + RVEC(2)*RDOTN*(1-COST)
      YVEC(3)=       RVEC(3)*RDOTN*(1-COST)+RVEC(1)*SINT

      RDOTN=RVEC(3)
      ZVEC(1)=       RVEC(1)*RDOTN*(1-COST)+RVEC(2)*SINT
      ZVEC(2)=       RVEC(2)*RDOTN*(1-COST)-RVEC(1)*SINT
      ZVEC(3)=COST + RVEC(3)*RDOTN*(1-COST)

!C     DO J1=1,3*NATOMS
!C        Q1(J1)=T1(J1)
!C     ENDDO

10    CONTINUE
!C     IF (DEBUG) THEN
!C        WRITE(4,'(I5)') NATOMS
!C        WRITE(4,'(A)') 'after first myorient rotation:'
!C        WRITE(4,'(A5,3F20.10)') 'LA   ',T1(1),T1(2),T1(3)
!C        DO J1=2,NATOMS
!C           WRITE(4,'(A5,3F20.10)') 'LB   ',T1(3*(J1-1)+1),T1(3*(J1-1)+2),T1(3*(J1-1)+3)
!C        ENDDO
!C     ENDIF
!C
!C  Find the atom with the largest distance from the z axis.
!C
      DMAX=-1.0D0
      DO J1=1,NATOMS
         DIST(J1)=SQRT(T1(3*(J1-1)+1)**2+T1(3*(J1-1)+2)**2)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1) 
      ENDDO
      DMAX2=-1.0D100

!C  PROJ is the sum of the x components. TVEC arguments are dummys here.
!C  Use T2 as a dummy in order not to change T1 until we have decided which 
!C  atom to put in the xz plane.
!C
      TXVEC(1:3)=0.0D0; TYVEC(1:3)=0.0D0;TZVEC(1:3)=0.0D0
      DO J1=1,NATOMS
!        PRINT '(A,I8,3F20.10)','z J1,DMAX,DIST(J1),ABS(DIST(J1)-DMAX)=',J1,DMAX,DIST(J1),ABS(DIST(J1)-DMAX)
!        PRINT '(A,2F20.10,L5)','ABS(DIST(J1)-DMAX),CUTOFF1,test=',ABS(DIST(J1)-DMAX),CUTOFF1,ABS(DIST(J1)-DMAX).LT.CUTOFF1
         IF (ABS(DIST(J1)-DMAX).LT.CUTOFF1) THEN
            T2(1:3*NATOMS)=T1(1:3*NATOMS)
            CALL ROTXZ(NATOMS,J1,T2,PROJ,DIST,TXVEC,TYVEC,TZVEC,CUTOFF1)
!           PRINT '(A,I8,3G20.10)','J1,DMAX2,PROJ,ABS(ABS(PROJ-DMAX2)=',J1,DMAX2,PROJ,ABS(PROJ-DMAX2)
            IF (ABS(PROJ-DMAX2).LT.CUTOFF1) THEN
               NORBIT2=NORBIT2+1
               IF (NORBIT2.EQ.NCHOOSE2) THEN
                  JMAX2=J1
               ENDIF
            ELSE IF (PROJ.GT.DMAX2) THEN
               NORBIT2=1
               DMAX2=PROJ
               JMAX2=J1
            ENDIF
!           PRINT '(A,2G20.10)','PROJ,DMAX2=',PROJ,DMAX2
         ENDIF
      ENDDO
!     IF (DEBUG) PRINT '(A,6I6)','myorient> NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX1,JMAX2=',
!    &                            NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX1,JMAX2
!C
!C  and rotate it into the xz plane.
!C
!     IF (DEBUG) PRINT '(A,I6,A,G20.10,A)','atom ',JMAX2,' is now the furthest from the  z axis, distance=', 
!    &                                        DIST(JMAX2),' rotate into xz plane'
!
! The correct distance from the z axis is used for each atom - they are saved in the vector DIST
! DMAX2 is a dummy here.
!
      CALL ROTXZ(NATOMS,JMAX2,T1,DMAX2,DIST,XVEC,YVEC,ZVEC,CUTOFF1)

!     IF (DEBUG) PRINT*,'after myorient'
!     IF (DEBUG) WRITE(*,'(3F20.10)') (T1(J1),J1=1,3*NATOMS)

      ROT(1:3,1)=XVEC(1:3); ROT(1:3,2)=YVEC(1:3); ROT(1:3,3)=ZVEC(1:3)
      ROTINV(1,1:3)=XVEC(1:3); ROTINV(2,1:3)=YVEC(1:3); ROTINV(3,1:3)=ZVEC(1:3)
!C
!C  Check orthogonality:
!C
      DO J1=1,3
         DO J2=1,3
            DUMMY=0.0D0
            DO J3=1,3
               DUMMY=DUMMY+ROT(J1,J3)*ROTINV(J3,J2)
            ENDDO
            TEMP(J1,J2)=DUMMY
         ENDDO
      ENDDO
!C     PRINT '(A)','myorient> transformation product:'
!C     PRINT '(3F20.10)',TEMP(1:3,1:3)
!C
!C  Compare Q1, ROT * T1, ROTINV * T1
!C
!C     PRINT '(A)','myorient> original coordinates:'
!C     PRINT '(6F20.10)',Q1(1:12)
!C     DO J1=1,NATOMS
!C        DO J2=1,3
!C           DUMMY=0.0D0
!C           DO J3=1,3
!C              DUMMY=DUMMY+ROT(J2,J3)*T1(3*(J1-1)+J3)
!C           ENDDO
!C           TVEC(3*(J1-1)+J2)=DUMMY
!C        ENDDO
!C     ENDDO
!C     PRINT '(A)','myorient> ROT * T1:'
!C     PRINT '(6F20.10)',TVEC(1:12)
!C     DO J1=1,NATOMS
!C        DO J2=1,3
!C           DUMMY=0.0D0
!C           DO J3=1,3
!C              DUMMY=DUMMY+ROTINV(J2,J3)*T1(3*(J1-1)+J3)
!C           ENDDO
!C           TVEC(3*(J1-1)+J2)=DUMMY
!C        ENDDO
!C     ENDDO
!C     PRINT '(A)','myorient> ROTINV * T1:'
!C     PRINT '(6F20.10)',TVEC(1:12)

      RETURN
      END

      SUBROUTINE ROTXZ(NATOMS,JDO,T1,PROJ,DIST,XVEC,YVEC,ZVEC,CUTOFF1)
      IMPLICIT NONE
      INTEGER NATOMS, JDO, J1, J3
      DOUBLE PRECISION T1(3*NATOMS), PROJ, DIST(NATOMS), RVEC(3), COST, SINT, RDOTN, TX, TY, TZ, CUTOFF1
      DOUBLE PRECISION XVEC(3), YVEC(3), ZVEC(3)

!C     PRINT '(I6)',NATOMS
!C     PRINT '(A,I6,G20.10)','before rotation'
!C     WRITE(*,'(A2,2X,3G20.10)') ('LA',T1(3*(J3-1)+1),T1(3*(J3-1)+2),T1(3*(J3-1)+3),J3=1,NATOMS)

      IF (ABS(T1(3*(JDO-1)+2)).LT.1.0D-8) THEN
         IF (T1(3*(JDO-1)+1).LT.0.0D0) THEN ! rotate about the z axis DO NOT INVERT!!
            DO J1=1,NATOMS
               T1(3*(J1-1)+1)=-T1(3*(J1-1)+1)
               T1(3*(J1-1)+2)=-T1(3*(J1-1)+2)
            ENDDO
            XVEC(1)=-XVEC(1); XVEC(2)=-XVEC(2)
            YVEC(1)=-YVEC(1); YVEC(2)=-YVEC(2)
            ZVEC(1)=-ZVEC(1); ZVEC(2)=-ZVEC(2)
         ENDIF
         GOTO 20
      ENDIF

      COST=T1(3*(JDO-1)+1)/DIST(JDO)
      SINT=T1(3*(JDO-1)+2)/DIST(JDO)
!C     PRINT '(A,4G20.10)','T1(3*(JDO-1)+2),COST,SINT,1=',T1(3*(JDO-1)+2),COST,SINT,COST**2+SINT**2
      IF (ABS(COST**2+SINT**2-1.0D0).GT.2.0D-3) THEN
         WRITE(*, '(A,G20.10)') 'ERROR - in ROTXZ cos**2+sin**2=',COST**2+SINT**2
         STOP
      ENDIF
      RVEC(1)=0.0D0
      RVEC(2)=0.0D0
      RVEC(3)=1.0D0
      DO J1=1,NATOMS
         IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=T1(3*(J1-1)+1)*RVEC(1)+T1(3*(J1-1)+2)*RVEC(2)+T1(3*(J1-1)+3)*RVEC(3)
            TX=T1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)+(T1(3*(J1-1)+2)*RVEC(3)-T1(3*(J1-1)+3)*RVEC(2))*SINT
            TY=T1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)+(T1(3*(J1-1)+3)*RVEC(1)-T1(3*(J1-1)+1)*RVEC(3))*SINT
            TZ=T1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)+(T1(3*(J1-1)+1)*RVEC(2)-T1(3*(J1-1)+2)*RVEC(1))*SINT
            T1(3*(J1-1)+1)=TX
            T1(3*(J1-1)+2)=TY
            T1(3*(J1-1)+3)=TZ
         ENDIF
      ENDDO
      RDOTN=XVEC(1)*RVEC(1)+XVEC(2)*RVEC(2)+XVEC(3)*RVEC(3)
      TX=XVEC(1)*COST + RVEC(1)*RDOTN*(1-COST)+(XVEC(2)*RVEC(3)-XVEC(3)*RVEC(2))*SINT
      TY=XVEC(2)*COST + RVEC(2)*RDOTN*(1-COST)+(XVEC(3)*RVEC(1)-XVEC(1)*RVEC(3))*SINT
      TZ=XVEC(3)*COST + RVEC(3)*RDOTN*(1-COST)+(XVEC(1)*RVEC(2)-XVEC(2)*RVEC(1))*SINT
      XVEC(1)=TX; XVEC(2)=TY; XVEC(3)=TZ
      RDOTN=YVEC(1)*RVEC(1)+YVEC(2)*RVEC(2)+YVEC(3)*RVEC(3)
      TX=YVEC(1)*COST + RVEC(1)*RDOTN*(1-COST)+(YVEC(2)*RVEC(3)-YVEC(3)*RVEC(2))*SINT
      TY=YVEC(2)*COST + RVEC(2)*RDOTN*(1-COST)+(YVEC(3)*RVEC(1)-YVEC(1)*RVEC(3))*SINT
      TZ=YVEC(3)*COST + RVEC(3)*RDOTN*(1-COST)+(YVEC(1)*RVEC(2)-YVEC(2)*RVEC(1))*SINT
      YVEC(1)=TX; YVEC(2)=TY; YVEC(3)=TZ
      RDOTN=ZVEC(1)*RVEC(1)+ZVEC(2)*RVEC(2)+ZVEC(3)*RVEC(3)
      TX=ZVEC(1)*COST + RVEC(1)*RDOTN*(1-COST)+(ZVEC(2)*RVEC(3)-ZVEC(3)*RVEC(2))*SINT
      TY=ZVEC(2)*COST + RVEC(2)*RDOTN*(1-COST)+(ZVEC(3)*RVEC(1)-ZVEC(1)*RVEC(3))*SINT
      TZ=ZVEC(3)*COST + RVEC(3)*RDOTN*(1-COST)+(ZVEC(1)*RVEC(2)-ZVEC(2)*RVEC(1))*SINT
      ZVEC(1)=TX; ZVEC(2)=TY; ZVEC(3)=TZ

20    CONTINUE

      PROJ=0.0D0
!
! For possibly sloppy alignments with LOCALPERMDIST it may be important to
! increase the cutoff and use CUTOFF1 ! DJW 27/8/09
!
      DO J1=1,NATOMS
!C        IF (T1(3*(J1-1)+3).GT.1.0D-2) PROJ=PROJ+T1(3*(J1-1)+1)
         IF (T1(3*(J1-1)+3).GT.CUTOFF1) PROJ=PROJ+T1(3*(J1-1)+1)
      ENDDO
!C     PRINT '(I6)',NATOMS
!C     PRINT '(A,I6,G20.10)','after rotation atom,PROJ=',JDO,PROJ
!C     WRITE(*,'(A2,2X,3G20.10)') ('LA',T1(3*(J3-1)+1),T1(3*(J3-1)+2),T1(3*(J3-1)+3),J3=1,NATOMS)

      RETURN
      END
