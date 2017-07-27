!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
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
SUBROUTINE BULKMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,TWOD,DEBUG,BOXLX,BOXLY,BOXLZ,PITEST,RESETA)
USE COMMONS,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, MYUNIT, GEOMDIFFTOL

IMPLICIT NONE
INTEGER J1, NATOMS, NPMIN, NGMIN, J2, PERM(NATOMS), PBEST(NATOMS), NDUMMY, NMATCHED, PATOMS, J3, J4, NMBEST, ND1
DOUBLE PRECISION DUMMYB(3*NATOMS),DUMMYA(3*NATOMS),DISTANCE,BOXLX,BOXLY,BOXLZ,XSHIFT,YSHIFT,ZSHIFT,XTEMP(3*NATOMS)
DOUBLE PRECISION XBEST(3*NATOMS), DMIN, DTOTAL, DIST, GDSQ
LOGICAL TWOD,DEBUG,PITEST,SAMEMIN,RESETA
COMMON /BULKSHIFT/ XSHIFT,YSHIFT,ZSHIFT
!
! Find smallest group of permutable atoms.
! Translate first atom of group to all positions and then find nearest atom within
! the same group for every other atom.
! Keep the best translation/permutation, which corresponds to the smallest
! minimum image distance.
!
DISTANCE=1.0D100
PITEST=.FALSE.
SAMEMIN=.TRUE.
GDSQ=GEOMDIFFTOL**2/NATOMS
NPMIN=HUGE(1)
! PRINT *,'DUMMYA in bulkmindist:'
! PRINT '(3G20.10)',DUMMYA(1:3*NATOMS)
! PRINT *,'DUMMYB in bulkmindist:'
! PRINT '(3G20.10)',DUMMYB(1:3*NATOMS)
DO J1=1,NPERMGROUP
   IF (NPERMSIZE(J1).LT.NPMIN) THEN
      NPMIN=NPERMSIZE(J1)
      NGMIN=J1
   ENDIF
ENDDO
ND1=0
DO J1=1,NGMIN-1
   ND1=ND1+NPERMSIZE(J1)
ENDDO
IF (DEBUG) WRITE(MYUNIT,'(3(A,I6))') 'bulkmindist> Smallest group of permutable atoms is number ',NGMIN,' with ',NPMIN,' members'
outer: DO J1=ND1+1,ND1+NPMIN
   J2=PERMGROUP(J1)
   XSHIFT=DUMMYA(3*(J2-1)+1)-DUMMYB(3*(ND1)+1)-BOXLX*NINT((DUMMYA(3*(J2-1)+1)-DUMMYB(3*(ND1)+1))/BOXLX)
   YSHIFT=DUMMYA(3*(J2-1)+2)-DUMMYB(3*(ND1)+2)-BOXLY*NINT((DUMMYA(3*(J2-1)+2)-DUMMYB(3*(ND1)+2))/BOXLY)
   IF (.NOT.TWOD) ZSHIFT=DUMMYA(3*(J2-1)+3)-DUMMYB(3*(ND1)+3)-BOXLZ*NINT((DUMMYA(3*(J2-1)+3)-DUMMYB(3*(ND1)+3))/BOXLZ)
   DO J2=1,NATOMS
      XTEMP(3*(J2-1)+1)=DUMMYA(3*(J2-1)+1)-XSHIFT
      XTEMP(3*(J2-1)+2)=DUMMYA(3*(J2-1)+2)-YSHIFT
      IF (.NOT.TWOD) XTEMP(3*(J2-1)+3)=DUMMYA(3*(J2-1)+3)-ZSHIFT
   ENDDO
   NDUMMY=1
   PERM(1:NATOMS)=-1
   NMATCHED=0
   DTOTAL=0.0D0
   DO J2=1,NPERMGROUP
      PATOMS=NPERMSIZE(J2)
      loop1: DO J3=1,PATOMS    ! for each atom in fixed structure B in group J2
         DMIN=1.0D100
         DO J4=1,PATOMS ! which is the closest atom in the same group for the structure in XTEMP (shifted A)?
            DIST=(XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+1)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+1) &
 &  - BOXLX*NINT((XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+1)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+1))/BOXLX))**2 &
            &  + (XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+2)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+2) &
 &  - BOXLY*NINT((XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+2)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+2))/BOXLY))**2 
            IF (.NOT.TWOD) DIST=DIST+(XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+3)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+3) &
 &  - BOXLZ*NINT((XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+3)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+3))/BOXLZ))**2
            IF (DIST.LT.DMIN) THEN
               DMIN=DIST
               PERM(PERMGROUP(NDUMMY+J3-1))=PERMGROUP(NDUMMY+J4-1)
               IF (DIST.LT.GDSQ) THEN
!                 WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)') ' match found between atom ',PERMGROUP(NDUMMY+J3-1), &
! &                                            ' and ',PERMGROUP(NDUMMY+J4-1),' DIST=',DIST
                  NMATCHED=NMATCHED+1
                  DTOTAL=DTOTAL+DMIN
                  IF (PERM(PERMGROUP(NDUMMY+J3-1)).NE.PERMGROUP(NDUMMY+J3-1)) SAMEMIN=.FALSE.
                  CYCLE loop1
               ENDIF
            ENDIF
         ENDDO
         DTOTAL=DTOTAL+DMIN
!        WRITE(MYUNIT,'(A,I6,A,G20.10,A,I6)') ' match failed for atom ',PERMGROUP(NDUMMY+J3-1),' DMIN=',DMIN,' J1=',J1
         CYCLE outer ! If we reached here then we don't have a permutational isomer because
                     ! the atom specified in the J3 loop does not have a partner.
      ENDDO loop1
      NDUMMY=NDUMMY+NPERMSIZE(J2)
   ENDDO
   IF (SAMEMIN) THEN
      IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'bulkmindist> identical isomers identified, distance=',SQRT(DTOTAL)
   ELSE
      IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'bulkmindist> permutational isomers identified, distance=',SQRT(DTOTAL)
   ENDIF
   PITEST=.TRUE.
   DISTANCE=DTOTAL
   IF (RESETA) THEN
      DO J2=1,NATOMS
         DUMMYA(3*(J2-1)+1)=XTEMP(3*(PERM(J2)-1)+1)-BOXLX*NINT(XTEMP(3*(PERM(J2)-1)+1)/BOXLX)
         DUMMYA(3*(J2-1)+2)=XTEMP(3*(PERM(J2)-1)+2)-BOXLY*NINT(XTEMP(3*(PERM(J2)-1)+2)/BOXLY)
         IF (.NOT.TWOD) DUMMYA(3*(J2-1)+3)=XTEMP(3*(PERM(J2)-1)+3)-BOXLZ*NINT(XTEMP(3*(PERM(J2)-1)+3)/BOXLZ)
      ENDDO
   ENDIF

   RETURN
ENDDO outer

IF (DEBUG) WRITE(MYUNIT,'(A)') 'bulkmindist> structures are not permutational isomers'

RETURN

END SUBROUTINE BULKMINDIST
!
! Apply Oh point group operation number OPNUM to coordinates in
! vector X of dimension 3*NLOCAL, returning the result in 
! vector Y.
!
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
