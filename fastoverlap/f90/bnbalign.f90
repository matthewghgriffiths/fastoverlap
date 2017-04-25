

INCLUDE "commons.f90"

MODULE GOPERMDIST

! SAVECOORDSA(3*NATOMS,NSTRUCTS) stores the centred candidate structures
! SAVECOORDSB(3*NATOMS) stores the centred target structure

! PERMCOORDSB(3,NATOMS,NPERMGROUP) stores the structures for the k-d tree


USE COMMONS, ONLY : PERMGROUP, NPERMSIZE, NPERMGROUP, BESTPERM, MYUNIT, &
 & NSETS, SETS, OHCELLT, PERMINVOPT, PERMDIST, PERMOPT
USE KDTREE2_MODULE, ONLY : KDTREE2, KDTREE2PTR, KDTREE2_CREATE
USE PRIORITYQUEUE, ONLY: QUEUE

IMPLICIT NONE

INTEGER, SAVE :: NATOMS, NCALC, NLAP, NQUENCH, NBAD
INTEGER, SAVE :: PMAXNEI = 60 ! Number of nearest neighbours to store
DOUBLE PRECISION, PARAMETER :: PSCALE = 1.D6 ! Scale for linear assignment problem
DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
DOUBLE PRECISION, SAVE :: ATOL=1E-3, RTOL=1E-1
LOGICAL, SAVE :: DEBUG = .TRUE.


DOUBLE PRECISION, SAVE :: LVECS(3,0:8), FVECS(4,6)

DOUBLE PRECISION, SAVE :: CMAX,CMAY,CMAZ,CMBX,CMBY,CMBZ
DOUBLE PRECISION, SAVE :: DUMMYRMAT(3,3), TRMAT(3,3)
LOGICAL, SAVE :: PERMINVOPTSAVE, NOINVERSIONSAVE

! Module saves periodic conditions variables
LOGICAL, SAVE :: BULKT
LOGICAL, SAVE :: OHCELLTSAVE
DOUBLE PRECISION, SAVE :: BOXLX,BOXLY,BOXLZ,BOXVEC(3)

! Arrays to store target and candidate structures and best found structures
DOUBLE PRECISION, SAVE, ALLOCATABLE  :: SAVECOORDSA(:,:),PERMCOORDSB(:,:,:), &
 & SAVECOORDSB(:), BESTCOORDSA(:,:), BESTRMAT(:,:,:)
DOUBLE PRECISION, SAVE, ALLOCATABLE  :: SAVERA(:,:), SAVERB(:)
INTEGER, SAVE, ALLOCATABLE :: BESTITERS(:)
INTEGER, SAVE :: BESTID, BESTITER


! Used when calculating Boundsin CALCBOUNDS
DOUBLE PRECISION :: BRANCHVECS(3,8)
DOUBLE PRECISION, SAVE, ALLOCATABLE  :: DUMMYCOORDSA(:,:), PDUMMYND(:)
! Arrays of distances and nearest neighbour distances
DOUBLE PRECISION, SAVE, ALLOCATABLE :: DUMMYDISTS(:,:), DUMMYNEARDISTS(:)
DOUBLE PRECISION, SAVE, ALLOCATABLE :: DUMMYDISPS(:,:,:)
! Arrays of bounded distances and nearest neighbour distances
DOUBLE PRECISION, SAVE, ALLOCATABLE :: DUMMYLDISTS(:,:), DUMMYNEARLDISTS(:), &
 & DUMMYLDISTS2(:,:), DUMMYDOTDISP(:,:,:)

INTEGER, SAVE, ALLOCATABLE :: DUMMYIDX(:,:), DINVIDX(:,:), DUMMYNEARIDX(:)
INTEGER, SAVE, ALLOCATABLE :: INVPERMGROUP(:)

! Used when solving assignment problem
DOUBLE PRECISION, SAVE, ALLOCATABLE :: PDUMMYA(:), PDUMMYB(:), DUMMYA(:), &
    & DUMMYB(:), XBESTA(:), XBESTASAVE(:)
INTEGER, SAVE, ALLOCATABLE :: NEWPERM(:), LPERM(:)

TYPE(KDTREE2PTR), ALLOCATABLE :: KDTREES(:)
TYPE(QUEUE) :: Q

DATA LVECS / &
 & 0.0D0, 0.0D0, 0.0D0, &
 & 1.0D0, 1.0D0, 1.0D0, &
 & 1.0D0, 1.0D0, -1.0D0, &
 & 1.0D0, -1.0D0, 1.0D0, &
 & 1.0D0, -1.0D0, -1.0D0, &
 & -1.0D0, 1.0D0, 1.0D0, &
 & -1.0D0, 1.0D0, -1.0D0, &
 & -1.0D0, -1.0D0, 1.0D0, &
 & -1.0D0, -1.0D0, -1.0D0 /

DATA FVECS / &
 &  1.0D0,  1.0D0,  1.0D0,  1.0D0, &
 &  1.0D0,  1.0D0, -1.0D0, -1.0D0, &
 &  1.0D0, -1.0D0,  1.0D0, -1.0D0, &
 & -1.0D0, -1.0D0, -1.0D0, -1.0D0, &
 & -1.0D0, -1.0D0,  1.0D0,  1.0D0, &
 & -1.0D0,  1.0D0, -1.0D0,  1.0D0 /

! Private so that module works with f2py and static linking to kdtree2.f90 and
! priorityqueue.f90
PRIVATE :: KDTREES, Q

CONTAINS

SUBROUTINE RUNGROUP(NITER, FORCE, IPRINT, BESTUPPER, NSTRUCTS, UPDATE)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NITER, IPRINT, NSTRUCTS, UPDATE
LOGICAL, INTENT(IN) :: FORCE
DOUBLE PRECISION, INTENT(INOUT) :: BESTUPPER(NSTRUCTS)


DOUBLE PRECISION LOWERBOUND, UPPERBOUND, VECTOR(3), WIDTH
INTEGER I,IDNUM,NODEITER,NSUCCESS

DO I=1,NITER

    CALL QUEUEGET(LOWERBOUND, UPPERBOUND, VECTOR, WIDTH, NODEITER, IDNUM)

    CALL BRANCH(VECTOR,WIDTH,IDNUM,BESTUPPER(IDNUM),FORCE)

    IF(DEBUG.AND.(IPRINT.GT.0).AND.(MOD(I,IPRINT).EQ.0)) THEN
        WRITE(MYUNIT,'(A)') &
         & "gopermdist> -----------------STATUS UPDATE----------------"
        WRITE(MYUNIT,'(A,I16)') &
         & "gopermdist> iteration  number           = ", I
!        WRITE(MYUNIT,'(A,G20.6)') &
!         & "gopermdist> lowest upper bound so far   = ", BESTUPPER
        WRITE(MYUNIT,'(A,G20.6)') &
         & "gopermdist> highest lower bound so far  = ", LOWERBOUND
        WRITE(MYUNIT,'(A,I16)') &
         & "gopermdist> total calculations so far   = ", NCALC
        WRITE(MYUNIT,'(A,I16)') &
         & "gopermdist> queue length                = ", QUEUELEN()
        WRITE(MYUNIT,'(A)') &
         & "gopermdist> ----------------------------------------------"
    ENDIF


    IF(QUEUELEN().LE.0) THEN
        IF(DEBUG) WRITE(MYUNIT,'(A)') &
             & "gopermdist> priority queue empty, stopping"
    END IF

!    IF((QUEUELEN().LE.0).OR.((LOWERBOUND).GT.(BESTUPPER - RTOL*BESTUPPER - ATOL))) THEN
!        IF(DEBUG) THEN 
!            WRITE(MYUNIT,'(A)') &
!             & "gopermdist> -------------------SUCCESS--------------------"
!!            WRITE(MYUNIT,'(A,G20.6)') &
!!             & "gopermdist> converged on minimum RMSD   = ", BESTUPPER
!            WRITE(MYUNIT,'(A,I16)') &
!             & "gopermdist> total calculations          = ", NCALC
!            WRITE(MYUNIT,'(A,I16)') &
!             & "gopermdist> found best on iteration     = ", BESTITER
!            WRITE(MYUNIT,'(A,I16)') &
!             & "gopermdist> best structure              = ", BESTID
!            WRITE(MYUNIT,'(A)') &
!             & "gopermdist> -------------------SUCCESS--------------------"
!        END IF
!        EXIT
!    END IF

END DO

END SUBROUTINE RUNGROUP

SUBROUTINE RUN(NITER, FORCE, IPRINT, BESTUPPER)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NITER, IPRINT
LOGICAL, INTENT(IN) :: FORCE
DOUBLE PRECISION, INTENT(INOUT) :: BESTUPPER

DOUBLE PRECISION LOWERBOUND, UPPERBOUND, VECTOR(3), WIDTH
INTEGER I,IDNUM,NODEITER

DO I=1,NITER

    CALL QUEUEGET(LOWERBOUND, UPPERBOUND, VECTOR, WIDTH, NODEITER, IDNUM)

    IF(DEBUG.AND.(IPRINT.GT.0).AND.(MOD(I,IPRINT).EQ.0)) THEN
        WRITE(MYUNIT,'(A)') &
         & "gopermdist> -----------------STATUS UPDATE----------------"
        WRITE(MYUNIT,'(A,I16)') &
         & "gopermdist> iteration  number           = ", I
        WRITE(MYUNIT,'(A,G20.6)') &
         & "gopermdist> lowest upper bound so far   = ", BESTUPPER
        WRITE(MYUNIT,'(A,G20.6)') &
         & "gopermdist> highest lower bound so far  = ", LOWERBOUND
        WRITE(MYUNIT,'(A,I16)') &
         & "gopermdist> total calculations so far   = ", NCALC
        WRITE(MYUNIT,'(A,I16)') &
         & "gopermdist> queue length                = ", QUEUELEN()
        WRITE(MYUNIT,'(A)') &
         & "gopermdist> ----------------------------------------------"
    ENDIF

    CALL BRANCH(VECTOR,WIDTH,IDNUM,BESTUPPER,FORCE)

    IF(QUEUELEN().LE.0) THEN
        IF(DEBUG) WRITE(MYUNIT,'(A)') &
             & "gopermdist> priority queue empty, stopping"
    END IF

    IF((QUEUELEN().LE.0).OR.((LOWERBOUND).GT.(BESTUPPER - RTOL*BESTUPPER - ATOL))) THEN
        IF(DEBUG) THEN 
            WRITE(MYUNIT,'(A)') &
             & "gopermdist> -------------------SUCCESS--------------------"
            WRITE(MYUNIT,'(A,G20.6)') &
             & "gopermdist> converged on minimum RMSD   = ", BESTUPPER
            WRITE(MYUNIT,'(A,I16)') &
             & "gopermdist> total calculations          = ", NCALC
            WRITE(MYUNIT,'(A,I16)') &
             & "gopermdist> found best on iteration     = ", BESTITER
            WRITE(MYUNIT,'(A,I16)') &
             & "gopermdist> best structure              = ", BESTID
            WRITE(MYUNIT,'(A)') &
             & "gopermdist> -------------------SUCCESS--------------------"
        END IF
        EXIT
    END IF

END DO

END SUBROUTINE

SUBROUTINE ADDNODE(VECTOR,WIDTH,IDNUM,BESTUPPER,FORCE)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: VECTOR(3), WIDTH
DOUBLE PRECISION, INTENT(INOUT) :: BESTUPPER
INTEGER, INTENT(IN) :: IDNUM
LOGICAL, INTENT(IN) :: FORCE

DOUBLE PRECISION :: LOWERBOUND, UPPERBOUND, DIST2

LOGICAL :: PERMINVOPTSAVE, OHCELLTSAVE

CALL CALCBOUNDS(LOWERBOUND,UPPERBOUND,VECTOR,WIDTH,IDNUM,BESTUPPER,FORCE)

IF (UPPERBOUND.LE.BESTUPPER) THEN
    IF(DEBUG) WRITE(MYUNIT, "(A,G20.5)") &
 & "gopermdist> NEW lowest upper bound RMSD = ", UPPERBOUND

    ! Don't need to test for inversion isomers
    PERMINVOPTSAVE = PERMINVOPT; OHCELLTSAVE = OHCELLT
    OHCELLT = .FALSE.; PERMINVOPT = .FALSE.

    CALL MINPERMDIST(SAVECOORDSB,DUMMYA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT, &
 & .FALSE.,BESTUPPER,DIST2,.FALSE.,DUMMYRMAT)

    ! Resetting keywords
    PERMINVOPT = PERMINVOPTSAVE; OHCELLT = OHCELLTSAVE
    NQUENCH = NQUENCH + 1

    IF(DEBUG) WRITE(MYUNIT, "(A,G20.5)") &
 & "gopermdist> post quench new lowest RMSD = ", BESTUPPER

    IF (.NOT.BULKT) BESTRMAT(:,:,IDNUM) = MATMUL(TRMAT,DUMMYRMAT)

    BESTCOORDSA(:,IDNUM) = DUMMYA
    BESTID = IDNUM
    BESTITER = NCALC
    CALL QUEUEPUT(LOWERBOUND,UPPERBOUND,VECTOR,WIDTH,NCALC,IDNUM)

ELSE IF( (LOWERBOUND ).LT.(BESTUPPER - RTOL*BESTUPPER - ATOL) ) THEN
    CALL QUEUEPUT(LOWERBOUND,UPPERBOUND,VECTOR,WIDTH,NCALC,IDNUM)
END IF



END SUBROUTINE ADDNODE

SUBROUTINE BRANCH(VECTOR,WIDTH,IDNUM,BESTUPPER,FORCE)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: VECTOR(3), WIDTH
DOUBLE PRECISION, INTENT(INOUT) :: BESTUPPER
INTEGER, INTENT(IN) :: IDNUM
LOGICAL, INTENT(IN) :: FORCE

DOUBLE PRECISION :: LOWERBOUND, UPPERBOUND, NEWVECT(3)

INTEGER I

DO I=1,8
    NEWVECT(:) = VECTOR + LVECS(:,I)*WIDTH*0.25D0
    CALL ADDNODE(NEWVECT,WIDTH*0.5D0,IDNUM,BESTUPPER,FORCE)
END DO

END SUBROUTINE BRANCH

SUBROUTINE CALCBOUNDS(LOWERBOUND,UPPERBOUND,VECTOR,WIDTH,IDNUM,BESTUPPER,FORCE)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: VECTOR(3), WIDTH, BESTUPPER
INTEGER, INTENT(IN) :: IDNUM
LOGICAL, INTENT(IN) :: FORCE

DOUBLE PRECISION, INTENT(OUT) :: LOWERBOUND, UPPERBOUND

DOUBLE PRECISION W,SINW,COSW,RA,RB,ESTLOWER,ESTUPPER,D
INTEGER I,J,J1,M,K,K1,IND,NDUMMY,NPERM,INFO,IA,IB
LOGICAL RECALC

DOUBLE PRECISION PERMDIST

W = SQRT(3.D0) * WIDTH * 0.5D0
COSW = COS(W)
SINW = SQRT(1.D0 - COSW**2)

IF(DEBUG) THEN 
    IF(BULKT) WRITE(MYUNIT, "(A,3F16.5)") &
 & "gopermdist> testing displacement vector = ", VECTOR
    IF(.NOT.BULKT) WRITE(MYUNIT, "(A,3F16.5)") &
 & "gopermdist> testing angle-axis vector   = ", VECTOR
    WRITE(MYUNIT, "(A,G20.5,A,I4)") &
 & "gopermdist> with width                  = ", WIDTH, &
 & "     on IDNUM    =", IDNUM
END IF

CALL TRANSFORM(DUMMYA, NATOMS, VECTOR, IDNUM)

! Find distance matrix
CALL PERMPAIRDISTS(SAVECOORDSB,DUMMYA,NATOMS,PMAXNEI,DUMMYDISTS,DUMMYIDX,NPERMGROUP)

!write(*,*) (dummyidx)

! Find bounded distanace matrix
IF(BULKT) THEN
    NDUMMY=0
    DO J1=1,NPERMGROUP
        NPERM=NPERMSIZE(J1)
        M = MIN(NPERM,PMAXNEI)
        DUMMYLDISTS(:NPERM*M,J1) = MAX(SQRT(DUMMYDISTS(:NPERM*M,J1)) - W,0.D0)**2
    ENDDO
ELSE
    NDUMMY=0
    DO J1=1,NPERMGROUP
        NPERM = NPERMSIZE(J1)
        M = MIN(NPERM,PMAXNEI)
        DO J=1,NPERM
            K=M*(J-1)
            RB = SAVERB(PERMGROUP(J+NDUMMY))
            DO I=1,M
                IND = K+I  
                RA = SAVERA(PERMGROUP(DUMMYIDX(IND,J1)+NDUMMY),IDNUM)
                DUMMYLDISTS(IND,J1) = BOUNDROTDISTANCE( &
                     & DUMMYDISTS(IND,J1),COSW,SINW,RB,RA)
            END DO
        ENDDO
    NDUMMY = NDUMMY + NPERMSIZE(J1)
    ENDDO
END IF

! Estimating upperbound by finding nearest neighbours
IF((.NOT.FORCE).OR.DEBUG) THEN
    CALL PERMNEARESTNEIGHBOURDISTS(DUMMYDISTS,DUMMYIDX,NATOMS,PMAXNEI, &
     & DUMMYNEARIDX,DUMMYNEARDISTS,NPERMGROUP)

    UPPERBOUND = SUM(DUMMYNEARDISTS)**0.5
    IF(DEBUG) WRITE(MYUNIT, "(A,G20.5,A)") &
 & "gopermdist> estimate for upper bound    = ", UPPERBOUND

    ! Check if permutation has been found anyway
    IF(UPPERBOUND.LT.BESTUPPER) THEN
        LPERM = 0
        DO J1=1,NATOMS
            LPERM(DUMMYNEARIDX(J1)) = 1
        END DO
        IF(ALL(LPERM.EQ.1)) THEN
            RECALC = .FALSE.
            IF(DEBUG) WRITE(MYUNIT, "(A)") &
 & "gopermdist> nearest neighbours are best permutation"
        ELSE
            RECALC = .TRUE.
        END IF
    ELSE
        RECALC = .FALSE.
    END IF
    ESTUPPER = UPPERBOUND
END IF


! Estimating Lower Bound by finding nearest neighbours
IF(DEBUG.OR.(.NOT.(FORCE.OR.RECALC))) THEN
    IF(BULKT) THEN
!        DO J1=1,NPERMGROUP
!            NPERM=NPERMSIZE(J1)
!            M = MIN(NPERM,PMAXNEI)
!            DUMMYLDISTS(:NPERM*M,J1) = MAX(SQRT(DUMMYDISTS(:NPERM*M,J1)) - W,0.D0)**2
!        ENDDO

        ! Find relative displacements
        DO J1=1,NPERMGROUP
            NPERM=NPERMSIZE(J1)
            M = MIN(NPERM,PMAXNEI)
            DO I=1,NPERM
                IB = PERMGROUP(I+NDUMMY)
                K = M*(I-1)
                DO J=1,M
                    IA = PERMGROUP(DUMMYIDX(K+J,J1)+NDUMMY)
                    DUMMYDISPS(:,K+J,J1) = SAVECOORDSB(3*IB-2:3*IB) - DUMMYA(3*IA-2:3*IA)
                    DUMMYDISPS(1,K+J,J1) = DUMMYDISPS(1,K+J,J1) - &
                     & NINT(DUMMYDISPS(1,K+J,J1)/BOXLX) * BOXLX
                    DUMMYDISPS(2,K+J,J1) = DUMMYDISPS(2,K+J,J1) - &
                     & NINT(DUMMYDISPS(2,K+J,J1)/BOXLY) * BOXLY
                    DUMMYDISPS(3,K+J,J1) = DUMMYDISPS(3,K+J,J1) - &
                     & NINT(DUMMYDISPS(3,K+J,J1)/BOXLZ) * BOXLZ

                    DUMMYDOTDISP(:,K+J,J1) = MATMUL(DUMMYDISPS(:,K+J,J1),LVECS(:,1:4))
                END DO
            END DO
            NDUMMY = NDUMMY + NPERM
        END DO

        ESTLOWER = HUGE(1.D0)
        DO I=1,6
            DO J1=1,NPERMGROUP
                NPERM=NPERMSIZE(J1)
                M = MIN(NPERM,PMAXNEI)
                DUMMYLDISTS2(:M*NPERM,J1) = MERGE(DUMMYDISTS(:M*NPERM,J1), &
                                                & DUMMYLDISTS(:M*NPERM,J1), & 
                 & MATMUL(FVECS(:,I),DUMMYDOTDISP(:,:M*NPERM,J1)).GT.0.D0)
            END DO
            
            CALL PERMNEARESTNEIGHBOURDISTS(DUMMYLDISTS2,DUMMYIDX,NATOMS, & 
             & PMAXNEI,DUMMYNEARIDX,DUMMYNEARLDISTS,NPERMGROUP)
            
            D = SUM(DUMMYNEARLDISTS)
            ESTLOWER = MIN(D, ESTLOWER)

            IF(DEBUG) WRITE(MYUNIT, "(A,I16,A,G10.5)") &
     & "gopermdist> estimating for face         = ", I, &
     & "         lower bound = ", D**0.5
        END DO
        ESTLOWER = SQRT(ESTLOWER)

    ELSE
        CALL PERMNEARESTNEIGHBOURDISTS(DUMMYLDISTS,DUMMYIDX,NATOMS,PMAXNEI, &
         & DUMMYNEARIDX,DUMMYNEARLDISTS,NPERMGROUP)
    
        ESTLOWER = SUM(DUMMYNEARLDISTS)**0.5
    END IF

    LOWERBOUND = ESTLOWER

    IF(DEBUG) WRITE(MYUNIT, "(A,G20.5)") &
     & "gopermdist> estimate for lower bound    = ", ESTLOWER
    
END IF




! If estimate of upperbound is lower than best found upperbound we need to
! solve assignment problem to find bounds
IF (FORCE.OR.RECALC) THEN

    ! Need to calculate this matrix to get total distance from reduced distance
    ! matrix and total permutation
    CALL INVPAIRDISTIDX(DUMMYIDX, DINVIDX, NATOMS, PMAXNEI, NPERMGROUP)

!    DINVIDX = -1
!    DO J1=1,NPERMGROUP
!        NPERM = NPERMSIZE(J1)
!        M = MIN(NPERM,PMAXNEI)
!        DO J=1,NPERM
!            K=M*(J-1)
!            K1 = NPERM*(J-1)
!            DO I=1,M
!                DINVIDX(K1+DUMMYIDX(K+I,J1),J1) = I
!            END DO
!        END DO
!    END DO

    IF(BULKT) THEN
        DO J1=1,NPERMGROUP
            NPERM=NPERMSIZE(J1)
!            M = MERGE(NPERM,PMAXNEI,NPERM.LT.PMAXNEI)
            M = MIN(NPERM,PMAXNEI)
            DUMMYLDISTS(:NPERM*M,J1) = MAX(SQRT(DUMMYDISTS(:NPERM*M,J1)) - W,0.D0)**2
        ENDDO
    END IF

    CALL FINDBESTPERM(DUMMYLDISTS,DUMMYIDX,NATOMS,PMAXNEI,NEWPERM, & 
     & LOWERBOUND,NPERMGROUP, INFO)
    
    CALL FINDPERMVAL(NEWPERM,NATOMS,DUMMYLDISTS,DINVIDX,PMAXNEI,NPERMGROUP,LOWERBOUND)
!    LOWERBOUND = 0.D0
!    NDUMMY = 0
!    DO J1=1,NPERMGROUP
!        NPERM = NPERMSIZE(J1)
!        M = MIN(NPERM,PMAXNEI)
!        DO J=1,NPERM
!!            K = M*(J-1)
!!            K1 = NPERM*(J-1)
!            IA = INVPERMGROUP(NEWPERM(PERMGROUP(J+NDUMMY)))-NDUMMY
!            I = DINVIDX(NPERM*(J-1)+IA,J1)
!            LOWERBOUND = LOWERBOUND + DUMMYLDISTS(M*(J-1)+I,J1)
!        END DO
!        NDUMMY = NDUMMY + NPERM
!    END DO 

!    LOWERBOUND = 0.D0
!    ! Perhaps there's a better way of calculating lowerbound from FINDBESTPERM?
!    DO J=1,NATOMS
!        I = NEWPERM(J)
!        LOWERBOUND = (LOWERBOUND + MAX(SQRT(PERMDIST( & 
!         & SAVECOORDSB(3*J-2:3*J),DUMMYA(3*I-2:3*I),BOXVEC,BULKT))-W,0.D0)**2)
!    END DO

    ! Check output of assignment problem
    IF(INFO.GT.0) THEN
        LOWERBOUND = 0.D0
        IF(DEBUG) WRITE(MYUNIT, "(A,I3)") &
 & "gopermdist> WARNING LAP algorithm failed to align npoints= ", INFO
    ELSE
        LOWERBOUND = SQRT(LOWERBOUND)
        IF(DEBUG) WRITE(MYUNIT, "(A,G20.5)") &
 & "gopermdist> calculated lower bound RMSD = ", LOWERBOUND
    END IF
    ! Calculate upperbound if lowerbound lower than bestupper
    IF((LOWERBOUND.LT.BESTUPPER).OR.FORCE) THEN
        CALL FINDBESTPERM(DUMMYDISTS,DUMMYIDX,NATOMS,PMAXNEI,LPERM, & 
         & UPPERBOUND,NPERMGROUP, INFO)

        CALL FINDPERMVAL(LPERM,NATOMS,DUMMYDISTS,DINVIDX,PMAXNEI,NPERMGROUP,UPPERBOUND)

!        UPPERBOUND = 0.D0
!        NDUMMY = 0
!        DO J1=1,NPERMGROUP
!            NPERM = NPERMSIZE(J1)
!            M = MIN(NPERM,PMAXNEI)
!            DO J=1,NPERM
!!                K = M*(J-1)
!!                K1 = NPERM*(J-1)
!                IA = INVPERMGROUP(NEWPERM(PERMGROUP(J+NDUMMY)))-NDUMMY
!                I = DINVIDX(NPERM*(J-1)+IA,J1)
!                UPPERBOUND = UPPERBOUND + DUMMYDISTS(M*(J-1)+I,J1)
!            END DO
!            NDUMMY = NDUMMY + NPERM
!        END DO

!        UPPERBOUND = 0.D0
!        DO J=1,NATOMS
!            I = LPERM(J)
!            UPPERBOUND = (UPPERBOUND + PERMDIST( & 
!         & SAVECOORDSB(3*J-2:3*J),DUMMYA(3*I-2:3*I),BOXVEC,BULKT))
!        END DO

        ! Check output of assignment problem
        IF(INFO.GT.0) THEN
            UPPERBOUND = HUGE(1.D0)
            IF(DEBUG) WRITE(MYUNIT, "(A,I3)") &
 & "gopermdist> WARNING LAP algorithm failed to align npoints= ", INFO
        ELSE
            UPPERBOUND = SQRT(UPPERBOUND)
            IF(DEBUG) WRITE(MYUNIT, "(A,G20.5)") &
 & "gopermdist> calculated upper bound RMSD = ", UPPERBOUND
        END IF
    ELSE
        UPPERBOUND = HUGE(1.D0)
    END IF
END IF

IF (DEBUG.AND.((ESTUPPER.GT.UPPERBOUND).OR.(ESTLOWER.GT.LOWERBOUND))) THEN
    WRITE(MYUNIT,"(A)") "gopermdist>************WARNING*********************"
    WRITE(MYUNIT,"(A)") "EST UPPER GT UPPERBOUND OR EST LOWER GT LOWERBOUND"
    WRITE(MYUNIT,"(A)") "gopermdist>************WARNING*********************"
    NBAD = NBAD + 1
ENDIF

NCALC = NCALC + 1

END SUBROUTINE CALCBOUNDS

SUBROUTINE FINDPERMVAL(PERM, NATOMS, MATVALS, DINVIDX, MAXNEI, NPERMGROUP, BEST)

IMPLICIT NONE
INTEGER, INTENT(IN) :: PERM(NATOMS), NATOMS, DINVIDX(NATOMS*NATOMS,NPERMGROUP), &
 & MAXNEI, NPERMGROUP
DOUBLE PRECISION, INTENT(IN) :: MATVALS(NATOMS*MAXNEI,NPERMGROUP)
DOUBLE PRECISION, INTENT(OUT) :: BEST

INTEGER J1,M,J,I,IA,NPERM,NDUMMY

BEST = 0.D0
NDUMMY = 0
DO J1=1,NPERMGROUP
    NPERM = NPERMSIZE(J1)
    M = MIN(NPERM,MAXNEI)
    DO J=1,NPERM
        IA = INVPERMGROUP(PERM(PERMGROUP(J+NDUMMY)))-NDUMMY
        I = DINVIDX(NPERM*(J-1)+IA,J1)
        BEST = BEST + MATVALS(M*(J-1)+I,J1)
    END DO
    NDUMMY = NDUMMY + NPERM
END DO 

END SUBROUTINE FINDPERMVAL

SUBROUTINE INVPAIRDISTIDX(DUMMYIDX, DINVIDX, NATOMS, MAXNEI, NPERMGROUP)

IMPLICIT NONE
INTEGER, INTENT(IN) :: DUMMYIDX(NATOMS*MAXNEI,NPERMGROUP), NATOMS, MAXNEI, NPERMGROUP
INTEGER, INTENT(OUT) :: DINVIDX(NATOMS*NATOMS,NPERMGROUP)
INTEGER J1,NPERM,I,J,M

DINVIDX = -1
DO J1=1,NPERMGROUP
    NPERM = NPERMSIZE(J1)
    M = MIN(NPERM,MAXNEI)
    DO J=1,NPERM
        DO I=1,M
            DINVIDX(NPERM*(J-1)+DUMMYIDX(M*(J-1)+I,J1),J1) = I
        END DO
    END DO
END DO

END SUBROUTINE INVPAIRDISTIDX

SUBROUTINE PERMNEARESTNEIGHBOURDISTS(NDISTS,NIDX,NATOMS,MAXNEI,NEARI,NEARD,NPERMGROUP)

IMPLICIT NONE
INTEGER, INTENT(IN) :: NATOMS,MAXNEI,NPERMGROUP,NIDX(MAXNEI*NATOMS,NPERMGROUP)
DOUBLE PRECISION, INTENT(IN) :: NDISTS(MAXNEI*NATOMS,NPERMGROUP)

INTEGER, INTENT(OUT) :: NEARI(NATOMS)
DOUBLE PRECISION, INTENT(OUT) :: NEARD(NATOMS)

INTEGER I, J1, J2, IND, NPERM, NDUMMY, M

NDUMMY = 0
DO J1=1,NPERMGROUP
    NPERM=NPERMSIZE(J1)
!    M = MERGE(NPERM,MAXNEI,NPERM.LT.MAXNEI)
    M = MIN(NPERM,PMAXNEI)
    CALL NEARESTNEIGHBOURDISTS(NDISTS(1:NPERM*M,J1),NIDX(1:NPERM*M,J1), & 
 & NPERM,M,LPERM(1:NPERM),PDUMMYND(1:NPERM))

    DO J2=1,NPERM
        IND = LPERM(J2)
        NEARI(PERMGROUP(NDUMMY+J2)) = PERMGROUP(NDUMMY + IND)
        NEARD(PERMGROUP(NDUMMY+J2)) = PDUMMYND(J2)
    END DO
    NDUMMY = NDUMMY + NPERM
END DO

END SUBROUTINE PERMNEARESTNEIGHBOURDISTS

SUBROUTINE NEARESTNEIGHBOURDISTS(CC, KK, N, MAXNEI, IDX, DISTS)

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, MAXNEI, KK(MAXNEI*N)
DOUBLE PRECISION, INTENT(IN) :: CC(MAXNEI*N)

INTEGER, INTENT(OUT) :: IDX(N)
DOUBLE PRECISION, INTENT(OUT) :: DISTS(N)

INTEGER I,J,K,M

M=MAXNEI
IF(N.LT.MAXNEI) M=N

DO I=1,N
    J = MINLOC(CC(M*(I-1)+1:M*I),1)
    DISTS(I) = CC(M*(I-1) + J)
    IDX(I)   = KK(M*(I-1) + J)
END DO

END SUBROUTINE NEARESTNEIGHBOURDISTS

SUBROUTINE FINDBESTPERM(NDISTS,NIDX,NATOMS,MAXNEI,PERM,DIST,NPERMGROUP,INFO)
! DISTANCE RETURN INACCURATE
IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS,NPERMGROUP,MAXNEI,NIDX(MAXNEI*NATOMS,NPERMGROUP)
DOUBLE PRECISION, INTENT(IN) :: NDISTS(MAXNEI*NATOMS,NPERMGROUP)

DOUBLE PRECISION, INTENT(OUT) :: DIST
INTEGER, INTENT(OUT) :: PERM(NATOMS), INFO

! COULD SET THESE AS MODULE VARIABLES
INTEGER*8 :: KK(NATOMS*MAXNEI), CC(NATOMS*MAXNEI)
INTEGER*8 :: FIRST(NATOMS+1), X(NATOMS), Y(NATOMS)
INTEGER*8 :: U(NATOMS), V(NATOMS), N8, SZ8, H
INTEGER N,M,I,J,K,K1,I1,J1,NDUMMY

DIST = 0.D0
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

    CALL JOVOSAP(N8, SZ8, CC(:M*N), KK(:M*N), FIRST(:N+1), Y(:N), X(:N), U(:N), V(:N), H)
    NLAP = NLAP + 1

    DO I=1,N
        IF (Y(I).GT.N) THEN
            Y(I)=N
            INFO = INFO + 1
        END IF
        IF (Y(I).LT.1) THEN
            Y(I)=1
            INFO = INFO + 1
        END IF
        PERM(PERMGROUP(NDUMMY+I)) = PERMGROUP(NDUMMY+Y(I))
    ENDDO
    DIST = DIST + H/PSCALE

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

SUBROUTINE PERMPAIRDISTS(COORDSB,COORDSA,NATOMS,MAXNEI,NDISTS,NIDX,NPERMGROUP)

! Uses module variables BOXLX, BOXLY, BOXLZ, BULKT when calculating periodic distances

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, NPERMGROUP, MAXNEI
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS)

INTEGER, INTENT(OUT) :: NIDX(MAXNEI*NATOMS,NPERMGROUP)
DOUBLE PRECISION, INTENT(OUT) :: NDISTS(MAXNEI*NATOMS,NPERMGROUP)

INTEGER NDUMMY,J1,J2,NPERM

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

FUNCTION BOUNDROTDISTANCE(D2,COSW,SINW,RA,RB) RESULT(LDIST)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: D2,COSW,SINW,RA,RB
DOUBLE PRECISION LDIST

DOUBLE PRECISION RARB,RA2RB2,COSAB,SINAB,MCOSAB

! Precalculate these?
RARB = 2*RA*RB
RA2RB2 = RA**2 + RB**2

COSAB = (RA2RB2 - D2)/RARB
SINAB = SQRT(1.D0-MIN(COSAB**2,1.D0)) ! Making sure sqrt is of positive number
MCOSAB = MERGE(1.D0, COSAB*COSW + SINAB*SINW, COSAB.GT.COSW)

LDIST = MAX(RA2RB2 - RARB*MCOSAB,0.D0)

END FUNCTION

FUNCTION QUEUELEN() RESULT(LENGTH)

IMPLICIT NONE
INTEGER LENGTH

LENGTH = Q%N

END FUNCTION

SUBROUTINE QUEUEGET(LOWERBOUND, UPPERBOUND, VECTOR, WIDTH, NITER, IDNUM)
USE PRIORITYQUEUE, ONLY: NODE, TOP

IMPLICIT NONE
DOUBLE PRECISION, INTENT(OUT) :: lowerbound, upperbound, vector(3), width
INTEGER, INTENT(OUT) :: niter, IDNUM

TYPE(NODE) RES

IF(Q%N.GT.0) THEN
    RES = TOP(Q)
    VECTOR = RES%VECTOR
    UPPERBOUND = RES%UPPERBOUND
    LOWERBOUND = RES%LOWERBOUND
    WIDTH = RES%WIDTH
    NITER = RES%NITER
    IDNUM = RES%IDNUM
ELSE IF(DEBUG) THEN
    WRITE(MYUNIT,"(A)") "gopermdist> warning, trying to read empty list"
ENDIF

END SUBROUTINE QUEUEGET

SUBROUTINE QUEUEPUT(LOWERBOUND, UPPERBOUND, VECTOR, WIDTH, NITER, IDNUM)
USE PRIORITYQUEUE, ONLY: ENQUEUE

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: lowerbound, upperbound, vector(3), width
INTEGER, INTENT(IN) :: niter, IDNUM

CALL ENQUEUE(Q, LOWERBOUND, UPPERBOUND, VECTOR, WIDTH, NITER, IDNUM)

END SUBROUTINE QUEUEPUT

SUBROUTINE QUEUECLEAR()
USE PRIORITYQUEUE, ONLY: NODE, TOP

IMPLICIT NONE
TYPE(NODE) RES

DO WHILE(Q%N.GT.0)
    RES = TOP(Q)
END DO

END SUBROUTINE QUEUECLEAR

SUBROUTINE INITIALISE(COORDSB,COORDSA,NATOMS,NBOXLX,NBOXLY,NBOXLZ,NBULKT)

!USE COMMONS, ONLY: PERMINVOPT, OHCELLT
IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS
DOUBLE PRECISION, INTENT(IN) :: COORDSB(3*NATOMS), COORDSA(3*NATOMS), &
 & NBOXLX, NBOXLY, NBOXLZ
LOGICAL, INTENT(IN) :: NBULKT

DOUBLE PRECISION BVEC(3)
INTEGER I, J, K, IND, NDUMMY, NUMSTRUCTS

BOXLX = NBOXLX
BOXLY = NBOXLY
BOXLZ = NBOXLZ
BOXVEC = (/BOXLX,BOXLY,BOXLZ/)
BULKT = NBULKT

NCALC   = 0
NLAP    = 0
NQUENCH = 0
NBAD = 0

! --------------------------------------------------------------------------- !
!    allocating memory to arrays
! --------------------------------------------------------------------------- !

NUMSTRUCTS = 1
IF (PERMINVOPT) THEN
    NUMSTRUCTS = 2
ELSE IF (BULKT.AND.OHCELLT) THEN
    NUMSTRUCTS = 48
ENDIF

CALL REALLOCATEARRAYS(NATOMS, NUMSTRUCTS, BULKT)

! --------------------------------------------------------------------------- !
!    calculate inverse permutation group
! --------------------------------------------------------------------------- !

DO I=1,NATOMS
    INVPERMGROUP(PERMGROUP(I)) = I
END DO

! --------------------------------------------------------------------------- !
!    storing coordinates to module
! --------------------------------------------------------------------------- !

NDUMMY = 0
IF(BULKT) THEN
!    Needed for k-d trees stuff
!    DO I=1,NPERMGROUP
!        DO J=1, NPERMSIZE(I)
!            IND = PERMGROUP(NDUMMY+J)
!            SAVECOORDSB(3*IND-2) = COORDSB(3*IND-2) - BOXLX*ANINT(COORDSB(3*IND-2)/BOXLX)
!            SAVECOORDSB(3*IND-1) = COORDSB(3*IND-1) - BOXLY*ANINT(COORDSB(3*IND-1)/BOXLY)
!            SAVECOORDSB(3 * IND) = COORDSB(3 * IND) - BOXLZ*ANINT(COORDSB(3 * IND)/BOXLZ)
!        ENDDO
!    NDUMMY = NDUMMY + NPERMSIZE(I)
!    ENDDO
    SAVECOORDSB = COORDSB
    IF(OHCELLT) THEN
        DO I=1,48
            CALL OHOPS(COORDSA,SAVECOORDSA(:,I),I,NATOMS)
        END DO
    ELSE
        SAVECOORDSA(:,1) = COORDSA
    END IF
ELSE
    ! Calculate COM
    DO J=1,NATOMS
        CMAX=CMAX+COORDSA(3*(J-1)+1)
        CMAY=CMAY+COORDSA(3*(J-1)+2)
        CMAZ=CMAZ+COORDSA(3*(J-1)+3)
    ENDDO
    CMAX=CMAX/NATOMS; CMAY=CMAY/NATOMS; CMAZ=CMAZ/NATOMS
    CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
    DO J=1,NATOMS
        CMBX=CMBX+COORDSB(3*(J-1)+1)
        CMBY=CMBY+COORDSB(3*(J-1)+2)
        CMBZ=CMBZ+COORDSB(3*(J-1)+3)
    ENDDO
    CMBX=CMBX/NATOMS; CMBY=CMBY/NATOMS; CMBZ=CMBZ/NATOMS

    ! Save COM centred coordinates
    DO I=1,NATOMS
        SAVECOORDSB(3*I-2) = COORDSB(3*I-2) - CMBX
        SAVECOORDSB(3*I-1) = COORDSB(3*I-1) - CMBY
        SAVECOORDSB(3 * I) = COORDSB(3 * I) - CMBZ
        SAVERB(I) = SQRT(SAVECOORDSB(3*I-2)**2+SAVECOORDSB(3*I-1)**2+ &
                       & SAVECOORDSB(3 * I)**2)
    ENDDO
    DO I=1,NATOMS
        SAVECOORDSA(3*I-2,1) = COORDSA(3*I-2) - CMAX
        SAVECOORDSA(3*I-1,1) = COORDSA(3*I-1) - CMAY
        SAVECOORDSA(3 * I,1) = COORDSA(3 * I) - CMAZ
        SAVERA(I,1) = SQRT(SAVECOORDSA(3*I-2,1)**2+SAVECOORDSA(3*I-1,1)**2+ &
                         & SAVECOORDSA(3 * I,1)**2)
    ENDDO
    ! Store inverted configuration
    IF (PERMINVOPT) THEN
        SAVECOORDSA(:,2) = -SAVECOORDSA(:,1)
        SAVERA(:,2) = SAVERA(:,1)
    END IF
END IF

! --------------------------------------------------------------------------- !
! Allocate and populate k-d trees, should be a faster way of finding nearest
! neighbours, currently isn't...
! --------------------------------------------------------------------------- !

IF(ALLOCATED(KDTREES)) DEALLOCATE(KDTREES)
ALLOCATE(KDTREES(NPERMGROUP))
NDUMMY = 0
DO I=1,NPERMGROUP
    IF(BULKT) THEN
    DO K=0,8
        BVEC(1) = BOXLX*LVECS(1,K)
        BVEC(2) = BOXLY*LVECS(2,K)
        BVEC(3) = BOXLZ*LVECS(3,K)
        DO J=1, NPERMSIZE(I)
            IND = PERMGROUP(NDUMMY+J)
            PERMCOORDSB(1,J+K*NPERMSIZE(I),I) = SAVECOORDSB(3*IND-2) + BVEC(1)
            PERMCOORDSB(2,J+K*NPERMSIZE(I),I) = SAVECOORDSB(3*IND-1) + BVEC(2)
            PERMCOORDSB(3,J+K*NPERMSIZE(I),I) = SAVECOORDSB(3*IND) + BVEC(3)
        ENDDO
    ENDDO
    KDTREES(I)%TREE => KDTREE2_CREATE(PERMCOORDSB(:,:NPERMSIZE(I)*9,I),NPERMSIZE(I)*9,.true.,.true.)
    ELSE
        DO J=1, NPERMSIZE(I)
            IND = PERMGROUP(NDUMMY+J)
            PERMCOORDSB(:,J,I)=COORDSB(3*IND-2:3*IND) - (/CMBX,CMBY,CMBZ/)
        ENDDO
        KDTREES(I)%TREE => KDTREE2_CREATE(PERMCOORDSB(:,:NPERMSIZE(I),I),NPERMSIZE(I),.true.,.true.)
    END IF
    NDUMMY = NDUMMY + NPERMSIZE(I)
ENDDO

CALL QUEUECLEAR()

END SUBROUTINE INITIALISE

SUBROUTINE SETNATOMS(NEWNATOMS)
! Checks if arrays need to be (re)allocated
IMPLICIT NONE

INTEGER, INTENT(IN) :: NEWNATOMS

IF(.NOT.(SIZE(PDUMMYA).EQ.(3*NEWNATOMS))) THEN
    IF(ALLOCATED(PDUMMYA)) THEN
        DEALLOCATE(PDUMMYA,PDUMMYB,DUMMYA,DUMMYB,XBESTA,XBESTASAVE)
        DEALLOCATE(NEWPERM, LPERM)
    ENDIF
    ALLOCATE(PDUMMYA(3*NEWNATOMS),PDUMMYB(3*NEWNATOMS),DUMMYA(3*NEWNATOMS), &
    &   DUMMYB(3*NEWNATOMS), XBESTA(3*NEWNATOMS), XBESTASAVE(3*NEWNATOMS))
    ALLOCATE(NEWPERM(NEWNATOMS), LPERM(NEWNATOMS))
ENDIF

END SUBROUTINE SETNATOMS

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

CALL SETNATOMS(NEWNATOMS)

NATOMS = NEWNATOMS
PERMGROUP = NEWPERMGROUP
NPERMSIZE = NEWNPERMSIZE
NSETS = 0

END SUBROUTINE SETPERM

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
!       maxnei: Maximum number of closest neighbours
      double precision scale, d, h

      parameter (scale = 1.0d6   )
!      parameter (maxnei = 60     )

      integer*8 first(n+1)!, x(n), y(n)
!      integer*8 u(n), v(n)
      integer   m, i, j, k, l, l2, t, a
      integer*8 n8, sz8
      integer J1

!     Distance function
      double precision permdist

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
               cc(k+j) = permdist(p(3*i-2), q(3*j-2), s, pbc)
               kk(k+j) = j
!              write(*,*) i, j, '-->', cc(k+j)
            enddo
         enddo
      else
!     We need to store the distances of the maxnei closeest neighbors
!     of each particle. The following builds a heap to keep track of
!     the maxnei closest neighbours seen so far. It might be more
!     efficient to use quick-select instead... (This is definately
!     true in the limit of infinite systems.)
        do i=1,n
           k = first(i)-1
           do j=1,m
              cc(k+j) = permdist(p(3*i-2), q(3*j-2), s, pbc)
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
              d = permdist(p(3*i-2), q(3*j-2), s, pbc)
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

SUBROUTINE TRANSFORM(NEWCOORDSA, NATOMS, VECTOR, IDNUM)

IMPLICIT NONE
INTEGER, INTENT(IN) :: NATOMS, IDNUM
DOUBLE PRECISION, INTENT(IN) :: VECTOR(3)

DOUBLE PRECISION, INTENT(OUT) :: NEWCOORDSA(3*NATOMS)

INTEGER I

IF(BULKT) THEN
    DO I=1,NATOMS
    NEWCOORDSA(3*I-2) = SAVECOORDSA(3*I-2,IDNUM) - VECTOR(1)
    NEWCOORDSA(3*I-1) = SAVECOORDSA(3*I-1,IDNUM) - VECTOR(2)
    NEWCOORDSA(3*I  ) = SAVECOORDSA(3*I  ,IDNUM) - VECTOR(3)
    ENDDO
ELSE
    CALL ANGLEAXIS2MAT(VECTOR, TRMAT)
    DO I=1,NATOMS
        NEWCOORDSA(3*I-2:3*I) = MATMUL(TRMAT,SAVECOORDSA(3*I-2:3*I,IDNUM))
    ENDDO
ENDIF

END SUBROUTINE TRANSFORM

SUBROUTINE ANGLEAXIS2MAT(VECTOR,RMAT)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: VECTOR(3)
DOUBLE PRECISION, INTENT(OUT) :: RMAT(3,3)

DOUBLE PRECISION THETA,X,Y,Z,S,C,C1,XS,YS,ZS,XC,YC,ZC,XYC,YZC,ZXC

THETA = SUM((VECTOR**2))**0.5

IF(THETA.EQ.0.D0) THEN
    RMAT = RESHAPE((/&
     & 1.00000000000D0,  0.0D0,  0.0D0,   & 
     & 0.0D0,  1.00000000000D0,  0.0D0,   & 
     & 0.0D0,  0.0D0,  1.00000000000D0/), (/3,3/))
ELSE
    X = VECTOR(1)/THETA; Y = VECTOR(2)/THETA; Z = VECTOR(3)/THETA
    S = SIN(THETA); C = COS(THETA); C1 = 1.D0 - C
    XS = X*S; YS = Y*S; ZS = Z*S
    XC = X*C1; YC = Y*C1; ZC = Z*C1
    XYC = X*YC; YZC = Y*ZC; ZXC = Z*XC
    
    RMAT = RESHAPE((/&
     & x * xC + c, xyC + zs, zxC - ys, &
     & xyC - zs, y * yC + c, yzC + xs, &
     & zxC + ys, yzC - xs, z * zC + c/), (/3,3/))
END IF

END SUBROUTINE ANGLEAXIS2MAT

SUBROUTINE MAT2ANGLEAXIS(VECTOR, RMAT)

IMPLICIT NONE
DOUBLE PRECISION, INTENT(OUT) :: VECTOR(3)
DOUBLE PRECISION, INTENT(IN) :: RMAT(0:2,0:2)

DOUBLE PRECISION TRACE, THETA

TRACE = RMAT(0,0)+RMAT(1,1)+RMAT(2,2)
THETA = ACOS(0.5D0*TRACE-0.5D0)
VECTOR = (/RMAT(2,1)-RMAT(1,2),RMAT(0,2)-RMAT(2,0),RMAT(1,0)-RMAT(0,1)/)
VECTOR = VECTOR * 0.5D0 * THETA / SIN(THETA)

END SUBROUTINE MAT2ANGLEAXIS

SUBROUTINE NEARESTNEIGHBOURSTREE(COORDSA, NATOMS, NDISTS, NIDX, NN, BULKT)
! Returns nearest neighbour distances and indexes of closest points in 
! PERMCOORDSB
! NDISTS returns the distance squared between each point.
! NIDX returns the sorted indexes of the closest points in PERMCOORDSB by
! permutation


USE KDTREE2_MODULE, ONLY : KDTREE2_RESULT, KDTREE2_N_NEAREST

IMPLICIT NONE
INTEGER, INTENT(IN) :: NATOMS, NN
DOUBLE PRECISION, INTENT(IN) :: COORDSA(3*NATOMS)
LOGICAL, INTENT(IN) :: BULKT

DOUBLE PRECISION, INTENT(OUT) :: NDISTS(NN,NATOMS)
INTEGER, INTENT(OUT) :: NIDX(NN,NATOMS)

INTEGER NDUMMY, I, J, K, IND, M
TYPE(KDTREE2_RESULT), ALLOCATABLE :: RESULTS(:)

M = MIN(NATOMS, NN)
ALLOCATE(RESULTS(M))

NIDX = -1
NDISTS = HUGE(1.D0)

NDUMMY = 0
DO I=1,NPERMGROUP
    DO J=1,NPERMSIZE(I)
        IND = PERMGROUP(NDUMMY+J)
        CALL KDTREE2_N_NEAREST(KDTREES(I)%TREE, & 
        & COORDSA(3*IND-2:3*IND), M, RESULTS)

        NIDX(:M,IND) = RESULTS(1:M)%IDX
        !For some reason RESULTS%DIS doesn't get the right answer...
        !NDISTS(:,IND) = RESULTS(1:NN)%DIS
        DO K=1,M
            NDISTS(K,IND) = SUM((PERMCOORDSB(:,NIDX(K,IND),I) - &
                                & COORDSA(3*IND-2:3*IND))**2)
        ENDDO
        IF(BULKT) THEN
            NIDX(:M,IND) = MODULO(NIDX(:M,IND)-1,NPERMSIZE(I)) + 1
        END IF
    ENDDO
    NDUMMY = NDUMMY + NPERMSIZE(I)
ENDDO

END SUBROUTINE NEARESTNEIGHBOURSTREE

SUBROUTINE FINDBESTPERMTREE(NDISTS, NIDX, NATOMS, NN, PERM, INFO)

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, NN, NIDX(NN, NATOMS)
DOUBLE PRECISION :: NDISTS(NN, NATOMS)
INTEGER, INTENT(OUT) :: PERM(NATOMS), INFO

INTEGER*8 :: KK(NATOMS*NN), CC(NATOMS*NN)
INTEGER*8 :: FIRST(NATOMS+1), X(NATOMS), Y(NATOMS)
INTEGER*8 :: U(NATOMS), V(NATOMS), H, N8, SZ8, D
INTEGER N,M,I,J,K,I1,NPERM,NDUMMY

INFO=0

NDUMMY=0
DO NPERM=1,NPERMGROUP

    N = NPERMSIZE(NPERM)
    M = NN
    IF(N.LE.NN) M=N
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
            KK(I+K) = NIDX(I,J+NDUMMY)
            CC(I+K) = INT(NDISTS(I,J+NDUMMY)*PSCALE, 8)
        ENDDO
    ENDDO

    CALL JOVOSAP(N8, SZ8, CC(:M*N), KK(:M*N), FIRST(:N+1), X(:N), Y(:N), U(:N), V(:N), H)
    NLAP = NLAP + 1

    DO I=1,N
        IF (Y(I).GT.N) THEN
            Y(I)=N
            INFO = INFO + 1
        END IF
        IF (Y(I).LT.1) THEN
            Y(I)=1
            INFO = INFO + 1
        END IF
        PERM(PERMGROUP(NDUMMY+I))=PERMGROUP(NDUMMY+Y(I))
    ENDDO

    ! untested!!
    IF (NSETS(NPERM).GT.0) THEN
        DO I=1,N
            DO K=1,NSETS(NPERM)
                PERM(SETS(PERMGROUP(NDUMMY+I),K))=SETS(PERM(PERMGROUP(NDUMMY+Y(I))),K)
            ENDDO
        ENDDO
    ENDIF

    NDUMMY = NDUMMY + NPERMSIZE(NPERM)
ENDDO

END SUBROUTINE

SUBROUTINE REALLOCATEARRAYS(NATOMS, NUMSTRUCTS, BULKT)

IMPLICIT NONE

INTEGER, INTENT(IN) :: NATOMS, NUMSTRUCTS
LOGICAL, INTENT(IN) :: BULKT

IF(ALLOCATED(PERMCOORDSB))  DEALLOCATE(PERMCOORDSB)
IF(BULKT) THEN
    ALLOCATE(PERMCOORDSB(3,9*NATOMS,NPERMGROUP))
ELSE
    ALLOCATE(PERMCOORDSB(3,NATOMS,NPERMGROUP))
END IF

IF(ALLOCATED(SAVECOORDSB))  DEALLOCATE(SAVECOORDSB,SAVECOORDSA)
IF(ALLOCATED(SAVERA)) DEALLOCATE(SAVERA,SAVERB,BESTCOORDSA,BESTRMAT,BESTITERS)
ALLOCATE(SAVECOORDSB(3*NATOMS),SAVECOORDSA(3*NATOMS,NUMSTRUCTS), &
 & SAVERB(NATOMS),SAVERA(NATOMS,NUMSTRUCTS),BESTCOORDSA(3*NATOMS,NUMSTRUCTS), & 
 & BESTRMAT(3,3,NUMSTRUCTS),BESTITERS(NUMSTRUCTS))

IF(ALLOCATED(PDUMMYA)) DEALLOCATE(PDUMMYA,PDUMMYB,DUMMYA,DUMMYB,NEWPERM,LPERM)
IF(ALLOCATED(INVPERMGROUP)) DEALLOCATE(INVPERMGROUP)
ALLOCATE(PDUMMYA(3*NATOMS),PDUMMYB(3*NATOMS),DUMMYA(3*NATOMS), & 
 & DUMMYB(3*NATOMS),NEWPERM(NATOMS),LPERM(NATOMS),INVPERMGROUP(NATOMS))

IF(ALLOCATED(DUMMYDISTS)) DEALLOCATE(DUMMYDISTS,DUMMYNEARDISTS,PDUMMYND, &
 & DUMMYDISPS,DUMMYIDX,DINVIDX,DUMMYNEARIDX,DUMMYLDISTS,DUMMYNEARLDISTS, &
 & DUMMYLDISTS2,DUMMYDOTDISP)
ALLOCATE(DUMMYDISTS(PMAXNEI*NATOMS,NPERMGROUP),DUMMYNEARDISTS(NATOMS), &
 & PDUMMYND(NATOMS),DUMMYIDX(PMAXNEI*NATOMS,NPERMGROUP),DUMMYNEARIDX(NATOMS), &
 & DINVIDX(NATOMS*NATOMS,NPERMGROUP),DUMMYLDISTS(PMAXNEI*NATOMS,NPERMGROUP), &
 & DUMMYNEARLDISTS(NATOMS),DUMMYLDISTS2(PMAXNEI*NATOMS,NPERMGROUP), &
 & DUMMYDISPS(3,NATOMS*PMAXNEI,NPERMGROUP),DUMMYDOTDISP(4,NATOMS*PMAXNEI,NPERMGROUP))

END SUBROUTINE REALLOCATEARRAYS

SUBROUTINE SETCLUSTER(INVERT)

USE COMMONS, ONLY : MYUNIT,NFREEZE,GEOMDIFFTOL,ORBITTOL,FREEZE,PULLT,TWOD,  &
    &   EFIELDT,AMBERT,QCIAMBERT,AMBER12T,CHRMMT,STOCKT,CSMT,PERMDIST,      &
    &   LOCALPERMDIST,LPERMDIST,OHCELLT,QCIPERMCHECK,PERMOPT,PERMINVOPT,    &
    &   NOINVERSION,GTHOMSONT,MKTRAPT,MULLERBROWNT,RIGID,OHCELLT

IMPLICIT NONE

LOGICAL, INTENT(IN) :: INVERT

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
PERMINVOPT = MERGE(.TRUE.,.FALSE.,INVERT)
NOINVERSION = .FALSE.
GTHOMSONT = .FALSE.
MKTRAPT = .FALSE.
MULLERBROWNT = .FALSE.
RIGID = .FALSE.
OHCELLT = .FALSE.

END SUBROUTINE SETCLUSTER

SUBROUTINE SETBULK(INVERT)

USE COMMONS, ONLY : MYUNIT,NFREEZE,GEOMDIFFTOL,ORBITTOL,FREEZE,PULLT,TWOD,  &
    &   EFIELDT,AMBERT,QCIAMBERT,AMBER12T,CHRMMT,STOCKT,CSMT,PERMDIST,      &
    &   LOCALPERMDIST,LPERMDIST,OHCELLT,QCIPERMCHECK,PERMOPT,PERMINVOPT,    &
    &   NOINVERSION,GTHOMSONT,MKTRAPT,MULLERBROWNT,RIGID,OHCELLT

IMPLICIT NONE

LOGICAL, INTENT(IN) :: INVERT

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
PERMDIST = .FALSE.
LOCALPERMDIST = .FALSE.
LPERMDIST = .FALSE.
QCIPERMCHECK = .FALSE.
PERMOPT = .TRUE.
PERMINVOPT = .FALSE.
NOINVERSION = .FALSE.
GTHOMSONT = .FALSE.
MKTRAPT = .FALSE.
MULLERBROWNT = .FALSE.
RIGID = .FALSE.
OHCELLT = MERGE(.TRUE.,.FALSE.,INVERT)

END SUBROUTINE SETBULK

END MODULE

INCLUDE "bulkmindist.f90"
INCLUDE "minpermdist.f90"
INCLUDE "minperm.f90"
INCLUDE "newmindist.f90"
INCLUDE "orient.f90"