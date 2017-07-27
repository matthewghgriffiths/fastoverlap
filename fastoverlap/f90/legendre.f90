! FCNPAK
!
!   DRIVER FOR TESTING CMLIB ROUTINES
!     SINGLE PRECISION LEGENDRE ROUTINES IN FCNPAK
!     DOUBLE PRECISION LEGENDRE ROUTINES IN FCNPAK
!
!    ONE INPUT DATA CARD IS REQUIRED
!         READ(LIN,1) KPRINT,TIMES
!    1    FORMAT(I1,E10.0)
!
!     KPRINT = 0   NO PRINTING
!              1   NO PRINTING FOR PASSED TESTS, SHORT MESSAGE
!                  FOR FAILED TESTS
!              2   PRINT SHORT MESSAGE FOR PASSED TESTS, FULLER
!                  INFORMATION FOR FAILED TESTS
!              3   PRINT COMPLETE QUICK-CHECK RESULTS
!
!                ***IMPORTANT NOTE***
!         ALL QUICK CHECKS USE ROUTINES R2MACH AND D2MACH
!         TO SET THE ERROR TOLERANCES.
!     TIMES IS A CONSTANT MULTIPLIER THAT CAN BE USED TO SCALE THE
!     VALUES OF R1MACH AND D1MACH SO THAT
!               R2MACH(I) = R1MACH(I) * TIMES   FOR I=3,4,5
!               D2MACH(I) = D1MACH(I) * TIMES   FOR I=3,4,5
!     THIS MAKES IT EASILY POSSIBLE TO CHANGE THE ERROR TOLERANCES
!     USED IN THE QUICK CHECKS.
!     IF TIMES .LE. 0.0 THEN TIMES IS DEFAULTED TO 1.0
!
!              ***END NOTE***
!


!FCNPAK   CMLIB
!========================================================================
!                           F C N P A K
!========================================================================
!  
!  
!Introduction
!============
!  
!FCNPAK is a package of mathematical function subroutines. Its purpose is
!to provide functions that are not readily available  elsewhere  or  that
!are available only under restrictive licensing agreements. All  programs
!in FCNPAK are coded  in  standard  Fortran.  They  are  designed  to  be
!installed with little difficulty on  most  conventional  computers.  The
!developers of FCNPAK - D. W. Lozier (NBS 711, lozier@nist.gov) and  J. M.
!Smith (formerly of NBS 715)- would appreciate receiving feedback  from
!users about the usefulness and applicability of the software.
!  
!  
!Legendre Functions
!==================
!  
!The package contains subroutines for computing the  associated  Legendre
!functions, or Ferrers functions,
!  
!         P-subNU-superMU(cosTHETA), Q-subNU-superMU(cosTHETA)
!  
!as well as the normalized Legendre polynomial
!  
!                     PBAR-subNU-superMU(cosTHETA)
!  
!in the ranges
!                        0 .LT. THETA .LE. PI/2
!                           MU = 0, 1, 2, ...
!                    -(1/2) .LE. NU .LT. INFINITY .
!Negative integral values of  MU  may  be  specified  also  for  P-subNU-
!superMU. NU is restricted to integers for PBAR-subNU-superMU.
!  
!An unusual feature of the FCNPAK subroutines for Legendre  functions  is
!the use of extended-range arithmetic, a software extension  of  ordinary
!floating-point arithmetic that greatly increases the exponent  range  of
!the representable numbers. In consequence, we avoid the need for scaling
!the solutions to lie within the exponent range of the  most  restrictive
!manufacturer's hardware. The increased exponent  range  is  achieved  by
!allocating an integer storage location together with each floating-point
!storage location. The increase in  the  time  required  to  execute  the
!algorithms in extended-range arithmetic depends on the  application.  In
!the case of the normalized Legendre polynomials, testing shows it to  be
!about a factor of two.  This is compensated in part by the lack  of  any
!need for scaling operations in the algorithms  for  the  functions.  The
!resulting function  values  are  supplied  in  ordinary  floating  point
!whenever possible.
!  
!P and Q are solutions of the associated Legendre  equation.  Definitions
!and properties of these and other solutions are supplied in
!      Bateman Manuscript Project, "Higher Transcendental Functions,"
!(A. Erdelyi, Ed.), v. 1, McGraw-Hill, New York, 1953
!      National Bureau of Standards, "Handbook of Mathematical
!Functions," AMS 55, (M. Abramowitz and I. A. Stegun, Eds.), U. S. GPO,
!Washington, D. C., 1964
!      F. W. J. Olver, "Asymptotics and Special Functions,"  Academic
!Press, New York, 1974.
!Algorithmic details for the subroutines in FCNPAK are supplied in
!      J. M. Smith, F. W. J. Olver, and D. W. Lozier, "Extended-Range
!Arithmetic and Normalized Legendre Polynomials," ACM Trans. on Math.
!Softw., v. 7, no. 1, March 1981
!      F. W. J. Olver and J. M. Smith, "Associated Legendre Functions On
!The Cut," J. Comp. Physics, v. 51, no. 3, September 1983.
!  
!The names of the subroutines used for Legendre functions are
!
!     XDLEGF     XSLEGF     XDNRMP     XSNRMP
! 
!where the LEGF subroutines compute the double and single precision
!associated Legendre functions and the NRMP subroutines compute the
!double and single precision normalized Legendre polynomials.   The
!algorithms used  in the LEGF  subroutines are  described in Smith,
!Olver and  Lozier (1983), and  those used in the  NRMP subroutines
!are described in  Olver and Smith (1983).  These subroutines  have
!been incorporated in the CMLIB library and (with slight changes of
!name) SLATEC library; see http://gams.nist.gov/serve.cgi/Class/C9.
!
!Two subroutines
!
!     XDCSRT     XSCSRT
!
!are provided for testing purposes. They can be used to construct tests
!of computed results based on Casoratian relations.
!
!The names flagged by * in the list
!
!      XDADD*    XDC210*     XDPNRM     XDQMU      XDSET*      
!      XDADJ*     XDPMU      XDPQNU     XDQNU
!      XDCON*    XDPMUP       XDPSI     XDRED*
!                                       
!      XSADD*    XSC210*     XSPNRM     XSQMU      XSSET*      
!      XSADJ*     XSPMU      XSPQNU     XSQNU
!      XSCON*    XSPMUP       XSPSI     XSRED*
!
!support extended-range arithmetic, and the others are subsidiary to the
!Legendre subroutines. The names that begin with XD are double-precision
!subroutines, and those which begin with XS are single-precision subroutines.
!XDSET  (and XSSET)  are initialization subroutines that must be called 
!before any other extended-range subroutine is called. The Legendre function
!subroutines do this for the user.
! 
!(updated July 21, 2010)



SUBROUTINE XDNRMP(NU,MU1,MU2,DARG,MODE,DPN,IPN,ISIG)
!***BEGIN PROLOGUE  XDNRMP
!***DATE WRITTEN   820712   (YYMMDD)
!***REVISION DATE  871110   (YYMMDD)
!***CATEGORY NO.  C3a2,C9
!***KEYWORDS  LEGENDRE FUNCTIONS
!***AUTHOR  LOZIER, DANIEL W. (NATIONAL BUREAU OF STANDARDS)
!           SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
!***PURPOSE  TO COMPUTE THE NORMALIZED LEGENDRE POLYNOMIAL
!***DESCRIPTION
!
!        SUBROUTINE TO CALCULATE NORMALIZED LEGENDRE POLYNOMIALS
!        (XSNRMP is single-precision version)
!        XDNRMP calculates normalized Legendre polynomials of varying
!        order and fixed argument and degree. The order MU and degree
!        NU are nonegative integers and the argument is real. Because
!        the algorithm requires the use of numbers outside the normal
!        machine range, this subroutine employs a special arithmetic
!        called extended-range arithmetic. See J.M. Smith, F.W.J. Olver,
!        and D.W. Lozier, Extended-Range Arithmetic and Normalized
!        Legendre Polynomials, ACM Transactions on Mathematical Soft-
!        ware, 93-105, March 1981, for a complete description of the
!        algorithm and special arithmetic. Also see program comments
!        in XDSET.
!
!        The normalized Legendre polynomials are multiples of the
!        associated Legendre polynomials of the first kind where the
!        normalizing coefficients are chosen so as to make the integral
!        from -1 to 1 of the square of each function equal to 1. See
!        E. Jahnke, F. Emde and F. Losch, Tables of Higher Functions,
!        McGraw-Hill, New York, 1960, p. 121.
!
!        The input values to XDNRMP are NU, MU1, MU2, DARG, and MODE.
!        These must satisfy
!          1. NU .GE. 0 specifies the degree of the normalized Legendre
!             polynomial that is wanted.
!          2. MU1 .GE. 0 specifies the lowest-order normalized Legendre
!             polynomial that is wanted.
!          3. MU2 .GE. MU1 specifies the highest-order normalized Leg-
!             endre polynomial that is wanted.
!         4a. MODE = 1 and -1.0D0 .LE. DARG .LE. 1.0D0 specifies that
!             Normalized Legendre(NU, MU, DARG) is wanted for MU = MU1,
!             MU1 + 1, ..., MU2.
!         4b. MODE = 2 and -3.14159... .LT. DARG .LT. 3.14159... spec-
!             ifies that Normalized Legendre(NU, MU, DCOS(DARG)) is want-
!             ed for MU = MU1, MU1 + 1, ..., MU2.
!
!        The output of XDNRMP consists of the two vectors DPN and IPN
!        and the error estimate ISIG. The computed values are stored as
!        extended-range numbers such that
!             (DPN(1),IPN(1))=NORMALIZED LEGENDRE(NU,MU1,DX)
!             (DPN(2),IPN(2))=NORMALIZED LEGENDRE(NU,MU1+1,DX)
!                .
!                .
!             (DPN(K),IPN(K))=NORMALIZED LEGENDRE(NU,MU2,DX)
!        where K = MU2 - MU1 + 1 and DX = DARG or DCOS(DARG) according
!        to whether MODE = 1 or 2. Finally, ISIG is an estimate of the
!        number of decimal digits lost through rounding errors in the
!        computation. For example if DARG is accurate to 12 significant
!        decimals, then the computed function values are accurate to
!        12 - ISIG significant decimals (except in neighborhoods of
!        zeros).
!
!        The interpretation of (DPN(I),IPN(I)) is DPN(I)*(IR**IPN(I))
!        where IR is the internal radix of the computer arithmetic. When
!        IPN(I) = 0 the value of the normalized Legendre polynomial is
!        contained entirely in DPN(I) and subsequent double-precision
!        computations can be performed without further consideration of
!        extended-range arithmetic. However, if IPN(I) .NE. 0 the corre-
!        sponding value of the normalized Legendre polynomial cannot be
!        represented in double-precision because of overflow or under-
!        flow. THE USER MUST TEST IPN(I) IN HIS/HER PROGRAM. In the event
!        that IPN(I) is nonzero, the user could rewrite his/her program
!        to use extended range arithmetic.
!
!
!
!        The interpretation of (DPN(I),IPN(I)) can be changed to
!        DPN(I)*(10**IPN(I)) by calling the extended-range subroutine
!        XDCON. This should be done before printing the computed values.
!        As an example of usage, the Fortran coding
!              J = K
!              DO 20 I = 1, K
!              CALL XDCON(DPN(I), IPN(I))
!              PRINT 10, DPN(I), IPN(I)
!           10 FORMAT(1X, D30.18 , I15)
!              IF ((IPN(I) .EQ. 0) .OR. (J .LT. K)) GO TO 20
!              J = I - 1
!           20 CONTINUE
!        will print all computed values and determine the largest J
!        such that IPN(1) = IPN(2) = ... = IPN(J) = 0. Because of the
!        change of representation caused by calling XDCON, (DPN(I),
!        IPN(I)) for I = J+1, J+2, ... cannot be used in subsequent
!        extended-range computations.
!
!***REFERENCES  (SEE DESCRIPTION ABOVE)
!***ROUTINES CALLED  XDADD, XDADJ, XDRED, XDSET
!***END PROLOGUE  XDNRMP
      INTEGER NU, MU1, MU2, MODE, IPN, ISIG
      DOUBLE PRECISION DARG, DPN
      DIMENSION DPN(*), IPN(*)
      DOUBLE PRECISION C1,C2,P,P1,P2,P3,S,SX,T,TX,X,DK
! CALL XDSET TO INITIALIZE EXTENDED-RANGE ARITHMETIC (SEE XDSET
! LISTING FOR DETAILS)
!***FIRST EXECUTABLE STATEMENT  XDNRMP
      CALL XDSET (0, 0, 0.0D0, 0)
!
!        TEST FOR PROPER INPUT VALUES.
!
      IF (NU.LT.0) GO TO 110
      IF (MU1.LT.0) GO TO 110
      IF (MU1.GT.MU2) GO TO 110
      IF (NU.EQ.0) GO TO 90
      IF (MODE.LT.1 .OR. MODE.GT.2) GO TO 110
      GO TO (10, 20), MODE
   10 IF (ABS(DARG).GT.1.0D0) GO TO 120
      IF (ABS(DARG).EQ.1.0D0) GO TO 90
      X = DARG
      SX = SQRT((1.0D0+ABS(X))*((0.5D0-ABS(X))+0.5D0))
      TX = X/SX
      ISIG = LOG10(2.0D0*DBLE(FLOAT(NU))*(5.0D0+TX**2))
      GO TO 30
   20 IF (ABS(DARG).GT.4.0D0*DATAN(1.0D0)) GO TO 120
      IF (DARG.EQ.0.0D0) GO TO 90
      X = COS(DARG)
      SX = ABS(SIN(DARG))
      TX = X/SX
      ISIG = LOG10(2.0D0*DBLE(FLOAT(NU))*(5.0D0+ABS(DARG*TX)))
!
!        BEGIN CALCULATION
!
   30 MU = MU2
      I = MU2 - MU1 + 1
!
!        IF MU.GT.NU, NORMALIZED LEGENDRE(NU,MU,X)=0.
!
   40 IF (MU.LE.NU) GO TO 50
      DPN(I) = 0.0D0
      IPN(I) = 0
      I = I - 1
      MU = MU - 1
      IF (I .GT. 0) GO TO 40
      ISIG = 0
      GO TO 160
   50 MU = NU
!
!        P1 = 0. = NORMALIZED LEGENDRE(NU,NU+1,X)
!
      P1 = 0.0D0
      IP1 = 0
!
!        CALCULATE P2 = NORMALIZED LEGENDRE(NU,NU,X)
!
      P2 = 1.0D0
      IP2 = 0
      P3 = 0.5D0
      DK = 2.0D0
      DO 60 J=1,NU
        P3 = ((DK+1.0D0)/DK)*P3
        P2 = P2*SX
        CALL XDADJ(P2, IP2)
        DK = DK + 2.0D0
   60 CONTINUE
      P2 = P2*SQRT(P3)
      CALL XDADJ(P2, IP2)
      S = 2.0D0*TX
      T = 1.0D0/DBLE(FLOAT(NU))
      IF (MU2.LT.NU) GO TO 70
      DPN(I) = P2
      IPN(I) = IP2
      I = I - 1
      IF (I .EQ. 0) GO TO 140
!
!        RECURRENCE PROCESS
!
   70 P = DBLE(FLOAT(MU))*T
      C1 = 1.0D0/DSQRT((1.0D0-P+T)*(1.0D0+P))
      C2 = S*P*C1*P2
      C1 = -DSQRT((1.0D0+P+T)*(1.0D0-P))*C1*P1
      CALL XDADD(C2, IP2, C1, IP1, P, IP)
      MU = MU - 1
      IF (MU.GT.MU2) GO TO 80
!
!        STORE IN ARRAY DPN FOR RETURN TO CALLING ROUTINE.
!
      DPN(I) = P
      IPN(I) = IP
      I = I - 1
      IF (I .EQ. 0) GO TO 140
   80 P1 = P2
      IP1 = IP2
      P2 = P
      IP2 = IP
      IF (MU.LE.MU1) GO TO 140
      GO TO 70
!
!        SPECIAL CASE WHEN X=-1 OR +1, OR NU=0.
!
   90 K = MU2 - MU1 + 1
      DO 100 I=1,K
        DPN(I) = 0.0D0
        IPN(I) = 0
  100 CONTINUE
      ISIG = 0
      IF (MU1.GT.0) GO TO 160
      ISIG = 1
      DPN(1) = DSQRT(DBLE(FLOAT(NU))+0.5D0)
      IPN(1) = 0
      IF (MOD(NU,2).EQ.0) GO TO 160
      IF (MODE.EQ.1 .AND. DARG.EQ.1.0D0) GO TO 160
      IF (MODE.EQ.2) GO TO 160
      DPN(1) = -DPN(1)
      GO TO 160
!
!          ERROR PRINTOUTS AND TERMINATION.
!
  110 WRITE(*,*) 'Err in XDNRMP...NU,MU1,MU2 or MODE not valid'
      GO TO 130
  120 WRITE(*,*) 'Err in XDNRMP...DARG out of range'
  130 RETURN
!
!        RETURN TO CALLING PROGRAM
!
  140 K = MU2 - MU1 + 1
      DO 150 I=1,K
        CALL XDRED(DPN(I),IPN(I))
  150 CONTINUE
  160 RETURN
      END
      SUBROUTINE XDADD(X,IX,Y,IY,Z,IZ)
!***BEGIN PROLOGUE  XDADD
!***DATE WRITTEN   820712   (YYMMDD)
!***REVISION DATE  831027   (YYMMDD)
!***CATEGORY NO.  A3d
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  LOZIER, DANIEL W. (NATIONAL BUREAU OF STANDARDS)
!           SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
!***PURPOSE  TO PROVIDE DOUBLE-PRECISION FLOATING-POINT ARITHMETIC
!            WITH AN EXTENDED EXPONENT RANGE
!***DESCRIPTION
!     DOUBLE PRECISION X, Y, Z
!     INTEGER IX, IY, IZ
!
!                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
!                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
!                  BEFORE RETURNING. THE INPUT OPERANDS
!                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
!                  PRINCIPAL PARTS MUST SATISFY
!                  RADIX**(-2L).LE.ABS(X).LE.RADIX**(2L),
!                  RADIX**(-2L).LE.ABS(Y).LE.RADIX**(2L).
!
!***REFERENCES  (PROGRAM LISTING FOR XDSET)
!***ROUTINES CALLED  XDADJ
!***COMMON BLOCKS    XDBLK2
!***END PROLOGUE  XDADD
      DOUBLE PRECISION X, Y, Z
      INTEGER IX, IY, IZ
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /XDBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE/XDBLK2/
      DOUBLE PRECISION S, T
!
!   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
! ARE
!     (1) 1 .LT. L .LE. 0.5D0*LOGR(0.5D0*DZERO)
!
!     (2) NRADPL .LT. L .LE. KMAX/6
!
!     (3) KMAX .LE. (2**NBITS - 4*L - 1)/2
!
! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE XDSET.
!
!***FIRST EXECUTABLE STATEMENT  XDADD
      IF (X.NE.0.0D0) GO TO 10
      Z = Y
      IZ = IY
      GO TO 220
   10 IF (Y.NE.0.0D0) GO TO 20
      Z = X
      IZ = IX
      GO TO 220
   20 CONTINUE
      IF (IX.GE.0 .AND. IY.GE.0) GO TO 40
      IF (IX.LT.0 .AND. IY.LT.0) GO TO 40
      IF (IABS(IX).LE.6*L .AND. IABS(IY).LE.6*L) GO TO 40
      IF (IX.GE.0) GO TO 30
      Z = Y
      IZ = IY
      GO TO 220
   30 CONTINUE
      Z = X
      IZ = IX
      GO TO 220
   40 I = IX - IY
      IF (I) 80, 50, 90
   50 IF (DABS(X).GT.1.0D0 .AND. DABS(Y).GT.1.0D0) GO TO 60
      IF (DABS(X).LT.1.0D0 .AND. DABS(Y).LT.1.0D0) GO TO 70
      Z = X + Y
      IZ = IX
      GO TO 220
   60 S = X/RADIXL
      T = Y/RADIXL
      Z = S + T
      IZ = IX + L
      GO TO 220
   70 S = X*RADIXL
      T = Y*RADIXL
      Z = S + T
      IZ = IX - L
      GO TO 220
   80 S = Y
      IS = IY
      T = X
      GO TO 100
   90 S = X
      IS = IX
      T = Y
  100 CONTINUE
!
!  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE
! LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL
! PART OF THE OTHER INPUT IS STORED IN T.
!
      I1 = IABS(I)/L
      I2 = MOD(IABS(I),L)
      IF (DABS(T).GE.RADIXL) GO TO 130
      IF (DABS(T).GE.1.0D0) GO TO 120
      IF (RADIXL*DABS(T).GE.1.0D0) GO TO 110
      J = I1 + 1
      T = T*RADIX**(L-I2)
      GO TO 140
  110 J = I1
      T = T*RADIX**(-I2)
      GO TO 140
  120 J = I1 - 1
      IF (J.LT.0) GO TO 110
      T = T*RADIX**(-I2)/RADIXL
      GO TO 140
  130 J = I1 - 2
      IF (J.LT.0) GO TO 120
      T = T*RADIX**(-I2)/RAD2L
  140 CONTINUE
!
!  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE
! AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT
! OF T.  THE SHIFTED VALUE OF T SATISFIES
!
!       RADIX**(-2*L) .LE. DABS(T) .LE. 1.0D0
!
! AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE.
!
      IF (J.EQ.0) GO TO 190
      IF (DABS(S).GE.RADIXL .OR. J.GT.3) GO TO 150
      IF (DABS(S).GE.1.0D0) GO TO (180, 150, 150), J
      IF (RADIXL*DABS(S).GE.1.0D0) GO TO (180, 170, 150), J
      GO TO (180, 170, 160), J
  150 Z = S
      IZ = IS
      GO TO 220
  160 S = S*RADIXL
  170 S = S*RADIXL
  180 S = S*RADIXL
  190 CONTINUE
!
!   AT THIS POINT, THE REMAINING DIFFERENCE IN THE
! AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT
! OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED
! RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE
! SUM.
!
      IF (DABS(S).GT.1.0D0 .AND. DABS(T).GT.1.0D0) GO TO 200
      IF (DABS(S).LT.1.0D0 .AND. DABS(T).LT.1.0D0) GO TO 210
      Z = S + T
      IZ = IS - J*L
      GO TO 220
  200 S = S/RADIXL
      T = T/RADIXL
      Z = S + T
      IZ = IS - J*L + L
      GO TO 220
  210 S = S*RADIXL
      T = T*RADIXL
      Z = S + T
      IZ = IS - J*L - L
  220 CALL XDADJ(Z, IZ)
      RETURN
      END
      SUBROUTINE XDADJ(X,IX)
!***BEGIN PROLOGUE  XDADJ
!***DATE WRITTEN   820712   (YYMMDD)
!***REVISION DATE  831027   (YYMMDD)
!***CATEGORY NO.  A3d
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  LOZIER, DANIEL W. (NATIONAL BUREAU OF STANDARDS)
!           SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
!***PURPOSE  TO PROVIDE DOUBLE-PRECISION FLOATING-POINT ARITHMETIC
!            WITH AN EXTENDED EXPONENT RANGE
!***DESCRIPTION
!     DOUBLE PRECISION X
!     INTEGER IX
!
!                  TRANSFORMS (X,IX) SO THAT
!                  RADIX**(-L) .LE. DABS(X) .LT. RADIX**L.
!                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
!                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
!                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
!
!***REFERENCES  (PROGRAM LISTING FOR XDSET)
!***ROUTINES CALLED  XERROR
!***COMMON BLOCKS    XDBLK2
!***END PROLOGUE  XDADJ
      DOUBLE PRECISION X
      INTEGER IX
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /XDBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /XDBLK2/
!
!   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
! IS
!     2*L .LE. KMAX
!
! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE XDSET.
!
!***FIRST EXECUTABLE STATEMENT  XDADJ
      IF (X.EQ.0.0D0) GO TO 50
      IF (DABS(X).GE.1.0D0) GO TO 20
      IF (RADIXL*DABS(X).GE.1.0D0) GO TO 60
      X = X*RAD2L
      IF (IX.LT.0) GO TO 10
      IX = IX - L2
      GO TO 70
   10 IF (IX.LT.-KMAX+L2) GO TO 40
      IX = IX - L2
      GO TO 70
   20 IF (DABS(X).LT.RADIXL) GO TO 60
      X = X/RAD2L
      IF (IX.GT.0) GO TO 30
      IX = IX + L2
      GO TO 70
   30 IF (IX.GT.KMAX-L2) GO TO 40
      IX = IX + L2
      GO TO 70
   40 WRITE(*,*) 'Err in XDADJ...overflow in auxiliary index'
      RETURN
   50 IX = 0
   60 IF (IABS(IX).GT.KMAX) GO TO 40
   70 RETURN
      END

SUBROUTINE XDRED(X,IX)
!***BEGIN PROLOGUE  XDRED
!***DATE WRITTEN   820712   (YYMMDD)
!***REVISION DATE  831027   (YYMMDD)
!***CATEGORY NO.  A3d
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  LOZIER, DANIEL W. (NATIONAL BUREAU OF STANDARDS)
!           SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
!***PURPOSE  TO PROVIDE DOUBLE-PRECISION FLOATING-POINT ARITHMETIC
!            WITH AN EXTENDED EXPONENT RANGE
!***DESCRIPTION
!     DOUBLE PRECISION X
!     INTEGER IX
!
!                  IF
!                  RADIX**(-2L) .LE. (DABS(X),IX) .LE. RADIX**(2L)
!                  THEN XDRED TRANSFORMS (X,IX) SO THAT IX=0.
!                  IF (X,IX) IS OUTSIDE THE ABOVE RANGE,
!                  THEN XDRED TAKES NO ACTION.
!                  THIS SUBROUTINE IS USEFUL IF THE
!                  RESULTS OF EXTENDED-RANGE CALCULATIONS
!                  ARE TO BE USED IN SUBSEQUENT ORDINARY
!                  DOUBLE-PRECISION CALCULATIONS.
!
!***REFERENCES  (PROGRAM LISTING FOR XDSET)
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    XDBLK2
!***END PROLOGUE  XDRED
      DOUBLE PRECISION X
      INTEGER IX
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R, XA
      INTEGER L, L2, KMAX
      COMMON /XDBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /XDBLK2/
!
!***FIRST EXECUTABLE STATEMENT  XDRED
      IF (X.EQ.0.0D0) GO TO 90
      XA = DABS(X)
      IF (IX.EQ.0) GO TO 70
      IXA = IABS(IX)
      IXA1 = IXA/L2
      IXA2 = MOD(IXA,L2)
      IF (IX.GT.0) GO TO 40
   10 CONTINUE
      IF (XA.GT.1.0D0) GO TO 20
      XA = XA*RAD2L
      IXA1 = IXA1 + 1
      GO TO 10
   20 XA = XA/RADIX**IXA2
      IF (IXA1.EQ.0) GO TO 70
      DO 30 I=1,IXA1
        IF (XA.LT.1.0D0) GO TO 100
        XA = XA/RAD2L
   30 CONTINUE
      GO TO 70
!
   40 CONTINUE
      IF (XA.LT.1.0D0) GO TO 50
      XA = XA/RAD2L
      IXA1 = IXA1 + 1
      GO TO 40
   50 XA = XA*RADIX**IXA2
      IF (IXA1.EQ.0) GO TO 70
      DO 60 I=1,IXA1
        IF (XA.GT.1.0D0) GO TO 100
        XA = XA*RAD2L
   60 CONTINUE
   70 IF (XA.GT.RAD2L) GO TO 100
      IF (XA.GT.1.0D0) GO TO 80
      IF (RAD2L*XA.LT.1.0D0) GO TO 100
   80 X = DSIGN(XA,X)
   90 IX = 0
  100 RETURN
      END
SUBROUTINE XDSET(IRAD,NRADPL,DZERO,NBITS)
!***BEGIN PROLOGUE  XDSET
!***DATE WRITTEN   820712   (YYMMDD)
!***REVISION DATE  871110   (YYMMDD)
!***CATEGORY NO.  A3d
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  LOZIER, DANIEL W. (NATIONAL BUREAU OF STANDARDS)
!           SMITH, JOHN M. (NBS AND GEORGE MASON UNIVERSITY)
!***PURPOSE  TO PROVIDE DOUBLE-PRECISION FLOATING-POINT ARITHMETIC
!            WITH AN EXTENDED EXPONENT RANGE
!***DESCRIPTION
!
!   SUBROUTINE  XDSET  MUST BE CALLED PRIOR TO CALLING ANY OTHER
! EXTENDED-RANGE SUBROUTINE. IT CALCULATES AND STORES SEVERAL
! MACHINE-DEPENDENT CONSTANTS IN COMMON BLOCKS. THE USER MUST
! SUPPLY FOUR CONSTANTS THAT PERTAIN TO HIS PARTICULAR COMPUTER.
! THE CONSTANTS ARE
!
!          IRAD = THE INTERNAL BASE OF DOUBLE-PRECISION
!                 ARITHMETIC IN THE COMPUTER.
!        NRADPL = THE NUMBER OF RADIX PLACES CARRIED IN
!                 THE DOUBLE-PRECISION REPRESENTATION.
!         DZERO = THE SMALLEST OF 1/DMIN, DMAX, DMAXLN WHERE
!                 DMIN = THE SMALLEST POSITIVE DOUBLE-PRECISION
!                 NUMBER OR AN UPPER BOUND TO THIS NUMBER,
!                 DMAX = THE LARGEST DOUBLE-PRECISION NUMBER
!                 OR A LOWER BOUND TO THIS NUMBER,
!                 DMAXLN = THE LARGEST DOUBLE-PRECISION NUMBER
!                 SUCH THAT DLOG10(DMAXLN) CAN BE COMPUTED BY THE
!                 FORTRAN SYSTEM (ON MOST SYSTEMS DMAXLN = DMAX).
!         NBITS = THE NUMBER OF BITS (EXCLUSIVE OF SIGN) IN
!                 AN INTEGER COMPUTER WORD.
!
! ALTERNATIVELY, ANY OR ALL OF THE CONSTANTS CAN BE GIVEN
! THE VALUE 0 (0.0D0 FOR DZERO). IF A CONSTANT IS ZERO, XDSET TRIES
! TO ASSIGN AN APPROPRIATE VALUE BY CALLING I1MACH
! (SEE P.A.FOX, A.D.HALL, N.L.SCHRYER, ALGORITHM 528 FRAMEWORK
! FOR A PORTABLE LIBRARY, ACM TRANSACTIONS ON MATH SOFTWARE,
! V.4, NO.2, JUNE 1978, 177-188).
!
!   THIS IS THE SETTING-UP SUBROUTINE FOR A PACKAGE OF SUBROUTINES
! THAT FACILITATE THE USE OF EXTENDED-RANGE ARITHMETIC. EXTENDED-RANGE
! ARITHMETIC ON A PARTICULAR COMPUTER IS DEFINED ON THE SET OF NUMBERS
! OF THE FORM
!
!               (X,IX) = X*RADIX**IX
!
! WHERE X IS A DOUBLE-PRECISION NUMBER CALLED THE PRINCIPAL PART,
! IX IS AN INTEGER CALLED THE AUXILIARY INDEX, AND RADIX IS THE
! INTERNAL BASE OF THE DOUBLE-PRECISION ARITHMETIC.  OBVIOUSLY,
! EACH REAL NUMBER IS REPRESENTABLE WITHOUT ERROR BY MORE THAN ONE
! EXTENDED-RANGE FORM.  CONVERSIONS BETWEEN  DIFFERENT FORMS ARE
! ESSENTIAL IN CARRYING OUT ARITHMETIC OPERATIONS.  WITH THE CHOICE
! OF RADIX WE HAVE MADE, AND THE SUBROUTINES WE HAVE WRITTEN, THESE
! CONVERSIONS ARE PERFORMED WITHOUT ERROR (AT LEAST ON MOST COMPUTERS).
! (SEE SMITH, J.M., OLVER, F.W.J., AND LOZIER, D.W., EXTENDED-RANGE
! ARITHMETIC AND NORMALIZED LEGENDRE POLYNOMIALS, ACM TRANSACTIONS ON
! MATHEMATICAL SOFTWARE, MARCH 1981).
!
!   AN EXTENDED-RANGE NUMBER  (X,IX)  IS SAID TO BE IN ADJUSTED FORM IF
! X AND IX ARE ZERO OR
!
!           RADIX**(-L) .LE. DABS(X) .LT. RADIX**L
!
! IS SATISFIED, WHERE L IS A COMPUTER-DEPENDENT INTEGER DEFINED IN THIS
! SUBROUTINE. TWO EXTENDED-RANGE NUMBERS IN ADJUSTED FORM CAN BE ADDED,
! SUBTRACTED, MULTIPLIED OR DIVIDED (IF THE DIVISOR IS NONZERO) WITHOUT
! CAUSING OVERFLOW OR UNDERFLOW IN THE PRINCIPAL PART OF THE RESULT.
! WITH PROPER USE OF THE EXTENDED-RANGE SUBROUTINES, THE ONLY OVERFLOW
! THAT CAN OCCUR IS INTEGER OVERFLOW IN THE AUXILIARY INDEX. IF THIS
! IS DETECTED, THE SOFTWARE CALLS XERROR (A GENERAL ERROR-HANDLING
! FORTRAN SUBROUTINE PACKAGE).
!
!   MULTIPLICATION AND DIVISION IS PERFORMED BY SETTING
!
!                 (X,IX)*(Y,IY) = (X*Y,IX+IY)
! OR
!                 (X,IX)/(Y,IY) = (X/Y,IX-IY).
!
! PRE-ADJUSTMENT OF THE OPERANDS IS ESSENTIAL TO AVOID
! OVERFLOW OR  UNDERFLOW OF THE PRINCIPAL PART. SUBROUTINE
! XDADJ (SEE BELOW) MAY BE CALLED TO TRANSFORM ANY EXTENDED-
! RANGE NUMBER INTO ADJUSTED FORM.
!
!   ADDITION AND SUBTRACTION REQUIRE THE USE OF SUBROUTINE XDADD
! (SEE BELOW).  THE INPUT OPERANDS NEED NOT BE IN ADJUSTED FORM.
! HOWEVER, THE RESULT OF ADDITION OR SUBTRACTION IS RETURNED
! IN ADJUSTED FORM.  THUS, FOR EXAMPLE, IF (X,IX),(Y,IY),
! (U,IU),  AND (V,IV) ARE IN ADJUSTED FORM, THEN
!
!                 (X,IX)*(Y,IY) + (U,IU)*(V,IV)
!
! CAN BE COMPUTED AND STORED IN ADJUSTED FORM WITH NO EXPLICIT
! CALLS TO XDADJ.
!
!   WHEN AN EXTENDED-RANGE NUMBER IS TO BE PRINTED, IT MUST BE
! CONVERTED TO AN EXTENDED-RANGE FORM WITH DECIMAL RADIX.  SUBROUTINE
! XDCON IS PROVIDED FOR THIS PURPOSE.
!
!   THE SUBROUTINES CONTAINED IN THIS PACKAGE ARE
!
!     SUBROUTINE XDADD
! USAGE
!                  CALL XDADD(X,IX,Y,IY,Z,IZ)
! DESCRIPTION
!                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
!                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
!                  BEFORE RETURNING. THE INPUT OPERANDS
!                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
!                  PRINCIPAL PARTS MUST SATISFY
!                  RADIX**(-2L).LE.DABS(X).LE.RADIX**(2L),
!                  RADIX**(-2L).LE.DABS(Y).LE.RADIX**(2L).
!
!     SUBROUTINE XDADJ
! USAGE
!                  CALL XDADJ(X,IX)
! DESCRIPTION
!                  TRANSFORMS (X,IX) SO THAT
!                  RADIX**(-L) .LE. DABS(X) .LT. RADIX**L.
!                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
!                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
!                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
!
!     SUBROUTINE XDC210
! USAGE
!                  CALL XDC210(K,Z,J)
! DESCRIPTION
!                  GIVEN K THIS SUBROUTINE COMPUTES J AND Z
!                  SUCH THAT  RADIX**K = Z*10**J, WHERE Z IS IN
!                  THE RANGE 1/10 .LE. Z .LT. 1.
!                  THE VALUE OF Z WILL BE ACCURATE TO FULL
!                  DOUBLE-PRECISION PROVIDED THE NUMBER
!                  OF DECIMAL PLACES IN THE LARGEST
!                  INTEGER PLUS THE NUMBER OF DECIMAL
!                  PLACES CARRIED IN DOUBLE-PRECISION DOES NOT
!                  EXCEED 60. XDC210 IS CALLED BY SUBROUTINE
!                  XDCON WHEN NECESSARY. THE USER SHOULD
!                  NEVER NEED TO CALL XDC210 DIRECTLY.
!
!     SUBROUTINE XDCON
! USAGE
!                  CALL XDCON(X,IX)
! DESCRIPTION
!                  CONVERTS (X,IX) = X*RADIX**IX
!                  TO DECIMAL FORM IN PREPARATION FOR
!                  PRINTING, SO THAT (X,IX) = X*10**IX
!                  WHERE 1/10 .LE. DABS(X) .LT. 1
!                  IS RETURNED, EXCEPT THAT IF
!                  (DABS(X),IX) IS BETWEEN RADIX**(-2L)
!                  AND RADIX**(2L) THEN THE REDUCED
!                  FORM WITH IX = 0 IS RETURNED.
!
!     SUBROUTINE XDRED
! USAGE
!                  CALL XDRED(X,IX)
! DESCRIPTION
!                  IF
!                  RADIX**(-2L) .LE. (DABS(X),IX) .LE. RADIX**(2L)
!                  THEN XDRED TRANSFORMS (X,IX) SO THAT IX=0.
!                  IF (X,IX) IS OUTSIDE THE ABOVE RANGE,
!                  THEN XDRED TAKES NO ACTION.
!                  THIS SUBROUTINE IS USEFUL IF THE
!                  RESULTS OF EXTENDED-RANGE CALCULATIONS
!                  ARE TO BE USED IN SUBSEQUENT ORDINARY
!                  DOUBLE-PRECISION CALCULATIONS.
!
!***REFERENCES  (SEE DESCRIPTION ABOVE)
!***ROUTINES CALLED  I1MACH, XERROR
!***COMMON BLOCKS    XDBLK1, XDBLK2, XDBLK3
!***END PROLOGUE  XDSET
      INTEGER IRAD, NRADPL, NBITS
      DOUBLE PRECISION DZERO, DZEROX
      COMMON /XDBLK1/ NBITSF
      SAVE /XDBLK1/
      DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
      INTEGER L, L2, KMAX
      COMMON /XDBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
      SAVE /XDBLK2/
      INTEGER NLG102, MLG102, LG102
      COMMON /XDBLK3/ NLG102, MLG102, LG102(21)
      SAVE /XDBLK3/
      INTEGER IFLAG
      SAVE IFLAG
!
      DIMENSION LOG102(20), LGTEMP(20)
!
!   LOG102 CONTAINS THE FIRST 60 DIGITS OF LOG10(2) FOR USE IN
! CONVERSION OF EXTENDED-RANGE NUMBERS TO BASE 10 .
      DATA LOG102 /301,029,995,663,981,195,213,738,894,724,493,026,768, &
     & 189,881,462,108,541,310,428/
!
! FOLLOWING CODING PREVENTS XDSET FROM BEING EXECUTED MORE THAN ONCE.
! THIS IS IMPORTANT BECAUSE SOME SUBROUTINES (SUCH AS XDNRMP AND
! XDLEGF) CALL XDSET TO MAKE SURE EXTENDED-RANGE ARITHMETIC HAS
! BEEN INITIALIZED. THE USER MAY WANT TO PRE-EMPT THIS CALL, FOR
! EXAMPLE WHEN I1MACH IS NOT AVAILABLE. SEE CODING BELOW.
      DATA IFLAG /0/
!***FIRST EXECUTABLE STATEMENT  XDSET
      IF (IFLAG .NE. 0) RETURN
      IFLAG = 1
      IRADX = IRAD
      NRDPLC = NRADPL
      DZEROX = DZERO
      IMINEX = 0
      IMAXEX = 0
      NBITSX = NBITS
! FOLLOWING 6 STATEMENTS SHOULD BE DELETED IF I1MACH
! NOT AVAILABLE OR NOT CONFIGURED TO RETURN THE CORRECT
! MACHINE-DEPENDENT VALUES.
      IF (IRADX .EQ. 0) IRADX = I1MACH (10)
      IF (NRDPLC .EQ. 0) NRDPLC = I1MACH (14)
      IF (DZEROX .EQ. 0.0D0) IMINEX = I1MACH (15)
      IF (DZEROX .EQ. 0.0D0) IMAXEX = I1MACH (16)
      IF (NBITSX .EQ. 0) NBITSX = I1MACH (8)
      IF (IRADX.EQ.2) GO TO 10
      IF (IRADX.EQ.4) GO TO 10
      IF (IRADX.EQ.8) GO TO 10
      IF (IRADX.EQ.16) GO TO 10
      WRITE(*,*) 'ERR IN XDSET...IMPROPER VALUE OF IRAD'
      GO TO 100
   10 CONTINUE
      LOG2R=0
      IF (IRADX.EQ.2) LOG2R = 1
      IF (IRADX.EQ.4) LOG2R = 2
      IF (IRADX.EQ.8) LOG2R = 3
      IF (IRADX.EQ.16) LOG2R = 4
      NBITSF=LOG2R*NRDPLC
      RADIX = IRADX
      DLG10R = DLOG10(RADIX)
      IF (DZEROX .NE. 0.0D0) GO TO 14
      L = MIN0 ((1-IMINEX)/2, (IMAXEX-1)/2)
      GO TO 16
   14 L = 0.5D0*DLOG10(DZEROX)/DLG10R
! RADIX**(2*L) SHOULD NOT OVERFLOW, BUT REDUCE L BY 1 FOR FURTHER
! PROTECTION.
      L=L-1
   16 L2 = 2*L
      IF (L.GE.4) GO TO 20
      WRITE(*,*) 'ERR IN XDSET...IMPROPER VALUE OF DZERO'
      GO TO 100
   20 RADIXL = RADIX**L
      RAD2L = RADIXL**2
!    IT IS NECESSARY TO RESTRICT NBITS (OR NBITSX) TO BE LESS THAN SOME
! UPPER LIMIT BECAUSE OF BINARY-TO-DECIMAL CONVERSION. SUCH CONVERSION
! IS DONE BY XDC210 AND REQUIRES A CONSTANT THAT IS STORED TO SOME FIXED
! PRECISION. THE CONSTANT THAT IS STORED (LOG102 IN THIS ROUTINE) PROVIDES
! FOR CONVERSIONS TO BE ACCURATE TO THE LAST DECIMAL DIGIT WHEN THE INTEGER
! WORD LENGTH DOES NOT EXCEED 63. A LOWER LIMIT OF 15 BITS IS IMPOSED
! BECAUSE THE SOFTWARE IS DESIGNED TO RUN ON COMPUTERS WITH INTEGER WORD
! LENGTH OF AT LEAST 16 BITS.
      IF (15.LE.NBITSX .AND. NBITSX.LE.63) GO TO 30
      WRITE(*,*) 'ERR IN XDSET...IMPROPER VALUE OF NBITS'
      GO TO 100
   30 CONTINUE
      KMAX = 2**(NBITSX-1) - L2
      NB = (NBITSX-1)/2
      MLG102 = 2**NB
      IF (1.LE.NRDPLC*LOG2R .AND. NRDPLC*LOG2R.LE.120) GO TO 40
      WRITE(*,*) 'ERR IN XDSET...IMPROPER VALUE OF NRADPL'
      GO TO 100
   40 CONTINUE
      NLG102 = NRDPLC*LOG2R/NB + 3
      NP1 = NLG102 + 1
!
!   AFTER COMPLETION OF THE FOLLOWING LOOP, IC CONTAINS
! THE INTEGER PART AND LGTEMP CONTAINS THE FRACTIONAL PART
! OF LOG10(IRADX) IN RADIX 1000.
      IC = 0
      DO 50 II=1,20
        I = 21 - II
        IT = LOG2R*LOG102(I) + IC
        IC = IT/1000
        LGTEMP(I) = MOD(IT,1000)
   50 CONTINUE
!
!   AFTER COMPLETION OF THE FOLLOWING LOOP, LG102 CONTAINS
! LOG10(IRADX) IN RADIX MLG102. THE RADIX POINT IS
! BETWEEN LG102(1) AND LG102(2).
      LG102(1) = IC
      DO 80 I=2,NP1
        LG102(I) = 0
        DO 70 J=1,NB
          IC = 0
          DO 60 KK=1,20
            K = 21 - KK
            IT = 2*LGTEMP(K) + IC
            IC = IT/1000
            LGTEMP(K) = MOD(IT,1000)
   60     CONTINUE
          LG102(I) = 2*LG102(I) + IC
   70   CONTINUE
   80 CONTINUE
!
! CHECK SPECIAL CONDITIONS REQUIRED BY SUBROUTINES...
      IF (NRDPLC.LT.L) GO TO 90
      WRITE(*,*) 'ERR IN XDSET...NRADPL .GE. L'
      GO TO 100
   90 IF (6*L.LE.KMAX) GO TO 100
      WRITE(*,*) 'ERR IN XDSET...6*L .GT. KMAX'
      GO TO 100
  100 CONTINUE
      RETURN
      END

INTEGER FUNCTION I1MACH(I)
      INTEGER I
!
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
!    INTEGERS HAVE FORM SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!    I1MACH( 7) = A, THE BASE.
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!    FLOATS HAVE FORM  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!               WHERE  EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, THE BASE.
!  SINGLE-PRECISION
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!  DOUBLE-PRECISION
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
      INTEGER IMACH(16), OUTPUT, SC, SMALL(2)
      SAVE IMACH, SC
      REAL RMACH
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
      INTEGER I3, J, K, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
!  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
!  INCLUDING AUTO-DOUBLE COMPILERS.
!  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!  ON THE NEXT LINE
      DATA SC/0/
!  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!          mail netlib@research.bell-labs.com
!          send old1mach from blas
!  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /   43 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   63 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SC/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SC/987/
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
!     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
!     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    6 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /, SC/987/
!
      IF (SC .NE. 987) THEN
!        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
!           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532   &
     &           .AND. SMALL(2) .EQ. -448790528)  & 
     &       .OR.     (SMALL(2) .EQ. 1117925532   &
     &           .AND. SMALL(1) .EQ. -448790528)) THEN
!               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935 &
     &          .AND. SMALL(2) .EQ. 10752) THEN
!               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943 &
     &          .AND. SMALL(2) .EQ. 704643072) THEN
!               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
!               WRITE(*,9010)
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
!               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SC = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
!               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
!               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
!              *** CONVEX C-1 ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -128
               IMACH(13) = 127
               IMACH(14) = 53
               IMACH(15) = -1024
               IMACH(16) = 1023
            ELSE
               DO 10 I3 = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(I3)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
!              *** CRAY T3E ***
               IMACH( 1) = 5
               IMACH( 2) = 6
               IMACH( 3) = 0
               IMACH( 4) = 0
               IMACH( 5) = 64
               IMACH( 6) = 8
               IMACH( 7) = 2
               IMACH( 8) = 63
               CALL I1MCR1(IMACH(9), K, 32767, 16777215, 16777215)
               IMACH(10) = 2
               IMACH(11) = 53
               IMACH(12) = -1021
               IMACH(13) = 1024
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               GO TO 35
 20            CALL I1MCR1(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
!                  WRITE(*,9020)
                  STOP 777
                  END IF
!              *** CRAY 1, XMP, 2, AND 3 ***
               IMACH(1) = 5
               IMACH(2) = 6
               IMACH(3) = 102
               IMACH(4) = 6
               IMACH(5) = 46
               IMACH(6) = 8
               IMACH(7) = 2
               IMACH(8) = 45
               CALL I1MCR1(IMACH(9), K, 0, 4194303, 16777215)
               IMACH(10) = 2
               IMACH(11) = 47
               IMACH(12) = -8188
               IMACH(13) = 8189
               IMACH(14) = 94
               IMACH(15) = -8141
               IMACH(16) = 8189
               GO TO 35
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
 35      SC = 987
         END IF
! 9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/
!     * ' statements appropriate for your machine and setting'/
!     * ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
! 9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/
!     * ' appropriate for your machine.')
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 40
      I1MACH = IMACH(I)
      RETURN
 40   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
      STOP
      END