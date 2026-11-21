!*==ZXSSQ.spg  processed by SPAG 6.72Dc at 22:52 on 28 May 2017
!   IMSL ROUTINE NAME   - ZXSSQ
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   LATEST REVISION     - JUNE 1, 1982
!                         august 6th of 2012. v-1.0.2. Codename: good bye "calling characters as 'integers'...."
!
!   PURPOSE             - MINIMUM OF THE SUM OF SQUARES OF M FUNCTIONS
!                           IN N VARIABLES USING A FINITE DIFFERENCE
!                           LEVENBERG-MARQUARDT ALGORITHM
!
!   USAGE               - CALL ZXSSQ(FUNC,M,N,NSIG,EPS,DELTA,MAXFN,IOPT,
!                           PARM,X,SSQ,F,XJAC,IXJAC,XJTJ,WORK,INFER,IER)
!
!   ARGUMENTS    FUNC   - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
!                           THE RESIDUAL VECTOR F(1),F(2),...,F(M) FOR
!                           GIVEN PARAMETER VALUE3S X(1),X(2),...,X(N).
!                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
!                             CALL FUNC(X,M,N,F)
!                             WHERE X IS A VECTOR OF LENGTH N AND F IS
!                               A VECTOR OF LENGTH M.
!                             FUNC MUST APPEAR IN AN EXTERNAL STATEMENT
!                               IN THE CALLING PROGRAM.
!                             FUNC MUST NOT ALTER THE VALUES OF
!                               X(I),I=1,...,N, M, OR N.
!                M      - THE NUMBER OF RESIDUALS OR OBSERVATIONS
!                           (INPUT)
!                N      - THE NUMBER OF UNKNOWN PARAMETERS (INPUT).
!                NSIG   - FIRST CONVERGENCE CRITERION. (INPUT)
!                           CONVERGENCE CONDITION SATISFIED IF ON TWO
!                           SUCCESSIVE ITERATIONS, THE PARAMETER
!                           ESTIMATES AGREE, COMPONENT BY COMPONENT,
!                           TO NSIG DIGITS.
!                EPS    - SECOND CONVERGENCE CRITERION. (INPUT)
!                           CONVERGENCE CONDITION SATISFIED IF, ON TWO
!                           SUCCESSIVE ITERATIONS THE RESIDUAL SUM
!                           OF SQUARES ESTIMATES HAVE RELATIVE
!                           DIFFERENCE LESS THAN OR EQUAL TO EPS. EPS
!                           MAY BE SET TO ZERO.
!                DELTA  - THIRD CONVERGENCE CRITERION. (INPUT)
!                           CONVERGENCE CONDITION SATISFIED IF THE
!                           (EUCLIDEAN) NORM OF THE APPROXIMATE
!                           GRADIENT IS LESS THAN OR EQUAL TO DELTA.
!                           DELTA MAY BE SET TO ZERO.
!                             NOTE, THE ITERATION IS TERMINATED, AND
!                             CONVERGENCE IS CONSIDERED ACHIEVED, IF
!                             ANY ONE OF THE THREE CONDITIONS IS
!                             SATISFIED.
!                MAXFN  - INPUT MAXIMUM NUMBER OF FUNCTION EVALUATIONS
!                           (I.E., CALLS TO SUBROUTINE FUNC) ALLOWED.
!                           THE ACTUAL NUMBER OF CALLS TO FUNC MAY
!                           EXCEED MAXFN SLIGHTLY.
!                IOPT   - INPUT OPTIONS PARAMETER.
!                         IOPT=0 IMPLIES BROWN'S ALGORITHM WITHOUT
!                           STRICT DESCENT IS DESIRED.
!                         IOPT=1 IMPLIES STRICT DESCENT AND DEFAULT
!                           VALUES FOR INPUT VECTOR PARM ARE DESIRED.
!                         IOPT=2 IMPLIES STRICT DESCENT IS DESIRED WITH
!                           USER PARAMETER CHOICES IN INPUT VECTOR PARM.
!                PARM   - INPUT VECTOR OF LENGTH 4 USED ONLY FOR
!                         IOPT EQUAL TWO.  PARM(I) CONTAINS, WHEN
!                           I=1, THE INITIAL VALUE OF THE MARQUARDT
!                             PARAMETER USED TO SCALE THE DIAGONAL OF
!                             THE APPROXIMATE HESSIAN MATRIX, XJTJ,
!                             BY THE FACTOR (1.0 + PARM(1)).  A SMALL
!                             VALUE GIVES A NEWTON STEP, WHILE A LARGE
!                             VALUE GIVES A STEEPEST DESCENT STEP.
!                             THE DEFAULT VALUE FOR PARM(1) IS 0.01.
!                           I=2, THE SCALING FACTOR USED TO MODIFY THE
!                             MARQUARDT PARAMETER, WHICH IS DECREASED
!                             BY PARM(2) AFTER AN IMMEDIATELY SUCCESSFUL
!                             DESCENT DIRECTION, AND INCREASED BY THE
!                             SQUARE OF PARM(2) IF NOT.  PARM(2) MUST
!                             BE GREATER THAN ONE, AND TWO IS DEFAULT.
!                           I=3, AN UPPER BOUND FOR INCREASING THE
!                             MARQUARDT PARAMETER.  THE SEARCH FOR A
!                             DESCENT POINT IS ABANDONED IF PARM(3) IS
!                             EXCEEDED.  PARM(3) GREATER THAN 100.0 IS
!                             RECOMMENDED.  DEFAULT IS 120.0.
!                           I=4, VALUE FOR INDICATING WHEN CENTRAL
!                             RATHER THAN FORWARD DIFFERENCING IS TO BE
!                             USED FOR CALCULATING THE JACOBIAN.  THE
!                             SWITCH IS MADE WHEN THE NORM OF THE
!                             GRADIENT OF THE SUM OF SQUARES FUNCTION
!                             BECOMES SMALLER THAN PARM(4).  CENTRAL
!                             DIFFERENCING IS GOOD IN THE VICINITY
!                             OF THE SOLUTION, SO PARM(4) SHOULD BE
!                             SMALL.  THE DEFAULT VALUE IS 0.10.
!                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
!                           VALUES.
!                         ON INPUT, X SHOULD CONTAIN THE INITIAL
!                           ESTIMATE OF THE LOCATION OF THE MINIMUM.
!                         ON OUTPUT, X CONTAINS THE FINAL ESTIMATE
!                           OF THE LOCATION OF THE MINIMUM.
!                SSQ    - OUTPUT SCALAR WHICH IS SET TO THE RESIDUAL
!                           SUMS OF SQUARES, F(1)**2+...+F(M)**2, FOR
!                           THE FINAL PARAMETER ESTIMATES.
!                F      - OUTPUT VECTOR OF LENGTH M CONTAINING THE
!                           RESIDUALS FOR THE FINAL PARAMETER ESTIMATES.
!                XJAC   - OUTPUT M BY N MATRIX CONTAINING THE
!                           APPROXIMATE JACOBIAN AT THE OUTPUT VECTOR X.
!                IXJAC  - INPUT ROW DIMENSION OF MATRIX XJAC EXACTLY
!                           AS SPECIFIED IN THE DIMENSION STATEMENT
!                           IN THE CALLING PROGRAM.
!                XJTJ   - OUTPUT VECTOR OF LENGTH (N+1)*N/2 CONTAINING
!                           THE N BY N MATRIX (XJAC-TRANSPOSED) * (XJAC)
!                           IN SYMMETRIC STORAGE MODE.
!                WORK   - WORK VECTOR OF LENGTH 5*N + 2*M + (N+1)*N/2.
!                         ON OUTPUT, WORK(I) CONTAINS FOR
!                           I=1, THE NORM OF THE GRADIENT DESCRIBED
!                             UNDER INPUT PARAMETERS DELTA AND PARM(4).
!                           I=2, THE NUMBER OF FUNCTION EVALUATIONS
!                             REQUIRED DURING THE WORK(5) ITERATIONS.
!                           I=3, THE ESTIMATED NUMBER OF SIGNIFICANT
!                             DIGITS IN OUTPUT VECTOR X.
!                           I=4, THE FINAL VALUE OF THE MARQUARDT
!                             SCALING PARAMETER DESCRIBED UNDER PARM(1).
!                           I=5, THE NUMBER OF ITERATIONS (I.E., CHANGES
!                             TO THE X VECTOR) PERFORMED.
!                           SEE PROGRAMMING NOTES FOR DESCRIPTION OF
!                             THE LATTER ELEMENTS OF WORK.
!                INFER  - AN INTEGER THAT IS SET, ON OUTPUT, TO
!                           INDICATE WHICH CONVERGENCE CRITERION WAS
!                           SATISFIED.
!                         INFER = 0 INDICATES THAT CONVERGENCE FAILED.
!                           IER GIVES FURTHER EXPLANATION.
!                         INFER = 1 INDICATES THAT THE FIRST CRITERION
!                           WAS SATISFIED.
!                         INFER = 2 INDICATES THAT THE SECOND CRITERION
!                           WAS SATISFIED.
!                         INFER = 4 INDICATES THAT THE THIRD CRITERION
!                           WAS SATISFIED.
!                         IF MORE THAN ONE OF THE CONVERGENCE CRITERIA
!                           WERE SATISFIED ON THE FINAL ITERATION,
!                           INFER CONTAINS THE CORRESPONDING SUM.
!                           (E.G., INFER = 3 IMPLIES FIRST AND SECOND
!                           CRITERIA SATISFIED SIMULTANEOUSLY).
!                IER    - ERROR PARAMETER (OUTPUT)
!                         TERMINAL ERROR
!                           IER=129 IMPLIES A SINGULARITY WAS DETECTED
!                             IN THE JACOBIAN AND RECOVERY FAILED.
!                           IER=130 IMPLIES AT LEAST ONE OF M, N, IOPT,
!                             PARM(1), OR PARM(2) WAS SPECIFIED
!                             INCORRECTLY.
!                           IER=132 IMPLIES THAT AFTER A SUCCESSFUL
!                             RECOVERY FROM A SINGULAR JACOBIAN, THE
!                             VECTOR X HAS CYCLED BACK TO THE
!                             FIRST SINGULARITY.
!                           IER=133 IMPLIES THAT MAXFN WAS EXCEEDED.
!                         WARNING ERROR
!                           IER=38 IMPLIES THAT THE JACOBIAN IS ZERO.
!                             THE SOLUTION X IS A STATIONARY POINT.
!                           IER=39 IMPLIES THAT THE MARQUARDT
!                             PARAMETER EXCEEDED PARM(3).  THIS
!                             USUALLY MEANS THAT THE REQUESTED
!                             ACCURACY WAS NOT ACHIEVED.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - LEQT1P,LUDECP,LUELMP,UERSET,UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZXSSQ(FUNC,M,N,Nsig,Eps,Delta,Maxfn,Iopt,Parm,X,Ssq,F, &
                     & Xjac,Ixjac,Xjtj,Work,Infer,Ier)
      IMPLICIT NONE
!*--ZXSSQ181
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER M , N , Nsig , Maxfn , Iopt , Ixjac , Infer , Ier
      DOUBLE PRECISION Eps , Delta , Parm(4) , X(N) , Ssq , F(M) ,      &
                     & Xjac(*) , Xjtj(*) , Work(*)
!                                  XJAC USED INTERNALLY IN PACKED FORM
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER imjc , igrad1 , igradl , igradu , idelx1 , idelxl ,       &
            & idelxu , iscal1 , iscall , iscalu , ixnew1 , ixnewl ,     &
            & ixbad1 , ifpl1 , ifpl , ifpu , ifml1 , ifml , ieval ,     &
            & ibad , isw , iter , j , ijac , i , k , l , is , js , li , &
            & lj , icount , izero , level , levold
      DOUBLE PRECISION al , cons2 , dnorm , dsq , erl2 , erl2x , f0 ,   &
                     & f0sq , f0sqs4 , g , half , hh , one , onep10 ,   &
                     & onep5 , onesf0 , ax , prec , rel , rhh , sig ,   &
                     & sqdif , ssqold , sum , ten , tenth , xdif ,      &
                     & xhold , up , zero , xdabs , relcon , p01 , two , &
                     & huntw , delta2
      DATA sig/16.8D0/
      DATA ax/0.1D0/
      DATA p01 , tenth , half , zero , one , onep5 , two , ten , huntw ,&
         & onep10/0.01D0 , 0.1D0 , 0.5D0 , 0.0D0 , 1.D0 , 1.5D0 , 2.D0 ,&
         & 10.0D0 , 1.2D2 , 1.D10/
!                                  ERROR CHECKS
!                                  FIRST EXECUTABLE STATEMENT
      Ier = 0
      level = 0
      CALL UERSET(level,levold)
      IF ( M.LE.0 .OR. M.GT.Ixjac .OR. N.LE.0 .OR. Iopt.LT.0 .OR.       &
         & Iopt.GT.2 ) GOTO 1100
      imjc = Ixjac - M
      IF ( Iopt.EQ.2 ) THEN
         IF ( Parm(2).LE.one .OR. Parm(1).LE.zero ) GOTO 1100
      ENDIF
!                                  MACHINE DEPENDENT CONSTANTS
      prec = ten**(-sig-one)
      rel = ten**(-sig*half)
      relcon = ten**(-Nsig)
!                                  WORK VECTOR IS CONCATENATION OF
!                                  SCALED HESSIAN,GRADIENT,DELX,SCALE,
!                                  XNEW,XBAD,F(X+DEL),F(X-DEL)
      igrad1 = ((N+1)*N)/2
      igradl = igrad1 + 1
      igradu = igrad1 + N
      idelx1 = igradu
      idelxl = idelx1 + 1
      idelxu = idelx1 + N
      iscal1 = idelxu
      iscall = iscal1 + 1
      iscalu = iscal1 + N
      ixnew1 = iscalu
      ixnewl = ixnew1 + 1
      ixbad1 = ixnew1 + N
      ifpl1 = ixbad1 + N
      ifpl = ifpl1 + 1
      ifpu = ifpl1 + M
      ifml1 = ifpu
      ifml = ifml1 + 1
      imjc = Ixjac - M
!                                  INITIALIZE VARIABLES
      al = one
      cons2 = tenth
      IF ( Iopt.NE.0 ) THEN
         IF ( Iopt.EQ.1 ) THEN
            al = p01
            f0 = two
            up = huntw
         ELSE
            al = Parm(1)
            f0 = Parm(2)
            up = Parm(3)
            cons2 = Parm(4)*half
         ENDIF
         onesf0 = one/f0
         f0sq = f0*f0
         f0sqs4 = f0sq**4
      ENDIF
      ieval = 0
      delta2 = Delta*half
      erl2 = onep10
      ibad = -99
      isw = 1
      iter = -1
      Infer = 0
      Ier = 0
      DO j = idelxl , idelxu
         Work(j) = zero
      ENDDO
      GOTO 500
!                                  JACOBIAN BY INCREMENTING X
 100  ijac = 0
      k = -imjc
      DO j = 1 , N
         k = k + imjc
         xdabs = DABS(X(j))
         hh = rel*(DMAX1(xdabs,ax))
         xhold = X(j)
         X(j) = X(j) + hh
         CALL FUNC(X,M,N,Work(ifpl))
         ieval = ieval + 1
         X(j) = xhold
         IF ( isw.EQ.1 ) THEN
!                                  FORWARD DIFFERENCES
            rhh = one/hh
            DO i = 1 , M
               k = k + 1
               Xjac(k) = (Work(ifpl1+i)-F(i))*rhh
            ENDDO
         ELSE
!                                  CENTRAL DIFFERENCES
            X(j) = xhold - hh
            CALL FUNC(X,M,N,Work(ifml))
            ieval = ieval + 1
            X(j) = xhold
            rhh = half/hh
            DO i = ifpl , ifpu
               k = k + 1
               Xjac(k) = (Work(i)-Work(i+M))*rhh
            ENDDO
         ENDIF
      ENDDO
!                                  CALCULATE GRADIENT
 200  erl2x = erl2
      erl2 = zero
      k = -imjc
      DO j = igradl , igradu
         k = k + imjc
         sum = zero
         DO i = 1 , M
            k = k + 1
            sum = sum + Xjac(k)*F(i)
         ENDDO
         Work(j) = sum
         erl2 = erl2 + sum*sum
      ENDDO
      erl2 = DSQRT(erl2)
!                                  CONVERGENCE TEST FOR NORM OF GRADIENT
      IF ( ijac.LE.0 ) THEN
         IF ( erl2.LE.delta2 ) Infer = Infer + 4
         IF ( erl2.LE.cons2 ) isw = 2
      ENDIF
!                                  CALCULATE THE LOWER SUPER TRIANGE OF
!                                  JACOBIAN (TRANSPOSED) * JACOBIAN
      l = 0
      is = -Ixjac
      DO i = 1 , N
         is = is + Ixjac
         js = -Ixjac
         DO j = 1 , i
            js = js + Ixjac
            l = l + 1
            sum = zero
            DO k = 1 , M
               li = is + k
               lj = js + k
               sum = sum + Xjac(li)*Xjac(lj)
            ENDDO
            Xjtj(l) = sum
         ENDDO
      ENDDO
!                                  CONVERGENCE CHECKS
      IF ( Infer.GT.0 ) GOTO 1300
      IF ( ieval.GE.Maxfn ) THEN
!                                  TERMINAL ERROR
         Ier = Ier + 1
         GOTO 1000
      ELSE
!                                  COMPUTE SCALING VECTOR
         IF ( Iopt.EQ.0 ) THEN
!                                  COMPUTE SCALING VECTOR AND NORM
            dnorm = zero
            k = 0
            DO j = 1 , N
               k = k + j
               Work(iscal1+j) = DSQRT(Xjtj(k))
               dnorm = dnorm + Xjtj(k)*Xjtj(k)
            ENDDO
            dnorm = one/DSQRT(dnorm)
!                                  NORMALIZE SCALING VECTOR
            DO j = iscall , iscalu
               Work(j) = Work(j)*dnorm*erl2
            ENDDO
         ELSE
            k = 0
            DO j = 1 , N
               k = k + j
               Work(iscal1+j) = Xjtj(k)
            ENDDO
         ENDIF
!                                  ADD L-M FACTOR TO DIAGONAL
         icount = 0
      ENDIF
 300  k = 0
      DO i = 1 , N
         DO j = 1 , i
            k = k + 1
            Work(k) = Xjtj(k)
         ENDDO
         Work(k) = Work(k) + Work(iscal1+i)*al
         Work(idelx1+i) = Work(igrad1+i)
      ENDDO
!                                  CHOLESKY DECOMPOSITION
 400  CALL LEQT1P(Work,1,N,Work(idelxl),N,0,g,xhold,Ier)
      IF ( Ier.EQ.0 ) THEN
         IF ( ibad.NE.-99 ) ibad = 0
      ELSE
         Ier = 0
         IF ( ijac.GT.0 ) GOTO 100
         IF ( ibad.LE.0 ) THEN
!                                  SINGULAR DECOMPOSITION
            IF ( ibad.LT.0 ) GOTO 800
            IF ( ibad.EQ.0 ) THEN
!                                  CHECK TO SEE IF CURRENT
!                                  ITERATE HAS CYCLED BACK TO
!                                  THE LAST SINGULAR POINT
               DO j = 1 , N
                  xhold = Work(ixbad1+j)
                  IF ( DABS(X(j)-xhold).GT.relcon*DMAX1(ax,DABS(xhold)) &
                     & ) GOTO 800
               ENDDO
               GOTO 1000
            ELSE
               GOTO 900
            ENDIF
         ELSE
            IF ( ibad.LT.2 ) GOTO 600
            GOTO 1200
         ENDIF
      ENDIF
!                                  CALCULATE SUM OF SQUARES
 500  DO j = 1 , N
         Work(ixnew1+j) = X(j) - Work(idelx1+j)
      ENDDO
      CALL FUNC(Work(ixnewl),M,N,Work(ifpl))
      ieval = ieval + 1
      Ssq = zero
      DO i = ifpl , ifpu
         Ssq = Ssq + Work(i)*Work(i)
      ENDDO
      IF ( iter.GE.0 ) THEN
         IF ( Iopt.EQ.0 ) GOTO 700
!                                  CHECK DESCENT PROPERTY
         IF ( Ssq.LE.ssqold ) THEN
!                                  ADJUST MARQUARDT PARAMETER
            IF ( icount.EQ.0 ) al = al/f0
            IF ( erl2x.GT.zero ) THEN
               g = erl2/erl2x
               IF ( erl2.LT.erl2x ) al = al*DMAX1(onesf0,g)
               IF ( erl2.GT.erl2x ) al = al*DMIN1(f0,g)
            ENDIF
            al = DMAX1(al,prec)
            GOTO 700
         ENDIF
      ELSE
!                                  SSQ FOR INITIAL ESTIMATES OF X
         iter = 0
         ssqold = Ssq
         DO i = 1 , M
            F(i) = Work(ifpl1+i)
         ENDDO
         GOTO 100
      ENDIF
!                                  INCREASE PARAMETER AND TRY AGAIN
 600  icount = icount + 1
      al = al*f0sq
      IF ( ijac.NE.0 ) THEN
         IF ( icount.GE.4 .OR. al.GT.up ) THEN
            al = al/f0sqs4
            GOTO 100
         ENDIF
      ENDIF
      IF ( al.LE.up ) GOTO 300
      IF ( ibad.EQ.1 ) GOTO 1200
      Ier = 39
      GOTO 1300
!                                  ONE ITERATION CYCLE COMPLETED
 700  iter = iter + 1
      DO j = 1 , N
         X(j) = Work(ixnew1+j)
      ENDDO
      DO i = 1 , M
         Work(ifml1+i) = F(i)
         F(i) = Work(ifpl1+i)
      ENDDO
!                                  RELATIVE CONVERGENCE TEST FOR X
      IF ( al.LE.5.0D0 ) THEN
         DO j = 1 , N
            xdif = DABS(Work(idelx1+j))/DMAX1(DABS(X(j)),ax)
            IF ( xdif.GT.relcon ) GOTO 750
         ENDDO
         Infer = 1
!                                  RELATIVE CONVERGENCE TEST FOR SSQ
 750     sqdif = DABS(Ssq-ssqold)/DMAX1(ssqold,ax)
         IF ( sqdif.LE.Eps ) Infer = Infer + 2
      ENDIF
!                                  MAIN LOOP
      ssqold = Ssq
!                                  CALCULATE JACOBIAN
      IF ( Infer.GT.0 .OR. ijac.GE.N .OR. Iopt.EQ.0 .OR. icount.GT.0 )  &
         & GOTO 100
!                                  RANK ONE UPDATE TO JACOBIAN
      ijac = ijac + 1
      dsq = zero
      DO j = idelxl , idelxu
         dsq = dsq + Work(j)*Work(j)
      ENDDO
      IF ( dsq.LE.zero ) GOTO 100
      DO i = 1 , M
         g = F(i) - Work(ifml1+i)
         k = i
         DO j = idelxl , idelxu
            g = g + Xjac(k)*Work(j)
            k = k + Ixjac
         ENDDO
         g = g/dsq
         k = i
         DO j = idelxl , idelxu
            Xjac(k) = Xjac(k) - g*Work(j)
            k = k + Ixjac
         ENDDO
      ENDDO
      GOTO 200
!                                  UPDATE THE BAD X VALUES
 800  DO j = 1 , N
         Work(ixbad1+j) = X(j)
      ENDDO
      ibad = 1
!                                  INCREASE DIAGONAL OF HESSIAN
 900  IF ( Iopt.NE.0 ) THEN
!                                  REPLACE ZEROES ON HESSIAN DIAGONAL
         izero = 0
         DO j = iscall , iscalu
            IF ( Work(j).LE.zero ) THEN
               izero = izero + 1
               Work(j) = one
            ENDIF
         ENDDO
         IF ( izero.LT.N ) GOTO 300
         Ier = 38
         GOTO 1300
      ELSE
         k = 0
         DO i = 1 , N
            DO j = 1 , i
               k = k + 1
               Work(k) = Xjtj(k)
            ENDDO
            Work(k) = onep5*(Xjtj(k)+al*erl2*Work(iscal1+i)) + rel
         ENDDO
         ibad = 2
         GOTO 400
      ENDIF
 1000 Ier = Ier + 1
      Ier = Ier + 1
 1100 Ier = Ier + 1
 1200 Ier = Ier + 129
      IF ( Ier.EQ.130 ) GOTO 1400
!                                  OUTPUT ERL2,IEVAL,NSIG,AL, AND ITER
 1300 g = sig
      DO j = 1 , N
         xhold = DABS(Work(idelx1+j))
         IF ( xhold.GT.zero ) g = DMIN1(g,-DLOG10(xhold)+DLOG10(DMAX1(ax&
                                & ,DABS(X(j)))))
      ENDDO
      IF ( N.LE.2 ) THEN
         DO j = 1 , N
            Work(j+5) = Work(j+igrad1)
         ENDDO
      ENDIF
      Work(1) = erl2 + erl2
      Work(2) = ieval
      Ssq = ssqold
      Work(3) = g
      Work(4) = al
      Work(5) = iter
 1400 CALL UERSET(levold,levold)
      IF ( Ier.NE.0 ) CALL UERTS(Ier,"ZXSSQ")
      END
!*==LEQT1P.spg  processed by SPAG 6.72Dc at 22:52 on 28 May 2017
!********************************************************************************
!   IMSL ROUTINE NAME   - LEQT1P
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - LINEAR EQUATION SOLUTION - POSITIVE DEFINITE
!                           MATRIX - SYMMETRIC STORAGE MODE - SPACE
!                           ECONOMIZER SOLUTION
!
!   USAGE               - CALL LEQT1P (A,M,N,B,IB,IDGT,D1,D2,IER)
!
!   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING THE
!                           N BY N COEFFICIENT MATRIX OF THE EQUATION
!                           AX = B. A IS A POSITIVE DEFINITE SYMMETRIC
!                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
!                         ON OUTPUT, A IS REPLACED BY THE LOWER
!                           TRIANGULAR MATRIX L WHERE A = L*L-TRANSPOSE.
!                           L IS STORED IN SYMMETRIC STORAGE MODE WITH
!                           THE DIAGONAL ELEMENTS OF L IN RECIPROCAL
!                           FORM.
!                M      - NUMBER OF RIGHT HAND SIDES (COLUMNS IN B).
!                           (INPUT)
!                N      - ORDER OF A AND NUMBER OF ROWS IN B. (INPUT)
!                B      - INPUT MATRIX OF DIMENSION N BY M CONTAINING
!                           THE RIGHT-HAND SIDES OF THE EQUATION
!                           AX = B.
!                         ON OUTPUT, THE N BY M SOLUTION MATRIX X
!                           REPLACES B.
!                IB     - ROW DIMENSION OF B EXACTLY AS SPECIFIED IN
!                           THE DIMENSION STATEMENT IN THE CALLING
!                           PROGRAM. (INPUT)
!                IDGT   - THE ELEMENTS OF A ARE ASSUMED TO BE CORRECT
!                           TO IDGT DECIMAL DIGITS. (CURRENTLY NOT USED)
!                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
!                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
!                IER    - ERROR PARAMETER. (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129 INDICATES THAT THE INPUT MATRIX
!                             A IS ALGORITHMICALLY NOT POSITIVE
!                             DEFINITE. (SEE THE CHAPTER L PRELUDE).
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - LUDECP,LUELMP,UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE LEQT1P(A,M,N,B,Ib,Idgt,D1,D2,Ier)
      IMPLICIT NONE
!*--LEQT1P624
!*** Start of declarations inserted by SPAG
      INTEGER i , Ib , Idgt , Ier , M , N
!*** End of declarations inserted by SPAG
!
      DIMENSION A(*) , B(Ib,*)
      DOUBLE PRECISION A , B , D1 , D2
!                                  FIRST EXECUTABLE STATEMENT
!                                  INITIALIZE IER
      Ier = 0
!                                  DECOMPOSE A
      CALL LUDECP(A,A,N,D1,D2,Ier)
      IF ( Ier.NE.0 ) THEN
!      CALL UERTS(IER,6HLEQT1P)
         CALL UERTS(Ier,"LEQT1P")
      ELSE
!                                  PERFORM ELIMINATION
         DO i = 1 , M
            CALL LUELMP(A,B(1,i),N,B(1,i))
         ENDDO
      ENDIF
      END
!*==LUDECP.spg  processed by SPAG 6.72Dc at 22:52 on 28 May 2017
!*****************************************************************
!   IMSL ROUTINE NAME   - LUDECP
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - DECOMPOSITION OF A POSITIVE DEFINITE MATRIX -
!                           SYMMETRIC STORAGE MODE
!
!   USAGE               - CALL LUDECP (A,UL,N,D1,D2,IER)
!
!   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
!                           THE N BY N POSITIVE DEFINITE SYMMETRIC
!                           MATRIX STORED IN SYMMETRIC STORAGE MODE.
!                UL     - OUTPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
!                           THE DECOMPOSED MATRIX L SUCH THAT A = L*
!                           L-TRANSPOSE. L IS STORED IN SYMMETRIC
!                           STORAGE MODE. THE DIAGONAL OF L CONTAINS THE
!                           RECIPROCALS OF THE ACTUAL DIAGONAL ELEMENTS.
!                N      - ORDER OF A. (INPUT)
!                D1,D2  - COMPONENTS OF THE DETERMINANT OF A.
!                           DETERMINANT(A) = D1*2.**D2. (OUTPUT)
!                IER    - ERROR PARAMETER. (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129 INDICATES THAT MATRIX A IS
!                             ALGORITHMICALLY NOT POSITIVE DEFINITE.
!                             (SEE THE CHAPTER L PRELUDE).
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE LUDECP(A,Ul,N,D1,D2,Ier)
      IMPLICIT NONE
!*--LUDECP697
!*** Start of declarations inserted by SPAG
      INTEGER i , Ier , ip , ip1 , iq , ir , j , k , N
!*** End of declarations inserted by SPAG
!
      DIMENSION A(*) , Ul(*)
      DOUBLE PRECISION A , Ul , D1 , D2 , zero , one , four , sixtn ,   &
                     & sixth , x , rn
      DATA zero , one , four , sixtn , sixth/0.0D0 , 1.D0 , 4.D0 ,      &
         & 16.D0 , .0625D0/
!                                  FIRST EXECUTABLE STATEMENT
      D1 = one
      D2 = zero
      rn = one/(N*sixtn)
      ip = 1
      Ier = 0
      DO i = 1 , N
         iq = ip
         ir = 1
         DO j = 1 , i
            x = A(ip)
            IF ( j.NE.1 ) THEN
               DO k = iq , ip1
                  x = x - Ul(k)*Ul(ir)
                  ir = ir + 1
               ENDDO
            ENDIF
            IF ( i.NE.j ) THEN
               Ul(ip) = x*Ul(ir)
            ELSE
               D1 = D1*x
               IF ( A(ip)+x*rn.LE.A(ip) ) GOTO 100
 10            IF ( DABS(D1).LE.one ) THEN
 15               IF ( DABS(D1).GE.sixth ) THEN
                     Ul(ip) = one/DSQRT(x)
                  ELSE
                     D1 = D1*sixtn
                     D2 = D2 - four
                     GOTO 15
                  ENDIF
               ELSE
                  D1 = D1*sixth
                  D2 = D2 + four
                  GOTO 10
               ENDIF
            ENDIF
            ip1 = ip
            ip = ip + 1
            ir = ir + 1
         ENDDO
      ENDDO
      GOTO 99999
 100  Ier = 129
!      CALL UERTS(IER,6HLUDECP)
      CALL UERTS(Ier,"LUDECP")
99999 END
!*==LUELMP.spg  processed by SPAG 6.72Dc at 22:52 on 28 May 2017
!************************************************************************
!   IMSL ROUTINE NAME   - LUELMP
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - VAX/DOUBLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - ELIMINATION PART OF THE SOLUTION OF AX=B -
!                           POSITIVE DEFINITE MATRIX - SYMMETRIC
!                           STORAGE MODE
!
!   USAGE               - CALL LUELMP (A,B,N,X)
!
!   ARGUMENTS    A      - INPUT VECTOR OF LENGTH N(N+1)/2 CONTAINING
!                           THE N BY N MATRIX L WHERE A = L*L-TRANSPOSE.
!                           L IS A LOWER TRIANGULAR MATRIX STORED IN
!                           SYMMETRIC STORAGE MODE. THE MAIN DIAGONAL
!                           ELEMENTS OF L ARE STORED IN RECIPROCAL
!                           FORM. MATRIX L MAY BE OBTAINED FROM IMSL
!                           ROUTINE LUDECP.
!                B      - VECTOR OF LENGTH N CONTAINING THE RIGHT HAND
!                           SIDE OF THE EQUATION AX = B. (INPUT)
!                N      - ORDER OF A AND THE LENGTH OF B AND X. (INPUT)
!                X      - VECTOR OF LENGTH N CONTAINING THE SOLUTION TO
!                           THE EQUATION AX = B. (OUTPUT)
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE LUELMP(A,B,N,X)
      IMPLICIT NONE
!*--LUELMP801
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , im1 , ip , iq , is , iw , k , kk , N , n1
!*** End of declarations inserted by SPAG
!
      DIMENSION A(*) , B(*) , X(*)
      DOUBLE PRECISION A , B , X , t , zero
      DATA zero/0.0D0/
!                                  FIRST EXECUTABLE STATEMENT
!                                  SOLUTION OF LY = B
      ip = 1
      iw = 0
      DO i = 1 , N
         t = B(i)
         im1 = i - 1
         IF ( iw.EQ.0 ) THEN
            IF ( t.NE.zero ) iw = i
            ip = ip + im1
         ELSE
            ip = ip + iw - 1
            DO k = iw , im1
               t = t - A(ip)*X(k)
               ip = ip + 1
            ENDDO
         ENDIF
         X(i) = t*A(ip)
         ip = ip + 1
      ENDDO
!                                  SOLUTION OF UX = Y
      n1 = N + 1
      DO i = 1 , N
         ii = n1 - i
         ip = ip - 1
         is = ip
         iq = ii + 1
         t = X(ii)
         IF ( N.GE.iq ) THEN
            kk = N
            DO k = iq , N
               t = t - A(is)*X(kk)
               kk = kk - 1
               is = is - kk
            ENDDO
         ENDIF
         X(ii) = t*A(is)
      ENDDO
      END
!*==UERSET.spg  processed by SPAG 6.72Dc at 22:52 on 28 May 2017
!*********************************************************************
!   IMSL ROUTINE NAME   - UERSET
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - VAX/SINGLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!                         august 6th of 2012...
!
!   PURPOSE             - SET MESSAGE LEVEL FOR IMSL ROUTINE UERTST
!
!   USAGE               - CALL UERSET (LEVEL,LEVOLD)
!
!   ARGUMENTS    LEVEL  - NEW VALUE FOR MESSAGE LEVEL. (INPUT)
!                           OUTPUT FROM IMSL ROUTINE UERTST IS
!                           CONTROLLED SELECTIVELY AS FOLLOWS,
!                             LEVEL = 4 CAUSES ALL MESSAGES TO BE
!                                       PRINTED,
!                             LEVEL = 3 MESSAGES ARE PRINTED IF IER IS
!                                       GREATER THAN 32,
!                             LEVEL = 2 MESSAGES ARE PRINTED IF IER IS
!                                       GREATER THAN 64,
!                             LEVEL = 1 MESSAGES ARE PRINTED IF IER IS
!                                       GREATER THAN 128,
!                             LEVEL = 0 ALL MESSAGE PRINTING IS
!                                       SUPPRESSED.
!                LEVOLD - PREVIOUS MESSAGE LEVEL. (OUTPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - UERTST,UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE UERSET(Level,Levold)
      IMPLICIT NONE
!*--UERSET897
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER Level , Levold
!                                  FIRST EXECUTABLE STATEMENT
      Levold = Level
!      CALL UERTS (LEVOLD,6HUERSET)
      CALL UERTS(Levold,"UERSET")
      END
!-----------------------------------------------------------------------
!   UERTS: error messages of the IMSL routines (replaces UERTST).
!
!   CALL UERTS(IER, NAME)
!     IER    - error code of the calling routine. IER = I + J, with
!                I = 128: terminal error
!                I =  64: warning with fix
!                I =  32: warning
!              and J the routine's own code.
!     NAME   - name of the calling routine ('ZXSSQ', 'LUDECP', ...).
!
!   Special case: UERSET calls with NAME = 'UERSET' to change the message
!   level. Then IER brings the new level and returns the previous one.
!     level 4: all errors are reported
!     level 3: only IER > 32;  level 2: only IER > 64;  level 1: only IER > 128
!     level 0: nothing is reported
!
!   Messages are appended to PARA.HIS. With IER = 0 nothing is written
!   (a false "UNDEFINED ERROR (IER = 0)" used to be written).
!-----------------------------------------------------------------------
      SUBROUTINE UERTS(Ier,Name)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: Ier
      CHARACTER(*), INTENT(IN) :: Name
      INTEGER, SAVE :: level = 4
      INTEGER :: previous_level, unit
      CHARACTER(30) :: kind
      LOGICAL :: report

      IF ( Name=='UERSET' ) THEN
         previous_level = level
         level = Ier
         IF ( level<0 .OR. level>4 ) level = 4
         Ier = previous_level
         RETURN
      ENDIF

      IF ( Ier==0 ) RETURN

      IF ( Ier>128 .AND. Ier<=999 ) THEN
         kind = 'TERMINAL ERROR'
         report = level>=1
      ELSEIF ( Ier>64 .AND. Ier<=128 ) THEN
         kind = 'WARNING WITH FIX ERROR'
         report = level>=2
      ELSEIF ( Ier>32 .AND. Ier<=64 ) THEN
         kind = 'WARNING ERROR'
         report = level>=3
      ELSE
         kind = 'UNDEFINED ERROR'
         report = level>=4
      ENDIF
      IF ( .NOT.report ) RETURN

!     NEWUNIT picks a free unit number, so no unit the program has open is
!     overwritten (unit 2 was used before, the same as tables.csv).
      OPEN (NEWUNIT=unit,FILE='PARA.HIS',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE (unit,'(" *** ",A," (IER = ",I0,") FROM IMSL ROUTINE ",A)')  &
            TRIM(kind), Ier, TRIM(Name)
      CLOSE (unit)
      END SUBROUTINE UERTS
