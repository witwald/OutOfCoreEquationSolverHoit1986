C=======================================================================
C
C HOIT, M. (1986). "AN UPDATED OUT-OF-CORE EQUATION SOLVER",
C MICROCOMPUTERS IN CIVIL ENGINEERING, VOLUME 1, ISSUE 3,
C PAGES 247-254.
C
C APPENDIX: LISTING OF ADD PROGRAM FROM SIMPAL (SUBROUTINES).
C
C TRANSCRIBED FROM THE PRINTED LISTING (PP. 250-253):
C
C    PROFIL  - FORM PROFILE POINTER
C    MAKNA   - FORM PROFILE CONNECTIVITY
C    MPOINT  - CONVERT NA TO POINTER ARRAY
C    BLKSET  - SET INFORMATION ABOUT THE BLOCKS (KA, KSIZB)
C    BLKPNT  - SET THE ELEMENT TO BLOCK POINTER (NELP, NBKP)
C    ASEMBL  - ASSEMBLE STIFFNESS MATRIX DIRECTLY INTO BLOCKS
C    ASEMLD  - ASSEMBLE THE LOAD VECTORS (BLOCKED) AND MASS
C    MLOAD   - ADD GLOBAL (NODAL) LOADS TO THE LOAD BLOCKS
C
C THE MAIN PROGRAM ADD IS REPLACED BY THE TEST DRIVER, WHICH
C REPRODUCES ITS MEMORY-SIZING LOGIC (SEE ADDDRV IN hoit_drv.f).
C
C THE CODE IS AS PRINTED EXCEPT WHERE MARKED "*** CHANGED".
C SIMPAL UTILITIES CALLED HERE (SOPEN, FCLOSE, ELREAD, RZERO,
C SAVBLK, GETBLK, FREE, FREEI, FREER, FREEPT) WERE NOT PUBLISHED;
C REPLACEMENTS ARE IN simpal_io.f.
C
C=======================================================================
C
C---------------------------------------------------------------- PROFIL
      SUBROUTINE PROFIL(NA,LM,MD,NQ)
C-----------------------------------------------------------------------
C      ROUTINE TO FORM PROFILE POINTER AND MASS
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 NAM,ELNAM
      DIMENSION NA(NQ),LM(MD),DUM(1)
      COMMON /CLT/ ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,
     * NUMEL,NEL(20),LENEL(20),NKBLK,NLBLK,KSIZ,KMEM,
     * KYDYN,KYSUB,KYLIN
      COMMON /PARC/ NAM(12),ELNAM(4,20)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C
      WRITE(NTM,2000)
C-----INITIALIZE POINTER ARRAY -----
      DO 10 I=1,NEQ
   10 NA(I) = I
C-----LOOP ON ALL ELEMENT TYPES -----
        DO 400 I=1,NUMEL
        NUML = NEL(I)
        CALL SOPEN(NT1,ELNAM(1,I),-LENEL(I))
C-----LOOP ON ALL ELEMENTS OF GROUP -----
          DO 300 J=1,NUML
          CALL ELREAD(J,LM,DUM,DUM,DUM,ND,LOAD,MAXD,1)
          CALL MAKNA(LM,NA,ND,NEQ,MD)
  300     CONTINUE
        CALL FCLOSE(NT1)
  400   CONTINUE
C
C-----CONVERT NA FROM CONNECTIVITY TO POINTER -----
  800 CALL MPOINT(NA,NEQ,KSIZ)
C
  900 RETURN
 2000 FORMAT(' SETTING UP STIFFNESS PROFILE')
      END
C
C----------------------------------------------------------------- MAKNA
      SUBROUTINE MAKNA(LM,NA,ND,NEQ,MD)
C-----------------------------------------------------------------------
C     ROUTINE TO FORM PROFILE CONNECTIVITY
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LM(MD),NA(NEQ)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C
C-----FORM CONNECTIVITY -----
   50 DO 200 I=1,ND
        LCR = LM(I)
        IF(LCR.GT.NEQ) GOTO 200
        IF(LCR.LE.0) GOTO 200
        NAC = NA(LCR)
        IF(NAC.EQ.1) GOTO 200
C-----FIND MINIMUM DOF CONNECTED TO CURRENT -----
        DO 100 J=1,ND
          LCK = LM(J)
          IF(LCK.GT.NEQ) GOTO 100
          IF(LCK.LE.0) GOTO 100
          IF(LCK.LT.NAC) NAC = LCK
  100   CONTINUE
        NA(LCR) = NAC
  200 CONTINUE
C
  900 RETURN
      END
C
C---------------------------------------------------------------- MPOINT
      SUBROUTINE MPOINT(NA,NEQ,KSIZ)
C-----------------------------------------------------------------------
C     ROUTINE TO CONVERT NA TO POINTER ARRAY
C-----------------------------------------------------------------------
      DIMENSION NA(NEQ)
C
      KSIZ = 1
      DO 100 I=2,NEQ
      KSIZ = KSIZ + I - NA(I) + 1
      NA(I) = KSIZ
  100 CONTINUE
C
      RETURN
      END
C
C---------------------------------------------------------------- BLKSET
      SUBROUTINE BLKSET(KA,KSIZB,NA,KSZKM,KSIZ,NKBLK,NEQ,MAXNA)
C-----------------------------------------------------------------------
C     ROUTINE TO SET INFORMATION ABOUT THE BLOCKS
C-----------------------------------------------------------------------
C     KA(4,NKBLK) = BLOCK POINTER ARRAY
C          1      = BLOCK NUMBER
C          2      = LOWEST CONNECTED BLOCK
C          3      = FIRST EQUATION NUMBER IN BLOCK
C          4      = LAST EQUATION NUMBER IN BLOCK
C    KSIZB(NKBLK) = ARRAY OF STIFFNESS BLOCK SIZES
C      NKBLK      = NUMBER OF STIFFNESS BLOCKS
C      KSZKM      = LARGEST SIZE OF A STIFFNESS BLOCK
C      MAXNA      = LARGEST NUMBER OF DIAGONAL TERMS IN A BLOCK
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NA(NEQ),KA(4,1),KSIZB(1)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C
      WRITE(NTM,1000)
C
C----- INITIALIZATION
      KSIZT=0
      KSZKM=0
      MAXNA=0
      NNA=0
      NKBLK=1
      NI=1
      KDIF=NEQ
      NDIF=1
   50 CONTINUE
      DO 200 I=NI+1,NEQ
        NCOLT=NA(I)-NA(I-1)
C----- SET LOWEST ROW NUMBER IN BLOCK (FOR BLOCK CONNECTION) -----
        KD=I-NCOLT+1
        IF(KD.LT.KDIF) KDIF=KD
C----- FIND NUMBER OF TERMS IF COLUMN ADDED TO BLOCK -----
        NDIFO=NDIF
        NDIF=NDIF+NCOLT
        IF(NDIF.GE.KSIZ.OR.I.EQ.NEQ) THEN
C----- SET SIZE OF CURRENT STIFFNESS BLOCK
          IF(NDIF.GT.KSIZ) THEN
            IT=I-1
            KSIZT=NDIFO
            ELSE
C
            IT=I
            KSIZT=NDIF
          ENDIF
          KSIZB(NKBLK)=KSIZT
C
C----- SET MAXIMUM SIZE OF A STIFFNESS BLOCK
          IF(KSIZT.GT.KSZKM) KSZKM=KSIZT
C----- CALCULATE NUMBER OF EQUATIONS IN CURRENT BLOCK
          NNA=IT-NI
          IF(NKBLK.EQ.1) NNA=NNA+1
C----- SET MAXIMUM NUMBER OF EQUATIONS IN A BLOCK
          IF(NNA.GT.MAXNA) MAXNA=NNA
C----- SET KA ARRAY FOR CURRENT BLOCK
          KA(1,NKBLK)=NKBLK
          NIS=NI+1
          IF(NKBLK.EQ.1) NIS=NI
          KA(3,NKBLK)=NIS
          KA(4,NKBLK)=IT
          DO 100 K=1,NKBLK
            IF(KDIF.LE.KA(4,K)) THEN
              KA(2,NKBLK)=K
              IF(IT.EQ.NEQ) GOTO 200
C----- RESET COUNTERS FOR ANOTHER BLOCK -----
C
                NDIF=0
                NKBLK=NKBLK+1
                KDIF=NEQ
                NI=IT
                GOTO 50
C
            ENDIF
C
  100     CONTINUE
C
        ENDIF
C
  200 CONTINUE
C
      RETURN
 1000 FORMAT(' SETTING STIFFNESS BLOCKING INFORMATION')
      END
C
C---------------------------------------------------------------- BLKPNT
      SUBROUTINE BLKPNT(LM,KA,NELP,NBKP,NBKT,NBP,NTNE)
C-----------------------------------------------------------------------
C     ROUTINE TO SET THE ELEMENT TO BLOCK POINTER
C-----------------------------------------------------------------------
C     NBKP() = ARRAY CONTAINING ALL BLOCK CONNECTED TO EACH ELEMENT
C    NELP(I) = LAST POSITION IN NBKP RELEVANT TO ELEMENT I.
C     NBKT   = DIMENSION OF NBKP ARRAY
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 NAM,ELNAM
      DIMENSION NBKP(1),NELP(1),KA(4,1),LM(1),DUM(1)
      COMMON /PARC/ NAM(12),ELNAM(4,20)
      COMMON /CLT/ ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,
     * NUMEL,NEL(20),LENEL(20),NKBLK,NLBLK,KSIZ,KMEM,
     * KYDYN,KYSUB,KYLIN
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C
      WRITE(NTM,1000)
C
C----- INITIALIZATION
      NBKT=0
      NCOUNT=0
C
C-----LOOP ON NUMBER OF ELEMENT GROUPS
      IF(NBP.EQ.1) THEN
C----- DEFAULT POINTER ARRAYS
        DO 50 I=1,NTNE
          NELP(I)=I
   50     NBKP(I)=1
          NBKT=NTNE
        ELSE
C
      DO 400 I=1,NUMEL
        NUML=NEL(I)
       CALL SOPEN(NT1,ELNAM(1,I),-LENEL(I))
C----- LOOP ON ELEMENTS
        DO 300 J=1,NUML
          NCOUNT=NCOUNT+1
         CALL ELREAD(J,LM,DUM,DUM,DUM,ND,LOAD,MAXD,1)
C----- LOOP ON NUMBER OF BLOCKS
          DO 200 K=1,NKBLK
            KL=KA(3,K)
            KG=KA(4,K)
C----- LOOP ON NUMBER OF DOF'S
            DO 100 L=1,ND
              LL=LM(L)
              IF(LL.GE.KL.AND.LL.LE.KG) THEN
C----- CURRENT ELEMENT STIFFNESS TERM IN CURRENT BLOCK
                NBKT=NBKT+1
                NBKP(NBKT)=K
                GOTO 200
              ENDIF
C
  100       CONTINUE
C
  200     CONTINUE
C
C-----SET POINTER FOR CURRENT ELEMENT (LAST TERM) -----
          NELP(NCOUNT)=NBKT
  300   CONTINUE
C
        CALL FCLOSE(NT1)
  400 CONTINUE
      ENDIF
C
      RETURN
 1000 FORMAT(' SETTING BLOCK TO ELEMENT POINTER INFORMATION')
      END
C
C---------------------------------------------------------------- ASEMBL
      SUBROUTINE ASEMBL(STIF,NA,ESTF,R,EMAS,KSIZB,LM,KA,
     * NELP,NBKP,MA,MD,NQ,ND,LD,KSZ)
C-----------------------------------------------------------------------
C      ROUTINE TO ASSEMBLE STIFFNESS MATRIX
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 NAM,ELNAM
      DIMENSION STIF(KSZ),NA(NQ),ESTF(MD*(MD+1)/2),LM(MD),NBKP(1),
     * R(1),EMAS(1),NELP(1),KA(4,1),KSIZB(1),MA(1)
      COMMON /PARC/ NAM(12),ELNAM(4,20)
      COMMON /CLT/ ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,
     * NUMEL,NEL(20),LENEL(20),NKBLK,NLBLK,KSIZ,KMEM,
     * KYDYN,KYSUB,KYLIN
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C
      WRITE(NTM,2000) NKBLK
C
C-----INITIALIZATION -----
      NOFSET=0
      DO 500 IB=1,NKBLK
        NDIP=1
        KL=KA(3,IB)
        KG=KA(4,IB)
        LBK=KA(2,IB)
        NCOUNT=0
        WRITE(NTM,2100) IB,KL,KG,LBK
C-----ZERO STIFFNESS -----
        CALL RZERO(STIF,KSIZ)
C-----LOOP ON ALL ELEMENT GROUPS -----
        DO 400 I=1,NUMEL
          NUML = NEL(I)
          IF(NUML.LE.0) GOTO 400
          CALL SOPEN(NT1,ELNAM(1,I),-LENEL(I))
C-----LOOP ON ALL ELEMENTS OF GROUP -----
          DO 300 J=1,NUML
            NCOUNT=NCOUNT+1
            NDIN=NELP(NCOUNT)
            DO 50 IN=NDIP,NDIN
              IF(NBKP(IN).EQ.IB) GOTO 75
   50       CONTINUE
            GOTO 250
   75       CALL ELREAD(J,LM,ESTF,R,EMAS,ND,LOAD,MAXD,0)
            NDO=(ND+1)*ND/2
            IEL = 1
            DO 200 K=1,ND
              JJ = LM(K)
              DO 100 L=1,K
                IF(JJ.GT.NEQ.OR.JJ.LE.0) GOTO 100
                II = LM(L)
                IF(II.GT.NEQ.OR.II.LE.0) GOTO 100
                IF(JJ.GE.II) THEN
C-----TERM IS UPPER HALF (II,JJ) -----
                  IF(JJ.GT.KG.OR.JJ.LT.KL) GOTO 100
                  M = NA(JJ) + II - JJ - NOFSET
C
                  ELSE
C-----TERM IS LOWER HALF - SWITCH TO UPPER HALF -----
                  IF(II.GT.KG.OR.II.LT.KL) GOTO 100
                  M = NA(II) + JJ - II - NOFSET
C
                ENDIF
C
                STIF(M) = STIF(M) + ESTF(IEL)
C-----ACCOUNT FOR CONSTRAINTS (ELEM HAS 2 DOF SAME) -----
                IF(II.EQ.JJ.AND.K.NE.L)
     *               STIF(M) = STIF(M) + ESTF(IEL)
  100         IEL = IEL + 1
  200       CONTINUE
  250       NDIP=NDIN+1
  300     CONTINUE
          CALL FCLOSE(NT1)
  400   CONTINUE
C
C-----SET UP BLOCK DIAGONAL POINTER -----
        NTS=KG-KL+1
        NAL=KL
        DO 450 N=1,NTS
        MA(N)=NA(NAL)-NOFSET
  450   NAL=NAL+1
C
C-----SAVE STIFFNESS -----
        CALL SAVBLK(STIF,KA(1,IB),MA,KSIZB(IB),IB,0)
C
  500 NOFSET=NA(KG)
C
  900 RETURN
 2000 FORMAT(/' *** ASSEMBLING GLOBAL STIFFNESS ***'/
     * '   THERE ARE A TOTAL OF ',I3,' BLOCKS'/)
 2100 FORMAT(/' FORMING STIFFNESS BLOCK NUMBER    ',I5/
     *        ' FIRST EQUATION NUMBER IN BLOCK IS ',I5/
     *        ' LAST  EQUATION NUMBER IN BLOCK IS ',I5/
     *        ' LOWEST CONNECTED BLOCK IS         ',I5)
      END
C
C---------------------------------------------------------------- ASEMLD
      SUBROUTINE ASEMLD(STIF,ESTF,KA,EMAS,TMAS,LM,R,ND,MD,
     *          KSZ)
C-----------------------------------------------------------------------
C     ROUTINE TO ASSEMBLE THE LOAD VECTORS, AND MASS
C-----------------------------------------------------------------------
C     NLBLK = NUMBER OF LOAD BLOCKS
C      NLDB = NUMBER OF LOAD CASES IN A BLOCK
C      NLEN = SIZE OF A BLOCK
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*1 NAM,ELNAM
      DIMENSION  STIF(KSZ),LM(MD),EMAS(MD),TMAS(1),
     * R(MD,1),KA(4,1),ESTF(1)
      COMMON /PARC/ NAM(12),ELNAM(4,20)
      COMMON /CLT/ ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,
     * NUMEL,NEL(20),LENEL(20),NKBLK,NLBLK,KSIZ,KMEM,
     * KYDYN,KYSUB,KYLIN
C
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C
C-----CALCULATE THE NUMBER OF BLOCKS NECESSARY FOR THE LOAD -----
      NLBLK=1+(LOAD*NEQ-1)/KSIZ
C-----CALCULATE THE NUMBER OF LOAD CASES IN A BLOCK -----
      NLDB=KSIZ/NEQ
C-----CALCULATE THE SIZE OF A BLOCK -----
      NLEN=NLDB*NEQ
C
      WRITE(NTM,1000) NLBLK
C-----INITIALIZATION -----
      NDIP=1
      DO 500 IB=1,NLBLK
        IF(IB.EQ.NLBLK) THEN
          NDIN=LOAD
          NLEN=(LOAD-NDIP+1)*NEQ
          ELSE
C
          NDIN=NDIP+NLDB-1
        ENDIF
C
        WRITE(NTM,1100) IB,NDIP,NDIN
C-----ZERO BLOCK SPACE -----
        CALL RZERO(STIF,NLEN)
C-----LOOP ON ALL ELEMENT GROUPS -----
        DO 400 I=1,NUMEL
          NUML=NEL(I)
          CALL SOPEN(NT1,ELNAM(1,I),-LENEL(I))
C-----LOOP ON ALL ELEMENT OF GROUP -----
          DO 300 J=1,NUML
            CALL ELREAD(J,LM,ESTF,R,EMAS,ND,LOAD,MAXD,0)
            KK=0
            DO 200 K=NDIP,NDIN
              DO 100 L=1,ND
                LL=LM(L)
                IF(LL.GT.NEQ.OR.LL.LE.0) GOTO 100
                IF(IB.EQ.1.AND.K.EQ.1) THEN
C-----ASSEMBLE GLOBAL MASS VECTOR -----
                  TMAS(LL)=TMAS(LL)+EMAS(L)
                ENDIF
C
C-----ASSEMBLE GLOBAL LOAD VECTORS -----
                M=KK+LL
                STIF(M)=STIF(M)+R(L,K)
  100         CONTINUE
C
              KK=KK+NEQ
  200       CONTINUE
C
  300     CONTINUE
C
          CALL FCLOSE(NT1)
  400   CONTINUE
C
C-----SAVE BLOCK AND BLOCK PARAMETERS -----
        NUMBLK=NKBLK+IB
        KA(1,NUMBLK)=NUMBLK
        KA(2,NUMBLK)=1
        KA(3,NUMBLK)=NDIP
        KA(4,NUMBLK)=NDIN
        CALL SAVBLK(STIF,KA(1,NUMBLK),KA,NLEN,NUMBLK,1)
        NDIP=NDIN+1
  500 CONTINUE
C
      RETURN
 1000 FORMAT(/'    ASSEMBLING LOAD BLOCKS'/
     * ' THERE ARE ',I3,' LOAD BLOCKS'/)
 1100 FORMAT(' ASSEMBLING LOAD BLOCK NUMBER ',I3/
     *       ' FIRST LOAD CASE IN BLOCK IS  ',I3/
     *       ' LAST LOAD CASE IN BLOCK IS   ',I3)
      END
C
C----------------------------------------------------------------- MLOAD
      SUBROUTINE MLOAD(STIF,RT,ID,NN,KA,NA,NEQ,LOAD,
     *  NUMNP,NKBLK,NLBLK,KSIZ)
C-----------------------------------------------------------------------
C     ROUTINE TO ADD GLOBAL LOADS TO LOAD VECTOR
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RT(6),ID(NUMNP,6),NN(3),KA(4,1),
     * STIF(KSIZ),NA(1),LOD(1)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
C *** CHANGED: ADDED. GFORTRAN HAS AN INTRINSIC NAMED FREE (CRAY
C *** POINTERS) THAT WOULD OTHERWISE HIDE SIMPAL'S FREE INPUT ROUTINE.
      EXTERNAL FREE
C
C-----INITIALIZATION -----
      LBLKRD=0
      NFLB=NKBLK+1
      NLLB=NKBLK+NLBLK
C
      WRITE(NTM,1000)
      WRITE(NOT,1000)
C-----READ LOAD LINE -----
      DO 400 I=1,9999
        CALL FREE
        CALL FREEI(' ',NN,3)
        N1 = NN(1)
        IF(N1.LE.0) GOTO 500
        CALL FREEI('E',LOD,1)
        CALL FREER('D',RT,6)
C-----CHECK LOAD CASE NUMBER -----
        IF(LOD(1).LE.0) LOD(1) = 1
        IF(LOD(1).LE.LOAD) GOTO 50
        WRITE(NTM,1100)
        CALL FREEPT
        GOTO 400
C----- FIND APPROPRIATE LOAD BLOCK -----
   50   DO 75 IH=NFLB,NLLB
          KF=KA(3,IH)
          KL=KA(4,IH)
          NLEN=(KL-KF+1)*NEQ
          IF(LOD(1).LE.KL.AND.LOD(1).GE.KF) THEN
            IF(IH.NE.LBLKRD) THEN
C *** CHANGED: PRINTED AS "IF(I.NE.1) CALL SAVBLK(...)". THAT SAVES
C *** BLOCK KA(1,0) IF THE FIRST LOAD LINE WAS SKIPPED. TEST LBLKRD.
              IF(LBLKRD.GT.0) CALL SAVBLK(STIF,KA(1,LBLKRD),
     *           KA,LNLEN,LBLKRD,1)
              CALL GETBLK(STIF,KA(1,IH),NA,IH)
              LBLKRD=IH
              LNLEN=NLEN
            ENDIF
C
            GOTO 80
          ENDIF
C
   75   CONTINUE
C-----SET GENERATION PARAMETERS -----
   80   N2 = NN(2)
        N3 = NN(3)
        IF(N2.LE.N1) N2 = N1
        IF(N3.LE.0) N3 = 1
C-----FORM LOADS FOR CURRENT LINE -----
        MT=(LOD(1)-KF)*NEQ
          DO 200 J=N1,N2,N3
            WRITE(NTM,1200) J,LOD(1),(RT(K),K=1,6)
            WRITE(NOT,1200) J,LOD(1),(RT(K),K=1,6)
            DO 100 K=1,6
              IDOF = ID(J,K)
              IF(IDOF.GT.NEQ.OR.IDOF.LE.0) GOTO 100
              M=MT+IDOF
              STIF(M)=STIF(M) + RT(K)
  100       CONTINUE
C
  200     CONTINUE
C
  400 CONTINUE
  500 WRITE(NTM,1300)
      WRITE(NOT,1300)
C
C-----LAST LOAD BLOCK SAVED -----
C *** CHANGED: GUARDED. AS PRINTED, NO LOAD LINES -> SAVES BLOCK 0.
      IF(LBLKRD.GT.0)
     *CALL SAVBLK(STIF,KA(1,LBLKRD),KA,LNLEN,LBLKRD,1)
C
  900 RETURN
 1000 FORMAT(/,32X,'GLOBAL  LOAD  DATA',//,7X,'LOAD',32X,
     * 'LOADS',/,1X,'NODE  CASE',5X,'F(X)',7X,'F(Y)',7X,
     * 'F(Z)',7X,'M(X)',7X,'M(Y)',7X,'M(Z)',/)
 1100 FORMAT(' LOAD CASE IS LARGER THAN THE TOTAL',
     * ' NUMBER OF LOAD CASES ',/,' THE FOLLOWING',
     * ' LOAD LINE WILL BE SKIPPED ',/)
 1200 FORMAT(1X,I4,3X,I2,1X,6(2X,F9.2) )
 1300 FORMAT(//)
      END
