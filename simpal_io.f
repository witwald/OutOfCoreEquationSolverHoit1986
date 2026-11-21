C=======================================================================
C
C REPLACEMENTS FOR THE SIMPAL UTILITY ROUTINES CALLED BY THE HOIT (1986)
C ASSEMBLY LISTING. THE PAPER DOES NOT PUBLISH THESE; ONLY THEIR CALLING
C SEQUENCES ARE KNOWN FROM THE LISTING. THESE VERSIONS ARE WRITTEN TO
C THOSE CALLING SEQUENCES.
C
C SOPEN  (NT,NAME,LEN)    OPEN FILE NAME(4 CHAR) ON UNIT NT.
C                         LEN<0 => DIRECT ACCESS, RECL=-LEN BYTES
C FCLOSE (NT)             CLOSE UNIT NT
C ELREAD (J,LM,ESTF,R,EMAS,ND,LOAD,MAXD,IFLAG)
C                         READ ELEMENT J OF THE GROUP OPEN ON NT1.
C                         IFLAG=1: CONNECTIVITY (ND,LM) ONLY - IT IS
C                         AT THE FRONT OF THE RECORD, AS THE PAPER
C                         DESCRIBES. IFLAG=0: WHOLE RECORD.
C ELWRIT (...)            WRITE AN ELEMENT RECORD (TEST HELPER)
C RZERO  (X,N)            ZERO N REAL*8 WORDS
C SAVBLK (X,KA,MA,LEN,IREC,IFLAG)
C                         WRITE BLOCK IREC TO THE DIRECT ACCESS
C                         STIFFNESS/LOAD FILE NT3. IFLAG=0 STIFFNESS
C                         (MA = DIAGONAL POINTERS), IFLAG=1 LOAD
C                         (MA NOT MEANINGFUL - WRITTEN AS ZEROS)
C GETBLK (X,KA,MA,IREC)   READ BLOCK IREC FROM NT3 (X AND KA ONLY;
C                         MLOAD PASSES NA AS MA, SO MA IS NOT
C                         OVERWRITTEN)
C FREE, FREEI, FREER, FREEPT
C                         FREE-FORMAT INPUT USED BY MLOAD. ONE LOAD
C                         LINE PER RECORD ON UNIT NIN, LIST-DIRECTED:
C                           N1 N2 N3 LOD FX FY FZ MX MY MZ
C                         N1 <= 0 (E.G. "0 /") ENDS THE LOAD DATA.
C
C ELEMENT RECORD LAYOUT (DIRECT ACCESS, ONE RECORD PER ELEMENT):
C    ND, LM(MAXD), ESTF(MAXD*(MAXD+1)/2), R(MAXD,LOAD), EMAS(MAXD)
C ESTF IS THE UPPER TRIANGLE STORED BY COLUMNS, K=1...ND, L=1...K,
C WHICH IS THE ORDER ASEMBL CONSUMES IT IN.
C
C=======================================================================
C
      SUBROUTINE SOPEN(NT,NAME,LEN)
      CHARACTER*1 NAME(4)
      CHARACTER*16 FN
      INTEGER I,N,NT,LEN
      LOGICAL OP
      FN='simpal_'
      N=7
      DO 10 I=1,4
        IF(NAME(I).EQ.' ') GOTO 10
        N=N+1
        FN(N:N)=NAME(I)
   10 CONTINUE
      FN(N+1:N+3)='.da'
      INQUIRE(UNIT=NT,OPENED=OP)
      IF(OP) CLOSE(NT)
      IF(LEN.LT.0) THEN
        OPEN(NT,FILE=FN,ACCESS='DIRECT',FORM='UNFORMATTED',
     *       RECL=-LEN,STATUS='UNKNOWN')
      ELSE
        OPEN(NT,FILE=FN,FORM='UNFORMATTED',STATUS='UNKNOWN')
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE FCLOSE(NT)
      CLOSE(NT)
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE ELREAD(J,LM,ESTF,R,EMAS,ND,LOAD,MAXD,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LM(MAXD),ESTF(*),R(*),EMAS(*)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
      KMAX=MAXD*(MAXD+1)/2
      IF(IFLAG.EQ.1) THEN
        READ(NT1,REC=J) ND,(LM(I),I=1,MAXD)
      ELSE
        READ(NT1,REC=J) ND,(LM(I),I=1,MAXD),(ESTF(I),I=1,KMAX),
     *       (R(I),I=1,MAXD*LOAD),(EMAS(I),I=1,MAXD)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE ELWRIT(NT,J,LM,ESTF,R,EMAS,ND,LOAD,MAXD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION LM(MAXD),ESTF(*),R(*),EMAS(*)
      KMAX=MAXD*(MAXD+1)/2
      WRITE(NT,REC=J) ND,(LM(I),I=1,MAXD),(ESTF(I),I=1,KMAX),
     *     (R(I),I=1,MAXD*LOAD),(EMAS(I),I=1,MAXD)
      RETURN
      END
C
C-----------------------------------------------------------------------
      INTEGER FUNCTION LENELR(MAXD,LOAD)
C     RECORD LENGTH IN BYTES OF AN ELEMENT RECORD (4-BYTE INTEGERS,
C     8-BYTE REALS; GFORTRAN DEFAULT RECL UNIT IS THE BYTE)
      KMAX=MAXD*(MAXD+1)/2
      LENELR=4*(1+MAXD)+8*(KMAX+MAXD*LOAD+MAXD)
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE RZERO(X,N)
      REAL*8 X(*)
      DO 10 I=1,N
   10 X(I)=0.0D0
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE SAVBLK(X,KA,MA,LEN,IREC,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),KA(4),MA(*)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
      COMMON /DAPAR/ MAXCDA
      KK=2
      IF(IFLAG.NE.0) KK=3
      CALL DABLK(NT3,IREC,X,LEN,LEN,MA,MAXCDA,KA,KK)
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE GETBLK(X,KA,MA,IREC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),KA(4),MA(*)
      COMMON /CLT/ ISTAG,NUMNP,NEQ,NEQS,LOAD,MAXD,MAXNA,
     * NUMEL,NEL(20),LENEL(20),NKBLK,NLBLK,KSIZ,KMEM,
     * KYDYN,KYSUB,KYLIN
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
      COMMON /DAPAR/ MAXCDA
      CALL DABLK(NT3,IREC,X,NW,KSIZ,MA,MAXCDA,KA,4)
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE FREE
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /IOLIST/ NTM,NTR,NIN,NOT,NT1,NFL,NT2,NT3
      COMMON /FREEC/ RFV(6),IFV(4)
      DO 10 I=1,4
   10 IFV(I)=0
      DO 20 I=1,6
   20 RFV(I)=0.0D0
      READ(NIN,*,END=30) (IFV(I),I=1,4),(RFV(I),I=1,6)
   30 RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE FREEI(KEY,IA,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 KEY
      DIMENSION IA(N)
      COMMON /FREEC/ RFV(6),IFV(4)
      IF(KEY.EQ.'E') THEN
        IA(1)=IFV(4)
      ELSE
        DO 10 I=1,N
   10   IA(I)=IFV(I)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE FREER(KEY,RA,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 KEY
      DIMENSION RA(N)
      COMMON /FREEC/ RFV(6),IFV(4)
      DO 10 I=1,N
   10 RA(I)=RFV(I)
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE FREEPT
      RETURN
      END
