************************************************************
      SUBROUTINE STRUCI(X,Q2,ITAR,UP,UPB,DO,DOB,ST,CH,BO,GL)
************************************************************
*
*   RETURNS LO/NLO PARTON UNPOLARIZED PARTON DISTRIBUTIONS
*       (!! NOTE: ALWAYS X TIMES THE DISTRIBUTION !!) 
*
*   CURRENTLY AVAILABLE SETS: GRV-98, CTEQ5, CTEQ6
*   NNPDF31@NLO ADDED Jul 2017 BY E.R. NOCERA (NNPDF COLLABORATION)
*
      IMPLICIT DOUBLE PRECISION (A-H,J-Z)
*
      INTEGER IIP,IIF,IRT,IORD
*
      DOUBLE PRECISION XPDFLH(-6:6)
*      
      COMMON /IPART/ IIP,IIF
      COMMON /ORDER/ IORD
*
****************************************************
* PARTON DISTRIBUTIONS FOR PROTON TARGET (IH1=IH2=0)
****************************************************
*
      IF (IIP.EQ.0) THEN   
*
****************
* GRV-98 (IIP=0)
****************
         STOP
*
*
      ELSE IF (IIP.EQ.1) THEN
*
***************
* CTEQ5 (IIP=1)
***************
*         
         IF(IORD.EQ.0) THEN     !LO
            UP=Ctq5Pd(3,1,X,DSQRT(Q2),Irt)*X
            DO=Ctq5Pd(3,2,X,DSQRT(Q2),Irt)*X 
            UPB=Ctq5Pd(3,-1,X,DSQRT(Q2),Irt)*X
            DOB=Ctq5Pd(3,-2,X,DSQRT(Q2),Irt)*X
            ST=Ctq5Pd(3,3,X,DSQRT(Q2),Irt)*X
            CH=Ctq5Pd(3,4,X,DSQRT(Q2),Irt)*X
            BO=Ctq5Pd(3,5,X,DSQRT(Q2),Irt)*X
            GL=Ctq5Pd(3,0,X,DSQRT(Q2),Irt)*X 
         ELSE                   !NLO
            UP=Ctq5Pd(1,1,X,DSQRT(Q2),Irt)*X
            DO=Ctq5Pd(1,2,X,DSQRT(Q2),Irt)*X 
            UPB=Ctq5Pd(1,-1,X,DSQRT(Q2),Irt)*X
            DOB=Ctq5Pd(1,-2,X,DSQRT(Q2),Irt)*X
            ST=Ctq5Pd(1,3,X,DSQRT(Q2),Irt)*X
            CH=Ctq5Pd(1,4,X,DSQRT(Q2),Irt)*X
            BO=Ctq5Pd(1,5,X,DSQRT(Q2),Irt)*X
            GL=Ctq5Pd(1,0,X,DSQRT(Q2),Irt)*X         
         ENDIF
*
      ELSE IF (IIP.EQ.2) THEN
*
***************
* CTEQ6 (IIP=2)
***************
*         
*   (!! NOTE: LO/NLO AUTOMATICALLY SELECTED HERE !!)
*
         UP=Ctq6Pdf(1,X,DSQRT(Q2))*X
         DO=Ctq6Pdf(2,X,DSQRT(Q2))*X 
         UPB=Ctq6Pdf(-1,X,DSQRT(Q2))*X
         DOB=Ctq6Pdf(-2,X,DSQRT(Q2))*X
         ST=Ctq6Pdf(3,X,DSQRT(Q2))*X
         CH=Ctq6Pdf(4,X,DSQRT(Q2))*X
         BO=Ctq6Pdf(5,X,DSQRT(Q2))*X
         GL=Ctq6Pdf(0,X,DSQRT(Q2))*X

c         UP=0.D0
c         DO=0.D0
c         UPB=0.D0
c         DOB=0.D0
c         ST=0.D0
         CH=0.D0
         BO=0.D0
c         GL=0.D0
*
      ELSE IF (IIP.EQ.3.OR.IIP.EQ.4) THEN
*
***************
* NNPDF31@NLO (IIP=3)
***************
*
         CALL EVOLVEPDFM(4,X,DSQRT(Q2),XPDFLH)
         UP=XPDFLH(2)
         DO=XPDFLH(1)
         UPB=XPDFLH(-2)
         DOB=XPDFLH(-1)
         ST=XPDFLH(3)
         CH=XPDFLH(4)
         BO=XPDFLH(5)
         GL=XPDFLH(0)
*
      ENDIF
*
**************************************
* SWITCH TO ANTIPROTON TARGET (ITAR=1)
**************************************
      IF(ITAR.EQ.1) THEN
         UHELP=UP
         DHELP=DO
         UP=UPB
         UPB=UHELP
         DO=DOB
         DOB=DHELP
      ENDIF
*
*************************************
* SWITCH TO ISOSCALAR TARGET (ITAR=2)
*************************************
      IF(ITAR.EQ.2) THEN
         UDN=0.5D0*(UP+DO)
         UDBN=0.5D0*(UPB+DOB)
         UP=UDN
         DO=UDN
         UPB=UDBN
         DOB=UDBN
      ENDIF
*
      RETURN
      END
*
*********************************************************************
      SUBROUTINE STRUCF(Z,Q2,DPU,DPUB,DPD,DPDB,DPS,DPSB,DPC,DPCB,DPG)
*********************************************************************
*
*   RETURNS LO/NLO UNPOLARIZED FRAGMENTATION FUNCTIONS
*      (!! NOTE: N O T  Z TIMES THE DISTRIBUTION !!) 
*
*   CURRENTLY AVAILABLE SETS: KRETZER, KKP, DSS, NNFF
*      (!! NOTE: CURRENTLY ONLY PI^0 PRODUCTION !!)
*
      IMPLICIT DOUBLE PRECISION (A-H,J-Z)
      INTEGER IIP,IIF,IORD,IHAD
      INTEGER IH,IC,IO
      DOUBLE PRECISION XFFLHP(-6:6),XFFLHM(-6:6),XFFLH(-6:6)
*
      DIMENSION UFF(2), DFF(2), SFF(2), CFF(2), BFF(2), DH(0:10)
*
      COMMON /IPART/ IIP,IIF
      COMMON /ORDER/ IORD
      COMMON /FRAGCUT/ZCUT

*
      IF (IIF.EQ.0) THEN   
*
*************
* DSS (IIF=0)
*************
*     
        IH=1 !1: PION 2: KAON
        IC=0  !CHARGE  (+1, -1, 0)
        if(q2.lt.0.45) q2=0.45d0
       call FDSS (IH,IC,1, Z, Q2, U, UB, D, DB, S, SB, C, B, GL)
*
        DPU=0.D0
        DPUB=0.D0
        DPD=0.D0
        DPDB=0.D0
        DPS=0.D0
        DPSB=0.D0
        DPC=0.D0
        DPCB=0.D0
        DPG=0.D0

         dpu = U/Z
         dpub = UB/Z
         dpd = D/Z 
         dpdb = DB/Z
         dps = S/Z 
         dpsb = SB/Z
         dpg = GL /Z 

C...K0-SHORT  u<->d
c         DPU  = D/Z
c         DPUB = DB/Z
c         DPD  = U/Z
c         DPDB = UB/Z
        IF(Z.LT.ZCUT) THEN
            DPU=0.D0
            DPUB=0.D0
            DPD=0.D0
            DPDB=0.D0
            DPS=0.D0
            DPSB=0.D0
            DPG=0.D0
         ENDIF

*
      ELSEIF (IIF.EQ.1) THEN
*
*****************
* KRETZER (IIF=1)
*****************
*

*$$$$$$$$OPTION DISABLED
*
cc        CALL PKHFF(IORD + 1,1,Z,Q2,uff,dff,sff,cff,bff,gff)  ! 1: PI+ 2: PI- FRAGMENATION
cc            CALL PKHFF(IORD + 3,1,Z,Q2,uff,dff,sff,cff,bff,gff)  ! 1: K+ 2: K- FRAGMENATION
*
C...Pi
c         DPU = uff(1)/2.
c         DPUB = uff(2)/2.
c         DPD = dff(1)/2.
c         DPDB = dff(2)/2
C...K0-SHORT  u<->d
c         DPU = dff(1)/2.
c         DPUB = dff(2)/2.
c         DPD = uff(1)/2.
c         DPDB = uff(2)/2
c...
c         DPS = sff(1)/2.
c         DPSB = sff(2)/2.
c         DPC = cff(1)/2.
c         DPCB = cff(2)/2.
 
         DPC=0.D0
         DPCB=0.D0

c         DPG = gff/2.

         DPC=0.D0
         DPCB=0.D0        
C...Pi-
         DPU=uff(1)
         DPUB=uff(2)
         DPD=dff(1)
         DPDB=dff(2)
         DPS=sff(1)
         DPSB=sff(2)
         DPG=gff

c            DPU=0.D0
c            DPUB=0.D0
c            DPD=0.D0
c            DPDB=0.D0
c            DPS=0.D0
c            DPSB=0.D0
c            DPG=0.D0

C ...   
         IF(Z.LT.ZCUT) THEN
            DPU=0.D0
            DPUB=0.D0
            DPD=0.D0
            DPDB=0.D0
            DPS=0.D0
            DPSB=0.D0
            DPG=0.D0
         ENDIF

      ELSEIF (IIF.EQ.2.OR.IIF.EQ.3) THEN

         CALL EVOLVEPDFM(1,Z,DSQRT(Q2),XFFLHP)
         CALL EVOLVEPDFM(2,Z,DSQRT(Q2),XFFLHM)

         DPU  = ( XFFLHP(2)  + XFFLHM(2)  ) / 2D0 / Z
         DPUB = ( XFFLHP(-2) + XFFLHM(-2) ) / 2D0 / Z
         DPD  = ( XFFLHP(1)  + XFFLHM(1)  ) / 2D0 / Z
         DPDB = ( XFFLHP(-1) + XFFLHM(-1) ) / 2D0 / Z
         DPS  = ( XFFLHP(3)  + XFFLHM(3)  ) / 2D0 / Z
         DPSB = ( XFFLHP(-3) + XFFLHM(-3) ) / 2D0 / Z
         DPC  = ( XFFLHP(4)  + XFFLHM(4)  ) / 2D0 / Z
         DPCB = ( XFFLHP(5)  + XFFLHM(5)  ) / 2D0 / Z
         DPG  = ( XFFLHP(0)  + XFFLHM(0)  ) / 2D0 / Z

         IF(Z.LT.ZCUT) THEN
            DPU  = 0.D0
            DPUB = 0.D0
            DPD  = 0.D0
            DPDB = 0.D0
            DPS  = 0.D0
            DPSB = 0.D0
            DPG  = 0.D0
         ENDIF

      ELSEIF (IIF.EQ.4.OR.IIF.EQ.5.OR.IIF.EQ.6) THEN

         CALL EVOLVEPDF(Z,DSQRT(Q2),XFFLH)

         DPU  = ( XFFLH(+2) ) / 4d0 / Z
         DPUB = ( XFFLH(-2) ) / 4d0 / Z
         DPD  = ( XFFLH(+1) ) / 4d0 / Z
         DPDB = ( XFFLH(-1) ) / 4d0 / Z
         DPS  = ( XFFLH(+3) ) / 4d0 / Z
         DPSB = ( XFFLH(-3) ) / 4d0 / Z
         DPC  = ( XFFLH(+4) ) / 4d0 / Z
         DPCB = ( XFFLH(+5) ) / 4d0 / Z
         DPG  = ( XFFLH(0) )  / 4d0 / Z

         IF(Z.LT.ZCUT) THEN
            DPU  = 0.D0
            DPUB = 0.D0
            DPD  = 0.D0
            DPDB = 0.D0
            DPS  = 0.D0
            DPSB = 0.D0
            DPG  = 0.D0
         ENDIF

      ENDIF
*
C!!!!!!!!!!!!!
ccc      DPG=0.
C!!!!!!!!!!!!!

      RETURN
      END

