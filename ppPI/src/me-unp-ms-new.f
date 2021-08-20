*****************************************
      DOUBLE PRECISION FUNCTION DPLUS(XX)
*****************************************
*
*   FUNCTION TO BE INTEGRATED BY VEGAS:
*   CALLS PARTON DENSITIES, FRAGMENTATION FUNCTIONS, AND
*   ALL PARTONIC CROSS SECTIONS
* 
      IMPLICIT REAL*8 (A-H,L-Z)
*
      INTEGER III
      DIMENSION XX(10),IA1(10),IA2(6),IA3(9),CC(16),F01(16),F02(16),
     #     FHODEL1(16),FHODEL2(16),FHOREST1(16),FHOREST2(16),
     #     FDELMU1(16),FDELMU2(16),BORN(16),DMU(16),
     #     GRRT(16),GRRC(16),GPPT(16),GPPC(16)
*
      COMMON /PARAMS/ CF,NC,NF,PI,HC2
      COMMON /CONS/ PPI,S,GV,GW,NNC,GTR,CCF,VC
      COMMON / PREFAC / CC
cms      COMMON /ORDE/ AL,CQ,V1,V2,V3,V4
      COMMON /KINVAR/ SQS,PT0,ETA0
      COMMON /PTREAL/ PPT
      COMMON /EDF/ J0
      COMMON /SCALESF/ SCFAC,SCMU,SCFRAG
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      COMMON /HADR/ IH1,IH2
      COMMON /PTINT/ PTDO,PTUP
      COMMON /IPTFL/ IPT 
      COMMON /YINT/ YDO,YUP
      COMMON /IYFL/ IY 
      COMMON /ORDER/ IORD
      COMMON /XSECT/ ISIGM
      COMMON /CONTROL/ZMINIMUM
*
      DATA IA1 /1,3,5,6,11,12,13,14,15,16/
      DATA IA2 /1,3,5,11,13,14/
      DATA IA3 /1,3,5,8,9,10,11,13,14/
C
      PPI=PI
      NNC=NC
      N=NNC
      CCF=CF
C
C *** THESE PARAMETERS ARE NOT TO BE CHANGED ! ***
cms      AL=1.D0                                 
cms      CQ=0.D0
C ************************************************
      S=SQS**2
C
      VC=(NNC)**2-1.
cms      V1=VC**2/NNC
cms      V2=VC/NNC
cms      V3=(NNC**4-1.)/2./NNC**2
cms      V4=VC**2/2./NNC**2
C
      GTR=NF/2.
C

C*******************************************************************
C
***********************************
* FIVE-DIMENSIONAL MC - INTEGRATION
************************** ********
      IF(IY.EQ.0) THEN
         ETA=ETA0
         XJETA = 1.D0
         IF(IPT.EQ.0) THEN 
            PT=PT0
            XJPT = 1.D0
         ELSE
C     ***  PT-INTEGRATION
            PT=PTDO+XX(4)*(PTUP-PTDO)
            XJPT = PTUP - PTDO
         ENDIF           
C ****************************************
       ELSE
C *** Y-INTEGRATION 
          ETA=YDO+XX(4)*(YUP-YDO)
          IF(IPT.EQ.0) THEN 
             PT=PT0
             XJPT = 1.D0    
         ELSE
C ***  PT-INTEGRATION
            PTUP2=SQS/2.D0/COSH(ETA)
            PTUP1=PTUP
            PTUP=MIN(PTUP1,PTUP2)
            PT=PTDO+XX(5)*(PTUP-PTDO)
            XJPT = PTUP - PTDO
         ENDIF           
cms         ETA=YDO+XX(4)*(YUP-YDO)
         XJETA = YUP-YDO
      ENDIF
C ****************************************
      XT=2.D0*PT/SQS
      PPT=PT
      Q2FAC=SCFAC*PT**2
      Q2MU=SCMU*PT**2
      Q2FRAG=SCFRAG*PT**2
c
c      write(6,*) pt**2,q2fac,q2mu,q2frag
c      IF(Q2FRAG.LE.1.D0) Q2FRAG=1.D0
c      IF(Q2FAC.LT.1.D0) Q2FAC=1.D0
c      IF(Q2MU.LT.1.D0) Q2MU=1.D0     
c      write(6,*) pt**2,q2fac,q2mu,q2frag
c
      ALPORD=4.D0*PI*ALPHAS(Q2MU)
      ALPASHO=ALPORD*DBLE(IORD)
*************
*   GV AND GW
*************
      GV=1.D0-PT/SQS*DEXP(-ETA)
      GW=PT**2/S/GV/(1.-GV)
********************
*   X3 - INTEGRATION
********************
      X3MIN=1.D0-GV+GV*GW
      X3MAX=1.D0
      X3=X3MIN+(X3MAX-X3MIN)*XX(3)
      IF(ZMINIMUM.GT.X3MIN) ZMINIMUM=X3MIN
*******************
*   V - INTEGRATION
*******************
      VMIN=GV*GW/X3
      VMAX=1.D0-(1.D0-GV)/X3
      V=VMIN+(VMAX-VMIN)*XX(2)
*******************
*   W - INTEGRATION
*******************
      WMIN=GV*GW/X3/V
      WMAX=1.D0
      W=WMIN+(WMAX-WMIN)*XX(1)
      UN=1.D0
******************************************
*   PARTON MOMENTUM FRACTIONS X1, X2 (NLO)
******************************************
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
*******************************************
*   PARTON MOMENTUM FRACTIONS X1, X2 (BORN)
*******************************************
      BX1=GV*GW/V/X3
      BX2=X2
*********************
*   PARTONIC S (BORN)
*********************
      SHD=BX1*BX2*S
**************************************
*   JACOBIANS FOR VEGAS (BORN AND NLO)
**************************************
      BXJAC=(X3MAX-X3MIN)*(VMAX-VMIN)/BX1/BX2/X3**2
      XJAC=(X3MAX-X3MIN)*(VMAX-VMIN)/X1/X2/X3**2
*******************************************
*   CALL FRAGMENTATIONS FUNCTIONS AND PDF's
*******************************************
      CALL STRUCF(X3,Q2FRAG,DPU,DPUB,DPD,DPDB,DPS,DPSB,DPC,DPCB,DPG)
*
      CALL STRUCI(BX1,Q2FAC,IH1,UP1,UPB1,DO1,DOB1,ST1,CH1,BO1,GL1)
      CALL STRUCI(BX2,Q2FAC,IH2,UP2,UPB2,DO2,DOB2,ST2,CH2,BO2,GL2)
*********************************************************
*   CALCULATE ALL RELEVANT COMBINATIONS OF PDF's AND FF's
*********************************************************
      CALL STRU(UP1,UPB1,DO1,DOB1,ST1,CH1,GL1,           
     &          UP2,UPB2,DO2,DOB2,ST2,CH2,GL2,           
     &          DPU,DPUB,DPD,DPDB,DPS,DPSB,DPC,DPCB,DPG,GRRT,GRRC) 
************************************************
*   CALCULATE AND SUM UP ALL BORN CROSS SECTIONS
************************************************
      CALL FBOR(V,SHD,F01)                          
      CALL FBOR(1.D0-V,SHD,F02)      
*                    
      DO 12 J0=1,16
         BORN(J0) = (F01(J0)*GRRT(J0) + F02(J0)*GRRC(J0)) * BXJAC 
 12   CONTINUE
*
***************************************************
*   LOOP OVER ALL CONTRIBUTION SUBPROCESSES IN NLO:
*     J0=1  : QI QK  ---> QI
*     J0=2  : QI QK  ---> G
*     J0=3  : QI QBK ---> QI
*     J0=4  : QI QBK ---> G
*     J0=5  : QI QBI ---> QK
*     J0=6  : QI QI  ---> QI 
*     J0=7  : QI QI  ---> G
*     J0=8  : QI G   ---> QK
*     J0=9  : QI G   ---> QBK
*     J0=10 : QI G   ---> QBI
*     J0=11 : QI QBI ---> QI
*     J0=12 : QI QBI ---> G
*     J0=13 : QI G   ---> QI
*     J0=14 : QI G   ---> G
*     J0=15 : G  G   ---> G 
*     J0=16 : G  G   ---> QJ
***************************************************
*
**************
*   INITIALIZE
**************
      ONE = 1.D0
      GHD = 0.D0
      GHE = 0.D0
*
      DO 1 J0=1,16
         FHODEL1(J0) = 0.D0
 1    CONTINUE
********************
*   PARTONIC S (NLO)
********************
         SHD = X1*X2*S
*************************************
*   CALL PDF FOR HADRON 1 WITH NLO X1
*************************************
         CALL STRUCI(X1,Q2FAC,IH1,UP1,UPB1,DO1,DOB1,ST1,CH1,BO1,GL1)
*********************************************************
*   CALCULATE ALL RELEVANT COMBINATIONS OF PDF's AND FF's
*********************************************************
         CALL STRU(UP1,UPB1,DO1,DOB1,ST1,CH1,GL1,           
     #        UP2,UPB2,DO2,DOB2,ST2,CH2,GL2,           
     #        DPU,DPUB,DPD,DPDB,DPS,DPSB,DPC,DPCB,DPG,GPPT,GPPC) 

C
      DO 11 III=1,10
      J0=IA1(III)
      FHODEL1(J0)=(FDEL1(V,X3)+FVWPL1(UN,V,X3)*DLOG(1.-WMIN)+
     #             FVLO1(UN,V,X3)*(DLOG(1.-WMIN))**2/2.D0)/(1.D0-V) -
     #             (WMAX-WMIN)/(1.D0-V) *
     #             (FVWPL1(UN,V,X3)+FVLO1(UN,V,X3)*DLOG(1.-W))/(1.-W)
      FHODEL1(J0)=FHODEL1(J0)/CC(J0)*BXJAC
 11   CONTINUE
C 
      DO 2 J0=1,16
        FHODEL2(J0)=FHODEL1(J0)
  2   CONTINUE
C
      DO 22 III=1,6
      J0=IA2(III)
      FHODEL2(J0)=(FDEL2(V,X3)+FVWPL2(UN,V,X3)*DLOG(1.-WMIN)+
     #             FVLO2(UN,V,X3)*(DLOG(1.-WMIN))**2/2.D0)/(1.D0-V) -
     #             (WMAX-WMIN)/(1.D0-V) *
     #             (FVWPL2(UN,V,X3)+FVLO2(UN,V,X3)*DLOG(1.-W))/(1.-W)
      FHODEL2(J0)=FHODEL2(J0)/CC(J0)*BXJAC
 22   CONTINUE
C
************************************************
*   STORE V**2, W**2, V**3, ETC. IN COMMON BLOCK
************************************************
      CALL PRECALC(V,W,SHD)
*
      DO 3 J0=1,16
      FHOREST1(J0)=(WMAX-WMIN)/(1.D0-V) * (FVWPL1(W,V,X3)/W/(1.-W)+
     #             FVLO1(W,V,X3)/W*DLOG(1.-W)/(1.-W)
     #            +FRESC1(W,V,X3)/W)
C
      FHOREST1(J0)=FHOREST1(J0)/CC(J0)*XJAC
 3    CONTINUE
C
      DO 4 J0=1,16
        FHOREST2(J0)=FHOREST1(J0)
 4    CONTINUE
C
************************************************
*   STORE V**2, W**2, V**3, ETC. IN COMMON BLOCK
************************************************
      VX=1.-V*W
      WX=(1.-V)/(1.-V*W)
      CALL PRECALC(VX,WX,SHD)
*
      DO 44 III=1,9
      J0=IA3(III)
      FHOREST2(J0)=(WMAX-WMIN)/(1.D0-V) * (FVWPL2(W,V,X3)/W/(1.-W)+
     #              FVLO2(W,V,X3)/W*DLOG(1.-W)/(1.-W)
     #             +FRESC2(W,V,X3)/W)
C
      FHOREST2(J0)=FHOREST2(J0)/CC(J0)*XJAC
 44   CONTINUE
C            
      DO J0=1,16
C ***   BORN ONLY ************
C       GHD=BORN(J0) + GHD
C       GHE=0.D0+GHE 
C ****************************
CMS        GHD=( FHODEL1(J0)*GRRT(J0) + FHODEL2(J0)*GRRC(J0) + DMU(J0) )
         GHD=( FHODEL1(J0)*GRRT(J0) + FHODEL2(J0)*GRRC(J0))
     #        * ALPASHO 
     #     +  BORN(J0) +  GHD
        GHE=( FHOREST1(J0)*GPPT(J0) + FHOREST2(J0)*GPPC(J0) )
     #        * ALPASHO + GHE
      ENDDO
      DPLUS=(GHD+GHE)*ALPORD**2
C
        IF (ISIGM.EQ.1) THEN
C ***   FOR DSIGMA/DY/DPT2
          FACIN=PI
        ELSE IF (ISIGM.EQ.2) THEN
C ***   FOR EDSIGMA/D3P
          FACIN=1.D0
        ELSE IF (ISIGM.EQ.3) THEN
C ***   FOR DSIGMA/DPT/DY
          FACIN=PI*2.D0*PT
        ELSE IF (ISIGM.EQ.4) THEN
C ***   FOR PT**3*DSIGMA/DPT/DY
          FACIN=PI*2.D0*PT**4
        ENDIF
C
      DPLUS=FACIN*DPLUS * HC2/PI/S  * XJETA * XJPT
C * PT
      RETURN
      END
C
*
*******************************
      SUBROUTINE PRECALC(V,W,S)
*******************************
*
      IMPLICIT DOUBLE PRECISION (A-Z)
*
      COMMON / SCALES / Q2FAC,Q2MU,Q2FRAG
      COMMON / PARAMS / CF,CA,NF,PI,HC2
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
      L1V  = DLOG(1.D0-V)
      LV   = DLOG(V)
      L1W  = DLOG(1.D0-W)
      LW   = DLOG(W)
      LVW  = DLOG((1.D0-V)/(1.D0-V*W))/(1.D0-W)
      L1VW = DLOG(1.D0-V+V*W)/(1.D0-W)
*
      LMU  = DLOG(Q2MU/S)
      LMS  = DLOG(Q2FAC/S)
      LMSS = DLOG(Q2FRAG/S)
C...
      CACF = CA*CF
      CA2  = CA**2
      CA4  = CA**4
C...
      V2  = V**2
      V3  = V**3
      V4  = V**4
      V5  = V**5
      V6  = V**6
      V7  = V**7
      V8  = V**8
      V9  = V**9
      V10 = V**10
      V11 = V**11
      V12 = V**12
c...
      W2  = W**2
      W3  = W**3
      W4  = W**4
      W5  = W**5
      W6  = W**6
      W7  = W**7
      W8  = W**8      
      W9  = W**9
      W10 = W**10
      W11 = W**11
      W12 = W**12
*
      RETURN
      END
*
*******************************
      SUBROUTINE FBOR(V,SHD,F0)
*******************************
*
*   BORN CROSS SECTIONS [INCL. PHASE SPACE FACTOR 1/v/(1-v)]
*
      IMPLICIT REAL*8(A-H,L-Z)
*
      DIMENSION F0(16) 
*
      COMMON / PARAMS / CF,NC,NF,PI,HC2
cms   COMMON /CONS/ PI,GS,GV,GW,N,GTR,CF,VC
*
      V2  = V**2
      V3  = V**3
      V4  = V**4
      VC  = (NC**2-1.d0)
      NC2 = NC**2
      VM  = (1.d0-V)
      VM2 = VM**2
      PRELO = PI/SHD/V/VM
*
      F0(1)=CF/NC*PRELO*(V2+1.)/VM2
*
      F0(2)=0.D0
*
      F0(3)=CF/NC*PRELO*(V2+1.)/VM2
*
      F0(4)=0.D0
*
      F0(5)=CF/NC*PRELO*(2.*V2-2.*V+1.)
*
      F0(6)=2.*CF/(NC2)*PRELO*
     &  (NC*V4-2.*NC*V3+4.*NC*V2+V2-(3.*NC+1.)*V+NC)/V2/VM2
*    
      F0(7)=0.D0
*
      F0(8)=0.D0
*
      F0(9)=0.D0
*
      F0(10)=0.D0
*
      F0(11)=2.*CF/(NC2)*PRELO*(NC*V4-(3.*NC+1.)*V3+
     &  (4.*NC+1.)*V2-2.*NC*V+NC)/VM2
*
      F0(12)=CF/(NC2)*PRELO*(2.*V2-2.*V+1.)*
     &  (2.*NC2*V2-2.*NC2*V+NC2-1.)/V/VM
*
      F0(13)=1.D0/(2.*NC2)*PRELO*(V2+1.)*((NC2-1.)*V2+
     &  2.*V+(NC2-1.))/V/VM2
*
      F0(14)=1.D0/(2.*NC2)*PRELO*(V2-2.*V+2.)*((NC2-1.)*V2-
     &  2.*NC2*V+2.*NC2)/V2/VM
*
      F0(15)=4.*NC2/VC*PRELO*(3.-V*VM+V/VM2+
     &  VM/V2)
*
      F0(16)=1.D0/(2.*NC)/VC*PRELO*(V2+VM2)*
     &  (2.*NC2*(V2-V)+NC2-1.)/V/VM
*
      RETURN
      END
C
*******************************************
      DOUBLE PRECISION FUNCTION FDEL1(V,X3)
*******************************************
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      BX1=GV*GW/V/X3
      BX2=(1.-GV)/(1.-V)/X3
      SHD=BX1*BX2*GS
      FKEL=AVDEL(V,SHD)
      FDEL1=FKEL/SHD
      RETURN
      END
C
*******************************************
      DOUBLE PRECISION FUNCTION FDEL2(V,X3)
*******************************************
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      BX1=GV*GW/V/X3
      BX2=(1.-GV)/(1.-V)/X3
      SHD=BX1*BX2*GS
      UN=1.D0
      FKELC=(AVDEL(1.-V,SHD)+DLOG(V/(1.-V))*AVWPL(UN,1.-V,SHD)
     # +.5*AVLO(UN,1.-V,SHD)*(DLOG((1.-V)/V))**2)*(1.-V)/V
      FDEL2=FKELC/SHD
      RETURN
      END
C
**********************************************
      DOUBLE PRECISION FUNCTION FVWPL1(W,V,X3)
**********************************************
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
      SH=X1*X2*GS
      RVWPL=AVWPL(W,V,SH)
      FVWPL1=RVWPL/SH
      RETURN
      END
C
**********************************************
      DOUBLE PRECISION FUNCTION FVWPL2(W,V,X3)
**********************************************
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
      VX=1.-V*W
      WX=(1.-V)/(1.-V*W)
      SH=X1*X2*GS
      RVWPLC=(AVWPL(WX,VX,SH)+AVLO(WX,VX,SH)*DLOG(V/VX))*VX/V
      FVWPL2=RVWPLC/SH
      RETURN
      END
C
*********************************************
      DOUBLE PRECISION FUNCTION FVLO1(W,V,X3)
*********************************************
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
      SH=X1*X2*GS
      RVWLO=AVLO(W,V,SH)
      FVLO1=RVWLO/SH
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION FVLO2(W,V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
      VX=1.-V*W
      WX=(1.-V)/(1.-V*W)
      SH=X1*X2*GS
      RVWLOC=AVLO(WX,VX,SH)*VX/V
      FVLO2=RVWLOC/SH
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION FRESC1(W,V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
      SH=X1*X2*GS
      RRESC=(STRUV(W,V,X3,SH)+AVGO(W,V))
      FRESC1=RRESC/SH
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION FRESC2(W,V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,VC
      X1=GV*GW/V/W/X3
      X2=(1.-GV)/(1.-V)/X3
      VX=1.-V*W
      WX=(1.-V)/(1.-V*W)
      SH=X1*X2*GS
      RRESCC=STRUV(WX,VX,X3,SH)+AVGO(WX,VX)
      FRESC2=RRESCC/SH
      RETURN
      END
C
********************************************
      DOUBLE PRECISION FUNCTION AVWPL(W,V,S)
********************************************
*
*   CONTAINS ALL 1/(1-W)+ PIECES
*
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      COMMON/ORDE/AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      M=DSQRT(Q2FAC)
      MP=DSQRT(Q2FRAG)  
      WPLUS=1.

      l1v = dlog(1.-v)
      lv = dlog(v)
      lms = dlog(M**2/S)
      lmss = dlog(MP**2/S)
      Nf = 2.*GTR 

      IF (J0.EQ.1) THEN
      AVWPL =          ((-6*CA*CF**2*(1 + v**2))/((1 - v)**2*v) -
     -     (8*(1 + CA**2)*CF*l1v*(1 + v**2))/((1 - v)**2*v) -
     -     (16*CA*CF**2*lms*(1 + v**2))/((1 - v)**2*v) -
     -     (8*CA*CF**2*lmss*(1 + v**2))/((1 - v)**2*v) -
     -     (4*(11 - 7*CA**2)*CF*lv*(1 + v**2))/((1 - v)**2*v))

      ELSE IF (J0.EQ.2) THEN
      AVWPL = 0.D0

      ELSE IF (J0.EQ.3) THEN
      AVWPL = ((-6*CA*CF**2*(1 + v**2))/((1 - v)**2*v) -
     -     (8*(1 + CA**2)*CF*l1v*(1 + v**2))/((1 - v)**2*v) -
     -     (16*CA*CF**2*lms*(1 + v**2))/((1 - v)**2*v) -
     -     (8*CA*CF**2*lmss*(1 + v**2))/((1 - v)**2*v) +
     -     (4*(5 + 3*CA**2)*CF*lv*(1 + v**2))/((1 - v)**2*v))

      ELSE IF (J0.EQ.4) THEN
      AVWPL = 0.D0

      ELSE IF (J0.EQ.5) THEN
      AVWPL = ((-6*CA*CF**2*(1 - 2*v + 2*v**2))/v -
     -     (8*(3 - CA**2)*CF*l1v*(1 - 2*v + 2*v**2))/v -
     -     (16*CA*CF**2*lms*(1 - 2*v + 2*v**2))/v -
     -     (8*CA*CF**2*lmss*(1 - 2*v + 2*v**2))/v +
     -     (4*(5 + 3*CA**2)*CF*lv*(1 - 2*v + 2*v**2))/v)

      ELSE IF (J0.EQ.6) THEN
      AVWPL = ((-12*CF**2*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v**3) -
     -     (32*CF**2*lms*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v**3) -
     -     (16*CF**2*lmss*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v**3) -
     -     (16*CF*l1v*(3*CA - CA**3 - v - 9*CA*v - 
     -      CA**2*v + 3*CA**3*v +
     -          v**2 + 11*CA*v**2 + CA**2*v**2 - 
     -      3*CA**3*v**2 - 6*CA*v**3 +
     -          2*CA**3*v**3 + 2*CA*v**4))/(CA*(1 - v)**2*v**3) -
     -     (8*CF*lv*(7*CA - 3*CA**3 - 7*v - 21*CA*v + 
     -      3*CA**2*v + 9*CA**3*v +
     -          7*v**2 + 30*CA*v**2 - 3*CA**2*v**2 - 
     -      14*CA**3*v**2 -
     -          14*CA*v**3 + 6*CA**3*v**3 + 9*CA*v**4 - 
     -      5*CA**3*v**4))/
     -      (CA*(1 - v)**2*v**3))

      ELSE IF (J0.EQ.7) THEN
      AVWPL = 0.D0

      ELSE IF (J0.EQ.8) THEN
      AVWPL = 0.D0

      ELSE IF (J0.EQ.9) THEN
      AVWPL = 0.D0

      ELSE IF (J0.EQ.10) THEN
      AVWPL = 0.D0

      ELSE IF (J0.EQ.11) THEN
      AVWPL = ((-12*CF**2*(CA - 2*CA*v + v**2 + 4*CA*v**2 - 
     -      v**3 - 3*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v) -
     -     (32*CF**2*lms*(CA - 2*CA*v + v**2 + 4*CA*v**2 - 
     -      v**3 - 3*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v) -
     -     (16*CF**2*lmss*(CA - 2*CA*v + v**2 + 4*CA*v**2 - 
     -      v**3 -
     -          3*CA*v**3 + CA*v**4))/((1 - v)**2*v) -
     -     (16*CF*l1v*(2*CA - 6*CA*v + 2*CA**3*v + v**2 + 11*CA*v**2 +
     -          CA**2*v**2 - 3*CA**3*v**2 - v**3 - 9*CA*v**3 - 
     -      CA**2*v**3 +
     -          3*CA**3*v**3 + 3*CA*v**4 - CA**3*v**4))/
     -     (CA*(1 - v)**2*v) +
     -     (8*CF*lv*(5*CA + 3*CA**3 - 10*CA*v - 6*CA**3*v + v**2 +
     -          20*CA*v**2 + 7*CA**2*v**2 + 12*CA**3*v**2 - v**3 -
     -          15*CA*v**3 - 7*CA**2*v**3 - 9*CA**3*v**3 + 5*CA*v**4 +
     -          3*CA**3*v**4))/(CA*(1 - v)**2*v))

      ELSE IF (J0.EQ.12) THEN
      AVWPL = ((-16*CA**2*CF*lmss*(1 - 2*v + 2*v**2)*
     -     (CF - CA*v + CA*v**2))/
     -      ((1 - v)*v**2) + (22*CA*CF*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/
     -     (3.*(1 - v)*v**2) +
     -     (16*CF**2*lms*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/((1 - v)*v**2) -
     -     (4*CF*Nf*(1 - 2*v + 2*v**2)*(-2*CA*CF + 2*CA**2*v - 
     -      2*CA**2*v**2))/
     -      (3.*(1 - v)*v**2) -
     -     (8*CF*l1v*(1 - 2*v + 2*v**2)*
     -        (1 + CA**4 + 2*CA**2*v - 2*CA**4*v - 2*CA**2*v**2))/
     -      (CA*(1 - v)*v**2) +
     -     (8*CF*lv*(1 - 2*v + 2*v**2)*
     -        (1 - 5*CA**2 + 4*CA**4 + 2*CA**2*v - 
     -      8*CA**4*v - 2*CA**2*v**2 +
     -          6*CA**4*v**2))/(CA*(1 - v)*v**2))

      ELSE IF (J0.EQ.13) THEN
      AVWPL = ((22*CA*CF*(1 + v**2)*(-2*CA*CF - 2*v + v**2 - 
     -      CA**2*v**2))/
     -      (3.*(1 - v)**2*v**2) -
     -     (4*(1 - 3*CA**2)*CF*lms*(1 + v**2)*
     -        (-2*CA*CF - 2*v + v**2 - CA**2*v**2))/
     -     (CA*(1 - v)**2*v**2) -
     -     (4*CF*Nf*(1 + v**2)*(-2*CA*CF - 2*v + v**2 - CA**2*v**2))/
     -      (3.*(1 - v)**2*v**2) -
     -     (8*CF**2*lmss*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -      CA**2*v**2))/
     -      ((1 - v)**2*v**2) +
     -     (8*CF*l1v*(1 + v**2)*
     -        (1 + 2*CA**2 - CA**4 - 2*v - 6*CA**2*v + v**2 + 
     -      2*CA**2*v**2 -
     -          CA**4*v**2))/(CA*(1 - v)**2*v**2) +
     -     (8*CF*lv*(1 + v**2)*(1 - 5*CA**2 + 4*CA**4 - 2*v + 
     -      8*CA**2*v +
     -          v**2 - 5*CA**2*v**2 + 2*CA**4*v**2))/
     -     (CA*(1 - v)**2*v**2))

      ELSE IF (J0.EQ.14) THEN
      AVWPL = ((-16*CA**3*CF*l1v*(1 - v)*(2 - 2*v + v**2))/v**3 -
     -     (6*CF**2*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     ((1 - v)*v**3) +
     -     (4*(1 - 3*CA**2)*CF*lms*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (CA*(1 - v)*v**3) -
     -     (8*CA*CF*lmss*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     ((1 - v)*v**3) -
     -     (4*CF*lv*(2 - 2*v + v**2)*
     -        (2*CA**2 - 18*CA**4 - 2*CA**2*v + 18*CA**4*v - v**2 +
     -          6*CA**2*v**2 - 9*CA**4*v**2))/(CA*(1 - v)*v**3))

      ELSE IF (J0.EQ.15) THEN
      AVWPL =  ((-256*CA**3*l1v*(1 - v + v**2)**2)/v**3 -
     -     (704*CA**3*(1 - v + v**2)**3)/(3.*(1 - v)**2*v**3) -
     -     (512*CA**3*lms*(1 - v + v**2)**3)/((1 - v)**2*v**3) -
     -     (256*CA**3*lmss*(1 - v + v**2)**3)/((1 - v)**2*v**3) +
     -     (128*CA**3*Nf*(1 - v + v**2)**3)/(9.*(1 - v)**2*v**3) +
     -     (256*CA**3*lv*(1 - v + v**2)**2*(5 - 5*v + 4*v**2))/
     -      ((1 - v)**2*v**3))

      ELSE IF (J0.EQ.16) THEN
      AVWPL =  ((-16*CA**3*CF*l1v*(1 - v)*(1 - 2*v + 2*v**2))/v**2 -
     -     (16*CA*CF**2*lmss*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      ((1 - v)*v**2) + (6*CF**2*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/((1 - v)*v**2) +
     -     (4*CF*lv*(1 - 2*v + 2*v**2)*
     -        (1 - 10*CA**2 + 9*CA**4 + 2*CA**2*v - 18*CA**4*v -
     -          2*CA**2*v**2 + 14*CA**4*v**2))/(CA*(1 - v)*v**2) +
     -     (8*CA*CF*lms*(-4*CA*CF - 6*v + 10*CA**2*v + 8*v**2 -
     -          24*CA**2*v**2 - 4*v**3 + 32*CA**2*v**3 -24*CA**2*v**4 +
     -          8*CA**2*v**5))/((1 - v)**2*v**2))

      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION AVDEL(V,S)
C     TERME EN DELTA(1-W)
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      COMMON/ORDE/AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      M=DSQRT(Q2FAC)
      MP=DSQRT(Q2FRAG)  
      DELW=1.
C    TERME EN DELTA

      l1v = dlog(1.-v)
      lv = dlog(v)
      lmu = dlog(Q2MU/S)
      lms = dlog(M**2/S)
      lmss = dlog(MP**2/S)
      Nf = 2.*GTR

      IF (J0.EQ.1) THEN
      AVDEL = ((-12*CA*CF**2*lms*(1 + v**2))/((1 - v)**2*v) +
     -     (8*CA*CF**2*l1v*lms*(1 + v**2))/((1 - v)**2*v) -
     -     (6*CA*CF**2*lmss*(1 + v**2))/((1 - v)**2*v) +
     -     (44*CA**2*CF*lmu*(1 + v**2))/(3.*(1 - v)**2*v) -
     -     (8*CA*CF**2*lms*lv*(1 + v**2))/((1 - v)**2*v) -
     -     (8*CA*CF**2*lmss*lv*(1 + v**2))/((1 - v)**2*v) -
     -     (40*CA*CF*Nf*(1 + v**2))/(9.*(1 - v)**2*v) +
     -     (8*CA*CF*l1v*Nf*(1 + v**2))/(3.*(1 - v)**2*v) -
     -     (8*CA*CF*lmu*Nf*(1 + v**2))/(3.*(1 - v)**2*v) -
     -     (2*CF*lv**2*(16 - 9*CA**2 + 20*v**2 - 11*CA**2*v**2))/
     -      ((1 - v)**2*v) + (4*CF*l1v*lv*
     -        (5 - 4*CA**2 + 9*v**2 - 6*CA**2*v**2))/((1 - v)**2*v) -
     -     (4*CF*l1v*(3 + 5*CA**2 - 3*CA**2*v + 15*v**2 + 
     -      2*CA**2*v**2))/
     -      (3.*(1 - v)**2*v) -
     -     (CF*lv*(5 - CA**2 - 8*v + 4*CA**2*v - 3*v**2 + 
     -      3*CA**2*v**2))/
     -      ((1 - v)**2*v) + (2*CF*l1v**2*
     -        (5 + 2*CA**2 - 3*v**2 + 4*CA**2*v**2))/((1 - v)**2*v) +
     -     (CF*(225 + 115*CA**2 + 42*Pi**2 + 12*CA**2*Pi**2 + 
     -      225*v**2 +
     -      115*CA**2*v**2 - 30*Pi**2*v**2 + 48*CA**2*Pi**2*v**2))/
     -      (9.*(1 - v)**2*v))

      ELSE IF (J0.EQ.2) THEN
      AVDEL = 0.D0

      ELSE IF (J0.EQ.3) THEN
      AVDEL = ((-12*CA*CF**2*lms*(1 + v**2))/((1 - v)**2*v) +
     -     (8*CA*CF**2*l1v*lms*(1 + v**2))/((1 - v)**2*v) -
     -     (6*CA*CF**2*lmss*(1 + v**2))/((1 - v)**2*v) +
     -     (44*CA**2*CF*lmu*(1 + v**2))/(3.*(1 - v)**2*v) -
     -     (8*CA*CF**2*lms*lv*(1 + v**2))/((1 - v)**2*v) -
     -     (8*CA*CF**2*lmss*lv*(1 + v**2))/((1 - v)**2*v) -
     -     (40*CA*CF*Nf*(1 + v**2))/(9.*(1 - v)**2*v) +
     -     (8*CA*CF*l1v*Nf*(1 + v**2))/(3.*(1 - v)**2*v) -
     -     (8*CA*CF*lmu*Nf*(1 + v**2))/(3.*(1 - v)**2*v) +
     -     (CF*lv*(11 - 3*CA**2 - 8*v + 3*v**2 - 3*CA**2*v**2))/
     -      ((1 - v)**2*v) - (2*CF*l1v**2*
     -        (3 - 4*CA**2 - 5*v**2 - 2*CA**2*v**2))/((1 - v)**2*v) +
     -     (4*CF*lv**2*(6 + CA**2 + 8*v**2 + CA**2*v**2))/
     -     ((1 - v)**2*v) -
     -     (4*CF*l1v*lv*(7 + CA**2 + 11*v**2 + CA**2*v**2))/
     -     ((1 - v)**2*v) -
     -     (4*CF*l1v*(15 + 2*CA**2 - 3*CA**2*v + 3*v**2 + 
     -      5*CA**2*v**2))/
     -      (3.*(1 - v)**2*v) +
     -     (CF*(225 + 115*CA**2 - 30*Pi**2 + 30*CA**2*Pi**2 + 
     -      225*v**2 +
     -      115*CA**2*v**2 + 42*Pi**2*v**2 + 30*CA**2*Pi**2*v**2))/
     -      (9.*(1 - v)**2*v))

      ELSE IF (J0.EQ.4) THEN
      AVDEL = 0.D0

      ELSE IF (J0.EQ.5) THEN
      AVDEL = (-4*(2 - CA**2)*CF*l1v - (12*CA*CF**2*lms*
     -     (1 - 2*v + 2*v**2))/v +
     -     (8*CA*CF**2*l1v*lms*(1 - 2*v + 2*v**2))/v -
     -     (6*CA*CF**2*lmss*(1 - 2*v + 2*v**2))/v +
     -     (44*CA**2*CF*lmu*(1 - 2*v + 2*v**2))/(3.*v) -
     -     (4*(7 - CA**2)*CF*l1v*lv*(1 - 2*v + 2*v**2))/v -
     -     (8*CA*CF**2*lms*lv*(1 - 2*v + 2*v**2))/v -
     -     (8*CA*CF**2*lmss*lv*(1 - 2*v + 2*v**2))/v -
     -     (40*CA*CF*Nf*(1 - 2*v + 2*v**2))/(9.*v) -
     -     (8*CA*CF*lmu*Nf*(1 - 2*v + 2*v**2))/(3.*v) +
     -     (CF*(225 + 115*CA**2 - 30*Pi**2 - 6*CA**2*Pi**2)*
     -        (1 - 2*v + 2*v**2))/(9.*v) +
     -     (CF*lv*(11 - 3*CA**2 - 14*v + 6*CA**2*v + 6*v**2 - 
     -      6*CA**2*v**2))/
     -      v - (2*CF*l1v**2*(3 - 2*CA**2 - 6*v + 4*CA**2*v + 2*v**2 -
     -          2*CA**2*v**2))/v +
     -     (4*CF*lv**2*(6 + CA**2 - 12*v - 2*CA**2*v + 14*v**2 +
     -          2*CA**2*v**2))/v)

      ELSE IF (J0.EQ.6) THEN
      AVDEL = ((-8*CF*l1v*Nf*(1 - v - CA*v - CA*v**3))/
     -     (3.*(1 - v)**2*v**2) +
     -     (8*CF*lv*Nf*(2*CA - v - 4*CA*v + 3*CA*v**2 - CA*v**3))/
     -      (3.*(1 - v)*v**3) +
     -     (4*CF*l1v*(6 - 6*CA + 2*CA**2 + 3*CA**3 + 
     -      9*CA*v - 2*CA**2*v -
     -          11*CA**3*v - 6*v**2 - 6*CA*v**2 + 
     -      6*CA**3*v**2 - 15*CA*v**3 -
     -          2*CA**3*v**3))/(3.*CA*(1 - v)**2*v**2) -
     -     (24*CF**2*lms*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v**3) +
     -     (16*CF**2*l1v*lms*(CA - v - 3*CA*v + v**2 + 4*CA*v**2 -
     -          2*CA*v**3 + CA*v**4))/((1 - v)**2*v**3) -
     -     (12*CF**2*lmss*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v**3) +
     -     (88*CA*CF*lmu*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/(3.*(1 - v)**2*v**3) -
     -     (16*CF**2*lms*lv*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/((1 - v)**2*v**3) -
     -     (16*CF**2*lmss*lv*(CA - v - 3*CA*v + v**2 + 4*CA*v**2 -
     -          2*CA*v**3 + CA*v**4))/((1 - v)**2*v**3) -
     -     (80*CF*Nf*(CA - v - 3*CA*v + v**2 + 4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/(9.*(1 - v)**2*v**3) -
     -     (16*CF*lmu*Nf*(CA - v - 3*CA*v + v**2 + 
     -      4*CA*v**2 - 2*CA*v**3 +
     -          CA*v**4))/(3.*(1 - v)**2*v**3) -
     -     (2*CF*lv**2*(16*CA - 8*CA**3 - 17*v - 56*CA*v + 6*CA**2*v +
     -          26*CA**3*v + 19*v**2 + 92*CA*v**2 - 6*CA**2*v**2 -
     -          42*CA**3*v**2 - 4*v**3 - 48*CA*v**3 + 20*CA**3*v**3 +
     -          2*v**4 + 32*CA*v**4 - 16*CA**3*v**4))/
     -     (CA*(1 - v)**2*v**3) +
     -     (4*CF*l1v*lv*(6*CA - 2*CA**3 - 5*v - 22*CA*v + 2*CA**2*v +
     -          8*CA**3*v + 7*v**2 + 36*CA*v**2 - 2*CA**2*v**2 -
     -          16*CA**3*v**2 - 4*v**3 - 20*CA*v**3 + 
     -      8*CA**3*v**3 + 2*v**4 +
     -          14*CA*v**4 - 8*CA**3*v**4))/(CA*(1 - v)**2*v**3) -
     -     (2*CF*l1v**2*(2*CA - 2*CA**3 + v - 10*CA*v + 4*CA**2*v +
     -          8*CA**3*v + v**2 + 12*CA*v**2 - 4*CA**2*v**2 -
     -          14*CA**3*v**2 - 4*v**3 - 12*CA*v**3 + 
     -      8*CA**3*v**3 + 2*v**4 +
     -          6*CA*v**4 - 6*CA**3*v**4))/(CA*(1 - v)**2*v**3) -
     -     (2*CF*lv*(27*CA + 17*CA**3 - 15*v - 105*CA*v - 13*CA**2*v -
     -          45*CA**3*v + 27*v**2 + 162*CA*v**2 + 13*CA**2*v**2 +
     -          46*CA**3*v**2 - 12*v**3 - 114*CA*v**3 - 22*CA**3*v**3 +
     -          21*CA*v**4 + 13*CA**3*v**4))/(3.*CA*(1 - v)**2*v**3) +
     -     (2*CF*(225*CA + 115*CA**3 + 6*CA*Pi**2 + 
     -      30*CA**3*Pi**2 - 225*v -
     -          675*CA*v - 115*CA**2*v - 345*CA**3*v - 15*Pi**2*v +
     -          18*CA*Pi**2*v - 30*CA**2*Pi**2*v - 108*CA**3*Pi**2*v +
     -          225*v**2 + 900*CA*v**2 + 115*CA**2*v**2 + 
     -      460*CA**3*v**2 -
     -          3*Pi**2*v**2 - 48*CA*Pi**2*v**2 + 30*CA**2*Pi**2*v**2 +
     -          156*CA**3*Pi**2*v**2 - 450*CA*v**3 - 230*CA**3*v**3 +
     -          36*Pi**2*v**3 + 60*CA*Pi**2*v**3 - 
     -      96*CA**3*Pi**2*v**3 +
     -          225*CA*v**4 + 115*CA**3*v**4 - 18*Pi**2*v**4 -
     -          30*CA*Pi**2*v**4 + 48*CA**3*Pi**2*v**4))/
     -      (9.*CA*(1 - v)**2*v**3))


      ELSE IF (J0.EQ.7) THEN
      AVDEL = 0.D0

      ELSE IF (J0.EQ.8) THEN
      AVDEL = 0.D0

      ELSE IF (J0.EQ.9) THEN
      AVDEL = 0.D0

      ELSE IF (J0.EQ.10) THEN
      AVDEL = 0.D0

      ELSE IF (J0.EQ.11) THEN
      AVDEL = ((8*CF*l1v*Nf*(CA + v**2 + CA*v**2 - v**3))/
     -     (3.*(1 - v)**2*v) -
     -     (4*CF*l1v*(15*CA + 2*CA**3 + 6*v + 6*CA*v - 6*CA**3*v -
     -          9*CA*v**2 + 2*CA**2*v**2 + 11*CA**3*v**2 - 6*v**3 +
     -          6*CA*v**3 - 2*CA**2*v**3 - 3*CA**3*v**3))/
     -     (3.*CA*(1 - v)**2*v)
     -       - (24*CF**2*lms*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 -
     -          3*CA*v**3 + CA*v**4))/((1 - v)**2*v) +
     -     (16*CF**2*l1v*lms*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 -
     -          3*CA*v**3 + CA*v**4))/((1 - v)**2*v) -
     -     (12*CF**2*lmss*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 -
     -          3*CA*v**3 + CA*v**4))/((1 - v)**2*v) +
     -     (88*CA*CF*lmu*(CA - 2*CA*v + v**2 + 4*CA*v**2 - 
     -      v**3 - 3*CA*v**3 +
     -          CA*v**4))/(3.*(1 - v)**2*v) -
     -     (16*CF**2*lms*lv*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 -
     -          3*CA*v**3 + CA*v**4))/((1 - v)**2*v) -
     -     (16*CF**2*lmss*lv*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 -
     -          3*CA*v**3 + CA*v**4))/((1 - v)**2*v) -
     -     (80*CF*Nf*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 - 
     -      3*CA*v**3 +
     -          CA*v**4))/(9.*(1 - v)**2*v) -
     -     (16*CF*lmu*Nf*(CA - 2*CA*v + v**2 + 4*CA*v**2 - 
     -      v**3 - 3*CA*v**3 +
     -          CA*v**4))/(3.*(1 - v)**2*v) +
     -     (2*CF*lv*(11*CA - 3*CA**3 - 22*CA*v + 6*CA**3*v + 3*v**2 +
     -          24*CA*v**2 - 3*CA**2*v**2 - 12*CA**3*v**2 - 3*v**3 -
     -          13*CA*v**3 + 3*CA**2*v**3 + 9*CA**3*v**3 + 3*CA*v**4 -
     -          3*CA**3*v**4))/(CA*(1 - v)**2*v) -
     -     (2*CF*l1v**2*(2 + 6*CA - 6*CA**3 - 4*v - 12*CA*v + 
     -      8*CA**3*v +
     -          v**2 + 12*CA*v**2 - 4*CA**2*v**2 - 
     -      14*CA**3*v**2 + v**3 -
     -          10*CA*v**3 + 4*CA**2*v**3 + 8*CA**3*v**3 + 2*CA*v**4 -
     -          2*CA**3*v**4))/(CA*(1 - v)**2*v) -
     -     (8*CF*l1v*lv*(7*CA - 14*CA*v + 2*CA**3*v + 
     -      4*v**2 + 30*CA*v**2 +
     -          4*CA**2*v**2 - 3*CA**3*v**2 - 4*v**3 - 21*CA*v**3 -
     -         4*CA**2*v**3 + 3*CA**3*v**3 + 7*CA*v**4 - CA**3*v**4))/
     -      (CA*(1 - v)**2*v) +
     -     (8*CF*lv**2*(6*CA + CA**3 - 12*CA*v - 2*CA**3*v + 3*v**2 +
     -          26*CA*v**2 + 5*CA**2*v**2 + 4*CA**3*v**2 - 3*v**3 -
     -          20*CA*v**3 - 5*CA**2*v**3 - 3*CA**3*v**3 + 7*CA*v**4 +
     -          CA**3*v**4))/(CA*(1 - v)**2*v) +
     -     (2*CF*(225*CA + 115*CA**3 - 30*CA*Pi**2 + 12*CA**3*Pi**2 -
     -          450*CA*v - 230*CA**3*v + 60*CA*Pi**2*v + 
     -      12*CA**3*Pi**2*v +
     -          225*v**2 + 900*CA*v**2 + 115*CA**2*v**2 + 
     -      460*CA**3*v**2 -
     -          12*Pi**2*v**2 - 84*CA*Pi**2*v**2 + 
     -      12*CA**2*Pi**2*v**2 -
     -          6*CA**3*Pi**2*v**2 - 225*v**3 - 675*CA*v**3 -
     -          115*CA**2*v**3 - 345*CA**3*v**3 + 12*Pi**2*v**3 +
     -          90*CA*Pi**2*v**3 - 12*CA**2*Pi**2*v**3 +
     -          18*CA**3*Pi**2*v**3 + 225*CA*v**4 + 115*CA**3*v**4 -
     -          30*CA*Pi**2*v**4 - 6*CA**3*Pi**2*v**4))/
     -      (9.*CA*(1 - v)**2*v))

      ELSE IF (J0.EQ.12) THEN
      AVDEL = ((2*CF*l1v*(2 + 2*CA**2 + v - 5*CA**2*v)*(1 - CA**2*v))/
     -      (CA*(1 - v)*v) + (16*CA*CF**2*l1v*lms*(1 - 2*v + 2*v**2)*
     -        (CF - CA*v + CA*v**2))/((1 - v)*v**2) -
     -     (44*CA**2*CF*lmss*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      (3.*(1 - v)*v**2) +
     -     (88*CA**2*CF*lmu*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      (3.*(1 - v)*v**2) -
     -     (16*CA*CF**2*lms*lv*(1 - 2*v + 2*v**2)*
     -     (CF - CA*v + CA*v**2))/
     -      ((1 - v)*v**2) - (16*CA**2*CF*lmss*lv*(1 - 2*v + 2*v**2)*
     -        (CF - CA*v + CA*v**2))/((1 - v)*v**2) +
     -     (8*CA*CF*lmss*Nf*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      (3.*(1 - v)*v**2) -
     -     (16*CA*CF*lmu*Nf*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      (3.*(1 - v)*v**2) +
     -     (12*CF**2*lms*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/((1 - v)*v**2) +
     -     (20*CF*Nf*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/
     -     (9.*(1 - v)*v**2) -
     -     (4*CF*lv*Nf*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/
     -     (3.*(1 - v)*v**2) -
     -     (4*CF*l1v*lv*(1 - 2*v + 2*v**2)*
     -        (1 - 4*CA**2 + 3*CA**4 + 2*CA**2*v - 6*CA**4*v - 
     -     2*CA**2*v**2 +
     -          2*CA**4*v**2))/(CA*(1 - v)*v**2) +
     -     (2*CF*lv*(9 - 7*CA**2 - 2*CA**4 - 12*v + 23*CA**2*v + 
     -     11*CA**4*v +
     -          3*v**2 - 8*CA**2*v**2 - 49*CA**4*v**2 + 3*CA**2*v**3 +
     -          73*CA**4*v**3 - 44*CA**4*v**4))/(3.*CA*(1 - v)*v**2) +
     -     (2*CF*l1v**2*(3 + CA**4 - 4*v + 2*CA**2*v - 6*CA**4*v + 
     -     3*v**2 -
     -          7*CA**2*v**2 + 13*CA**4*v**2 + 7*CA**2*v**3 - 
     -     12*CA**4*v**3 -
     -          4*CA**2*v**4 + 4*CA**4*v**4))/(CA*(1 - v)*v**2) +
     -     (2*CF*lv**2*(2 - 12*CA**2 + 10*CA**4 - 2*v + 27*CA**2*v -
     -          40*CA**4*v + 3*v**2 - 30*CA**2*v**2 + 73*CA**4*v**2 +
     -          9*CA**2*v**3 - 68*CA**4*v**3 - 4*CA**2*v**4 + 
     -     28*CA**4*v**4))/
     -      (CA*(1 - v)*v**2) -
     -     (2*CF*(63 - 59*CA**2 - 4*CA**4 - 9*Pi**2 + 12*CA**2*Pi**2 -
     -          3*CA**4*Pi**2 - 126*v + 253*CA**2*v + 25*CA**4*v +
     -          18*Pi**2*v - 42*CA**2*Pi**2*v + 12*CA**4*Pi**2*v + 
     -     126*v**2 -
     -          541*CA**2*v**2 - 77*CA**4*v**2 - 18*Pi**2*v**2 +
     -          78*CA**2*Pi**2*v**2 - 24*CA**4*Pi**2*v**2 + 
     -     576*CA**2*v**3 +
     -          104*CA**4*v**3 - 72*CA**2*Pi**2*v**3 + 
     -     24*CA**4*Pi**2*v**3 -
     -          288*CA**2*v**4 - 52*CA**4*v**4 + 36*CA**2*Pi**2*v**4 -
     -          12*CA**4*Pi**2*v**4))/(9.*CA*(1 - v)*v**2))

      ELSE IF (J0.EQ.13) THEN
      AVDEL = (-((9 - 31*CA**2)*CF*lms*(1 + v**2)*
     -         (-2*CA*CF - 2*v + v**2 - CA**2*v**2))/
     -     (3.*CA*(1 - v)**2*v**2)
     -      + (20*CF*Nf*(1 + v**2)*(-2*CA*CF - 2*v + v**2 - 
     -     CA**2*v**2))/
     -      (9.*(1 - v)**2*v**2) -
     -     (4*CF*lms*Nf*(1 + v**2)*(-2*CA*CF - 2*v + v**2 - 
     -     CA**2*v**2))/
     -      (3.*(1 - v)**2*v**2) -
     -     (4*CF*lv*Nf*(1 + v**2)*(-2*CA*CF - 2*v + v**2 - 
     -     CA**2*v**2))/
     -      (3.*(1 - v)**2*v**2) +
     -     (8*CA*CF*l1v*lms*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -     CA**2*v**2))/
     -      ((1 - v)**2*v**2) -
     -     (6*CF**2*lmss*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -     CA**2*v**2))/
     -      ((1 - v)**2*v**2) +
     -     (44*CA*CF*lmu*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -     CA**2*v**2))/
     -      (3.*(1 - v)**2*v**2) -
     -     (8*CA*CF*lms*lv*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -     CA**2*v**2))/
     -      ((1 - v)**2*v**2) -
     -     (8*CF**2*lmss*lv*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -     CA**2*v**2))/
     -      ((1 - v)**2*v**2) -
     -     (8*CF*lmu*Nf*(1 + v**2)*(2*CA*CF + 2*v - v**2 + 
     -     CA**2*v**2))/
     -      (3.*(1 - v)**2*v**2) +
     -     (2*CF*l1v*(4 - CA**2 + CA**4 - 8*v - 10*CA**2*v + 
     -     10*CA**4*v +
     -          4*v**2 - CA**2*v**2 + CA**4*v**2))/(CA*(1 - v)**2*v) +
     -     (2*CF*lv*(9 - 7*CA**2 - 2*CA**4 - 24*v + 5*CA**2*v - 
     -     3*CA**4*v +
     -          21*v**2 + 19*CA**2*v**2 - 28*CA**4*v**2 - 6*v**3 -
     -          28*CA**2*v**3 + 11*CA**2*v**4 - 11*CA**4*v**4))/
     -      (3.*CA*(1 - v)**2*v**2) +
     -     (4*CF*l1v*lv*(4*CA**2 - 4*CA**4 + 2*v - 9*CA**2*v - 5*v**2 +
     -          5*CA**2*v**2 - 5*CA**4*v**2 + 4*v**3 - 6*CA**2*v**3 -
     -          2*CA**4*v**3 - v**4 + 2*CA**2*v**4 - CA**4*v**4))/
     -      (CA*(1 - v)**2*v**2) +
     -     (2*CF*l1v**2*(3 + CA**4 - 10*v - CA**2*v + 2*CA**4*v + 
     -     14*v**2 +
     -          2*CA**2*v**2 + 2*CA**4*v**2 - 10*v**3 - CA**2*v**3 +
     -          2*CA**4*v**3 + 3*v**4 + CA**4*v**4))/
     -     (CA*(1 - v)**2*v**2) +
     -     (2*CF*lv**2*(2 - 12*CA**2 + 10*CA**4 - 6*v + 
     -     21*CA**2*v + 9*v**2 -
     -          21*CA**2*v**2 + 13*CA**4*v**2 - 8*v**3 + 
     -     18*CA**2*v**3 +
     -         2*CA**4*v**3 + 3*v**4 - 10*CA**2*v**4 + 3*CA**4*v**4))/
     -      (CA*(1 - v)**2*v**2) -
     -     (2*CF*(63 - 59*CA**2 - 4*CA**4 - 9*Pi**2 + 12*CA**2*Pi**2 -
     -          3*CA**4*Pi**2 - 126*v - 17*CA**2*v - 
     -     9*CA**4*v + 36*Pi**2*v -
     -          15*CA**2*Pi**2*v + 126*v**2 - 136*CA**2*v**2 -
     -          26*CA**4*v**2 - 63*Pi**2*v**2 - 3*CA**2*Pi**2*v**2 -
     -          15*CA**4*Pi**2*v**2 - 126*v**3 - 17*CA**2*v**3 -
     -          9*CA**4*v**3 + 54*Pi**2*v**3 + 12*CA**2*Pi**2*v**3 -
     -          18*CA**4*Pi**2*v**3 + 63*v**4 - 59*CA**2*v**4 -
     -          4*CA**4*v**4 - 18*Pi**2*v**4 - 6*CA**2*Pi**2*v**4 -
     -          12*CA**4*Pi**2*v**4))/(9.*CA*(1 - v)**2*v**2))

      ELSE IF (J0.EQ.14) THEN
      AVDEL = ((-2*CF*l1v*(CA**2 - v)*(1 - 5*CA**2 + 2*v + 2*CA**2*v))/
     -      (CA*(1 - v)*v**2) +
     -     ((9 - 31*CA**2)*CF*lms*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (3.*CA*(1 - v)*v**3)
     -      + (8*CA*CF*l1v*lms*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     ((1 - v)*v**3) -
     -     (22*CA*CF*lmss*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (3.*(1 - v)*v**3) +
     -     (44*CA*CF*lmu*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (3.*(1 - v)*v**3) -
     -     (8*CA*CF*lms*lv*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     ((1 - v)*v**3) -
     -     (8*CA*CF*lmss*lv*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     ((1 - v)*v**3) +
     -     (4*CF*lms*Nf*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (3.*(1 - v)*v**3) +
     -     (4*CF*lmss*Nf*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (3.*(1 - v)*v**3) -
     -     (8*CF*lmu*Nf*(2 - 2*v + v**2)*
     -        (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (3.*(1 - v)*v**3) +
     -     (2*CF*l1v**2*(8*CA**4 - CA**2*v - 20*CA**4*v + v**2 - 
     -     CA**2*v**2 +
     -          21*CA**4*v**2 - 2*v**3 - 10*CA**4*v**3 + 
     -     2*v**4 + 2*CA**4*v**4
     -          ))/(CA*(1 - v)*v**3) -
     -     (CF*lv*(12*CA**2 - 12*CA**4 - 24*CA**2*v + 
     -     24*CA**4*v - 2*v**2 -
     -          4*CA**2*v**2 - 2*CA**4*v**2 + 2*v**3 + 16*CA**2*v**3 -
     -          10*CA**4*v**3 + 3*v**4 - 6*CA**2*v**4 + 3*CA**4*v**4))/
     -      (CA*(1 - v)*v**3) -
     -     (4*CF*l1v*lv*(16*CA**4 - CA**2*v - 36*CA**4*v + v**2 +
     -          3*CA**2*v**2 + 37*CA**4*v**2 - 2*v**3 - 4*CA**2*v**3 -
     -          18*CA**4*v**3 + 2*v**4 + 2*CA**2*v**4 + 4*CA**4*v**4))/
     -      (CA*(1 - v)*v**3) +
     -     (2*CF*lv**2*(48*CA**4 - 96*CA**4*v + 2*v**2 - 5*CA**2*v**2 +
     -          94*CA**4*v**2 - 2*v**3 + 5*CA**2*v**3 - 46*CA**4*v**3 +
     -          3*v**4 - 2*CA**2*v**4 + 11*CA**4*v**4))/
     -     (CA*(1 - v)*v**3) +
     -     (CF*(108*CA**2 - 60*CA**4 + 8*CA**2*Pi**2 + 40*CA**4*Pi**2 -
     -          216*CA**2*v + 120*CA**4*v - 22*CA**2*Pi**2*v -
     -          104*CA**4*Pi**2*v - 42*v**2 + 240*CA**2*v**2 -
     -          138*CA**4*v**2 + 2*Pi**2*v**2 + 14*CA**2*Pi**2*v**2 +
     -          110*CA**4*Pi**2*v**2 + 42*v**3 - 132*CA**2*v**3 +
     -          78*CA**4*v**3 - 8*Pi**2*v**3 - 12*CA**2*Pi**2*v**3 -
     -          52*CA**4*Pi**2*v**3 - 21*v**4 + 42*CA**2*v**4 -
     -          21*CA**4*v**4 + 10*Pi**2*v**4 + 4*CA**2*Pi**2*v**4 +
     -          10*CA**4*Pi**2*v**4))/(3.*CA*(1 - v)*v**3))

      ELSE IF (J0.EQ.15) THEN
      AVDEL = ((-16*CA**3*lv**2*Nf*(1 - v - v**2))/(3.*(1 - v)*v) +
     -     (16*CA**3*l1v**2*Nf*(1 - 3*v + v**2))/(3.*v**2) -
     -     (1408*CA**3*lms*(1 - v + v**2)**3)/(3.*(1 - v)**2*v**3) +
     -     (256*CA**3*l1v*lms*(1 - v + v**2)**3)/((1 - v)**2*v**3) -
     -     (704*CA**3*lmss*(1 - v + v**2)**3)/(3.*(1 - v)**2*v**3) +
     -     (1408*CA**3*lmu*(1 - v + v**2)**3)/(3.*(1 - v)**2*v**3) -
     -     (256*CA**3*lms*lv*(1 - v + v**2)**3)/((1 - v)**2*v**3) -
     -     (256*CA**3*lmss*lv*(1 - v + v**2)**3)/((1 - v)**2*v**3) +
     -     (256*CA**3*lms*Nf*(1 - v + v**2)**3)/(9.*(1 - v)**2*v**3) +
     -     (128*CA**3*lmss*Nf*(1 - v + v**2)**3)/(9.*(1 - v)**2*v**3) -
     -     (256*CA**3*lmu*Nf*(1 - v + v**2)**3)/(9.*(1 - v)**2*v**3) -
     -     (32*CA**3*l1v*lv*Nf*(1 - 2*v + 2*v**2))/(3.*(1 - v)*v**2) -
     -     (32*CA**3*l1v*Nf*(1 - v + v**2)*(5 - 2*v + 5*v**2))/
     -      (9.*(1 - v)**2*v**2) +
     -     (64*CA**3*l1v*(1 - v + v**2)*(7 + 8*v + 7*v**2))/
     -      (3.*(1 - v)**2*v**2) +
     -     (64*CA**3*lv*(1 - v + v**2)*
     -    (11 - 22*v - 4*v**2 + 15*v**3 - 11*v**4))/
     -     (3.*(1 - v)**2*v**3)
     -      - (32*CA**3*lv*Nf*(1 - v + v**2)*
     -     (4 - 8*v + v**2 + 3*v**3 - 4*v**4))/(9.*(1 - v)**2*v**3) +
     -(64*CA**3*l1v**2*(2 - 7*v + 14*v**2 - 16*v**3 + 14*v**4 - 7*v**5+
     -          2*v**6))/((1 - v)**2*v**3) -
     -     (64*CA**3*l1v*lv*(8 - 26*v + 47*v**2 - 50*v**3 + 37*v**4 -
     -          16*v**5 + 4*v**6))/((1 - v)**2*v**3) +
     -     (64*CA**3*lv**2*(12 - 36*v + 66*v**2 - 72*v**3 + 57*v**4 -
     -          27*v**5 + 8*v**6))/((1 - v)**2*v**3) +
     -  (16*CA**3*Nf*(40 - 120*v + 9*Pi**2*v + 294*v**2 -27*Pi**2*v**2-
     -          388*v**3 + 36*Pi**2*v**3 + 294*v**4 - 18*Pi**2*v**4 -
     -          120*v**5 + 40*v**6))/(27.*(1 - v)**2*v**3) -
     -     (32*CA**3*(134 - 24*Pi**2 - 402*v + 90*Pi**2*v + 831*v**2 -
     -          171*Pi**2*v**2 - 992*v**3 + 186*Pi**2*v**3 + 831*v**4 -
     -          153*Pi**2*v**4 - 402*v**5 + 72*Pi**2*v**5 + 134*v**6 -
     -          24*Pi**2*v**6))/(9.*(1 - v)**2*v**3))

      ELSE IF (J0.EQ.16) THEN
      AVDEL = ((2*CF*l1v*(2 + 2*CA**2 + v - 5*CA**2*v)*(1 - CA**2*v))/
     -(CA*(1 - v)*v) + (16*CA*CF*l1v*lv*(1 + CA - CA*v)*
     -     (1 - CA + CA*v)*
     -        (1 - 2*v + 2*v**2))/((1 - v)*v**2) +
     -  (16*CA**2*CF*l1v*lms*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      ((1 - v)*v**2) - (12*CA*CF**2*lmss*(1 - 2*v + 2*v**2)*
     -        (CF - CA*v + CA*v**2))/((1 - v)*v**2) +
     -     (88*CA**2*CF*lmu*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      (3.*(1 - v)*v**2) -
     -  (16*CA**2*CF*lms*lv*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      ((1 - v)*v**2) - (16*CA*CF**2*lmss*lv*(1 - 2*v + 2*v**2)*
     -        (CF - CA*v + CA*v**2))/((1 - v)*v**2) -
     -     (16*CA*CF*lmu*Nf*(1 - 2*v + 2*v**2)*(CF - CA*v + CA*v**2))/
     -      (3.*(1 - v)*v**2) +
     -     (44*CA*CF*lms*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/
     -     (3.*(1 - v)*v**2) -
     -     (8*CF*lms*Nf*(1 - 2*v + 2*v**2)*
     -        (-2*CA*CF + 2*CA**2*v - 2*CA**2*v**2))/
     -     (3.*(1 - v)*v**2) +
     -  (CF*lv*(3 - 6*CA**2 + 3*CA**4 - 2*v + 12*CA**2*v - 10*CA**4*v -
     -          4*v**2 + 10*CA**2*v**2 + 2*CA**4*v**2 - 22*CA**2*v**3 +
     -          14*CA**4*v**3 + 12*CA**2*v**4 - 12*CA**4*v**4))/
     -      (CA*(1 - v)*v**2) +
     - (2*CF*l1v**2*(2 + 2*CA**4 - 2*v - 10*CA**4*v + 
     -     v**2 - CA**2*v**2+
     -  21*CA**4*v**2 - CA**2*v**3 - 20*CA**4*v**3 + 8*CA**4*v**4))/
     -      (CA*(1 - v)*v**2) +
     -     (2*CF*lv**2*(2 - 12*CA**2 + 10*CA**4 - 2*v + 27*CA**2*v -
     -          40*CA**4*v + 3*v**2 - 30*CA**2*v**2 + 73*CA**4*v**2 +
     -   9*CA**2*v**3 - 68*CA**4*v**3 - 4*CA**2*v**4 + 28*CA**4*v**4))/
     -      (CA*(1 - v)*v**2) -
     -     (CF*(21 - 42*CA**2 + 21*CA**4 - 4*Pi**2 + 8*CA**2*Pi**2 -
     -     4*CA**4*Pi**2 - 42*v + 132*CA**2*v - 78*CA**4*v + 
     -     8*Pi**2*v -
     -          24*CA**2*Pi**2*v + 16*CA**4*Pi**2*v + 42*v**2 -
     -          240*CA**2*v**2 + 138*CA**4*v**2 - 8*Pi**2*v**2 +
     -     40*CA**2*Pi**2*v**2 - 32*CA**4*Pi**2*v**2 + 216*CA**2*v**3 -
     -     120*CA**4*v**3 - 32*CA**2*Pi**2*v**3 + 32*CA**4*Pi**2*v**3 -
     -          108*CA**2*v**4 + 60*CA**4*v**4 + 16*CA**2*Pi**2*v**4 -
     -          16*CA**4*Pi**2*v**4))/(3.*CA*(1 - v)*v**2))

      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION AVLO(W,V,S)
C     TERME EN (LOG(1-W)/(1-W)+
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      COMMON/ORDE/AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      M=DSQRT(Q2FAC)
      MP=DSQRT(Q2FRAG)  
C    TERME EN LOPLUS


      l1v = dlog(1.-v)
      lv = dlog(v)
      lmu = dlog(Q2MU/S)
      lms = dlog(M**2/S)
      lmss = dlog(MP**2/S)

      IF (J0.EQ.1) THEN
      AVLO = (40*CA*CF**2*(1 + v**2))/((-1 + v)**2*v)

      ELSE IF (J0.EQ.2) THEN
      AVLO = 0.D0

      ELSE IF (J0.EQ.3) THEN
      AVLO = (40*CA*CF**2*(1 + v**2))/((-1 + v)**2*v)

      ELSE IF (J0.EQ.4) THEN
      AVLO = 0.D0

      ELSE IF (J0.EQ.5) THEN
      AVLO = (40*CA*CF**2*(1 - 2*v + 2*v**2))/v

      ELSE IF (J0.EQ.6) THEN
      AVLO = (80*CF**2*(CA - v - 3*CA*v + v**2 + 4*CA*v**2 - 
     -     2*CA*v**3 +
     -       CA*v**4))/((-1 + v)**2*v**3)

      ELSE IF (J0.EQ.7) THEN
      AVLO = 0.D0

      ELSE IF (J0.EQ.8) THEN
      AVLO = 0.D0

      ELSE IF (J0.EQ.9) THEN
      AVLO = 0.D0

      ELSE IF (J0.EQ.10) THEN
      AVLO = 0.D0

      ELSE IF (J0.EQ.11) THEN
      AVLO = (80*CF**2*(CA - 2*CA*v + v**2 + 4*CA*v**2 - v**3 - 
     -     3*CA*v**3 +
     -       CA*v**4))/((-1 + v)**2*v)

      ELSE IF (J0.EQ.12) THEN
      AVLO = (-16*(-2 + 3*CA**2)*CF*(1 - 2*v + 2*v**2)*
     -     (CF - CA*v + CA*v**2))/
     -   ((-1 + v)*v**2)

      ELSE IF (J0.EQ.13) THEN
      AVLO = (8*(-2 + 3*CA**2)*CF*(1 + v**2)*
     -     (2*CA*CF + 2*v - v**2 + CA**2*v**2))/(CA*(-1 + v)**2*v**2)

      ELSE IF (J0.EQ.14) THEN
      AVLO = (-4*(-1 + 3*CA)*(1 + 3*CA)*CF*(2 - 2*v + v**2)*
     -     (2*CA**2 - 2*CA**2*v - v**2 + CA**2*v**2))/
     -     (CA*(-1 + v)*v**3)

      ELSE IF (J0.EQ.15) THEN
      AVLO = (1280*CA**3*(1 - v + v**2)**3)/((-1 + v)**2*v**3)

      ELSE IF (J0.EQ.16) THEN
      AVLO =  (-8*(-1 + 3*CA)*(1 + 3*CA)*CF*(1 - 2*v + 2*v**2)*
     -     (CF - CA*v + CA*v**2))/((-1 + v)*v**2)

      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION AVGO(W,V)
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      COMMON/ORDE/AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0

      AVGO = 0.D0

      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/EDF/J0
*
cc      CALL PRECALC(V,W,S)
*
      IF (J0.EQ.1) THEN
      STRUV=STRUV1(W,V,X3,S)
      ELSE IF (J0.EQ.2) THEN
      STRUV=STRUV2(W,V,X3,S)
      ELSE IF (J0.EQ.3) THEN
      STRUV=STRUV3(W,V,X3,S)
      ELSE IF (J0.EQ.4) THEN
      STRUV=STRUV4(W,V,X3,S)
      ELSE IF (J0.EQ.5) THEN
      STRUV=STRUV5(W,V,X3,S)
      ELSE IF (J0.EQ.6) THEN
      STRUV=STRUV6(W,V,X3,S)
      ELSE IF (J0.EQ.7) THEN
      STRUV=STRUV7(W,V,X3,S)
      ELSE IF (J0.EQ.8) THEN
      STRUV=STRUV8(W,V,X3,S)
      ELSE IF (J0.EQ.9) THEN
      STRUV=STRUV9(W,V,X3,S)
      ELSE IF (J0.EQ.10) THEN
      STRUV=STRUV10(W,V,X3,S)
      ELSE IF (J0.EQ.11) THEN
      STRUV=STRUV11(W,V,X3,S)
      ELSE IF (J0.EQ.12) THEN
      STRUV=STRUV12(W,V,X3,S)
      ELSE IF (J0.EQ.13) THEN
      STRUV=STRUV13(W,V,X3,S)
      ELSE IF (J0.EQ.14) THEN
      STRUV=STRUV14(W,V,X3,S)
      ELSE IF (J0.EQ.15) THEN
      STRUV=STRUV15(W,V,X3,S)
      ELSE IF (J0.EQ.16) THEN
      STRUV=STRUV16(W,V,X3,S)
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV1(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON /SCALES/ Q2FAC,Q2MU,Q2FRAG
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (-4*CF*l1v*(2 - 4*CA2 - 4*v + CA2*v - CA2*v2 + 4*v*w -
     &   3*CA2*v*w - 4*v2*w + 3*CA2*v2*w + 2*v2*w2))/
     &   ((1 - v)**2*(1 - v*w)) 
      part2 = (4*CA*CF**2*lmss*(4 - 5*v + 2*v2 - v3 + 9*v*w - 8*v2*w +
     &   3*v3*w + 6*v2*w2 - 4*v3*w2 + 2*v3*w3))/
     &   ((1 - v)**2*(1 - v + v*w)) 
      part3 = (4*CF*lvw*(2*CA*CF - v2 + CA2*v2 + 5*w + 2*CA2*w - v*w -
     &   4*CA2*v*w - v2*w + 3*CA2*v2*w + v3*w - CA2*v3*w + v*w2 + 
     &   4*CA2*v*w2 - 2*v2*w2 - 3*CA2*v2*w2 - 3*v3*w2 + 3*CA2*v3*w2 + 
     &   4*CA2*v2*w3 + 4*v3*w3 - 4*CA2*v3*w3 - 2*v3*w4 + 2*CA2*v3*w4))/
     &   ((1 - v)**2*v*w) 
      part4 = -(4*CF*l1vw*(2*CA2 - 4*CA*CF + 3*v - 6*CA2*v + 12*CA*CF*v-
     &   7*v2 + 7*CA2*v2 - 14*CA*CF*v2 + 5*v3 - 4*CA2*v3 + 8*CA*CF*v3 -
     &   v4 + CA2*v4 - 2*CA*CF*v4 - 3*v*w + 6*CA2*v*w - 12*CA*CF*v*w +
     &   10*v2*w - 13*CA2*v2*w + 28*CA*CF*v2*w - 11*v3*w + 
     &   11*CA2*v3*w - 24*CA*CF*v3*w + 4*v4*w - 4*CA2*v4*w + 
     &   8*CA*CF*v4*w - 7*v2*w2 + 8*CA2*v2*w2 - 18*CA*CF*v2*w2 +
     &   12*v3*w2 - 13*CA2*v3*w2 + 28*CA*CF*v3*w2 - 7*v4*w2 + 
     &   7*CA2*v4*w2 - 14*CA*CF*v4*w2 - 6*v3*w3 + 6*CA2*v3*w3 - 
     &   12*CA*CF*v3*w3 + 6*v4*w3 - 6*CA2*v4*w3 + 12*CA*CF*v4*w3 -
     &   2*v4*w4 + 2*CA2*v4*w4 - 4*CA*CF*v4*w4))/
     &   ((1 - v)**2*v*(1 - v + v*w)) 
      part5 = -(2*CF*lms*(4*CA*CF + 4*CA*CF*v2 - 2*CA*CF*w + 2*v*w -
     &   18*CA*CF*v*w - 4*CA*CF*v2*w + 2*v3*w - 8*CA*CF*v3*w +
     &   6*CA*CF*v*w2 - 5*v2*w2 + CA2*v2*w2 + 18*CA*CF*v2*w2 - 
     &   4*v3*w2 + 16*CA*CF*v3*w2 - v4*w2 + CA2*v4*w2 + 8*CA*CF*v4*w2 -
     &   6*CA*CF*v2*w3 + 8*v3*w3 - 2*CA2*v3*w3 - 6*CA*CF*v3*w3 + 
     &   2*v4*w3 - 2*CA2*v4*w3 - 16*CA*CF*v4*w3 + 2*v5*w3 - 
     &   4*CA*CF*v5*w3 + 2*CA*CF*v3*w4 - 7*v4*w4 + 3*CA2*v4*w4 +
     &   2*CA*CF*v4*w4 - 4*v5*w4 + 4*CA*CF*v5*w4 - v6*w4 + CA2*v6*w4 +
     &   6*v5*w5 - 2*CA2*v5*w5 + 2*v6*w5 - 2*CA2*v6*w5 - 2*v6*w6 +
     &   2*CA2*v6*w6))/((1 - v)**2*v*w*(1 - v*w)**3) -
     &   (2*CF*l1w*(2 - 2*CA2 - 2*v + 2*CA2*v + 2*v2 - 2*CA2*v2 -
     &   2*v3 + 2*CA2*v3 - w + CA2*w - 24*v*w + 18*CA2*v*w +
     &   34*v2*w - 22*CA2*v2*w - 12*v3*w + 6*CA2*v3*w + 7*v4*w - 
     &   7*CA2*v4*w + 2*v*w2 - 2*CA2*v*w2 + 7*v2*w2 - 3*CA2*v2*w2 - 
     &   27*v3*w2 + 17*CA2*v3*w2 - 2*v4*w2 + 4*CA2*v4*w2 - 8*v5*w2 + 
     &   12*CA2*v5*w2 + 39*v3*w3 - 31*CA2*v3*w3 - 16*v4*w3 + 
     &   14*CA2*v4*w3 + 22*v5*w3 - 28*CA2*v5*w3 + 7*v6*w3 - 
     &   7*CA2*v6*w3 - 2*v3*w4 + 2*CA2*v3*w4 - 19*v4*w4 + 
     &   11*CA2*v4*w4 + 11*v5*w4 + 3*CA2*v5*w4 - 24*v6*w4 + 
     &   18*CA2*v6*w4 - 2*v7*w4 + 2*CA2*v7*w4 + v4*w5 - CA2*v4*w5 -
     &   13*v5*w5 + 7*CA2*v5*w5 + 18*v6*w5 - 12*CA2*v6*w5 + 6*v7*w5 - 
     &   6*CA2*v7*w5 - 8*v7*w6 + 8*CA2*v7*w6 + 4*v7*w7 -
     &   4*CA2*v7*w7))/((1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w)) 
      part6 = -(2*CF*lv*(2 - 2*CA2 - 2*v + 2*CA2*v + 2*v2 - 2*CA2*v2 -
     &   2*v3 + 2*CA2*v3 - w + CA2*w - 44*v*w + 26*CA2*v*w + 54*v2*w -
     &   32*CA2*v2*w - 20*v3*w + 10*CA2*v3*w + 15*v4*w - 9*CA2*v4*w + 
     &   2*v*w2 - 2*CA2*v*w2 + 11*v2*w2 - 5*CA2*v2*w2 - 35*v3*w2 +
     &   23*CA2*v3*w2 - 10*v4*w2 + 4*CA2*v4*w2 - 24*v5*w2 + 
     &   16*CA2*v5*w2 + 71*v3*w3 - 45*CA2*v3*w3 - 40*v4*w3 + 
     &   26*CA2*v4*w3 + 62*v5*w3 - 40*CA2*v5*w3 + 15*v6*w3 -
     &   9*CA2*v6*w3 - 2*v3*w4 + 2*CA2*v3*w4 - 19*v4*w4 + 
     &   13*CA2*v4*w4 + 3*v5*w4 + CA2*v5*w4 - 48*v6*w4 + 26*CA2*v6*w4 -
     &   2*v7*w4 + 2*CA2*v7*w4 + v4*w5 - CA2*v4*w5 - 25*v5*w5 + 
     &   13*CA2*v5*w5 + 38*v6*w5 - 18*CA2*v6*w5 + 6*v7*w5 - 
     &   6*CA2*v7*w5 - 4*v6*w6 - 8*v7*w6 + 8*CA2*v7*w6 + 4*v7*w7 -
     &   4*CA2*v7*w7))/((1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w))  
      part7 = -(CF*(6 - 6*CA2 + 4*CA*CF - 14*v + 14*CA2*v - 12*CA*CF*v+
     &   14*v2 - 14*CA2*v2 + 12*CA*CF*v2 - 6*v3 + 6*CA2*v3 - 
     &   4*CA*CF*v3 + 4*w - 4*CA*CF*w - 25*v*w + 5*CA2*v*w + 
     &   4*CA*CF*v*w + 70*v2*w - 30*CA2*v2*w + 24*CA*CF*v2*w - 
     &   83*v3*w + 47*CA2*v3*w - 44*CA*CF*v3*w + 38*v4*w - 
     &   26*CA2*v4*w + 20*CA*CF*v4*w - 8*v*w2 + 8*CA*CF*v*w2 - 
     &   2*v2*w2 + 10*CA2*v2*w2 - 36*CA*CF*v2*w2 + v3*w2 -
     &   13*CA2*v3*w2 + 20*CA*CF*v3*w2 + 39*v4*w2 - 11*CA2*v4*w2 + 
     &   40*CA*CF*v4*w2 - 40*v5*w2 + 24*CA2*v5*w2 - 32*CA*CF*v5*w2 + 
     &   48*v3*w3 - 8*CA2*v3*w3 + 36*CA*CF*v3*w3 - 103*v4*w3 +
     &   15*CA2*v4*w3 - 64*CA*CF*v4*w3 + 65*v5*w3 - 5*CA2*v5*w3 + 
     &   8*CA*CF*v5*w3 + 6*v6*w3 - 18*CA2*v6*w3 + 20*CA*CF*v6*w3 + 
     &   8*v3*w4 - 8*CA*CF*v3*w4 - 20*v4*w4 + 12*CA2*v4*w4 +
     &   5*v5*w4 - 17*CA2*v5*w4 + 32*CA*CF*v5*w4 - 13*v6*w4 + 
     &   25*CA2*v6*w4 - 20*CA*CF*v6*w4 - 6*v7*w4 + 6*CA2*v7*w4 - 
     &   4*CA*CF*v7*w4 - 4*v4*w5 + 4*CA*CF*v4*w5 + 13*v5*w5 -
     &   CA2*v5*w5 - 8*CA*CF*v5*w5 + 9*v6*w5 - 17*CA2*v6*w5 + 
     &   10*v7*w5 - 10*CA2*v7*w5 + 4*CA*CF*v7*w5 - 8*v6*w6 + 
     &   8*CA2*v6*w6 - 8*v7*w6 + 8*CA2*v7*w6 + 4*v7*w7 -
     &   4*CA2*v7*w7))/((1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w)) 
      part8 = -(2*CF*lw*(3 - 3*CA2 + 2*CA*CF - 3*v + 3*CA2*v - 
     &   2*CA*CF*v + 3*v2 - 3*CA2*v2 + 2*CA*CF*v2 - 3*v3 + 3*CA2*v3 -
     &   2*CA*CF*v3 + 14*w - 8*CA2*w + 8*CA*CF*w - 48*v*w +
     &   30*CA2*v*w - 28*CA*CF*v*w + 57*v2*w - 37*CA2*v2*w +
     &   34*CA*CF*v2*w - 36*v3*w + 22*CA2*v3*w - 20*CA*CF*v3*w + 
     &   19*v4*w - 13*CA2*v4*w + 10*CA*CF*v4*w + w2 - CA2*w2 + 
     &   2*CA*CF*w2 - v*w2 + CA2*v*w2 - 2*CA*CF*v*w2 + 23*v2*w2 -
     &   13*CA2*v2*w2 + 10*CA*CF*v2*w2 - 50*v3*w2 +
     &   32*CA2*v3*w2 - 24*CA*CF*v3*w2 + 33*v4*w2 - 17*CA2*v4*w2 + 
     &   14*CA*CF*v4*w2 - 28*v5*w2 + 16*CA2*v5*w2 - 12*CA*CF*v5*w2 - 
     &   2*v*w3 + 2*CA2*v*w3 - 4*CA*CF*v*w3 - 3*v2*w3 +
     &   CA2*v2*w3 + 2*CA*CF*v2*w3 + 10*v3*w3 - 8*CA2*v3*w3 + 
     &   8*CA*CF*v3*w3 + 2*v4*w3 - 6*CA2*v4*w3 + 8*CA*CF*v4*w3 + 
     &   16*v5*w3 - 6*CA2*v5*w3 - 4*CA*CF*v5*w3 + 15*v6*w3 -
     &   9*CA2*v6*w3 + 10*CA*CF*v6*w3 - 19*v3*w4 + 13*CA2*v3*w4 - 
     &   14*CA*CF*v3*w4 + 30*v4*w4 - 10*CA2*v4*w4 - 4*CA*CF*v4*w4 - 
     &   26*v5*w4 + 8*CA2*v5*w4 + 12*CA*CF*v5*w4 - 30*v6*w4 +
     &   18*CA2*v6*w4 - 20*CA*CF*v6*w4 - v7*w4 + CA2*v7*w4 - 
     &   2*CA*CF*v7*w4 + 2*v3*w5 - 2*CA2*v3*w5 + 4*CA*CF*v3*w5 - 
     &   5*v4*w5 - 3*CA2*v4*w5 + 10*CA*CF*v4*w5 - 4*v5*w5 +
     &   8*CA2*v5*w5 - 16*CA*CF*v5*w5 + 45*v6*w5 - 25*CA2*v6*w5 + 
     &   22*CA*CF*v6*w5 + 4*v7*w5 - 4*CA2*v7*w5 + 8*CA*CF*v7*w5 - 
     &   v4*w6 + CA2*v4*w6 - 2*CA*CF*v4*w6 + 16*v5*w6 - 10*CA2*v5*w6 +
     &   8*CA*CF*v5*w6 - 34*v6*w6 + 16*CA2*v6*w6 - 12*CA*CF*v6*w6 - 
     &   7*v7*w6 + 7*CA2*v7*w6 - 14*CA*CF*v7*w6 + 4*v6*w7 +
     &   6*v7*w7 - 6*CA2*v7*w7 + 12*CA*CF*v7*w7 - 2*v7*w8 + 
     &   2*CA2*v7*w8 - 4*CA*CF*v7*w8))/
     &   ((1 - v)**2*v*(1 - w)*w*(1 - v*w)**3*(1 - v + v*w))
      struv1 = part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV2(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (16*CA*CF**2*l1vw*(1 - w)*(1 - v - v2 + v3 + v*w + 
     &   4*v2*w - 5*v3*w - 3*v2*w2 + 5*v3*w2 - v3*w3))/((1 - v)*v*w*
     &   (1 - v + v*w)**2) 
      part2 = -(4*CF*l1v*(1 - v - v*w)*(2*CA2 - 2*CA2*v + 2*CA2*v2 - 
     &   2*CA2*v3 - 2*v*w + CA2*v*w - 2*CA2*v2*w + 2*v3*w + 
     &   5*CA2*v3*w + 2*v2*w2 + CA2*v2*w2 - 2*v3*w2 - 5*CA2*v3*w2 + 
     &   2*CA2*v3*w3))/((1 - v)**2*v2*w2*(1 - v + v*w)) 
      part3 = -(8*CA*CF**2*lmss*(1 + v2 - 2*v2*w + v2*w2)*
     &   (1 - 4*v + 6*v2 - 4*v3 + v4 + v*w - 3*v2*w + 3*v3*w -
     &   v4*w + v2*w2 - 2*v3*w2 + v4*w2 + v3*w3 - v4*w3 + v4*w4))/
     &   ((1 - v)**2*v2*w2*(1 - v + v*w)**2) 
      part4 = -(4*CF*lw*(2 - 4*CA2 - 6*v + 6*CA2*v + 8*v2 - 8*CA2*v2 -
     &   8*v3 + 8*CA2*v3 + 6*v4 - 4*CA2*v4 - 2*v5 + 2*CA2*v5 + 2*v*w +
     &   2*CA2*v*w - 6*v2*w + 2*CA2*v2*w + 10*v3*w - 14*CA2*v3*w - 
     &   10*v4*w + 6*CA2*v4*w + 4*v5*w - 4*CA2*v5*w - 3*v2*w2 -
     &   CA2*v2*w2 - 5*v3*w2 + 13*CA2*v3*w2 + 11*v4*w2 - 3*CA2*v4*w2 -
     &   3*v5*w2 + 3*CA2*v5*w2 + 5*v3*w3 - 5*CA2*v3*w3 - 6*v4*w3 -
     &   2*CA2*v4*w3 + v5*w3 - CA2*v5*w3 + 2*CA2*v4*w4))/
     &   ((1 - v)**2*v2*w2*(1 - v + v*w)) 
      part5 = (4*CF*lvw*(1 - w)*(4*CA2 + 4*CA2*v2 + 2*v*w - 9*CA2*v*w -
     &   6*CA2*v2*w - 2*v3*w - CA2*v3*w - 5*v2*w2 + 12*CA2*v2*w2 +
     &   2*v3*w2 + 3*CA2*v3*w2 - v4*w2 + CA2*v4*w2 + 2*v3*w3 - 
     &   6*CA2*v3*w3 + 2*v4*w3 - 2*CA2*v4*w3 - 2*v4*w4 +
     &   2*CA2*v4*w4))/((1 - v)**2*v2*w2) 
      part6 = (2*CF*lms*(2 - 6*CA2 - 4*v + 4*v2 - 8*CA2*v2 - 4*v3 + 
     &   2*v4 - 2*CA2*v4 - 2*w + 2*CA2*w + 20*CA2*v*w + 4*v2*w + 
     &   12*CA2*v2*w - 4*v3*w + 16*CA2*v3*w + 6*v4*w + 2*CA2*v4*w - 
     &   4*v5*w + 4*CA2*v5*w + w2 - CA2*w2 + 2*v*w2 - 4*CA2*v*w2 - 
     &   2*v2*w2 - 36*CA2*v2*w2 + 2*v3*w2 - 24*CA2*v3*w2 - v4*w2 - 
     &   13*CA2*v4*w2 - 4*CA2*v5*w2 + 2*v6*w2 - 2*CA2*v6*w2 - 2*v*w3 + 
     &   2*CA2*v*w3 + 2*v2*w3 + 2*CA2*v2*w3 - 6*v3*w3 + 40*CA2*v3*w3 -
     &   4*v4*w3 + 20*CA2*v4*w3 + 6*CA2*v5*w3 - 2*v6*w3 + 2*CA2*v6*w3 +
     &   v2*w4 - CA2*v2*w4 - 2*v3*w4 + 11*v4*w4 - 27*CA2*v4*w4 + 
     &   2*v5*w4 - 8*CA2*v5*w4 + 2*v6*w4 - 2*CA2*v6*w4 - 6*v5*w5 +
     &   10*CA2*v5*w5 - 2*v6*w5 + 2*CA2*v6*w5 + 2*v6*w6 - 
     &   2*CA2*v6*w6))/((1 - v)**2*v2*w2*(1 - v*w)**2) 
      part7 = -(2*CF*l1w*(4 - 12*CA2 - 16*v + 32*CA2*v + 28*v2 - 
     &   44*CA2*v2 - 32*v3 + 48*CA2*v3 + 28*v4 - 36*CA2*v4 - 16*v5 + 
     &   16*CA2*v5 + 4*v6 - 4*CA2*v6 - 2*w + 2*CA2*w + 18*v*w + 
     &   8*CA2*v*w - 28*v2*w - 28*CA2*v2*w + 32*v3*w - 4*CA2*v3*w - 
     &   38*v4*w + 14*CA2*v4*w + 6*v5*w + 20*CA2*v5*w + 20*v6*w - 
     &   20*CA2*v6*w - 8*v7*w + 8*CA2*v7*w + w2 - CA2*w2 - 4*v*w2 + 
     &   2*CA2*v*w2 - 17*v2*w2 + CA2*v2*w2 + 28*v3*w2 + 46*CA2*v3*w2 -
     &   13*v4*w2 - 13*CA2*v4*w2 + 72*v5*w2 - 88*CA2*v5*w2 - 79*v6*w2 +
     &   49*CA2*v6*w2 + 8*v7*w2 - 8*CA2*v7*w2 + 4*v8*w2 - 4*CA2*v8*w2 +
     &   6*v2*w3 - 6*CA2*v2*w3 - 6*v3*w3 - 12*CA2*v3*w3 - 
     &   32*CA2*v4*w3 - 92*v5*w3 + 126*CA2*v5*w3 + 74*v6*w3 -
     &   22*CA2*v6*w3 + 30*v7*w3 - 18*CA2*v7*w3 - 12*v8*w3 + 
     &   12*CA2*v8*w3 - 2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 - 2*CA2*v3*w4 -
     &   2*v4*w4 + 20*CA2*v4*w4 + 64*v5*w4 - 66*CA2*v5*w4 - 8*v6*w4 -
     &   40*CA2*v6*w4 - 66*v7*w4 + 38*CA2*v7*w4 + 16*v8*w4 -
     &   16*CA2*v8*w4 - 4*v4*w5 + 4*CA2*v4*w5 - 16*v5*w5 + 
     &   8*CA2*v5*w5 - 38*v6*w5 + 50*CA2*v6*w5 + 50*v7*w5 - 
     &   22*CA2*v7*w5 - 16*v8*w5 + 16*CA2*v8*w5 + v4*w6 - CA2*v4*w6 -
     &   2*v5*w6 + 23*v6*w6 - 17*CA2*v6*w6 - 10*v7*w6 - 2*CA2*v7*w6 +
     &   16*v8*w6 - 16*CA2*v8*w6 - 4*v7*w7 + 4*CA2*v7*w7 - 12*v8*w7 +
     &   12*CA2*v8*w7 + 4*v8*w8 - 4*CA2*v8*w8))/
     &   ((1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part8 = -(2*CF*lv*(4 - 16*CA2 - 16*v + 44*CA2*v + 28*v2 -
     &   60*CA2*v2 - 32*v3 + 64*CA2*v3 + 28*v4 - 48*CA2*v4 - 16*v5 +
     &   20*CA2*v5 + 4*v6 - 4*CA2*v6 - 2*w + 2*CA2*w + 22*v*w + 
     &   14*CA2*v*w - 36*v2*w - 44*CA2*v2*w + 32*v3*w + 8*CA2*v3*w -
     &   30*v4*w + 6*CA2*v4*w + 2*v5*w + 34*CA2*v5*w + 20*v6*w - 
     &   28*CA2*v6*w - 8*v7*w + 8*CA2*v7*w + w2 - CA2*w2 - 4*v*w2 + 
     &   2*CA2*v*w2 - 29*v2*w2 + 3*CA2*v2*w2 + 56*v3*w2 + 
     &   52*CA2*v3*w2 - 25*v4*w2 - 7*CA2*v4*w2 + 60*v5*w2 -
     &   114*CA2*v5*w2 - 71*v6*w2 + 57*CA2*v6*w2 + 8*v7*w2 - 
     &   4*CA2*v7*w2 + 4*v8*w2 - 4*CA2*v8*w2 + 6*v2*w3 - 6*CA2*v2*w3 +
     &   2*v3*w3 - 20*CA2*v3*w3 - 32*v4*w3 - 40*CA2*v4*w3 - 64*v5*w3 +
     &   148*CA2*v5*w3 + 74*v6*w3 - 10*CA2*v6*w3 + 26*v7*w3 -
     &   28*CA2*v7*w3 - 12*v8*w3 + 12*CA2*v8*w3 - 2*v2*w4 + 
     &   2*CA2*v2*w4 + 6*v3*w4 - 2*CA2*v3*w4 + 6*v4*w4 + 28*CA2*v4*w4 +
     &   72*v5*w4 - 74*CA2*v5*w4 - 28*v6*w4 - 70*CA2*v6*w4 - 62*v7*w4 +
     &   44*CA2*v7*w4 + 16*v8*w4 - 16*CA2*v8*w4 - 4*v4*w5 +
     &   4*CA2*v4*w5 - 28*v5*w5 + 6*CA2*v5*w5 - 30*v6*w5 + 
     &   74*CA2*v6*w5 + 54*v7*w5 - 16*CA2*v7*w5 - 16*v8*w5 + 
     &   16*CA2*v8*w5 + v4*w6 - CA2*v4*w6 - 2*v5*w6 + 27*v6*w6 -
     &   23*CA2*v6*w6 - 14*v7*w6 - 12*CA2*v7*w6 + 16*v8*w6 - 
     &   16*CA2*v8*w6 - 4*v7*w7 + 8*CA2*v7*w7 - 12*v8*w7 + 
     &   12*CA2*v8*w7 + 4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part9 = -(2*CF*(2 - 6*CA2 - 12*v + 28*CA2*v + 30*v2 - 54*CA2*v2 -
     &   40*v3 + 56*CA2*v3 + 30*v4 - 34*CA2*v4 - 12*v5 + 12*CA2*v5 +
     &   2*v6 - 2*CA2*v6 + w - CA2*w - v*w + CA2*v*w - 2*v2*w - 
     &   6*CA2*v2*w - 4*v3*w + 28*CA2*v3*w + 23*v4*w - 47*CA2*v4*w -
     &   31*v5*w + 39*CA2*v5*w + 18*v6*w - 18*CA2*v6*w - 4*v7*w +
     &   4*CA2*v7*w - w2 + CA2*w2 + 9*v*w2 - 7*CA2*v*w2 - 32*v2*w2 +
     &   22*CA2*v2*w2 + 34*v3*w2 - 38*CA2*v3*w2 - 11*v4*w2 +
     &   39*CA2*v4*w2 + 13*v5*w2 - 27*CA2*v5*w2 - 14*v6*w2 + 
     &   12*CA2*v6*w2 + 2*v8*w2 - 2*CA2*v8*w2 - 4*v2*w3 + 4*CA2*v2*w3 +
     &   58*v3*w3 - 28*CA2*v3*w3 - 80*v4*w3 + 32*CA2*v4*w3 - 4*v5*w3 +
     &   32*v6*w3 - 8*CA2*v6*w3 + 4*v7*w3 - 6*CA2*v7*w3 - 6*v8*w3 +
     &   6*CA2*v8*w3 + 2*v2*w4 - 2*CA2*v2*w4 - 16*v3*w4 + 
     &   12*CA2*v3*w4 - 12*v4*w4 + 2*CA2*v4*w4 + 70*v5*w4 -
     &   32*CA2*v5*w4 - 44*v6*w4 + 18*CA2*v6*w4 - 10*v7*w4 + 
     &   8*CA2*v7*w4 + 10*v8*w4 - 10*CA2*v8*w4 + 3*v4*w5 - 
     &   3*CA2*v4*w5 - 11*v5*w5 + 13*CA2*v5*w5 + 10*v6*w5 -
     &   6*CA2*v6*w5 + 10*v7*w5 - 8*CA2*v7*w5 - 12*v8*w5 + 
     &   12*CA2*v8*w5 - v4*w6 + CA2*v4*w6 + 7*v5*w6 - 5*CA2*v5*w6 -
     &   12*v6*w6 + 4*CA2*v6*w6 - 2*v7*w6 + 4*CA2*v7*w6 + 10*v8*w6 -
     &   10*CA2*v8*w6 + 2*v7*w7 - 2*CA2*v7*w7 - 6*v8*w7 +
     &   6*CA2*v8*w7 + 2*v8*w8 - 2*CA2*v8*w8))/
     &   ((1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      struv2 = part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV3(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (4*CF*l1v*(6 + 2*CA2 - 4*v + CA2*v + 
     &     CA2*v2 + 4*v*w +
     &        CA2*v*w - 4*v2*w - CA2*v2*w - 2*v2*w2))/
     &    ((1 - v)**2*(1 - v*w)) 
      part2 = (4*CA*CF**2*lmss*(4 - 5*v + 2*v2 - v3 + 9*v*w - 8*v2*w +
     &        3*v3*w + 6*v2*w2 - 4*v3*w2 + 2*v3*w3))/
     &    ((1 - v)**2*(1 - v + v*w)) 
      part3 = (4*CF*lvw*(2*CA*CF - v2 + CA2*v2 - 3*w + 
     &     4*CA2*w + 7*v*w -
     &        6*CA2*v*w - v2*w + 3*CA2*v2*w + 
     &     v3*w - CA2*v3*w -
     &        7*v*w2 + 6*CA2*v*w2 + 6*v2*w2 - 
     &     5*CA2*v2*w2 -
     &        3*v3*w2 + 3*CA2*v3*w2 + 4*CA2*v2*w3 +
     &        4*v3*w3 - 4*CA2*v3*w3 - 2*v3*w4 +
     &        2*CA2*v3*w4))/((1 - v)**2*v*w) 
      part4 = (4*CF*l1vw*(4 - 3*CA2 + 4*CA*CF - 11*v + 
     &     8*CA2*v - 12*CA*CF*v +
     &        11*v2 - 8*CA2*v2 + 14*CA*CF*v2 - 5*v3 +
     &        4*CA2*v3 - 8*CA*CF*v3 + v4 - CA2*v4 +
     &        2*CA*CF*v4 + 11*v*w - 8*CA2*v*w + 12*CA*CF*v*w -
     &        22*v2*w + 16*CA2*v2*w - 28*CA*CF*v2*w + 
     &     15*v3*w -
     &        12*CA2*v3*w + 24*CA*CF*v3*w - 4*v4*w + 
     &     4*CA2*v4*w -
     &        8*CA*CF*v4*w + 11*v2*w2 - 9*CA2*v2*w2 +
     &        18*CA*CF*v2*w2 - 16*v3*w2 + 14*CA2*v3*w2 -
     &        28*CA*CF*v3*w2 + 7*v4*w2 - 7*CA2*v4*w2 +
     &        14*CA*CF*v4*w2 + 6*v3*w3 - 6*CA2*v3*w3 +
     &        12*CA*CF*v3*w3 - 6*v4*w3 + 6*CA2*v4*w3 -
     &        12*CA*CF*v4*w3 + 2*v4*w4 - 2*CA2*v4*w4 +
     &        4*CA*CF*v4*w4))/((1 - v)**2*v*(1 - v + v*w)) 
      part5 = -(2*CF*lms*(4*CA*CF + 4*CA*CF*v2 - 2*CA*CF*w + 2*v*w -
     &        18*CA*CF*v*w - 4*CA*CF*v2*w + 
     &     2*v3*w - 8*CA*CF*v3*w +
     &        6*CA*CF*v*w2 - 5*v2*w2 + CA2*v2*w2 +
     &        18*CA*CF*v2*w2 - 4*v3*w2 + 16*CA*CF*v3*w2 -
     &        v4*w2 + CA2*v4*w2 + 8*CA*CF*v4*w2 -
     &        6*CA*CF*v2*w3 + 8*v3*w3 - 2*CA2*v3*w3 -
     &        6*CA*CF*v3*w3 + 2*v4*w3 - 2*CA2*v4*w3 -
     &        16*CA*CF*v4*w3 + 2*v5*w3 - 4*CA*CF*v5*w3 +
     &        2*CA*CF*v3*w4 - 7*v4*w4 + 3*CA2*v4*w4 +
     &        2*CA*CF*v4*w4 - 4*v5*w4 + 4*CA*CF*v5*w4 -
     &        v6*w4 + CA2*v6*w4 + 
     &     6*v5*w5 - 2*CA2*v5*w5 +
     &        2*v6*w5 - 2*CA2*v6*w5 - 2*v6*w6 +
     &        2*CA2*v6*w6))/((1 - v)**2*v*w*(1 - v*w)**3) 
      part6 = -(2*CF*l1w*(2 - 2*CA2 - 2*v + 2*CA2*v + 
     &     2*v2 - 2*CA2*v2 -
     &        2*v3 + 2*CA2*v3 - w + CA2*w - 
     &     8*v*w + 14*CA2*v*w +
     &        2*v2*w - 14*CA2*v2*w + 4*v3*w + 2*CA2*v3*w +
     &        7*v4*w - 7*CA2*v4*w + 2*v*w2 - 2*CA2*v*w2 +
     &        7*v2*w2 - 3*CA2*v2*w2 - 11*v3*w2 +
     &        13*CA2*v3*w2 - 18*v4*w2 + 8*CA2*v4*w2 -
     &        8*v5*w2 + 12*CA2*v5*w2 + 7*v3*w3 -
     &        23*CA2*v3*w3 + 32*v4*w3 + 2*CA2*v4*w3 +
     &        6*v5*w3 - 24*CA2*v5*w3 + 7*v6*w3 -
     &        7*CA2*v6*w3 - 2*v3*w4 + 2*CA2*v3*w4 -
     &        19*v4*w4 + 11*CA2*v4*w4 - 5*v5*w4 +
     &        7*CA2*v5*w4 - 8*v6*w4 + 14*CA2*v6*w4 -
     &        2*v7*w4 + 2*CA2*v7*w4 + 
     &     v4*w5 - CA2*v4*w5 +
     &        3*v5*w5 + 3*CA2*v5*w5 + 2*v6*w5 -
     &        8*CA2*v6*w5 + 6*v7*w5 - 6*CA2*v7*w5 -
     &        8*v7*w6 + 8*CA2*v7*w6 + 4*v7*w7 -
     &        4*CA2*v7*w7))/((1 - v)**2*v*w*(1 - v*w)**3*
     &     (1 - v + v*w))
      part7 = -(2*CF*lv*(2 - 2*CA2 - 2*v + 2*CA2*v + 
     &     2*v2 - 2*CA2*v2 -
     &        2*v3 + 2*CA2*v3 - w + CA2*w + 
     &     20*v*w + 10*CA2*v*w -
     &        26*v2*w - 12*CA2*v2*w + 12*v3*w + 
     &     2*CA2*v3*w -
     &        v4*w - 5*CA2*v4*w + 2*v*w2 - 2*CA2*v*w2 -
     &        5*v2*w2 - CA2*v2*w2 + 13*v3*w2 +
     &        11*CA2*v3*w2 - 10*v4*w2 + 4*CA2*v4*w2 +
     &        8*v5*w2 + 8*CA2*v5*w2 - 41*v3*w3 -
     &        17*CA2*v3*w3 + 56*v4*w3 + 2*CA2*v4*w3 -
     &        34*v5*w3 - 16*CA2*v5*w3 - v6*w3 -
     &        5*CA2*v6*w3 - 2*v3*w4 + 2*CA2*v3*w4 -
     &        3*v4*w4 + 9*CA2*v4*w4 - 13*v5*w4 +
     &        5*CA2*v5*w4 + 16*v6*w4 + 10*CA2*v6*w4 -
     &        2*v7*w4 + 2*CA2*v7*w4 + 
     &     v4*w5 - CA2*v4*w5 +
     &        23*v5*w5 + CA2*v5*w5 - 10*v6*w5 -
     &        6*CA2*v6*w5 + 6*v7*w5 - 6*CA2*v7*w5 -
     &        4*v6*w6 - 8*v7*w6 + 
     &     8*CA2*v7*w6 + 4*v7*w7 -
     &        4*CA2*v7*w7))/((1 - v)**2*v*w*(1 - v*w)**3*
     &     (1 - v + v*w))
      part8 = -(CF*(6 - 6*CA2 + 4*CA*CF - 14*v + 
     &     14*CA2*v - 12*CA*CF*v +
     &        14*v2 - 14*CA2*v2 + 12*CA*CF*v2 - 6*v3 +
     &        6*CA2*v3 - 4*CA*CF*v3 - 12*w + 
     &     4*CA2*w - 4*CA*CF*w +
     &        39*v*w - 11*CA2*v*w + 4*CA*CF*v*w - 42*v2*w -
     &        2*CA2*v2*w + 24*CA*CF*v2*w + 13*v3*w +
     &        23*CA2*v3*w - 44*CA*CF*v3*w + 6*v4*w -
     &        18*CA2*v4*w + 20*CA*CF*v4*w + 24*v*w2 -
     &        8*CA2*v*w2 + 8*CA*CF*v*w2 - 34*v2*w2 +
     &        18*CA2*v2*w2 - 36*CA*CF*v2*w2 - 31*v3*w2 -
     &        5*CA2*v3*w2 + 20*CA*CF*v3*w2 + 71*v4*w2 -
     &        19*CA2*v4*w2 + 40*CA*CF*v4*w2 - 40*v5*w2 +
     &        24*CA2*v5*w2 - 32*CA*CF*v5*w2 + 16*v3*w3 +
     &        36*CA*CF*v3*w3 - 7*v4*w3 - 9*CA2*v4*w3 -
     &        64*CA*CF*v4*w3 - 31*v5*w3 + 19*CA2*v5*w3 +
     &        8*CA*CF*v5*w3 + 38*v6*w3 - 26*CA2*v6*w3 +
     &        20*CA*CF*v6*w3 - 24*v3*w4 + 8*CA2*v3*w4 -
     &        8*CA*CF*v3*w4 + 12*v4*w4 + 4*CA2*v4*w4 +
     &        37*v5*w4 - 25*CA2*v5*w4 + 32*CA*CF*v5*w4 -
     &        45*v6*w4 + 33*CA2*v6*w4 - 20*CA*CF*v6*w4 -
     &        6*v7*w4 + 6*CA2*v7*w4 - 4*CA*CF*v7*w4 +
     &        12*v4*w5 - 4*CA2*v4*w5 + 4*CA*CF*v4*w5 -
     &        19*v5*w5 + 7*CA2*v5*w5 - 8*CA*CF*v5*w5 +
     &        25*v6*w5 - 21*CA2*v6*w5 + 10*v7*w5 -
     &        10*CA2*v7*w5 + 4*CA*CF*v7*w5 - 8*v6*w6 +
     &        8*CA2*v6*w6 - 8*v7*w6 + 8*CA2*v7*w6 +
     &        4*v7*w7 - 4*CA2*v7*w7))/
     &    ((1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w)) 
      part9 = -(2*CF*lw*(3 - 3*CA2 + 2*CA*CF - 3*v + 
     &     3*CA2*v - 2*CA*CF*v +
     &        3*v2 - 3*CA2*v2 + 2*CA*CF*v2 - 
     &     3*v3 + 3*CA2*v3 -
     &        2*CA*CF*v3 - 10*w - 2*CA2*w + 8*CA*CF*w + 16*v*w +
     &        14*CA2*v*w - 28*CA*CF*v*w - 15*v2*w - 
     &     19*CA2*v2*w +
     &        34*CA*CF*v2*w + 12*v3*w + 10*CA2*v3*w -
     &        20*CA*CF*v3*w + 3*v4*w - 9*CA2*v4*w + 
     &     10*CA*CF*v4*w +
     &        w2 - CA2*w2 + 2*CA*CF*w2 + 7*v*w2 - 
     &     CA2*v*w2 -
     &        2*CA*CF*v*w2 - 17*v2*w2 - 3*CA2*v2*w2 +
     &        10*CA*CF*v2*w2 + 14*v3*w2 + 16*CA2*v3*w2 -
     &        24*CA*CF*v3*w2 - 31*v4*w2 - CA2*v4*w2 +
     &        14*CA*CF*v4*w2 + 4*v5*w2 + 8*CA2*v5*w2 -
     &        12*CA*CF*v5*w2 - 2*v*w3 + 2*CA2*v*w3 -
     &        4*CA*CF*v*w3 - 3*v2*w3 + CA2*v2*w3 +
     &        2*CA*CF*v2*w3 + 10*v3*w3 - 8*CA2*v3*w3 +
     &        8*CA*CF*v3*w3 + 34*v4*w3 - 14*CA2*v4*w3 +
     &        8*CA*CF*v4*w3 - 2*CA2*v5*w3 - 
     &     4*CA*CF*v5*w3 -
     &        v6*w3 - 5*CA2*v6*w3 + 10*CA*CF*v6*w3 +
     &        13*v3*w4 + 5*CA2*v3*w4 - 14*CA*CF*v3*w4 -
     &        66*v4*w4 + 14*CA2*v4*w4 - 4*CA*CF*v4*w4 +
     &        6*v5*w4 + 12*CA*CF*v5*w4 + 2*v6*w4 +
     &        10*CA2*v6*w4 - 20*CA*CF*v6*w4 - v7*w4 +
     &        CA2*v7*w4 - 2*CA*CF*v7*w4 + 2*v3*w5 -
     &        2*CA2*v3*w5 + 4*CA*CF*v3*w5 + 19*v4*w5 -
     &        9*CA2*v4*w5 + 10*CA*CF*v4*w5 + 28*v5*w5 -
     &        16*CA*CF*v5*w5 - 11*v6*w5 - 11*CA2*v6*w5 +
     &        22*CA*CF*v6*w5 + 4*v7*w5 - 4*CA2*v7*w5 +
     &        8*CA*CF*v7*w5 - v4*w6 + CA2*v4*w6 -
     &        2*CA*CF*v4*w6 - 24*v5*w6 + 8*CA*CF*v5*w6 +
     &        6*v6*w6 + 6*CA2*v6*w6 - 12*CA*CF*v6*w6 -
     &        7*v7*w6 + 7*CA2*v7*w6 - 14*CA*CF*v7*w6 +
     &        4*v6*w7 + 6*v7*w7 - 6*CA2*v7*w7 +
     &        12*CA*CF*v7*w7 - 2*v7*w8 + 2*CA2*v7*w8 -
     &        4*CA*CF*v7*w8))/
     &    ((1 - v)**2*v*(1 - w)*w*(1 - v*w)**3*(1 - v + v*w))
      struv3 = part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV4(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (8*CF*l1vw*(1 - w)*(3 - 3*v - 3*v2 + 
     &     3*v3 + 3*v*w + 4*v2*w +
     &        2*CA2*v2*w - 7*v3*w - 2*CA2*v3*w - v2*w2 -
     &        2*CA2*v2*w2 + 7*v3*w2 + 2*CA2*v3*w2 -
     &        3*v3*w3))/((1 - v)*v*w*(1 - v + v*w)**2) 
      part2 = -(8*CF*l1v*(1 - v - v*w)*
     &      (CA2 - CA2*v + CA2*v2 - CA2*v3 + v*w - 
     &     CA2*v2*w -
     &        v3*w + 3*CA2*v3*w - v2*w2 + CA2*v2*w2 +
     &        v3*w2 - 3*CA2*v3*w2 + CA2*v3*w3))/
     &    ((1 - v)**2*v2*w2*(1 - v + v*w)) 
      part3 = -(8*CA*CF**2*lmss*(1 + v2 - 2*v2*w + v2*w2)*
     &      (1 - 4*v + 6*v2 - 4*v3 + v4 + v*w - 
     &     3*v2*w + 3*v3*w -
     &        v4*w + v2*w2 - 2*v3*w2 + v4*w2 + 
     &     v3*w3 -
     &        v4*w3 + v4*w4))/((1 - v)**2*v2*w2*
     &     (1 - v + v*w)**2)
      part4 = -(4*CF*lw*(2 - 4*CA2 - 6*v + 6*CA2*v + 
     &     8*v2 - 8*CA2*v2 -
     &        8*v3 + 8*CA2*v3 + 6*v4 - 4*CA2*v4 - 2*v5 +
     &        2*CA2*v5 + 2*v*w + 2*CA2*v*w - 6*v2*w +
     &        2*CA2*v2*w + 10*v3*w - 14*CA2*v3*w - 
     &     10*v4*w +
     &        6*CA2*v4*w + 4*v5*w - 4*CA2*v5*w + 
     &     5*v2*w2 -
     &        3*CA2*v2*w2 - 5*v3*w2 + 13*CA2*v3*w2 +
     &        3*v4*w2 - CA2*v4*w2 - 3*v5*w2 +
     &        3*CA2*v5*w2 - 3*v3*w3 - 3*CA2*v3*w3 +
     &        2*v4*w3 - 4*CA2*v4*w3 + v5*w3 - 
     &     CA2*v5*w3 +
     &        2*CA2*v4*w4))/((1 - v)**2*v2*w2*
     &     (1 - v + v*w)) 
      part5 = (4*CF*lvw*(1 - w)*(4*CA2 + 4*CA2*v2 - 2*v*w - 
     &     8*CA2*v*w -
     &        6*CA2*v2*w + 2*v3*w - 2*CA2*v3*w - v2*w2 +
     &        11*CA2*v2*w2 - 2*v3*w2 + 4*CA2*v3*w2 -
     &        v4*w2 + CA2*v4*w2 + 2*v3*w3 - 
     &     6*CA2*v3*w3 +
     &        2*v4*w3 - 2*CA2*v4*w3 - 2*v4*w4 +
     &        2*CA2*v4*w4))/((1 - v)**2*v2*w2) 
      part6 = (2*CF*lms*(2 - 6*CA2 - 4*v + 4*v2 - 
     &     8*CA2*v2 - 4*v3 +
     &        2*v4 - 2*CA2*v4 - 2*w + 2*CA2*w + 20*CA2*v*w +
     &        4*v2*w + 12*CA2*v2*w - 4*v3*w + 16*CA2*v3*w +
     &        6*v4*w + 2*CA2*v4*w - 4*v5*w + 
     &     4*CA2*v5*w + w2 -
     &        CA2*w2 + 2*v*w2 - 4*CA2*v*w2 - 2*v2*w2 -
     &        36*CA2*v2*w2 + 2*v3*w2 - 24*CA2*v3*w2 -
     &        v4*w2 - 13*CA2*v4*w2 - 4*CA2*v5*w2 +
     &        2*v6*w2 - 2*CA2*v6*w2 - 
     &     2*v*w3 + 2*CA2*v*w3 +
     &        2*v2*w3 + 2*CA2*v2*w3 - 6*v3*w3 +
     &        40*CA2*v3*w3 - 4*v4*w3 + 20*CA2*v4*w3 +
     &        6*CA2*v5*w3 - 2*v6*w3 + 2*CA2*v6*w3 +
     &        v2*w4 - CA2*v2*w4 - 
     &     2*v3*w4 + 11*v4*w4 -
     &        27*CA2*v4*w4 + 2*v5*w4 - 8*CA2*v5*w4 +
     &        2*v6*w4 - 2*CA2*v6*w4 - 6*v5*w5 +
     &        10*CA2*v5*w5 - 2*v6*w5 + 2*CA2*v6*w5 +
     &        2*v6*w6 - 2*CA2*v6*w6))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2)  
      part7 = -(2*CF*l1w*(4 - 12*CA2 - 16*v + 32*CA2*v + 28*v2 -
     &        44*CA2*v2 - 32*v3 + 48*CA2*v3 + 28*v4 -
     &        36*CA2*v4 - 16*v5 + 16*CA2*v5 + 4*v6 -
     &        4*CA2*v6 - 2*w + 2*CA2*w - 6*v*w + 14*CA2*v*w +
     &        20*v2*w - 40*CA2*v2*w - 32*v3*w + 
     &     12*CA2*v3*w +
     &        42*v4*w - 6*CA2*v4*w - 34*v5*w + 
     &     30*CA2*v5*w +
     &        20*v6*w - 20*CA2*v6*w - 8*v7*w + 
     &     8*CA2*v7*w +
     &        w2 - CA2*w2 - 4*v*w2 + 
     &     2*CA2*v*w2 + 23*v2*w2 -
     &        9*CA2*v2*w2 - 12*v3*w2 + 56*CA2*v3*w2 -
     &        37*v4*w2 - 7*CA2*v4*w2 + 16*v5*w2 -
     &        74*CA2*v5*w2 + v6*w2 + 29*CA2*v6*w2 +
     &        8*v7*w2 - 8*CA2*v7*w2 + 4*v8*w2 -
     &        4*CA2*v8*w2 + 6*v2*w3 - 6*CA2*v2*w3 -
     &        38*v3*w3 - 4*CA2*v3*w3 + 64*v4*w3 -
     &        48*CA2*v4*w3 + 44*v5*w3 + 92*CA2*v5*w3 -
     &        54*v6*w3 + 10*CA2*v6*w3 - 10*v7*w3 -
     &        8*CA2*v7*w3 - 12*v8*w3 + 12*CA2*v8*w3 -
     &        2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 -
     &        2*CA2*v3*w4 - 2*v4*w4 + 20*CA2*v4*w4 -
     &        96*v5*w4 - 26*CA2*v5*w4 + 48*v6*w4 -
     &        54*CA2*v6*w4 + 38*v7*w4 + 12*CA2*v7*w4 +
     &        16*v8*w4 - 16*CA2*v8*w4 - 4*v4*w5 +
     &        4*CA2*v4*w5 + 40*v5*w5 - 6*CA2*v5*w5 +
     &        10*v6*w5 + 38*CA2*v6*w5 - 54*v7*w5 +
     &        4*CA2*v7*w5 - 16*v8*w5 + 16*CA2*v8*w5 +
     &        v4*w6 - CA2*v4*w6 - 2*v5*w6 - 
     &     17*v6*w6 -
     &        7*CA2*v6*w6 + 30*v7*w6 - 12*CA2*v7*w6 +
     &        16*v8*w6 - 16*CA2*v8*w6 - 4*v7*w7 +
     &        4*CA2*v7*w7 - 12*v8*w7 + 12*CA2*v8*w7 +
     &        4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part8 = -(2*CF*lv*(4 - 16*CA2 - 16*v + 44*CA2*v + 28*v2 -
     &        60*CA2*v2 - 32*v3 + 64*CA2*v3 + 28*v4 -
     &        48*CA2*v4 - 16*v5 + 20*CA2*v5 + 4*v6 -
     &        4*CA2*v6 - 2*w + 2*CA2*w - 10*v*w + 22*CA2*v*w +
     &        28*v2*w - 60*CA2*v2*w - 32*v3*w + 
     &     24*CA2*v3*w +
     &        34*v4*w - 10*CA2*v4*w - 30*v5*w + 
     &     42*CA2*v5*w +
     &        20*v6*w - 28*CA2*v6*w - 8*v7*w + 8*CA2*v7*w +
     &        w2 - CA2*w2 - 4*v*w2 + 
     &     2*CA2*v*w2 + 35*v2*w2 -
     &        13*CA2*v2*w2 - 40*v3*w2 + 76*CA2*v3*w2 -
     &        25*v4*w2 - 7*CA2*v4*w2 + 28*v5*w2 -
     &        106*CA2*v5*w2 - 7*v6*w2 + 41*CA2*v6*w2 +
     &        8*v7*w2 - 4*CA2*v7*w2 + 4*v8*w2 -
     &        4*CA2*v8*w2 + 6*v2*w3 - 6*CA2*v2*w3 -
     &        46*v3*w3 - 8*CA2*v3*w3 + 96*v4*w3 -
     &        72*CA2*v4*w3 + 16*v5*w3 + 128*CA2*v5*w3 -
     &        54*v6*w3 + 22*CA2*v6*w3 - 6*v7*w3 -
     &        20*CA2*v7*w3 - 12*v8*w3 + 12*CA2*v8*w3 -
     &        2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 -
     &        2*CA2*v3*w4 - 10*v4*w4 + 32*CA2*v4*w4 -
     &        104*v5*w4 - 30*CA2*v5*w4 + 68*v6*w4 -
     &        94*CA2*v6*w4 + 34*v7*w4 + 20*CA2*v7*w4 +
     &        16*v8*w4 - 16*CA2*v8*w4 - 4*v4*w5 +
     &        4*CA2*v4*w5 + 52*v5*w5 - 14*CA2*v5*w5 +
     &        2*v6*w5 + 66*CA2*v6*w5 - 58*v7*w5 +
     &        12*CA2*v7*w5 - 16*v8*w5 + 16*CA2*v8*w5 +
     &        v4*w6 - CA2*v4*w6 - 2*v5*w6 - 
     &     21*v6*w6 -
     &        11*CA2*v6*w6 + 34*v7*w6 - 24*CA2*v7*w6 +
     &        16*v8*w6 - 16*CA2*v8*w6 - 4*v7*w7 +
     &        8*CA2*v7*w7 - 12*v8*w7 + 12*CA2*v8*w7 +
     &        4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      part9 = -(2*CF*(2 - 6*CA2 - 12*v + 28*CA2*v + 30*v2 - 
     &     54*CA2*v2 -
     &        40*v3 + 56*CA2*v3 + 30*v4 - 
     &     34*CA2*v4 - 12*v5 +
     &        12*CA2*v5 + 2*v6 - 2*CA2*v6 + w - 
     &     CA2*w - v*w +
     &        CA2*v*w - 2*v2*w - 6*CA2*v2*w - 4*v3*w +
     &        28*CA2*v3*w + 23*v4*w - 47*CA2*v4*w - 
     &     31*v5*w +
     &        39*CA2*v5*w + 18*v6*w - 18*CA2*v6*w - 
     &     4*v7*w +
     &        4*CA2*v7*w - w2 + CA2*w2 + v*w2 - 
     &     5*CA2*v*w2 +
     &        8*v2*w2 + 12*CA2*v2*w2 + 2*v3*w2 -
     &        30*CA2*v3*w2 - 43*v4*w2 + 47*CA2*v4*w2 +
     &        53*v5*w2 - 37*CA2*v5*w2 - 22*v6*w2 +
     &        14*CA2*v6*w2 + 2*v8*w2 - 2*CA2*v8*w2 -
     &        4*v2*w3 + 4*CA2*v2*w3 - 14*v3*w3 -
     &        10*CA2*v3*w3 + 48*v4*w3 - 36*v5*w3 +
     &        8*CA2*v5*w3 + 12*v7*w3 - 8*CA2*v7*w3 -
     &        6*v8*w3 + 6*CA2*v8*w3 + 2*v2*w4 -
     &        2*CA2*v2*w4 + 8*CA2*v3*w4 - 4*v4*w4 -
     &        2*v5*w4 - 14*CA2*v5*w4 + 12*v6*w4 +
     &        4*CA2*v6*w4 - 18*v7*w4 + 10*CA2*v7*w4 +
     &        10*v8*w4 - 10*CA2*v8*w4 + 3*v4*w5 -
     &        3*CA2*v4*w5 - 3*v5*w5 + 11*CA2*v5*w5 -
     &        6*v6*w5 - 2*CA2*v6*w5 + 18*v7*w5 -
     &        10*CA2*v7*w5 - 12*v8*w5 + 12*CA2*v8*w5 -
     &        v4*w6 + CA2*v4*w6 - v5*w6 - 
     &     3*CA2*v5*w6 +
     &        4*v6*w6 - 10*v7*w6 + 6*CA2*v7*w6 + 
     &     10*v8*w6 -
     &        10*CA2*v8*w6 + 2*v7*w7 - 2*CA2*v7*w7 -
     &        6*v8*w7 + 6*CA2*v8*w7 + 2*v8*w8 -
     &        2*CA2*v8*w8))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      struv4 = part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV5(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (-4*CF*l1v*(10 - 2*CA2 - 20*v + 5*CA2*v + 14*v2 -
     &        4*CA2*v2 - 8*v*w + CA2*v*w + 2*v2*w2))/
     &     (1 - v + v*w)
      part2 = -(4*CF*lvw*(4 - 3*CA2 - 5*v + 4*CA2*v + 2*v2 -
     &        2*CA2*v2 - 7*v*w + 5*CA2*v*w + 6*v2*w -
     &        4*CA2*v2*w - 2*v3*w + 2*CA2*v3*w + 
     &     4*v2*w2 -
     &        4*CA2*v2*w2 - 2*v3*w3 + 2*CA2*v3*w3))/
     &    ((1 - v)*v) 
      part3 = -(4*CA*CF**2*lms*
     &      (1 - 3*v + 4*v2 - 2*v3 - w + 6*v*w - 
     &     8*v2*w + 2*v3*w +
     &        2*v4*w + v*w2 - 12*v2*w2 + 14*v3*w2 - 
     &     6*v4*w2 +
     &        6*v3*w3 - 2*v4*w3 - 2*v4*w4))/
     &    ((1 - v)*v*w*(1 - v*w)) 
      part4 = (4*CF*l1vw*(4 - 5*CA2 + 4*CA*CF - 14*v + 18*CA2*v - 
     &     12*CA*CF*v +
     &        20*v2 - 28*CA2*v2 + 14*CA*CF*v2 - 14*v3 +
     &        22*CA2*v3 - 8*CA*CF*v3 + 4*v4 - 7*CA2*v4 +
     &        2*CA*CF*v4 + 6*v*w - 8*CA2*v*w + 
     &     4*CA*CF*v*w - 16*v2*w +
     &        22*CA2*v2*w - 12*CA*CF*v2*w + 18*v3*w -
     &        24*CA2*v3*w + 12*CA*CF*v3*w - 8*v4*w +
     &        10*CA2*v4*w - 4*CA*CF*v4*w - 4*v2*w2 -
     &        2*CA2*v2*w2 + 6*CA*CF*v2*w2 + 6*v3*w2 +
     &        2*CA2*v3*w2 - 8*CA*CF*v3*w2 - 
     &     4*CA2*v4*w2 +
     &        4*CA*CF*v4*w2 - 10*v3*w3 + 4*CA*CF*v3*w3 +
     &        8*v4*w3 + 6*CA2*v4*w3 - 4*CA*CF*v4*w3 -
     &        4*v4*w4 - 5*CA2*v4*w4 + 2*CA*CF*v4*w4))/
     &    (v*(1 - v + v*w)**3) 
      part5 = -(2*CF*lmss*(2*CA*CF + 2*v - 2*CA2*v - 8*CA*CF*v - 9*v2 +
     &        9*CA2*v2 + 14*CA*CF*v2 + 18*v3 - 18*CA2*v3 -
     &        12*CA*CF*v3 - 19*v4 + 19*CA2*v4 + 4*CA*CF*v4 +
     &        10*v5 - 10*CA2*v5 - 2*v6 + 2*CA2*v6 - 2*v*w +
     &        12*CA*CF*v*w + 14*v2*w - 4*CA2*v2*w - 
     &     54*CA*CF*v2*w -
     &        40*v3*w + 18*CA2*v3*w + 94*CA*CF*v3*w + 
     &     54*v4*w -
     &        28*CA2*v4*w - 82*CA*CF*v4*w - 34*v5*w +
     &        18*CA2*v5*w + 38*CA*CF*v5*w + 
     &     8*v6*w - 4*CA2*v6*w -
     &        8*CA*CF*v6*w - 5*v2*w2 + CA2*v2*w2 +
     &        12*CA*CF*v2*w2 + 30*v3*w2 - 10*CA2*v3*w2 -
     &        46*CA*CF*v3*w2 - 58*v4*w2 + 22*CA2*v4*w2 +
     &        74*CA*CF*v4*w2 + 48*v5*w2 - 20*CA2*v5*w2 -
     &        54*CA*CF*v5*w2 - 14*v6*w2 + 6*CA2*v6*w2 +
     &        16*CA*CF*v6*w2 - 8*v3*w3 + 2*CA2*v3*w3 +
     &        4*CA*CF*v3*w3 + 30*v4*w3 - 12*CA2*v4*w3 -
     &        14*CA*CF*v4*w3 - 40*v5*w3 + 20*CA2*v5*w3 +
     &        14*CA*CF*v5*w3 + 16*v6*w3 - 8*CA2*v6*w3 -
     &        8*CA*CF*v6*w3 - 7*v4*w4 + 3*CA2*v4*w4 +
     &        2*CA*CF*v4*w4 + 22*v5*w4 - 10*CA2*v5*w4 +
     &        2*CA*CF*v5*w4 - 14*v6*w4 + 6*CA2*v6*w4 -
     &        6*v5*w5 + 2*CA2*v5*w5 + 8*v6*w5 -
     &        4*CA2*v6*w5 - 2*v6*w6 + 2*CA2*v6*w6))/
     &    ((1 - v)*v*w*(1 - v + v*w)**3) 
      part6 = -(CF*(4 - 4*CA2 + 4*CA*CF - 16*v + 
     &     16*CA2*v - 16*CA*CF*v +
     &        24*v2 - 24*CA2*v2 + 24*CA*CF*v2 - 16*v3 +
     &        16*CA2*v3 - 16*CA*CF*v3 + 4*v4 - 4*CA2*v4 +
     &        4*CA*CF*v4 - 12*w + 4*CA2*w - 4*CA*CF*w + 39*v*w -
     &        15*CA2*v*w + 24*CA*CF*v*w - 28*v2*w + 
     &     4*CA2*v2*w -
     &        36*CA*CF*v2*w - 36*v3*w + 36*CA2*v3*w +
     &        4*CA*CF*v3*w + 66*v4*w - 42*CA2*v4*w +
     &        24*CA*CF*v4*w - 35*v5*w + 11*CA2*v5*w -
     &        12*CA*CF*v5*w + 6*v6*w + 2*CA2*v6*w - 24*v*w2 +
     &        8*CA2*v*w2 - 8*CA*CF*v*w2 + 102*v2*w2 -
     &        26*CA2*v2*w2 + 12*CA*CF*v2*w2 - 131*v3*w2 +
     &        27*CA2*v3*w2 + 40*CA*CF*v3*w2 + 63*v4*w2 -
     &        23*CA2*v4*w2 - 80*CA*CF*v4*w2 - 27*v5*w2 +
     &        27*CA2*v5*w2 + 32*CA*CF*v5*w2 + 23*v6*w2 -
     &        11*CA2*v6*w2 + 4*CA*CF*v6*w2 - 6*v7*w2 -
     &        2*CA2*v7*w2 - 20*v3*w3 - 36*CA*CF*v3*w3 +
     &        53*v4*w3 - 9*CA2*v4*w3 + 80*CA*CF*v4*w3 -
     &        4*CA2*v5*w3 - 32*CA*CF*v5*w3 - 45*v6*w3 +
     &        9*CA2*v6*w3 - 12*CA*CF*v6*w3 + 12*v7*w3 +
     &        4*CA2*v7*w3 + 24*v3*w4 - 8*CA2*v3*w4 +
     &        8*CA*CF*v3*w4 - 82*v4*w4 + 38*CA2*v4*w4 -
     &        32*CA*CF*v4*w4 + 47*v5*w4 - 23*CA2*v5*w4 +
     &        16*CA*CF*v5*w4 + 21*v6*w4 - CA2*v6*w4 +
     &        12*CA*CF*v6*w4 - 6*v7*w4 - 10*CA2*v7*w4 +
     &        12*v4*w5 - 4*CA2*v4*w5 + 4*CA*CF*v4*w5 -
     &        11*v5*w5 - CA2*v5*w5 - 4*CA*CF*v5*w5 -
     &        5*v6*w5 - 7*CA2*v6*w5 - 4*CA*CF*v6*w5 +
     &        16*CA2*v7*w5 + 8*CA2*v6*w6 - 
     &     8*CA2*v7*w6))/
     &    ((1 - v)*v*w*(1 - v*w)*(1 - v + v*w)**3) 
      part7 = -(2*CF*l1w*(2 - 2*CA2 - 12*v + 12*CA2*v + 32*v2 -
     &        32*CA2*v2 - 48*v3 + 48*CA2*v3 + 42*v4 -
     &        42*CA2*v4 - 20*v5 + 20*CA2*v5 + 4*v6 -
     &        4*CA2*v6 - w + CA2*w + 15*v*w - 21*CA2*v*w -
     &        67*v2*w + 91*CA2*v2*w + 141*v3*w - 
     &     177*CA2*v3*w -
     &        152*v4*w + 176*CA2*v4*w + 76*v5*w - 
     &     82*CA2*v5*w -
     &        8*v6*w + 8*CA2*v6*w - 4*v7*w + 4*CA2*v7*w -
     &        2*v*w2 + 2*CA2*v*w2 + 19*v2*w2 - 
     &     15*CA2*v2*w2 -
     &        54*v3*w2 + 32*CA2*v3*w2 + 48*v4*w2 -
     &        10*CA2*v4*w2 + 28*v5*w2 - 54*CA2*v5*w2 -
     &        67*v6*w2 + 73*CA2*v6*w2 + 28*v7*w2 -
     &        28*CA2*v7*w2 - 7*v3*w3 + 23*CA2*v3*w3 +
     &        60*v4*w3 - 90*CA2*v4*w3 - 144*v5*w3 +
     &        156*CA2*v5*w3 + 143*v6*w3 - 
     &     141*CA2*v6*w3 -
     &        52*v7*w3 + 52*CA2*v7*w3 + 2*v3*w4 -
     &        2*CA2*v3*w4 - 27*v4*w4 + 19*CA2*v4*w4 +
     &        74*v5*w4 - 52*CA2*v5*w4 - 83*v6*w4 +
     &        69*CA2*v6*w4 + 36*v7*w4 - 36*CA2*v7*w4 +
     &        v4*w5 - CA2*v4*w5 - 6*v5*w5 + 
     &     11*v6*w5 -
     &        5*CA2*v6*w5 - 12*v7*w5 + 12*CA2*v7*w5 +
     &        8*v7*w6 - 8*CA2*v7*w6 - 4*v7*w7 +
     &        4*CA2*v7*w7))/((1 - v)*v*w*(1 - v*w)*
     &     (1 - v + v*w)**3) 
      part8 = -(2*CF*lv*(2 - 2*CA2 - 12*v + 12*CA2*v + 32*v2 -
     &        32*CA2*v2 - 48*v3 + 48*CA2*v3 + 42*v4 -
     &        42*CA2*v4 - 20*v5 + 20*CA2*v5 + 4*v6 -
     &        4*CA2*v6 - w + CA2*w - 13*v*w - 17*CA2*v*w +
     &        73*v2*w + 69*CA2*v2*w - 147*v3*w - 
     &     127*CA2*v3*w +
     &        152*v4*w + 118*CA2*v4*w - 88*v5*w - 
     &     48*CA2*v5*w +
     &        28*v6*w - 4*v7*w + 4*CA2*v7*w - 2*v*w2 +
     &        2*CA2*v*w2 + 7*v2*w2 - 13*CA2*v2*w2 -
     &        18*v3*w2 + 24*CA2*v3*w2 + 32*v4*w2 -
     &        2*CA2*v4*w2 - 36*v5*w2 - 46*CA2*v5*w2 +
     &        25*v6*w2 + 55*CA2*v6*w2 - 8*v7*w2 -
     &        20*CA2*v7*w2 + 41*v3*w3 + 17*CA2*v3*w3 -
     &        108*v4*w3 - 66*CA2*v4*w3 + 112*v5*w3 +
     &        112*CA2*v5*w3 - 65*v6*w3 - 99*CA2*v6*w3 +
     &        20*v7*w3 + 36*CA2*v7*w3 + 2*v3*w4 -
     &        2*CA2*v3*w4 - 11*v4*w4 + 17*CA2*v4*w4 +
     &        34*v5*w4 - 44*CA2*v5*w4 - 27*v6*w4 +
     &        55*CA2*v6*w4 + 4*v7*w4 - 28*CA2*v7*w4 +
     &        v4*w5 - CA2*v4*w5 - 26*v5*w5 +
     &        2*CA2*v5*w5 + 39*v6*w5 - 7*CA2*v6*w5 -
     &        20*v7*w5 + 12*CA2*v7*w5 - 4*v6*w6 +
     &        12*v7*w6 - 8*CA2*v7*w6 - 4*v7*w7 +
     &        4*CA2*v7*w7))/((1 - v)*v*w*(1 - v*w)*
     &     (1 - v + v*w)**3) 
      part9 = -(2*CF*lw*(3 - 3*CA2 + 2*CA*CF - 18*v + 
     &     18*CA2*v - 12*CA*CF*v +
     &        48*v2 - 48*CA2*v2 + 32*CA*CF*v2 - 72*v3 +
     &        72*CA2*v3 - 48*CA*CF*v3 + 63*v4 - 63*CA2*v4 +
     &        42*CA*CF*v4 - 30*v5 + 30*CA2*v5 - 20*CA*CF*v5 +
     &        6*v6 - 6*CA2*v6 + 4*CA*CF*v6 - 10*w - 2*CA2*w +
     &        8*CA*CF*w + 54*v*w - 28*CA*CF*v*w - 129*v2*w +
     &        23*CA2*v2*w + 34*CA*CF*v2*w + 173*v3*w -
     &        55*CA2*v3*w - 10*CA*CF*v3*w - 129*v4*w +
     &        51*CA2*v4*w - 10*CA*CF*v4*w + 39*v5*w -
     &        11*CA2*v5*w + 2*CA*CF*v5*w + 
     &     8*v6*w - 12*CA2*v6*w +
     &        8*CA*CF*v6*w - 6*v7*w + 6*CA2*v7*w - 
     &     4*CA*CF*v7*w +
     &        w2 - CA2*w2 + 2*CA*CF*w2 - 14*v*w2 + 
     &     8*CA2*v*w2 -
     &        12*CA*CF*v*w2 + 46*v2*w2 - 30*CA2*v2*w2 +
     &        40*CA*CF*v2*w2 - 69*v3*w2 + 49*CA2*v3*w2 -
     &        66*CA*CF*v3*w2 + 30*v4*w2 - 22*CA2*v4*w2 +
     &        48*CA*CF*v4*w2 + 49*v5*w2 - 35*CA2*v5*w2 +
     &        2*CA*CF*v5*w2 - 65*v6*w2 + 49*CA2*v6*w2 -
     &        26*CA*CF*v6*w2 + 22*v7*w2 - 18*CA2*v7*w2 +
     &        12*CA*CF*v7*w2 + 2*v*w3 - 2*CA2*v*w3 +
     &        4*CA*CF*v*w3 - 15*v2*w3 + 13*CA2*v2*w3 -
     &        22*CA*CF*v2*w3 + 35*v3*w3 - 27*CA2*v3*w3 +
     &        42*CA*CF*v3*w3 + 4*v4*w3 + 4*CA2*v4*w3 -
     &        20*CA*CF*v4*w3 - 102*v5*w3 + 52*CA2*v5*w3 -
     &        28*CA*CF*v5*w3 + 114*v6*w3 - 66*CA2*v6*w3 +
     &        44*CA*CF*v6*w3 - 38*v7*w3 + 26*CA2*v7*w3 -
     &        20*CA*CF*v7*w3 - 13*v3*w4 - 5*CA2*v3*w4 +
     &        14*CA*CF*v3*w4 - 14*v4*w4 + 34*CA2*v4*w4 -
     &        60*CA*CF*v4*w4 + 114*v5*w4 - 72*CA2*v5*w4 +
     &        84*CA*CF*v5*w4 - 132*v6*w4 + 72*CA2*v6*w4 -
     &        64*CA*CF*v6*w4 + 46*v7*w4 - 30*CA2*v7*w4 +
     &        28*CA*CF*v7*w4 - 2*v3*w5 + 2*CA2*v3*w5 -
     &        4*CA*CF*v3*w5 + 27*v4*w5 - 17*CA2*v4*w5 +
     &        26*CA*CF*v4*w5 - 97*v5*w5 + 39*CA2*v5*w5 -
     &        38*CA*CF*v5*w5 + 110*v6*w5 - 46*CA2*v6*w5 +
     &        36*CA*CF*v6*w5 - 42*v7*w5 + 26*CA2*v7*w5 -
     &        28*CA*CF*v7*w5 - v4*w6 + CA2*v4*w6 -
     &        2*CA*CF*v4*w6 + 27*v5*w6 - 3*CA2*v5*w6 -
     &        2*CA*CF*v5*w6 - 45*v6*w6 + 9*CA2*v6*w6 -
     &        2*CA*CF*v6*w6 + 26*v7*w6 - 14*CA2*v7*w6 +
     &        20*CA*CF*v7*w6 + 4*v6*w7 - 10*v7*w7 +
     &        6*CA2*v7*w7 - 12*CA*CF*v7*w7 + 2*v7*w8 -
     &        2*CA2*v7*w8 + 4*CA*CF*v7*w8))/
     &    ((1 - v)*v*(1 - w)*w*(1 - v*w)*(1 - v + v*w)**3)
      struv5= part1 + part2 + part3 + part4 +
     &        part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV6(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (4*CF*l1v*(4*CA - 12*CA*v + 14*CA*v**2 -
     &     8*CA*v3 + 2*CA*v4 -
     &        12*CA*w + 4*CA**3*w + 32*CA*v*w + CA2*v*w - 
     &     14*CA**3*v*w +
     &        4*v2*w - 32*CA*v2*w + 19*CA**3*v2*w - 7*v3*w +
     &        14*CA*v3*w - 2*CA2*v3*w - 12*CA**3*v3*w + 
     &     4*v4*w +
     &        2*CA2*v4*w + 3*CA**3*v4*w - v5*w - 2*CA*v5*w -
     &        CA2*v5*w + 12*CA*v*w2 - 4*CA**3*v*w2 - 
     &     5*v2*w2 -
     &        32*CA*v2*w2 - 4*CA2*v2*w2 + 
     &     13*CA**3*v2*w2 +
     &        8*v3*w2 + 34*CA*v3*w2 + 5*CA2*v3*w2 -
     &        13*CA**3*v3*w2 - 5*v4*w2 - 20*CA*v4*w2 -
     &        3*CA2*v4*w2 + 10*CA**3*v4*w2 + 3*v5*w2 +
     &        8*CA*v5*w2 + CA2*v5*w2 - 2*CA**3*v5*w2 -
     &        v6*w2 + CA2*v6*w2 + v3*w3 - 
     &     4*CA*v3*w3 -
     &        CA2*v3*w3 + CA**3*v3*w3 + 4*CA*v4*w3 +
     &        2*CA2*v4*w3 + CA**3*v4*w3 - 3*v5*w3 +
     &        CA2*v5*w3 - 2*CA**3*v5*w3 + 2*v6*w3 -
     &        2*CA2*v6*w3 + v5*w4 - 2*CA*v5*w4 -
     &        CA2*v5*w4 - v6*w4 + CA2*v6*w4))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v*w)) 
      part2 = -(4*CF*lvw*(4*CA - 4*CA**3 - 8*CA*v + 8*CA**3*v + 
     &   6*CA*v2 -
     &        6*CA**3*v2 - 2*CA*v3 + 2*CA**3*v3 - 
     &     v*w - 8*CA*v*w +
     &        6*CA2*v*w + 6*CA**3*v*w + 3*v2*w + 15*CA*v2*w -
     &        11*CA2*v2*w - 11*CA**3*v2*w - 
     &     3*v3*w - 8*CA*v3*w +
     &        8*CA2*v3*w + 6*CA**3*v3*w + v4*w + 3*CA*v4*w -
     &        3*CA2*v4*w - 3*CA**3*v4*w - CA*v2*w2 +
     &        5*CA2*v2*w2 - 5*CA**3*v2*w2 - 4*CA*v3*w2 -
     &        11*CA2*v3*w2 + 7*CA**3*v3*w2 + 2*CA*v4*w2 +
     &        6*CA2*v4*w2 - 3*CA**3*v4*w2 - CA*v5*w2 +
     &        CA**3*v5*w2 - v3*w3 - 2*CA*v3*w3 +
     &        3*CA2*v3*w3 - 3*CA**3*v3*w3 + v4*w3 +
     &        3*CA*v4*w3 - 3*CA2*v4*w3 + 2*CA**3*v4*w3 +
     &        3*CA*v5*w3 - 3*CA**3*v5*w3 - 4*CA**3*v4*w4 -
     &        4*CA*v5*w4 + 4*CA**3*v5*w4 + 2*CA*v5*w5 -
     &        2*CA**3*v5*w5))/(CA*(1 - v)**2*v3*w2) 
      part3 = -(8*CF**2*lmss*(2*CA - 10*CA*v + 21*CA*v2 - 24*CA*v3 +
     &        16*CA*v4 - 6*CA*v5 + CA*v6 + 2*CA*w - 
     &     2*v*w - 6*CA*v*w +
     &        8*v2*w + 3*CA*v2*w - 13*v3*w + 
     &     9*CA*v3*w + 11*v4*w -
     &        15*CA*v4*w - 5*v5*w + 9*CA*v5*w + 
     &     v6*w - 2*CA*v6*w +
     &        2*CA*v*w2 - 2*v2*w2 - 
     &     6*CA*v2*w2 + 9*v3*w2 +
     &        3*CA*v3*w2 - 15*v4*w2 + 4*CA*v4*w2 + 
     &     11*v5*w2 -
     &        5*CA*v5*w2 - 3*v6*w2 + 2*CA*v6*w2 + 
     &     4*v4*w3 -
     &        5*CA*v4*w3 - 7*v5*w3 + 5*CA*v5*w3 + 
     &     3*v6*w3 -
     &        2*CA*v6*w3 + v5*w4 - 3*CA*v5*w4 - 
     &     v6*w4 +
     &        2*CA*v6*w4 - CA*v6*w5))/
     &    ((1 - v)**2*v3*w2*(1 - v + v*w)) 
      part4 = (4*CF*l1vw*(4*CA - 4*CA**3 + 8*CA2*CF - 20*CA*v + 
     &     20*CA**3*v -
     &        40*CA2*CF*v + 42*CA*v2 - 42*CA**3*v2 + 
     &     84*CA2*CF*v2 -
     &        48*CA*v3 + 48*CA**3*v3 - 96*CA2*CF*v3 + 
     &     32*CA*v4 -
     &        32*CA**3*v4 + 64*CA2*CF*v4 - 12*CA*v5 + 
     &     12*CA**3*v5 -
     &        24*CA2*CF*v5 + 2*CA*v6 - 2*CA**3*v6 + 
     &     4*CA2*CF*v6 -
     &        4*v*w + 4*CA*v*w + 4*CA2*v*w - 8*CA**3*v*w - 
     &     8*CA*CF*v*w +
     &        16*CA2*CF*v*w + 16*v2*w - 22*CA*v2*w - 
     &     16*CA2*v2*w +
     &        35*CA**3*v2*w + 32*CA*CF*v2*w - 72*CA2*CF*v2*w -
     &        26*v3*w + 48*CA*v3*w + 26*CA2*v3*w - 
     &     63*CA**3*v3*w -
     &        52*CA*CF*v3*w + 132*CA2*CF*v3*w + 22*v4*w -
     &        52*CA*v4*w - 22*CA2*v4*w + 59*CA**3*v4*w +
     &        44*CA*CF*v4*w - 124*CA2*CF*v4*w - 10*v5*w +
     &        28*CA*v5*w + 10*CA2*v5*w - 29*CA**3*v5*w -
     &        20*CA*CF*v5*w + 60*CA2*CF*v5*w + 2*v6*w - 
     &     6*CA*v6*w -
     &        2*CA2*v6*w + 6*CA**3*v6*w + 4*CA*CF*v6*w -
     &        12*CA2*CF*v6*w - 12*v2*w2 + 4*CA*v2*w2 +
     &        12*CA2*v2*w2 - 9*CA**3*v2*w2 - 
     &     24*CA*CF*v2*w2 +
     &        20*CA2*CF*v2*w2 + 40*v3*w2 - 24*CA*v3*w2 -
     &        40*CA2*v3*w2 + 33*CA**3*v3*w2 + 
     &     80*CA*CF*v3*w2 -
     &        72*CA2*CF*v3*w2 - 52*v4*w2 + 44*CA*v4*w2 +
     &        52*CA2*v4*w2 - 47*CA**3*v4*w2 - 
     &     104*CA*CF*v4*w2 +
     &        100*CA2*CF*v4*w2 + 32*v5*w2 - 32*CA*v5*w2 -
     &        32*CA2*v5*w2 + 31*CA**3*v5*w2 + 
     &     64*CA*CF*v5*w2 -
     &        64*CA2*CF*v5*w2 - 8*v6*w2 + 8*CA*v6*w2 +
     &        8*CA2*v6*w2 - 8*CA**3*v6*w2 - 
     &     16*CA*CF*v6*w2 +
     &        16*CA2*CF*v6*w2 - 14*v3*w3 + 8*CA*v3*w3 +
     &        14*CA2*v3*w3 - 10*CA**3*v3*w3 - 
     &     28*CA*CF*v3*w3 +
     &        20*CA2*CF*v3*w3 + 38*v4*w3 - 24*CA*v4*w3 -
     &        38*CA2*v4*w3 + 25*CA**3*v4*w3 + 
     &     76*CA*CF*v4*w3 -
     &        52*CA2*CF*v4*w3 - 36*v5*w3 + 24*CA*v5*w3 +
     &        36*CA2*v5*w3 - 23*CA**3*v5*w3 - 
     &     72*CA*CF*v5*w3 +
     &        48*CA2*CF*v5*w3 + 12*v6*w3 - 8*CA*v6*w3 -
     &        12*CA2*v6*w3 + 8*CA**3*v6*w3 + 
     &     24*CA*CF*v6*w3 -
     &        16*CA2*CF*v6*w3 - 8*v4*w4 + 8*CA*v4*w4 +
     &        8*CA2*v4*w4 - 9*CA**3*v4*w4 - 
     &     16*CA*CF*v4*w4 +
     &        20*CA2*CF*v4*w4 + 16*v5*w4 - 14*CA*v5*w4 -
     &        16*CA2*v5*w4 + 15*CA**3*v5*w4 + 
     &     32*CA*CF*v5*w4 -
     &        32*CA2*CF*v5*w4 - 8*v6*w4 + 8*CA*v6*w4 +
     &        8*CA2*v6*w4 - 8*CA**3*v6*w4 - 
     &     16*CA*CF*v6*w4 +
     &        16*CA2*CF*v6*w4 - 2*v5*w5 + 6*CA*v5*w5 +
     &        2*CA2*v5*w5 - 6*CA**3*v5*w5 - 
     &     4*CA*CF*v5*w5 +
     &        12*CA2*CF*v5*w5 + 2*v6*w5 - 6*CA*v6*w5 -
     &        2*CA2*v6*w5 + 6*CA**3*v6*w5 + 
     &     4*CA*CF*v6*w5 -
     &        12*CA2*CF*v6*w5 + 2*CA*v6*w6 - 
     &     2*CA**3*v6*w6 +
     &        4*CA2*CF*v6*w6))/(CA*(1 - v)**2*v3*w2*
     &     (1 - v + v*w))
      part5 = -(2*CF*lms*(8*CA2 + 8*CA*CF - 24*CA2*v - 
     &     16*CA*CF*v - 4*v2 +
     &        32*CA2*v2 + 12*CA*CF*v2 + 8*v3 - 24*CA2*v3 -
     &        4*CA*CF*v3 - 6*v4 + 10*CA2*v4 + 
     &     2*v5 - 2*CA2*v5 -
     &        8*CA2*w + 12*CA*CF*w - 12*CF*v*w - 
     &     60*CA*CF*v*w + 3*v2*w +
     &        41*CA2*v2*w + 20*CF*v2*w + 
     &     84*CA*CF*v2*w + 4*v3*w -
     &        72*CA2*v3*w - 12*CF*v3*w - 
     &     48*CA*CF*v3*w - 19*v4*w +
     &        63*CA2*v4*w + 4*CF*v4*w + 
     &     16*CA*CF*v4*w + 16*v5*w -
     &        28*CA2*v5*w - 6*v6*w + 
     &     6*CA2*v6*w + 4*CA2*w2 -
     &        4*CA*CF*w2 + 12*CA2*v*w2 + 4*CF*v*w2 -
     &        24*CA*CF*v*w2 - 2*v2*w2 - 32*CA2*v2*w2 +
     &        32*CF*v2*w2 + 116*CA*CF*v2*w2 - 4*v3*w2 +
     &        10*CA2*v3*w2 - 48*CF*v3*w2 - 
     &     156*CA*CF*v3*w2 +
     &        11*v4*w2 + 27*CA2*v4*w2 + 16*CF*v4*w2 +
     &        72*CA*CF*v4*w2 + 9*v5*w2 - 43*CA2*v5*w2 -
     &        4*CF*v5*w2 - 28*CA*CF*v5*w2 - 12*v6*w2 +
     &        24*CA2*v6*w2 + 6*v7*w2 - 6*CA2*v7*w2 -
     &        12*CA2*v*w3 + 12*CA*CF*v*w3 + 12*CA2*v2*w3 -
     &        12*CF*v2*w3 + 6*v3*w3 + 16*CA2*v3*w3 -
     &        24*CF*v3*w3 - 68*CA*CF*v3*w3 - 8*v4*w3 -
     &        32*CA2*v4*w3 + 40*CF*v4*w3 + 
     &     104*CA*CF*v4*w3 -
     &        15*v5*w3 + 25*CA2*v5*w3 - 4*CF*v5*w3 -
     &        28*CA*CF*v5*w3 + 5*v6*w3 - CA2*v6*w3 +
     &        20*CA*CF*v6*w3 - 4*CA2*v7*w3 - 2*v8*w3 +
     &        2*CA2*v8*w3 + 12*CA2*v2*w4 - 
     &     12*CA*CF*v2*w4 -
     &        28*CA2*v3*w4 + 12*CF*v3*w4 + 
     &     24*CA*CF*v3*w4 -
     &        6*v4*w4 + 24*CA2*v4*w4 - 12*CA*CF*v4*w4 +
     &        18*v5*w4 - 8*CA2*v5*w4 - 12*CF*v5*w4 -
     &        8*CA*CF*v5*w4 + 3*v6*w4 - 13*CA2*v6*w4 -
     &        8*CA*CF*v6*w4 - v7*w4 + 7*CA2*v7*w4 -
     &        8*CA*CF*v7*w4 + 2*v8*w4 - 2*CA2*v8*w4 -
     &        4*CA2*v3*w5 + 4*CA*CF*v3*w5 + 
     &     12*CA2*v4*w5 -
     &        4*CF*v4*w5 - 12*CA*CF*v4*w5 + 2*v5*w5 -
     &        16*CA2*v5*w5 + 4*CF*v5*w5 + 
     &     16*CA*CF*v5*w5 -
     &        13*v6*w5 + 17*CA2*v6*w5 - 12*CA*CF*v6*w5 -
     &        v7*w5 - 5*CA2*v7*w5 + 8*CA*CF*v7*w5 -
     &        2*v8*w5 + 2*CA2*v8*w5 + 6*v7*w6 -
     &        2*CA2*v7*w6 + 2*v8*w6 - 2*CA2*v8*w6 -
     &        2*v8*w7 + 2*CA2*v8*w7))/
     &    ((1 - v)**2*v3*w2*(1 - v*w)**3) 
      part6 = -(2*CF*l1w*(16*CA - 16*CA**3 - 64*CA*v + 
     &     64*CA**3*v + 112*CA*v2 -
     &        112*CA**3*v2 - 112*CA*v3 + 112*CA**3*v3 + 
     &     68*CA*v4 -
     &        68*CA**3*v4 - 24*CA*v5 + 24*CA**3*v5 + 4*CA*v6 -
     &        4*CA**3*v6 + 18*CA*w - 10*CA**3*w - 18*v*w - 
     &     120*CA*v*w +
     &        12*CA2*v*w + 76*CA**3*v*w + 64*v2*w + 303*CA*v2*w -
     &        34*CA2*v2*w - 213*CA**3*v2*w - 92*v3*w -
     &        397*CA*v3*w + 38*CA2*v3*w + 311*CA**3*v3*w +
     &        68*v4*w + 323*CA*v4*w - 22*CA2*v4*w -
     &        285*CA**3*v4*w - 26*v5*w - 179*CA*v5*w + 
     &     6*CA2*v5*w +
     &        173*CA**3*v5*w + 4*v6*w + 
     &     64*CA*v6*w - 64*CA**3*v6*w -
     &        12*CA*v7*w + 12*CA**3*v7*w - 2*CA*w2 - 
     &     2*CA**3*w2 +
     &        2*v*w2 - 28*CA*v*w2 - 2*CA2*v*w2 + 
     &     28*CA**3*v*w2 +
     &        30*v2*w2 + 182*CA*v2*w2 - 26*CA2*v2*w2 -
     &        110*CA**3*v2*w2 - 114*v3*w2 - 451*CA*v3*w2 +
     &        82*CA2*v3*w2 + 245*CA**3*v3*w2 + 
     &     158*v4*w2 +
     &        515*CA*v4*w2 - 86*CA2*v4*w2 - 
     &     281*CA**3*v4*w2 -
     &        116*v5*w2 - 313*CA*v5*w2 + 48*CA2*v5*w2 +
     &        199*CA**3*v5*w2 + 48*v6*w2 + 137*CA*v6*w2 -
     &        16*CA2*v6*w2 - 119*CA**3*v6*w2 - 8*v7*w2 -
     &        48*CA*v7*w2 + 48*CA**3*v7*w2 + 12*CA*v8*w2 -
     &        12*CA**3*v8*w2 + 4*CA*v*w3 + 4*CA**3*v*w3 -
     &        4*v2*w3 - 18*CA*v2*w3 + 4*CA2*v2*w3 -
     &        18*CA**3*v2*w3 + 10*v3*w3 + 8*CA*v3*w3 -
     &        6*CA2*v3*w3 + 28*CA**3*v3*w3 + 20*v4*w3 +
     &        120*CA*v4*w3 - 32*CA2*v4*w3 - 
     &     52*CA**3*v4*w3 -
     &        42*v5*w3 - 212*CA*v5*w3 + 44*CA2*v5*w3 +
     &        52*CA**3*v5*w3 + 32*v6*w3 + 89*CA*v6*w3 -
     &        26*CA2*v6*w3 + 9*CA**3*v6*w3 - 20*v7*w3 -
     &        15*CA*v7*w3 + 16*CA2*v7*w3 + CA**3*v7*w3 +
     &        4*v8*w3 - 4*CA*v9*w3 + 4*CA**3*v9*w3 +
     &        42*CA*v3*w4 - 14*CA**3*v3*w4 - 42*v4*w4 -
     &        156*CA*v4*w4 + 38*CA2*v4*w4 + 
     &     56*CA**3*v4*w4 +
     &        58*v5*w4 + 210*CA*v5*w4 - 42*CA2*v5*w4 -
     &        94*CA**3*v5*w4 - 34*v6*w4 - 72*CA*v6*w4 +
     &        22*CA2*v6*w4 + 48*CA**3*v6*w4 + 20*v7*w4 +
     &        21*CA*v7*w4 - 12*CA2*v7*w4 - 
     &     47*CA**3*v7*w4 -
     &        2*v8*w4 - CA*v8*w4 - 6*CA2*v8*w4 +
     &        7*CA**3*v8*w4 + 8*CA*v9*w4 - 8*CA**3*v9*w4 -
     &        4*CA*v3*w5 - 4*CA**3*v3*w5 + 4*v4*w5 -
     &        8*CA*v4*w5 - 4*CA2*v4*w5 + 20*CA**3*v4*w5 +
     &        20*v5*w5 + 64*CA*v5*w5 - 14*CA2*v5*w5 -
     &        40*CA**3*v5*w5 - 32*v6*w5 - 113*CA*v6*w5 +
     &        18*CA2*v6*w5 + 47*CA**3*v6*w5 + 12*v7*w5 +
     &        57*CA*v7*w5 - 12*CA2*v7*w5 - 
     &     7*CA**3*v7*w5 -
     &        4*v8*w5 - 24*CA*v8*w5 + 12*CA2*v8*w5 +
     &        12*CA**3*v8*w5 - 8*CA*v9*w5 + 8*CA**3*v9*w5 +
     &        2*CA*v4*w6 + 2*CA**3*v4*w6 - 2*v5*w6 -
     &        6*CA*v5*w6 + 2*CA2*v5*w6 - 6*CA**3*v5*w6 +
     &        4*v6*w6 + 6*CA*v6*w6 - 4*CA2*v6*w6 +
     &        6*CA**3*v6*w6 - 13*CA*v7*w6 + 8*CA2*v7*w6 +
     &        3*CA**3*v7*w6 - 2*v8*w6 + 15*CA*v8*w6 -
     &        6*CA2*v8*w6 - 9*CA**3*v8*w6 + 8*CA*v9*w6 -
     &        8*CA**3*v9*w6 - 4*v7*w7 + 4*v8*w7 -
     &        8*CA*v9*w7 + 8*CA**3*v9*w7 + 4*CA*v9*w8 -
     &        4*CA**3*v9*w8))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v*w)**3*(1 - v + v*w)) 
      part7 = -(2*CF*lv*(24*CA - 16*CA**3 - 96*CA*v + 
     &     64*CA**3*v + 164*CA*v2 -
     &        112*CA**3*v2 - 156*CA*v3 + 
     &     112*CA**3*v3 + 88*CA*v4 -
     &        68*CA**3*v4 - 28*CA*v5 + 24*CA**3*v5 + 4*CA*v6 -
     &        4*CA**3*v6 + 26*CA*w - 2*CA**3*w - 
     &     26*v*w - 152*CA*v*w +
     &        10*CA2*v*w + 40*CA**3*v*w + 88*v2*w + 367*CA*v2*w -
     &        32*CA2*v2*w - 147*CA**3*v2*w - 122*v3*w -
     &        485*CA*v3*w + 42*CA2*v3*w + 249*CA**3*v3*w +
     &        90*v4*w + 407*CA*v4*w - 30*CA2*v4*w -
     &        255*CA**3*v4*w - 36*v5*w - 227*CA*v5*w +
     &        12*CA2*v5*w + 167*CA**3*v5*w + 
     &     6*v6*w + 76*CA*v6*w -
     &        2*CA2*v6*w - 64*CA**3*v6*w - 12*CA*v7*w +
     &        12*CA**3*v7*w - 2*CA*w2 - 2*CA**3*w2 + 2*v*w2 -
     &        44*CA*v*w2 - 2*CA2*v*w2 + 12*CA**3*v*w2 +
     &        48*v2*w2 + 222*CA*v2*w2 - 16*CA2*v2*w2 -
     &        32*CA**3*v2*w2 - 164*v3*w2 - 487*CA*v3*w2 +
     &        60*CA2*v3*w2 + 99*CA**3*v3*w2 + 210*v4*w2 +
     &        511*CA*v4*w2 - 74*CA2*v4*w2 - 
     &     135*CA**3*v4*w2 -
     &        152*v5*w2 - 317*CA*v5*w2 + 52*CA2*v5*w2 +
     &        121*CA**3*v5*w2 + 70*v6*w2 + 169*CA*v6*w2 -
     &        26*CA2*v6*w2 - 103*CA**3*v6*w2 - 14*v7*w2 -
     &        60*CA*v7*w2 + 6*CA2*v7*w2 + 
     &     48*CA**3*v7*w2 +
     &        12*CA*v8*w2 - 12*CA**3*v8*w2 + 4*CA*v*w3 +
     &        4*CA**3*v*w3 - 4*v2*w3 - 18*CA*v2*w3 +
     &        4*CA2*v2*w3 - 18*CA**3*v2*w3 + 10*v3*w3 -
     &        10*CA2*v3*w3 + 8*CA**3*v3*w3 + 30*v4*w3 +
     &        108*CA*v4*w3 - 10*CA2*v4*w3 + 
     &     16*CA**3*v4*w3 -
     &        40*v5*w3 - 140*CA*v5*w3 + 16*CA2*v5*w3 -
     &        40*CA**3*v5*w3 + 22*v6*w3 + 5*CA*v6*w3 -
     &        6*CA2*v6*w3 + 75*CA**3*v6*w3 - 28*v7*w3 -
     &        15*CA*v7*w3 + 12*CA2*v7*w3 - 
     &     13*CA**3*v7*w3 +
     &        10*v8*w3 + 4*CA*v8*w3 - 6*CA2*v8*w3 -
     &        4*CA*v9*w3 + 4*CA**3*v9*w3 + 58*CA*v3*w4 +
     &        2*CA**3*v3*w4 - 58*v4*w4 - 172*CA*v4*w4 +
     &        26*CA2*v4*w4 + 8*CA**3*v4*w4 + 74*v5*w4 +
     &        226*CA*v5*w4 - 30*CA2*v5*w4 - 
     &     54*CA**3*v5*w4 -
     &        56*v6*w4 - 92*CA*v6*w4 + 20*CA2*v6*w4 +
     &        38*CA**3*v6*w4 + 54*v7*w4 + 85*CA*v7*w4 -
     &        26*CA2*v7*w4 - 65*CA**3*v7*w4 - 12*v8*w4 -
     &        5*CA*v8*w4 + 8*CA2*v8*w4 + 11*CA**3*v8*w4 -
     &        2*v9*w4 + 8*CA*v9*w4 + 2*CA2*v9*w4 -
     &        8*CA**3*v9*w4 - 4*CA*v3*w5 - 4*CA**3*v3*w5 +
     &        4*v4*w5 - 16*CA*v4*w5 - 4*CA2*v4*w5 +
     &        12*CA**3*v4*w5 + 24*v5*w5 + 72*CA*v5*w5 -
     &        8*CA2*v5*w5 - 16*CA**3*v5*w5 - 20*v6*w5 -
     &        97*CA*v6*w5 + 8*CA2*v6*w5 + 
     &     17*CA**3*v6*w5 -
     &        6*v7*w5 + 25*CA*v7*w5 + 6*CA2*v7*w5 +
     &        15*CA**3*v7*w5 - 8*v8*w5 - 40*CA*v8*w5 +
     &        4*CA2*v8*w5 + 12*CA**3*v8*w5 + 6*v9*w5 -
     &        8*CA*v9*w5 - 6*CA2*v9*w5 + 8*CA**3*v9*w5 +
     &        2*CA*v4*w6 + 2*CA**3*v4*w6 - 2*v5*w6 -
     &        6*CA*v5*w6 + 2*CA2*v5*w6 - 6*CA**3*v5*w6 +
     &        2*v6*w6 + 6*CA*v6*w6 - 2*CA2*v6*w6 +
     &        8*CA**3*v6*w6 - 6*v7*w6 - 25*CA*v7*w6 +
     &        2*CA2*v7*w6 + 5*CA**3*v7*w6 + 12*v8*w6 +
     &        35*CA*v8*w6 - 8*CA2*v8*w6 - 
     &     13*CA**3*v8*w6 -
     &        6*v9*w6 + 8*CA*v9*w6 + 6*CA2*v9*w6 -
     &        8*CA**3*v9*w6 - 2*v8*w7 - 4*CA*v8*w7 +
     &        2*CA2*v8*w7 + 2*v9*w7 - 8*CA*v9*w7 -
     &        2*CA2*v9*w7 + 8*CA**3*v9*w7 + 4*CA*v9*w8 -
     &        4*CA**3*v9*w8))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v*w)**3*(1 - v + v*w)) 
      part8 = -(CF*(8*CA*v2 - 8*CA**3*v2 - 24*CA*v3 + 24*CA**3*v3 +
     &   28*CA*v4 - 28*CA**3*v4 - 16*CA*v5 + 16*CA**3*v5 +
     &   4*CA*v6 - 4*CA**3*v6 + 6*CA*w - 14*CA**3*w - 24*v*w -
     &   24*CA*v*w + 24*CA2*v*w + 56*CA**3*v*w - 40*CA*CF*v*w +
     &   76*v2*w + 39*CA*v2*w - 76*CA2*v2*w - 91*CA**3*v2*w + 
     &   136*CA*CF*v2*w + 8*CA2*CF*v2*w - 108*v3*w - 45*CA*v3*w + 
     &   108*CA2*v3*w + 89*CA**3*v3*w - 192*CA*CF*v3*w - 
     &   24*CA2*CF*v3*w + 92*v4*w + 63*CA*v4*w - 92*CA2*v4*w - 
     &   83*CA**3*v4*w + 144*CA*CF*v4*w + 24*CA2*CF*v4*w - 44*v5*w -
     &   67*CA*v5*w + 44*CA2*v5*w + 71*CA**3*v5*w - 56*CA*CF*v5*w - 
     &   8*CA2*CF*v5*w + 8*v6*w + 40*CA*v6*w - 8*CA2*v6*w - 
     &   40*CA**3*v6*w + 8*CA*CF*v6*w - 12*CA*v7*w + 12*CA**3*v7*w + 
     &   8*CA*w2 + 8*CA**3*w2 - 2*v*w2 - 76*CA*v*w2 - 10*CA2*v*w2 + 
     &   4*CA**3*v*w2 + 8*CA*CF*v*w2 + 64*v2*w2 + 220*CA*v2*w2 -
     &   12*CA2*v2*w2 - 108*CA**3*v2*w2 + 56*CA*CF*v2*w2 -
     &   8*CA2*CF*v2*w2 - 182*v3*w2 - 275*CA*v3*w2 + 102*CA2*v3*w2 + 
     &   223*CA**3*v3*w2 - 232*CA*CF*v3*w2 + 8*CA2*CF*v3*w2 + 
     &   220*v4*w2 + 179*CA*v4*w2 - 172*CA2*v4*w2 - 207*CA**3*v4*w2 +
     &   344*CA*CF*v4*w2 + 48*CA2*CF*v4*w2 - 180*v5*w2 -
     &   99*CA*v5*w2 + 176*CA2*v5*w2 + 127*CA**3*v5*w2 -
     &   296*CA*CF*v5*w2 - 88*CA2*CF*v5*w2 + 104*v6*w2 + 59*CA*v6*w2 -
     &   108*CA2*v6*w2 - 63*CA**3*v6*w2 + 144*CA*CF*v6*w2 + 
     &   40*CA2*CF*v6*w2 - 24*v7*w2 - 24*CA*v7*w2 + 24*CA2*v7*w2 + 
     &   24*CA**3*v7*w2 - 24*CA*CF*v7*w2 + 12*CA*v8*w2 - 
     &   12*CA**3*v8*w2 - 16*CA*v*w3 - 16*CA**3*v*w3 + 4*v2*w3 +
     &   136*CA*v2*w3 + 20*CA2*v2*w3 + 56*CA**3*v2*w3 -
     &   16*CA*CF*v2*w3 - 34*v3*w3 - 406*CA*v3*w3 - 82*CA2*v3*w3 - 
     &   10*CA**3*v3*w3 + 56*CA*CF*v3*w3 + 16*CA2*CF*v3*w3 + 
     &   102*v4*w3 + 482*CA*v4*w3 + 90*CA2*v4*w3 - 122*CA**3*v4*w3 - 
     &   32*CA*CF*v4*w3 - 72*CA2*CF*v4*w3 - 70*v5*w3 - 216*CA*v5*w3 -
     &   42*CA2*v5*w3 + 136*CA**3*v5*w3 - 24*CA*CF*v5*w3 +
     &   40*CA2*CF*v5*w3 + 28*v6*w3 + 43*CA*v6*w3 - 36*CA2*v6*w3 - 
     &   71*CA**3*v6*w3 + 104*CA*CF*v6*w3 + 80*CA2*CF*v6*w3 - 
     &   54*v7*w3 - 21*CA*v7*w3 + 74*CA2*v7*w3 + 25*CA**3*v7*w3 - 
     &   112*CA*CF*v7*w3 - 64*CA2*CF*v7*w3 + 24*v8*w3 - 8*CA*v8*w3 -
     &   24*CA2*v8*w3 + 8*CA**3*v8*w3 + 24*CA*CF*v8*w3 - 4*CA*v9*w3 +
     &   4*CA**3*v9*w3 - 12*CA*v3*w4 - 52*CA**3*v3*w4 - 42*v4*w4 + 
     &   138*CA*v4*w4 + 78*CA2*v4*w4 + 118*CA**3*v4*w4 -
     &   104*CA*CF*v4*w4 + 32*v5*w4 - 226*CA*v5*w4 - 156*CA2*v5*w4 - 
     &   70*CA**3*v5*w4 + 192*CA*CF*v5*w4 + 72*CA2*CF*v5*w4 - 
     &   44*v6*w4 + 64*CA*v6*w4 + 168*CA2*v6*w4 - 208*CA*CF*v6*w4 -
     &        128*CA2*CF*v6*w4 + 88*v7*w4 + 35*CA*v7*w4 -
     &        108*CA2*v7*w4 + 25*CA**3*v7*w4 +
     &        120*CA*CF*v7*w4 + 16*CA2*CF*v7*w4 - 
     &     26*v8*w4 +
     &        9*CA*v8*w4 + 10*CA2*v8*w4 - 
     &     29*CA**3*v8*w4 +
     &        8*CA*CF*v8*w4 + 40*CA2*CF*v8*w4 - 8*v9*w4 +
     &        8*CA*v9*w4 + 8*CA2*v9*w4 - 8*CA**3*v9*w4 -
     &        8*CA*CF*v9*w4 + 16*CA*v3*w5 + 
     &     16*CA**3*v3*w5 -
     &        4*v4*w5 - 110*CA*v4*w5 - 20*CA2*v4*w5 -
     &        10*CA**3*v4*w5 + 16*CA*CF*v4*w5 + 50*v5*w5 +
     &        182*CA*v5*w5 + 18*CA2*v5*w5 - 
     &     54*CA**3*v5*w5 +
     &        16*CA*CF*v5*w5 - 16*CA2*CF*v5*w5 - 
     &     34*v6*w5 -
     &        89*CA*v6*w5 + 2*CA2*v6*w5 + 
     &     85*CA**3*v6*w5 -
     &        32*CA*CF*v6*w5 - 24*v7*w5 - 7*CA*v7*w5 -
     &        20*CA2*v7*w5 - 53*CA**3*v7*w5 + 
     &     40*CA*CF*v7*w5 +
     &        64*CA2*CF*v7*w5 - 8*v8*w5 - 2*CA*v8*w5 +
     &        40*CA2*v8*w5 + 26*CA**3*v8*w5 - 
     &     56*CA*CF*v8*w5 -
     &        40*CA2*CF*v8*w5 + 20*v9*w5 - 16*CA*v9*w5 -
     &        20*CA2*v9*w5 + 16*CA**3*v9*w5 + 
     &     16*CA*CF*v9*w5 -
     &        8*CA2*CF*v9*w5 - 8*CA*v4*w6 - 
     &     8*CA**3*v4*w6 +
     &        2*v5*w6 + 56*CA*v5*w6 + 10*CA2*v5*w6 +
     &        16*CA**3*v5*w6 - 8*CA*CF*v5*w6 - 14*v6*w6 -
     &        110*CA*v6*w6 - 26*CA2*v6*w6 - 
     &     2*CA**3*v6*w6 +
     &        16*CA*CF*v6*w6 + 8*CA2*CF*v6*w6 + 
     &     14*v7*w6 +
     &        89*CA*v7*w6 + 30*CA2*v7*w6 - 
     &     13*CA**3*v7*w6 -
     &        24*CA*CF*v7*w6 - 16*CA2*CF*v7*w6 + 
     &     18*v8*w6 -
     &        15*CA*v8*w6 - 34*CA2*v8*w6 - 
     &     5*CA**3*v8*w6 +
     &        32*CA*CF*v8*w6 - 20*v9*w6 + 16*CA*v9*w6 +
     &        20*CA2*v9*w6 - 16*CA**3*v9*w6 - 
     &     16*CA*CF*v9*w6 +
     &        8*CA2*CF*v9*w6 - 8*v8*w7 - 8*CA*v8*w7 +
     &        8*CA2*v8*w7 + 8*CA**3*v8*w7 - 
     &     8*CA*CF*v8*w7 +
     &        8*v9*w7 - 8*CA*v9*w7 - 8*CA2*v9*w7 +
     &        8*CA**3*v9*w7 + 8*CA*CF*v9*w7 + 4*CA*v9*w8 -
     &        4*CA**3*v9*w8))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v*w)**3*(1 - v + v*w)) 
      part9 = -(4*CF*lw*(8*CA - 12*CA**3 + 8*CA2*CF - 32*CA*v + 
     &   48*CA**3*v - 32*CA2*CF*v + 58*CA*v2 - 84*CA**3*v2 + 
     &   56*CA2*CF*v2 - 62*CA*v3 + 84*CA**3*v3 - 56*CA2*CF*v3 + 
     &   41*CA*v4 - 51*CA**3*v4 + 34*CA2*CF*v4 - 16*CA*v5 + 
     &   18*CA**3*v5 - 12*CA2*CF*v5 + 3*CA*v6 - 3*CA**3*v6 + 
     &   2*CA2*CF*v6 + CA*w + 3*CA**3*w + 2*CA2*CF*w - 10*v*w - 
     &   20*CA*v*w + 10*CA2*v*w + 12*CA**3*v*w - 10*CA*CF*v*w - 
     &   24*CA2*CF*v*w + 32*v2*w + 64*CA*v2*w - 32*CA2*v2*w - 
     &   74*CA**3*v2*w + 32*CA*CF*v2*w + 74*CA2*CF*v2*w - 42*v3*w -
     &   99*CA*v3*w + 42*CA2*v3*w + 141*CA**3*v3*w -
     &   42*CA*CF*v3*w - 112*CA2*CF*v3*w + 30*v4*w + 108*CA*v4*w - 
     &   30*CA2*v4*w - 156*CA**3*v4*w + 30*CA*CF*v4*w + 
     &   110*CA2*CF*v4*w - 12*v5*w - 84*CA*v5*w + 12*CA2*v5*w + 
     &   110*CA**3*v5*w - 12*CA*CF*v5*w - 74*CA2*CF*v5*w + 2*v6*w +
     &   39*CA*v6*w - 2*CA2*v6*w - 45*CA**3*v6*w +
     &   2*CA*CF*v6*w + 30*CA2*CF*v6*w - 9*CA*v7*w + 9*CA**3*v7*w - 
     &   6*CA2*CF*v7*w - 10*CA*v*w2 - 3*CA2*v*w2 - 4*CA**3*v*w2 - 
     &   4*CA2*CF*v*w2 + 11*v2*w2 + 54*CA*v2*w2 + CA2*v2*w2 +
     &   5*CA**3*v2*w2 + 8*CA*CF*v2*w2 + 30*CA2*CF*v2*w2 -
     &   38*v3*w2 - 98*CA*v3*w2 + 22*CA2*v3*w2 + 13*CA**3*v3*w2 - 
     &   30*CA*CF*v3*w2 - 70*CA2*CF*v3*w2 + 52*v4*w2 + 58*CA*v4*w2 -
     &   44*CA2*v4*w2 - 14*CA**3*v4*w2 + 44*CA*CF*v4*w2 +
     &   62*CA2*CF*v4*w2 - 44*v5*w2 - 11*CA*v5*w2 + 43*CA2*v5*w2 +
     &   16*CA**3*v5*w2 - 40*CA*CF*v5*w2 - 32*CA2*CF*v5*w2 +
     &   25*v6*w2 + 22*CA*v6*w2 - 25*CA2*v6*w2 - 37*CA**3*v6*w2 + 
     &   24*CA*CF*v6*w2 + 28*CA2*CF*v6*w2 - 6*v7*w2 - 21*CA*v7*w2 +
     &   6*CA2*v7*w2 + 27*CA**3*v7*w2 - 6*CA*CF*v7*w2 -
     &   18*CA2*CF*v7*w2 + 9*CA*v8*w2 - 9*CA**3*v8*w2 +
     &   6*CA2*CF*v8*w2 + CA*w3 - CA**3*w3 + 2*CA2*CF*w3 -
     &   v*w3 - 4*CA*v*w3 + CA2*v*w3 + 4*CA**3*v*w3 -
     &   2*CA*CF*v*w3 - 8*CA2*CF*v*w3 + v2*w3 +
     &   23*CA*v2*w3 + CA2*v2*w3 - 10*CA**3*v2*w3 + 4*CA*CF*v2*w3 +
     &   14*CA2*CF*v2*w3 + 9*v3*w3 - 77*CA*v3*w3 - 23*CA2*v3*w3 + 
     &   8*CA**3*v3*w3 + 10*CA*CF*v3*w3 - 18*CA2*CF*v3*w3 - 
     &   18*v4*w3 + 142*CA*v4*w3 + 44*CA2*v4*w3 - 
     &   16*CA**3*v4*w3 - 28*CA*CF*v4*w3 + 40*CA2*CF*v4*w3 + 
     &   27*v5*w3 - 84*CA*v5*w3 - 43*CA2*v5*w3 - 18*CA**3*v5*w3 +
     &   36*CA*CF*v5*w3 - 22*CA2*CF*v5*w3 - 17*v6*w3 - 32*CA*v6*w3 +
     &   19*CA2*v6*w3 + 72*CA**3*v6*w3 - 20*CA*CF*v6*w3 - 
     &   30*CA2*CF*v6*w3 - 7*v7*w3 + 34*CA*v7*w3 + 7*CA2*v7*w3 - 
     &   42*CA**3*v7*w3 - 6*CA*CF*v7*w3 + 24*CA2*CF*v7*w3 + 6*v8*w3 -
     &   11*CA*v8*w3 - 6*CA2*v8*w3 + 9*CA**3*v8*w3 +
     &   6*CA*CF*v8*w3 - 6*CA2*CF*v8*w3 - 3*CA*v9*w3 +
     &   3*CA**3*v9*w3 - 2*CA2*CF*v9*w3 - 2*CA*v*w4 + 2*CA**3*v*w4 -
     &   4*CA2*CF*v*w4 + 2*v2*w4 + 9*CA*v2*w4 - 2*CA2*v2*w4 - 
     &   9*CA**3*v2*w4 + 4*CA*CF*v2*w4 + 18*CA2*CF*v2*w4 - 5*v3*w4 -
     &   15*CA*v3*w4 + 7*CA2*v3*w4 + 21*CA**3*v3*w4 -
     &   10*CA*CF*v3*w4 - 30*CA2*CF*v3*w4 - v4*w4 + 25*CA*v4*w4 +
     &   CA2*v4*w4 - 17*CA**3*v4*w4 + 14*CA2*CF*v4*w4 - 92*CA*v5*w4 - 
     &   13*CA2*v5*w4 + 47*CA**3*v5*w4 + 10*CA*CF*v5*w4 -
     &   22*CA2*CF*v5*w4 - 19*v6*w4 + 112*CA*v6*w4 + 31*CA2*v6*w4 -
     &   69*CA**3*v6*w4 - 26*CA*CF*v6*w4 + 44*CA2*CF*v6*w4 + 
     &   39*v7*w4 - 26*CA*v7*w4 - 40*CA2*v7*w4 + 13*CA**3*v7*w4 + 
     &   38*CA*CF*v7*w4 - 14*CA2*CF*v7*w4 - 14*v8*w4 - CA*v8*w4 +
     &   14*CA2*v8*w4 + 8*CA**3*v8*w4 - 14*CA*CF*v8*w4 -
     &   2*CA2*CF*v8*w4 - 2*v9*w4 + 9*CA*v9*w4 + 2*CA2*v9*w4 -
     &   9*CA**3*v9*w4 - 2*CA*CF*v9*w4 + 6*CA2*CF*v9*w4 - 3*CA*v3*w5 + 
     &   3*CA**3*v3*w5 - 6*CA2*CF*v3*w5 + 3*v4*w5 - 8*CA*v4*w5 -
     &   5*CA2*v4*w5 - 8*CA**3*v4*w5 + 6*CA*CF*v4*w5 +
     &   16*CA2*CF*v4*w5 - 6*v5*w5 + 33*CA*v5*w5 + 10*CA2*v5*w5 -
     &   5*CA**3*v5*w5 - 10*CA*CF*v5*w5 - 10*CA2*CF*v5*w5 + 33*v6*w5 -
     &   6*CA*v6*w5 - 27*CA2*v6*w5 + 2*CA**3*v6*w5 + 
     &   24*CA*CF*v6*w5 - 4*CA2*CF*v6*w5 - 42*v7*w5 - 38*CA*v7*w5 +
     &   34*CA2*v7*w5 + 28*CA**3*v7*w5 - 32*CA*CF*v7*w5 -
     &   6*CA2*CF*v7*w5 + 4*v8*w5 + 10*CA*v8*w5 - 4*CA2*v8*w5 - 
     &   16*CA**3*v8*w5 + 4*CA*CF*v8*w5 + 4*CA2*CF*v8*w5 + 8*v9*w5 -
     &   11*CA*v9*w5 - 8*CA2*v9*w5 + 11*CA**3*v9*w5 + 
     &   8*CA*CF*v9*w5 - 8*CA2*CF*v9*w5 + 2*CA*v3*w6 - 2*CA**3*v3*w6 +
     &   4*CA2*CF*v3*w6 - 2*v4*w6 - 5*CA*v4*w6 + 2*CA2*v4*w6 +
     &   5*CA**3*v4*w6 - 4*CA*CF*v4*w6 - 10*CA2*CF*v4*w6 + 3*v5*w6 +
     &   13*CA*v5*w6 - 5*CA**3*v5*w6 + 2*CA*CF*v5*w6 + 
     &   10*CA2*CF*v5*w6 - 14*v6*w6 - 32*CA*v6*w6 + 4*CA2*v6*w6 +
     &   11*CA**3*v6*w6 - 2*CA*CF*v6*w6 - 2*CA2*CF*v6*w6 +
     &   9*v7*w6 + 28*CA*v7*w6 - 6*CA2*v7*w6 - 16*CA**3*v7*w6 +
     &   4*CA*CF*v7*w6 + 16*v8*w6 + 7*CA*v8*w6 - 12*CA2*v8*w6 +
     &   2*CA**3*v8*w6 + 12*CA*CF*v8*w6 + 4*CA2*CF*v8*w6 - 
     &   12*v9*w6 + 8*CA*v9*w6 + 12*CA2*v9*w6 - 8*CA**3*v9*w6 -
     &   12*CA*CF*v9*w6 + 8*CA2*CF*v9*w6 - CA*v4*w7 + CA**3*v4*w7 -
     &   2*CA2*CF*v4*w7 + v5*w7 + 3*CA*v5*w7 - CA2*v5*w7 -
     &   3*CA**3*v5*w7 + 2*CA*CF*v5*w7 + 6*CA2*CF*v5*w7 - 4*CA*v6*w7 +
     &   3*CA**3*v6*w7 - 2*CA*CF*v6*w7 - 8*CA2*CF*v6*w7 + 7*v7*w7 +
     &   10*CA*v7*w7 - CA2*v7*w7 - 5*CA**3*v7*w7 + 2*CA*CF*v7*w7 + 
     &   8*CA2*CF*v7*w7 - 16*v8*w7 - 16*CA*v8*w7 + 10*CA2*v8*w7 +
     &   6*CA**3*v8*w7 - 10*CA*CF*v8*w7 - 6*CA2*CF*v8*w7 +
     &   8*v9*w7 - 5*CA*v9*w7 - 8*CA2*v9*w7 + 5*CA**3*v9*w7 + 
     &   8*CA*CF*v9*w7 - 8*CA2*CF*v9*w7 - 2*v7*w8 + 4*v8*w8 +
     &   2*CA*v8*w8 -
     &        2*CA2*v8*w8 + 2*CA*CF*v8*w8 - 2*v9*w8 +
     &        3*CA*v9*w8 + 2*CA2*v9*w8 - 3*CA**3*v9*w8 -
     &        2*CA*CF*v9*w8 + 6*CA2*CF*v9*w8 - CA*v9*w9 +
     &        CA**3*v9*w9 - 2*CA2*CF*v9*w9))/
     &    (CA*(1 - v)**2*v3*(1 - w)*w2*(1 - v*w)**3*(1 - v + v*w))
      struv6= part1 + part2 + part3 + part4 +
     &        part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV7(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (16*CA*CF**2*l1vw*(1 - w)*
     &      (1 - v - v2 + v3 + v*w + 4*v2*w - 
     &     5*v3*w - 3*v2*w2 +
     &        5*v3*w2 - v3*w3))/((1 - v)*v*w*
     &     (1 - v + v*w)**2) 
      part2 = -(8*CF**2*lmss*(1 + v2 - 2*v2*w + v2*w2)*
     &      (CA - 4*CA*v + 6*CA*v2 - 4*CA*v3 + CA*v4 - 
     &     v*w + CA*v*w +
     &        3*v2*w - 3*CA*v2*w - 3*v3*w + 3*CA*v3*w + 
     &     v4*w -
     &        CA*v4*w - 2*v2*w2 + CA*v2*w2 + 4*v3*w2 -
     &        2*CA*v3*w2 - 2*v4*w2 + CA*v4*w2 - v3*w3 +
     &        CA*v3*w3 + v4*w3 - CA*v4*w3 + CA*v4*w4))/
     &    ((1 - v)**2*v2*w2*(1 - v + v*w)**2) 
      part3 = (4*CF*lvw*(1 - w)*(4*CA2 + 4*CA2*v2 + 2*v*w - 
     &     9*CA2*v*w -
     &        6*CA2*v2*w - 2*v3*w - CA2*v3*w - 5*v2*w2 +
     &        12*CA2*v2*w2 + 2*v3*w2 + 3*CA2*v3*w2 -
     &        v4*w2 + CA2*v4*w2 + 2*v3*w3 - 
     &     6*CA2*v3*w3 +
     &        2*v4*w3 - 2*CA2*v4*w3 - 2*v4*w4 +
     &        2*CA2*v4*w4))/((1 - v)**2*v2*w2) 
      part4 = -(4*CF*lw*(2*CA - 4*CA**3 - 6*CA*v + 
     &        6*CA**3*v + 8*CA*v2 -
     &        8*CA**3*v2 - 8*CA*v3 + 8*CA**3*v3 + 6*CA*v4 -
     &        4*CA**3*v4 - 2*CA*v5 + 2*CA**3*v5 + 
     &     w - CA2*w - 3*v*w +
     &        3*CA2*v*w + 6*CA**3*v*w + 4*v2*w - 4*CA2*v2*w -
     &        4*CA**3*v2*w - 4*v3*w + 2*CA*v3*w + 
     &     4*CA2*v3*w -
     &        6*CA**3*v3*w + 3*v4*w - 2*CA*v4*w - 
     &     3*CA2*v4*w -
     &        2*CA**3*v4*w - v5*w - 2*CA*v5*w + CA2*v5*w +
     &        2*CA*v6*w - 2*CA**3*v6*w - v*w2 - CA2*v*w2 -
     &        v2*w2 - 5*CA*v2*w2 + CA2*v2*w2 -
     &        3*CA**3*v2*w2 + 4*v3*w2 + CA*v3*w2 -
     &        4*CA2*v3*w2 + 11*CA**3*v3*w2 - 4*v4*w2 +
     &        CA*v4*w2 + 4*CA2*v4*w2 + 11*CA**3*v4*w2 +
     &        v5*w2 + 7*CA*v5*w2 + CA2*v5*w2 -
     &        3*CA**3*v5*w2 + v6*w2 - 4*CA*v6*w2 -
     &        CA2*v6*w2 + 4*CA**3*v6*w2 + v2*w3 +
     &        CA2*v2*w3 + 8*CA*v3*w3 + 2*CA2*v3*w3 -
     &        4*CA**3*v3*w3 - CA*v4*w3 + 2*CA2*v4*w3 -
     &        15*CA**3*v4*w3 + 2*v5*w3 - 10*CA*v5*w3 -
     &        8*CA2*v5*w3 + 2*CA**3*v5*w3 - 3*v6*w3 +
     &        3*CA*v6*w3 + 3*CA2*v6*w3 - 3*CA**3*v6*w3 -
     &        v3*w4 - CA2*v3*w4 + v4*w4 - 
     &     5*CA*v4*w4 -
     &        5*CA2*v4*w4 + 7*CA**3*v4*w4 - 3*v5*w4 +
     &        6*CA*v5*w4 + 9*CA2*v5*w4 + 2*CA**3*v5*w4 +
     &        3*v6*w4 - CA*v6*w4 - 3*CA2*v6*w4 +
     &        CA**3*v6*w4 + 2*CA2*v4*w5 + v5*w5 -
     &        3*CA2*v5*w5 - 2*CA**3*v5*w5 - v6*w5 +
     &        CA2*v6*w5))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)*(1 - v + v*w))
      part5 = -(4*CF*l1v*(2*CA**3 - 4*CA**3*v + 4*CA**3*v2 - 
     &     4*CA**3*v3 +
     &        2*CA**3*v4 - w - CA2*w + v*w - 2*CA*v*w + CA2*v*w -
     &        3*CA**3*v*w + 2*CA*v2*w + 3*CA**3*v2*w + 
     &     2*CA*v3*w +
     &        CA**3*v3*w + v4*w - 2*CA*v4*w + CA2*v4*w +
     &        CA**3*v4*w - v5*w - CA2*v5*w - 2*CA**3*v5*w +
     &        v*w2 + CA2*v*w2 + v2*w2 + 6*CA*v2*w2 -
     &        CA2*v2*w2 + CA**3*v2*w2 - 2*v3*w2 -
     &        6*CA*v3*w2 + 2*CA2*v3*w2 - 3*CA**3*v3*w2 -
     &        2*v4*w2 - 2*CA*v4*w2 - 6*CA2*v4*w2 -
     &        5*CA**3*v4*w2 + v5*w2 + 2*CA*v5*w2 +
     &        5*CA2*v5*w2 + 3*CA**3*v5*w2 + v6*w2 -
     &        CA2*v6*w2 - v2*w3 - CA2*v2*w3 -
     &        6*CA*v3*w3 - 2*CA2*v3*w3 + CA**3*v3*w3 +
     &        2*v4*w3 + 6*CA*v4*w3 + 8*CA2*v4*w3 +
     &        7*CA**3*v4*w3 + 2*v5*w3 - 8*CA2*v5*w3 -
     &        3*v6*w3 + 3*CA2*v6*w3 + v3*w4 + 
     &     CA2*v3*w4 -
     &        v4*w4 + 2*CA*v4*w4 - 3*CA2*v4*w4 -
     &        3*CA**3*v4*w4 - 3*v5*w4 - 2*CA*v5*w4 +
     &        5*CA2*v5*w4 - 3*CA**3*v5*w4 + 3*v6*w4 -
     &        3*CA2*v6*w4 + v5*w5 - CA2*v5*w5 +
     &        2*CA**3*v5*w5 - v6*w5 + CA2*v6*w5))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)*(1 - v + v*w)) 
      part6 = (2*CF*lms*(2 - 6*CA2 - 4*v + 4*v2 - 8*CA2*v2 - 
     &     4*v3 +
     &        2*v4 - 2*CA2*v4 - 2*w + 2*CA2*w + 20*CA2*v*w +
     &        4*v2*w + 12*CA2*v2*w - 4*v3*w + 16*CA2*v3*w +
     &        6*v4*w + 2*CA2*v4*w - 4*v5*w + 
     &     4*CA2*v5*w + w2 -
     &        CA2*w2 + 2*v*w2 - 4*CA2*v*w2 - 2*v2*w2 -
     &        36*CA2*v2*w2 + 2*v3*w2 - 24*CA2*v3*w2 -
     &        v4*w2 - 13*CA2*v4*w2 - 4*CA2*v5*w2 +
     &        2*v6*w2 - 2*CA2*v6*w2 - 
     &     2*v*w3 + 2*CA2*v*w3 +
     &        2*v2*w3 + 2*CA2*v2*w3 - 6*v3*w3 +
     &        40*CA2*v3*w3 - 4*v4*w3 + 20*CA2*v4*w3 +
     &        6*CA2*v5*w3 - 2*v6*w3 + 2*CA2*v6*w3 +
     &        v2*w4 - CA2*v2*w4 - 
     &     2*v3*w4 + 11*v4*w4 -
     &        27*CA2*v4*w4 + 2*v5*w4 - 8*CA2*v5*w4 +
     &        2*v6*w4 - 2*CA2*v6*w4 - 6*v5*w5 +
     &        10*CA2*v5*w5 - 2*v6*w5 + 2*CA2*v6*w5 +
     &        2*v6*w6 - 2*CA2*v6*w6))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2) 
      part7 = -(2*CF*l1w*(4*CA - 12*CA**3 - 16*CA*v + 32*CA**3*v + 
     &     28*CA*v2 -
     &        44*CA**3*v2 - 32*CA*v3 + 48*CA**3*v3 + 28*CA*v4 -
     &        36*CA**3*v4 - 16*CA*v5 + 16*CA**3*v5 + 4*CA*v6 -
     &        4*CA**3*v6 - 2*w - 2*CA*w + 2*CA2*w + 
     &     2*CA**3*w + 6*v*w +
     &        18*CA*v*w - 2*CA2*v*w + 8*CA**3*v*w - 8*v2*w -
     &        28*CA*v2*w - 4*CA2*v2*w - 28*CA**3*v2*w + 
     &     8*v3*w +
     &        32*CA*v3*w + 8*CA2*v3*w - 4*CA**3*v3*w - 
     &     6*v4*w -
     &        38*CA*v4*w - 10*CA2*v4*w + 14*CA**3*v4*w + 
     &     2*v5*w +
     &        6*CA*v5*w + 10*CA2*v5*w + 20*CA**3*v5*w +
     &        20*CA*v6*w - 4*CA2*v6*w - 20*CA**3*v6*w - 
     &     8*CA*v7*w +
     &        8*CA**3*v7*w + CA*w2 - CA**3*w2 + 2*v*w2 -
     &        4*CA*v*w2 + 2*CA2*v*w2 + 2*CA**3*v*w2 - 
     &     6*v2*w2 -
     &        17*CA*v2*w2 - 2*CA2*v2*w2 + CA**3*v2*w2 +
     &        4*v3*w2 + 28*CA*v3*w2 + 46*CA**3*v3*w2 -
     &        4*v4*w2 - 13*CA*v4*w2 + 12*CA2*v4*w2 -
     &        13*CA**3*v4*w2 + 10*v5*w2 + 72*CA*v5*w2 -
     &        22*CA2*v5*w2 - 88*CA**3*v5*w2 - 6*v6*w2 -
     &        79*CA*v6*w2 + 6*CA2*v6*w2 + 
     &     49*CA**3*v6*w2 +
     &        8*CA*v7*w2 + 4*CA2*v7*w2 - 8*CA**3*v7*w2 +
     &        4*CA*v8*w2 - 4*CA**3*v8*w2 + 6*CA*v2*w3 -
     &        4*CA2*v2*w3 - 6*CA**3*v2*w3 - 6*CA*v3*w3 +
     &        4*CA2*v3*w3 - 12*CA**3*v3*w3 + 12*v4*w3 -
     &        12*CA2*v4*w3 - 32*CA**3*v4*w3 - 24*v5*w3 -
     &        92*CA*v5*w3 + 24*CA2*v5*w3 + 
     &     126*CA**3*v5*w3 +
     &        8*v6*w3 + 74*CA*v6*w3 + 4*CA2*v6*w3 -
     &        22*CA**3*v6*w3 + 4*v7*w3 + 30*CA*v7*w3 -
     &        16*CA2*v7*w3 - 18*CA**3*v7*w3 - 
     &     12*CA*v8*w3 +
     &        12*CA**3*v8*w3 - 2*CA*v2*w4 + 2*CA**3*v2*w4 +
     &        6*CA*v3*w4 - 2*CA**3*v3*w4 - 8*v4*w4 -
     &        2*CA*v4*w4 + 8*CA2*v4*w4 + 20*CA**3*v4*w4 +
     &        12*v5*w4 + 64*CA*v5*w4 - 20*CA2*v5*w4 -
     &        66*CA**3*v5*w4 + 8*v6*w4 - 8*CA*v6*w4 -
     &        16*CA2*v6*w4 - 40*CA**3*v6*w4 - 12*v7*w4 -
     &        66*CA*v7*w4 + 28*CA2*v7*w4 + 
     &     38*CA**3*v7*w4 +
     &        16*CA*v8*w4 - 16*CA**3*v8*w4 + 2*v4*w5 -
     &        4*CA*v4*w5 - 2*CA2*v4*w5 + 4*CA**3*v4*w5 +
     &        2*v5*w5 - 16*CA*v5*w5 + 10*CA2*v5*w5 +
     &        8*CA**3*v5*w5 - 16*v6*w5 - 38*CA*v6*w5 +
     &        20*CA2*v6*w5 + 50*CA**3*v6*w5 + 12*v7*w5 +
     &        50*CA*v7*w5 - 28*CA2*v7*w5 - 
     &     22*CA**3*v7*w5 -
     &        16*CA*v8*w5 + 16*CA**3*v8*w5 + CA*v4*w6 -
     &        CA**3*v4*w6 - 2*v5*w6 - 2*CA*v5*w6 -
     &        2*CA2*v5*w6 + 6*v6*w6 + 23*CA*v6*w6 -
     &        14*CA2*v6*w6 - 17*CA**3*v6*w6 - 4*v7*w6 -
     &        10*CA*v7*w6 + 16*CA2*v7*w6 - 
     &     2*CA**3*v7*w6 +
     &        16*CA*v8*w6 - 16*CA**3*v8*w6 + 
     &     4*CA2*v6*w7 -
     &        4*CA*v7*w7 - 4*CA2*v7*w7 + 4*CA**3*v7*w7 -
     &        12*CA*v8*w7 + 12*CA**3*v8*w7 + 4*CA*v8*w8 -
     &        4*CA**3*v8*w8))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part8 = -(2*CF*lv*(4*CA - 16*CA**3 - 16*CA*v + 44*CA**3*v + 
     &     28*CA*v2 -
     &        60*CA**3*v2 - 32*CA*v3 + 64*CA**3*v3 + 28*CA*v4 -
     &        48*CA**3*v4 - 16*CA*v5 + 20*CA**3*v5 + 4*CA*v6 -
     &        4*CA**3*v6 - 2*CA*w + 2*CA**3*w - 2*v*w + 22*CA*v*w +
     &        6*CA2*v*w + 14*CA**3*v*w + 6*v2*w - 36*CA*v2*w -
     &        18*CA2*v2*w - 44*CA**3*v2*w - 8*v3*w + 
     &     32*CA*v3*w +
     &        24*CA2*v3*w + 8*CA**3*v3*w + 8*v4*w - 
     &     30*CA*v4*w -
     &        24*CA2*v4*w + 6*CA**3*v4*w - 6*v5*w + 
     &     2*CA*v5*w +
     &        18*CA2*v5*w + 34*CA**3*v5*w + 2*v6*w + 
     &     20*CA*v6*w -
     &        6*CA2*v6*w - 28*CA**3*v6*w - 8*CA*v7*w +
     &        8*CA**3*v7*w + CA*w2 - CA**3*w2 - 4*CA*v*w2 +
     &        2*CA**3*v*w2 - 4*v2*w2 - 29*CA*v2*w2 +
     &        3*CA**3*v2*w2 + 8*v3*w2 + 56*CA*v3*w2 -
     &        4*CA2*v3*w2 + 52*CA**3*v3*w2 - 12*v4*w2 -
     &        25*CA*v4*w2 + 20*CA2*v4*w2 - 
     &     7*CA**3*v4*w2 +
     &        12*v5*w2 + 60*CA*v5*w2 - 20*CA2*v5*w2 -
     &        114*CA**3*v5*w2 - 71*CA*v6*w2 - 
     &     4*CA2*v6*w2 +
     &        57*CA**3*v6*w2 - 4*v7*w2 + 8*CA*v7*w2 +
     &        8*CA2*v7*w2 - 4*CA**3*v7*w2 + 4*CA*v8*w2 -
     &        4*CA**3*v8*w2 + 6*CA*v2*w3 - 6*CA**3*v2*w3 +
     &        2*v3*w3 + 2*CA*v3*w3 - 2*CA2*v3*w3 -
     &        20*CA**3*v3*w3 + 2*v4*w3 - 32*CA*v4*w3 -
     &        2*CA2*v4*w3 - 40*CA**3*v4*w3 - 4*v5*w3 -
     &        64*CA*v5*w3 - 12*CA2*v5*w3 + 
     &     148*CA**3*v5*w3 -
     &        16*v6*w3 + 74*CA*v6*w3 + 40*CA2*v6*w3 -
     &        10*CA**3*v6*w3 + 14*v7*w3 + 26*CA*v7*w3 -
     &        22*CA2*v7*w3 - 28*CA**3*v7*w3 + 2*v8*w3 -
     &        12*CA*v8*w3 - 2*CA2*v8*w3 + 
     &     12*CA**3*v8*w3 -
     &        2*CA*v2*w4 + 2*CA**3*v2*w4 + 6*CA*v3*w4 -
     &        2*CA**3*v3*w4 + 6*CA*v4*w4 + 28*CA**3*v4*w4 -
     &        4*v5*w4 + 72*CA*v5*w4 + 20*CA2*v5*w4 -
     &        74*CA**3*v5*w4 + 28*v6*w4 - 28*CA*v6*w4 -
     &        44*CA2*v6*w4 - 70*CA**3*v6*w4 - 16*v7*w4 -
     &        62*CA*v7*w4 + 16*CA2*v7*w4 + 
     &     44*CA**3*v7*w4 -
     &        8*v8*w4 + 16*CA*v8*w4 + 8*CA2*v8*w4 -
     &        16*CA**3*v8*w4 - 4*CA*v4*w5 + 4*CA**3*v4*w5 +
     &        2*v5*w5 - 28*CA*v5*w5 - 6*CA2*v5*w5 +
     &        6*CA**3*v5*w5 - 18*v6*w5 - 30*CA*v6*w5 +
     &        14*CA2*v6*w5 + 74*CA**3*v6*w5 + 4*v7*w5 +
     &        54*CA*v7*w5 + 4*CA2*v7*w5 - 
     &     16*CA**3*v7*w5 +
     &        12*v8*w5 - 16*CA*v8*w5 - 12*CA2*v8*w5 +
     &        16*CA**3*v8*w5 + CA*v4*w6 - CA**3*v4*w6 -
     &        2*CA*v5*w6 + 4*v6*w6 + 27*CA*v6*w6 -
     &        23*CA**3*v6*w6 + 4*v7*w6 - 14*CA*v7*w6 -
     &        8*CA2*v7*w6 - 12*CA**3*v7*w6 - 8*v8*w6 +
     &        16*CA*v8*w6 + 8*CA2*v8*w6 - 
     &     16*CA**3*v8*w6 -
     &        2*v7*w7 - 4*CA*v7*w7 + 2*CA2*v7*w7 +
     &        8*CA**3*v7*w7 + 2*v8*w7 - 12*CA*v8*w7 -
     &        2*CA2*v8*w7 + 12*CA**3*v8*w7 + 4*CA*v8*w8 -
     &        4*CA**3*v8*w8))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part9 = -(2*CF*(2*CA - 6*CA**3 - 12*CA*v + 28*CA**3*v + 
     &      30*CA*v2 -
     &        54*CA**3*v2 - 40*CA*v3 + 56*CA**3*v3 + 30*CA*v4 -
     &        34*CA**3*v4 - 12*CA*v5 + 12*CA**3*v5 + 2*CA*v6 -
     &        2*CA**3*v6 + CA*w - CA**3*w - 2*v*w - CA*v*w + 
     &     2*CA2*v*w +
     &        CA**3*v*w + 10*v2*w - 2*CA*v2*w - 10*CA2*v2*w -
     &        6*CA**3*v2*w - 20*v3*w - 4*CA*v3*w + 
     &     20*CA2*v3*w +
     &        28*CA**3*v3*w + 20*v4*w + 23*CA*v4*w - 
     &     20*CA2*v4*w -
     &        47*CA**3*v4*w - 10*v5*w - 31*CA*v5*w + 
     &     10*CA2*v5*w +
     &        39*CA**3*v5*w + 2*v6*w + 18*CA*v6*w - 
     &     2*CA2*v6*w -
     &        18*CA**3*v6*w - 4*CA*v7*w + 4*CA**3*v7*w - 
     &     CA*w2 +
     &        CA**3*w2 + 9*CA*v*w2 - 7*CA**3*v*w2 - 4*v2*w2 -
     &        32*CA*v2*w2 + 4*CA2*v2*w2 + 
     &     22*CA**3*v2*w2 +
     &        12*v3*w2 + 34*CA*v3*w2 - 12*CA2*v3*w2 -
     &        38*CA**3*v3*w2 - 8*v4*w2 - 11*CA*v4*w2 +
     &        8*CA2*v4*w2 + 39*CA**3*v4*w2 - 8*v5*w2 +
     &        13*CA*v5*w2 + 8*CA2*v5*w2 - 
     &     27*CA**3*v5*w2 +
     &        12*v6*w2 - 14*CA*v6*w2 - 12*CA2*v6*w2 +
     &        12*CA**3*v6*w2 - 4*v7*w2 + 4*CA2*v7*w2 +
     &        2*CA*v8*w2 - 2*CA**3*v8*w2 - 4*CA*v2*w3 +
     &        4*CA**3*v2*w3 + 2*v3*w3 + 58*CA*v3*w3 -
     &        2*CA2*v3*w3 - 28*CA**3*v3*w3 - 18*v4*w3 -
     &        80*CA*v4*w3 + 18*CA2*v4*w3 + 
     &     32*CA**3*v4*w3 +
     &        40*v5*w3 - 4*CA*v5*w3 - 40*CA2*v5*w3 -
     &        32*v6*w3 + 32*CA*v6*w3 + 32*CA2*v6*w3 -
     &        8*CA**3*v6*w3 + 6*v7*w3 + 4*CA*v7*w3 -
     &        6*CA2*v7*w3 - 6*CA**3*v7*w3 + 2*v8*w3 -
     &        6*CA*v8*w3 - 2*CA2*v8*w3 + 6*CA**3*v8*w3 +
     &        2*CA*v2*w4 - 2*CA**3*v2*w4 - 16*CA*v3*w4 +
     &        12*CA**3*v3*w4 + 8*v4*w4 - 12*CA*v4*w4 -
     &        8*CA2*v4*w4 + 2*CA**3*v4*w4 - 24*v5*w4 +
     &        70*CA*v5*w4 + 24*CA2*v5*w4 - 
     &     32*CA**3*v5*w4 +
     &        16*v6*w4 - 44*CA*v6*w4 - 16*CA2*v6*w4 +
     &        18*CA**3*v6*w4 + 8*v7*w4 - 10*CA*v7*w4 -
     &        8*CA2*v7*w4 + 8*CA**3*v7*w4 - 8*v8*w4 +
     &        10*CA*v8*w4 + 8*CA2*v8*w4 - 
     &     10*CA**3*v8*w4 +
     &        3*CA*v4*w5 - 3*CA**3*v4*w5 + 2*v5*w5 -
     &        11*CA*v5*w5 - 2*CA2*v5*w5 + 
     &     13*CA**3*v5*w5 +
     &        6*v6*w5 + 10*CA*v6*w5 - 6*CA2*v6*w5 -
     &        6*CA**3*v6*w5 - 20*v7*w5 + 10*CA*v7*w5 +
     &        20*CA2*v7*w5 - 8*CA**3*v7*w5 + 12*v8*w5 -
     &        12*CA*v8*w5 - 12*CA2*v8*w5 + 
     &     12*CA**3*v8*w5 -
     &        CA*v4*w6 + CA**3*v4*w6 + 7*CA*v5*w6 -
     &        5*CA**3*v5*w6 - 4*v6*w6 - 12*CA*v6*w6 +
     &        4*CA2*v6*w6 + 4*CA**3*v6*w6 + 12*v7*w6 -
     &        2*CA*v7*w6 - 12*CA2*v7*w6 + 4*CA**3*v7*w6 -
     &        8*v8*w6 + 10*CA*v8*w6 + 8*CA2*v8*w6 -
     &        10*CA**3*v8*w6 - 2*v7*w7 + 2*CA*v7*w7 +
     &        2*CA2*v7*w7 - 2*CA**3*v7*w7 + 2*v8*w7 -
     &        6*CA*v8*w7 - 2*CA2*v8*w7 + 6*CA**3*v8*w7 +
     &        2*CA*v8*w8 - 2*CA**3*v8*w8))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      struv7= part1 + part2 + part3 + part4 +
     &        part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV8(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (-2*CF*lmss*(2 - 4*v + 2*v2 + 
     &     2*v*w - 2*v2*w + v2*w2)*
     &      (1 - 2*v + 2*v2 + 2*v*w - 4*v2*w + 2*v2*w2)*
     &      (2*CA2 - 4*CA2*v + 2*CA2*v2 + 2*CA2*v*w -
     &        2*CA2*v2*w - v2*w2 + CA2*v2*w2))/
     &    ((1 - v)*v2*w2*(1 - v + v*w)**2) 
      part2 = -(4*CF*l1vw*(1 - w)*(2 - 8*v + 2*CA2*v + 10*v2 - 
     &     4*CA2*v2 -
     &        4*v3 + 2*CA2*v3 + 5*v*w - CA2*v*w - 13*v2*w +
     &        3*CA2*v2*w + 8*v3*w - 2*CA2*v3*w + 
     &     4*v2*w2 -
     &        6*v3*w2 + 2*v3*w3))/(v*w*(1 - v + v*w)**2) 
      part3 = -(8*CF*l1v*(CA2 - 2*CA2*v + 2*CA2*v2 - 
     &     3*v*w + CA2*v*w +
     &        6*v2*w - 3*CA2*v2*w - 8*v3*w - 2*v2*w2 +
     &        2*CA2*v2*w2 + 10*v3*w2 - CA2*v3*w2 +
     &        2*CA2*v4*w2 - 5*v3*w3 - 2*CA2*v4*w3 +
     &        CA2*v4*w4))/(v2*w2*(1 - v*w)) 
      part4 = (8*CF*lvw*(1 - w)*(2*CA*CF + 2*v - 2*CA2*v - 2*v2 +
     &        2*CA2*v2 + 2*v*w + CA2*v*w - 5*v2*w - 
     &     4*CA2*v2*w +
     &        4*v3*w + 2*CA2*v3*w + 3*CA2*v2*w2 + 
     &     v3*w2 -
     &        4*CA2*v3*w2 - 2*v4*w2 + 2*CA2*v4*w2 -
     &        v3*w3 + CA2*v3*w3 + 2*v4*w3 - 
     &     2*CA2*v4*w3 -
     &        v4*w4 + CA2*v4*w4))/((1 - v)*v2*w2) 
      part5 = -(4*CF*lw*(2 - 4*CA2 - 4*v + 14*CA2*v + 4*v2 - 
     &     24*CA2*v2 +
     &        20*CA2*v3 - 8*CA2*v4 - 2*v*w -
     &      2*CA2*v*w + 2*v2*w +
     &        10*CA2*v2*w - 4*v3*w - 4*CA2*v3*w - 
     &     8*CA2*v4*w +
     &        8*CA2*v5*w + 5*v2*w2 - 3*CA2*v2*w2 -
     &        10*v3*w2 - 4*CA2*v3*w2 + 8*v4*w2 +
     &        16*CA2*v4*w2 - 12*CA2*v5*w2 + 3*v3*w3 +
     &        3*CA2*v3*w3 - 4*v4*w3 - 10*CA2*v4*w3 +
     &        8*CA2*v5*w3 + 2*CA2*v4*w4 - 
     &     2*CA2*v5*w4))/
     &    ((1 - v)*v2*w2*(1 - v*w)) 
      part6 = (2*CF*lms*(4 - 4*CA2 - 8*v + 12*CA2*v + 8*v2 - 
     &     20*CA2*v2 +
     &        16*CA2*v3 - 8*CA2*v4 - 2*w + 2*CA2*w - 2*v*w -
     &        2*CA2*v*w + 4*v2*w - 12*v3*w + 20*CA2*v3*w -
     &        24*CA2*v4*w + 16*CA2*v5*w + w2 - CA2*w2 +
     &        2*v*w2 + 2*CA2*v2*w2 + 4*v3*w2 -
     &        16*CA2*v3*w2 + 8*v4*w2 + 8*CA2*v4*w2 -
     &        8*CA2*v6*w2 - 2*v*w3 + 2*CA2*v*w3 + 
     &     2*v2*w3 -
     &        6*CA2*v2*w3 - 4*v3*w3 + 12*CA2*v3*w3 -
     &        4*v4*w3 - 4*v5*w3 - 4*CA2*v5*w3 +
     &        8*CA2*v6*w3 + v2*w4 - CA2*v2*w4 - 
     &     2*v3*w4 +
     &        4*CA2*v3*w4 + 6*v4*w4 - 12*CA2*v4*w4 +
     &        8*CA2*v5*w4 + 4*v6*w4 - 8*CA2*v6*w4 -
     &        2*v5*w5 + 2*CA2*v5*w5 - 4*v6*w5 +
     &        4*CA2*v6*w5 + 2*v6*w6 - 2*CA2*v6*w6))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2) 
      part7 = -(2*CF*l1w*(4 - 12*CA2 - 16*v + 64*CA2*v + 28*v2 -
     &        156*CA2*v2 - 24*v3 + 216*CA2*v3 + 8*v4 -
     &        176*CA2*v4 + 80*CA2*v5 - 16*CA2*v6 - 2*w +
     &        2*CA2*w + 22*v*w - 30*CA2*v*w - 78*v2*w +
     &        114*CA2*v2*w + 150*v3*w - 178*CA2*v3*w - 
     &     168*v4*w +
     &        84*CA2*v4*w + 108*v5*w + 72*CA2*v5*w - 
     &     32*v6*w -
     &        96*CA2*v6*w + 32*CA2*v7*w + w2 - CA2*w2 -
     &        4*v*w2 + 6*CA2*v*w2 + 23*v2*w2 - 
     &     23*CA2*v2*w2 -
     &        98*v3*w2 + 12*CA2*v3*w2 + 178*v4*w2 +
     &        138*CA2*v4*w2 - 124*v5*w2 - 
     &     292*CA2*v5*w2 -
     &        4*v6*w2 + 204*CA2*v6*w2 + 32*v7*w2 -
     &        32*CA2*v7*w2 - 16*CA2*v8*w2 + 6*v2*w3 -
     &        6*CA2*v2*w3 + 2*v3*w3 + 40*CA2*v3*w3 -
     &        36*v4*w3 - 158*CA2*v4*w3 - 40*v5*w3 +
     &        260*CA2*v5*w3 + 172*v6*w3 - 
     &     132*CA2*v6*w3 -
     &        116*v7*w3 - 40*CA2*v7*w3 + 48*CA2*v8*w3 -
     &        2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 -
     &        10*CA2*v3*w4 - 2*v4*w4 + 40*CA2*v4*w4 +
     &        84*v5*w4 - 74*CA2*v5*w4 - 222*v6*w4 -
     &        2*CA2*v6*w4 + 144*v7*w4 + 92*CA2*v7*w4 +
     &        8*v8*w4 - 64*CA2*v8*w4 - 4*v4*w5 +
     &        4*CA2*v4*w5 - 24*v5*w5 - 10*CA2*v5*w5 +
     &        106*v6*w5 + 44*CA2*v6*w5 - 70*v7*w5 -
     &        78*CA2*v7*w5 - 24*v8*w5 + 56*CA2*v8*w5 +
     &        v4*w6 - CA2*v4*w6 - 2*v5*w6 + 
     &     4*CA2*v5*w6 -
     &        17*v6*w6 - 13*CA2*v6*w6 + 6*v7*w6 +
     &        30*CA2*v7*w6 + 28*v8*w6 - 36*CA2*v8*w6 +
     &        4*v7*w7 - 4*CA2*v7*w7 - 16*v8*w7 +
     &        16*CA2*v8*w7 + 4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part8 = -(2*CF*lv*(4 - 16*CA2 - 16*v + 84*CA2*v + 28*v2 -
     &        200*CA2*v2 - 24*v3 + 268*CA2*v3 + 8*v4 -
     &        208*CA2*v4 + 88*CA2*v5 - 16*CA2*v6 - 2*w +
     &        2*CA2*w + 26*v*w - 38*CA2*v*w - 98*v2*w +
     &        150*CA2*v2*w + 186*v3*w - 238*CA2*v3*w - 
     &     196*v4*w +
     &        120*CA2*v4*w + 116*v5*w + 76*CA2*v5*w - 
     &     32*v6*w -
     &        104*CA2*v6*w + 32*CA2*v7*w + w2 - CA2*w2 -
     &        4*v*w2 + 6*CA2*v*w2 + 35*v2*w2 - 
     &     27*CA2*v2*w2 -
     &        142*v3*w2 + 16*CA2*v3*w2 + 230*v4*w2 +
     &        178*CA2*v4*w2 - 144*v5*w2 - 
     &     380*CA2*v5*w2 -
     &        4*v6*w2 + 260*CA2*v6*w2 + 32*v7*w2 -
     &        40*CA2*v7*w2 - 16*CA2*v8*w2 + 6*v2*w3 -
     &        6*CA2*v2*w3 + 10*v3*w3 + 44*CA2*v3*w3 -
     &        44*v4*w3 - 202*CA2*v4*w3 - 60*v5*w3 +
     &        360*CA2*v5*w3 + 200*v6*w3 - 
     &     196*CA2*v6*w3 -
     &        124*v7*w3 - 44*CA2*v7*w3 + 56*CA2*v8*w3 -
     &        2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 -
     &        10*CA2*v3*w4 - 10*v4*w4 + 52*CA2*v4*w4 +
     &        124*v5*w4 - 118*CA2*v5*w4 - 274*v6*w4 +
     &        18*CA2*v6*w4 + 164*v7*w4 + 
     &     128*CA2*v7*w4 +
     &        8*v8*w4 - 88*CA2*v8*w4 - 4*v4*w5 +
     &        4*CA2*v4*w5 - 36*v5*w5 - 2*CA2*v5*w5 +
     &        134*v6*w5 + 48*CA2*v6*w5 - 86*v7*w5 -
     &        118*CA2*v7*w5 - 24*v8*w5 + 84*CA2*v8*w5 +
     &        v4*w6 - CA2*v4*w6 - 2*v5*w6 + 
     &     4*CA2*v5*w6 -
     &        21*v6*w6 - 17*CA2*v6*w6 + 10*v7*w6 +
     &        50*CA2*v7*w6 + 28*v8*w6 - 52*CA2*v8*w6 +
     &        4*v7*w7 - 8*CA2*v7*w7 - 16*v8*w7 +
     &        20*CA2*v8*w7 + 4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part9 = -(2*CF*(4*v - 12*CA2*v - 12*v2 + 52*CA2*v2 + 12*v3 -
     &        92*CA2*v3 - 4*v4 + 84*CA2*v4 - 40*CA2*v5 +
     &        8*CA2*v6 + w - CA2*w - 9*v*w + 17*CA2*v*w + 
     &     27*v2*w -
     &        75*CA2*v2*w - 29*v3*w + 125*CA2*v3*w + 
     &     8*v4*w -
     &        72*CA2*v4*w + 2*v5*w - 26*CA2*v5*w +
     &        48*CA2*v6*w - 16*CA2*v7*w - w2 + CA2*w2 +
     &        7*v*w2 - 3*CA2*v*w2 - 17*v2*w2 + 
     &     11*CA2*v2*w2 +
     &        v3*w2 + 7*CA2*v3*w2 + 34*v4*w2 -
     &        112*CA2*v4*w2 - 28*v5*w2 + 
     &     204*CA2*v5*w2 +
     &        4*v6*w2 - 140*CA2*v6*w2 + 24*CA2*v7*w2 +
     &        8*CA2*v8*w2 - 4*v2*w3 + 4*CA2*v2*w3 +
     &        34*v3*w3 - 24*CA2*v3*w3 - 76*v4*w3 +
     &        86*CA2*v4*w3 + 68*v5*w3 - 144*CA2*v5*w3 -
     &        32*v6*w3 + 80*CA2*v6*w3 + 10*v7*w3 +
     &        30*CA2*v7*w3 - 32*CA2*v8*w3 + 2*v2*w4 -
     &        2*CA2*v2*w4 - 12*v3*w4 + 4*CA2*v3*w4 +
     &        27*v4*w4 - 5*CA2*v4*w4 - 36*v5*w4 +
     &        20*CA2*v5*w4 + 36*v6*w4 + 10*CA2*v6*w4 -
     &        8*v7*w4 - 84*CA2*v7*w4 - 8*v8*w4 +
     &        56*CA2*v8*w4 + 3*v4*w5 - 3*CA2*v4*w5 -
     &        5*v5*w5 - CA2*v5*w5 - 3*v6*w5 -
     &        15*CA2*v6*w5 - 19*v7*w5 + 79*CA2*v7*w5 +
     &        24*v8*w5 - 60*CA2*v8*w5 - v4*w6 +
     &        CA2*v4*w6 + 5*v5*w6 - CA2*v5*w6 - 
     &     6*v6*w6 +
     &        6*CA2*v6*w6 + 23*v7*w6 - 43*CA2*v7*w6 -
     &        26*v8*w6 + 42*CA2*v8*w6 - 6*v7*w7 +
     &        10*CA2*v7*w7 + 12*v8*w7 - 16*CA2*v8*w7 -
     &        2*v8*w8 + 2*CA2*v8*w8))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      struv8= part1 + part2 + part3 + part4 +
     &        part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV9(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (-2*CF*lmss*(2 - 4*v + 2*v2 + 2*v*w - 
     &     2*v2*w + v2*w2)*
     &      (1 - 2*v + 2*v2 + 2*v*w - 4*v2*w + 2*v2*w2)*
     &      (2*CA2 - 4*CA2*v + 2*CA2*v2 + 2*CA2*v*w -
     &        2*CA2*v2*w - v2*w2 + CA2*v2*w2))/
     &    ((1 - v)*v2*w2*(1 - v + v*w)**2) 
      part2 = (4*CF*l1vw*(1 - w)*(2 - CA2 - 8*v + 2*CA2*v + 10*v2 -
     &        CA2*v2 - 4*v3 + 7*v*w - 2*CA2*v*w - 19*v2*w +
     &        5*CA2*v2*w + 12*v3*w - 3*CA2*v3*w + 
     &     8*v2*w2 -
     &        3*CA2*v2*w2 - 10*v3*w2 + 4*CA2*v3*w2 +
     &        2*v3*w3 - CA2*v3*w3))/(v*w*(1 - v + v*w)**2) 
      part3 = (16*CA*CF**2*lvw*(1 - w)*
     &      (1 - 2*v + 2*v2 + 2*v*w - 7*v2*w + 4*v3*w + 
     &     4*v2*w2 -
     &        5*v3*w2 + 2*v4*w2 + v3*w3 - 
     &     2*v4*w3 + v4*w4)
     &      )/((1 - v)*v2*w2) 
      part4 = -(4*CF*l1v*(2*CA2 - 4*CA2*v + 4*CA2*v2 + 
     &     6*v*w - CA2*v*w -
     &        12*v2*w + 16*v3*w - 8*CA2*v3*w + 4*v2*w2 +
     &        2*CA2*v2*w2 - 20*v3*w2 + 8*CA2*v3*w2 +
     &        4*CA2*v4*w2 + 10*v3*w3 - 5*CA2*v3*w3 -
     &        4*CA2*v4*w3 + 2*CA2*v4*w4))/
     &     (v2*w2*(1 - v*w)) 
      part5 = -(4*CF*lw*(2 - 4*CA2 - 4*v + 14*CA2*v + 4*v2 - 
     &     24*CA2*v2 +
     &        20*CA2*v3 - 8*CA2*v4 - 2*v*w - 
     &     2*CA2*v*w + 2*v2*w +
     &        10*CA2*v2*w - 4*v3*w - 4*CA2*v3*w - 
     &     8*CA2*v4*w +
     &        8*CA2*v5*w - 3*v2*w2 - CA2*v2*w2 + 
     &     14*v3*w2 -
     &        10*CA2*v3*w2 - 8*v4*w2 + 20*CA2*v4*w2 -
     &        12*CA2*v5*w2 - 5*v3*w3 + 5*CA2*v3*w3 +
     &        4*v4*w3 - 12*CA2*v4*w3 + 8*CA2*v5*w3 +
     &        2*CA2*v4*w4 - 2*CA2*v5*w4))/
     &    ((1 - v)*v2*w2*(1 - v*w)) 
      part6 = (2*CF*lms*(4 - 4*CA2 - 8*v + 12*CA2*v + 8*v2 - 
     &     20*CA2*v2 +
     &        16*CA2*v3 - 8*CA2*v4 - 2*w + 2*CA2*w - 2*v*w -
     &        2*CA2*v*w + 4*v2*w - 12*v3*w + 20*CA2*v3*w -
     &        24*CA2*v4*w + 16*CA2*v5*w + w2 - CA2*w2 +
     &        2*v*w2 + 2*CA2*v2*w2 + 4*v3*w2 -
     &        16*CA2*v3*w2 + 8*v4*w2 + 8*CA2*v4*w2 -
     &        8*CA2*v6*w2 - 2*v*w3 + 2*CA2*v*w3 + 
     &     2*v2*w3 -
     &        6*CA2*v2*w3 - 4*v3*w3 + 12*CA2*v3*w3 -
     &        4*v4*w3 - 4*v5*w3 - 4*CA2*v5*w3 +
     &        8*CA2*v6*w3 + v2*w4 - CA2*v2*w4 - 
     &     2*v3*w4 +
     &        4*CA2*v3*w4 + 6*v4*w4 - 12*CA2*v4*w4 +
     &        8*CA2*v5*w4 + 4*v6*w4 - 8*CA2*v6*w4 -
     &        2*v5*w5 + 2*CA2*v5*w5 - 4*v6*w5 +
     &        4*CA2*v6*w5 + 2*v6*w6 - 2*CA2*v6*w6))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2) 
      part7 = -(2*CF*l1w*(4 - 12*CA2 - 16*v + 64*CA2*v + 28*v2 -
     &        156*CA2*v2 - 24*v3 + 216*CA2*v3 + 8*v4 -
     &        176*CA2*v4 + 80*CA2*v5 - 16*CA2*v6 - 2*w +
     &        2*CA2*w - 2*v*w - 24*CA2*v*w + 42*v2*w +
     &        84*CA2*v2*w - 130*v3*w - 108*CA2*v3*w + 
     &     192*v4*w -
     &        6*CA2*v4*w - 132*v5*w + 132*CA2*v5*w + 
     &     32*v6*w -
     &        112*CA2*v6*w + 32*CA2*v7*w + w2 - CA2*w2 -
     &        4*v*w2 + 6*CA2*v*w2 - 17*v2*w2 - 
     &     13*CA2*v2*w2 +
     &        102*v3*w2 - 38*CA2*v3*w2 - 198*v4*w2 +
     &        232*CA2*v4*w2 + 124*v5*w2 - 
     &     354*CA2*v5*w2 +
     &        28*v6*w2 + 196*CA2*v6*w2 - 32*v7*w2 -
     &        16*CA2*v7*w2 - 16*CA2*v8*w2 + 6*v2*w3 -
     &        6*CA2*v2*w3 - 30*v3*w3 + 48*CA2*v3*w3 +
     &        60*v4*w3 - 182*CA2*v4*w3 + 32*v5*w3 +
     &        242*CA2*v5*w3 - 172*v6*w3 - 
     &     46*CA2*v6*w3 +
     &        92*v7*w3 - 92*CA2*v7*w3 + 48*CA2*v8*w3 -
     &        2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 -
     &        10*CA2*v3*w4 - 2*v4*w4 + 40*CA2*v4*w4 -
     &        76*v5*w4 - 34*CA2*v5*w4 + 202*v6*w4 -
     &        108*CA2*v6*w4 - 120*v7*w4 + 
     &     158*CA2*v7*w4 +
     &        8*v8*w4 - 64*CA2*v8*w4 - 4*v4*w5 +
     &        4*CA2*v4*w5 + 32*v5*w5 - 24*CA2*v5*w5 -
     &        110*v6*w5 + 98*CA2*v6*w5 + 90*v7*w5 -
     &        118*CA2*v7*w5 - 24*v8*w5 + 56*CA2*v8*w5 +
     &        v4*w6 - CA2*v4*w6 - 2*v5*w6 + 
     &     4*CA2*v5*w6 +
     &        23*v6*w6 - 23*CA2*v6*w6 - 34*v7*w6 +
     &        40*CA2*v7*w6 + 28*v8*w6 - 36*CA2*v8*w6 +
     &        4*v7*w7 - 4*CA2*v7*w7 - 16*v8*w7 +
     &        16*CA2*v8*w7 + 4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part8 = -(2*CF*lv*(4 - 16*CA2 - 16*v + 84*CA2*v + 28*v2 -
     &        200*CA2*v2 - 24*v3 + 268*CA2*v3 + 8*v4 -
     &        208*CA2*v4 + 88*CA2*v5 - 16*CA2*v6 - 2*w +
     &        2*CA2*w - 6*v*w - 30*CA2*v*w + 62*v2*w +
     &        110*CA2*v2*w - 166*v3*w - 
     &     150*CA2*v3*w + 220*v4*w +
     &        16*CA2*v4*w - 140*v5*w + 140*CA2*v5*w + 
     &     32*v6*w -
     &        120*CA2*v6*w + 32*CA2*v7*w + w2 - CA2*w2 -
     &        4*v*w2 + 6*CA2*v*w2 - 29*v2*w2 - 
     &     11*CA2*v2*w2 +
     &        146*v3*w2 - 56*CA2*v3*w2 - 250*v4*w2 +
     &        298*CA2*v4*w2 + 144*v5*w2 - 
     &     452*CA2*v5*w2 +
     &        28*v6*w2 + 252*CA2*v6*w2 - 32*v7*w2 -
     &        24*CA2*v7*w2 - 16*CA2*v8*w2 + 6*v2*w3 -
     &        6*CA2*v2*w3 - 38*v3*w3 + 56*CA2*v3*w3 +
     &        68*v4*w3 - 230*CA2*v4*w3 + 52*v5*w3 +
     &        332*CA2*v5*w3 - 200*v6*w3 - 
     &     96*CA2*v6*w3 +
     &        100*v7*w3 - 100*CA2*v7*w3 + 
     &     56*CA2*v8*w3 -
     &        2*v2*w4 + 2*CA2*v2*w4 + 6*v3*w4 -
     &        10*CA2*v3*w4 + 6*v4*w4 + 48*CA2*v4*w4 -
     &        116*v5*w4 - 58*CA2*v5*w4 + 254*v6*w4 -
     &        114*CA2*v6*w4 - 140*v7*w4 + 
     &     204*CA2*v7*w4 +
     &        8*v8*w4 - 88*CA2*v8*w4 - 4*v4*w5 +
     &        4*CA2*v4*w5 + 44*v5*w5 - 22*CA2*v5*w5 -
     &        138*v6*w5 + 116*CA2*v6*w5 + 106*v7*w5 -
     &        166*CA2*v7*w5 - 24*v8*w5 + 
     &     84*CA2*v8*w5 +
     &        v4*w6 - CA2*v4*w6 - 
     &     2*v5*w6 + 4*CA2*v5*w6 +
     &        27*v6*w6 - 29*CA2*v6*w6 - 38*v7*w6 +
     &        62*CA2*v7*w6 + 28*v8*w6 - 52*CA2*v8*w6 +
     &        4*v7*w7 - 8*CA2*v7*w7 - 16*v8*w7 +
     &        20*CA2*v8*w7 + 4*v8*w8 - 4*CA2*v8*w8))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part9 = -(2*CF*(4*v - 12*CA2*v - 12*v2 + 52*CA2*v2 + 12*v3 -
     &        92*CA2*v3 - 4*v4 + 84*CA2*v4 - 40*CA2*v5 +
     &        8*CA2*v6 + w - CA2*w - 9*v*w + 
     &     17*CA2*v*w + 27*v2*w -
     &        75*CA2*v2*w - 29*v3*w +
     &      125*CA2*v3*w + 8*v4*w -
     &        72*CA2*v4*w + 2*v5*w - 26*CA2*v5*w +
     &        48*CA2*v6*w - 16*CA2*v7*w - w2 + CA2*w2 -
     &        v*w2 - CA2*v*w2 - v2*w2 + 7*CA2*v2*w2 +
     &        41*v3*w2 - 3*CA2*v3*w2 - 94*v4*w2 -
     &        80*CA2*v4*w2 + 84*v5*w2 + 176*CA2*v5*w2 -
     &        28*v6*w2 - 132*CA2*v6*w2 + 24*CA2*v7*w2 +
     &        8*CA2*v8*w2 - 4*v2*w3 + 4*CA2*v2*w3 -
     &        38*v3*w3 - 6*CA2*v3*w3 + 156*v4*w3 +
     &        28*CA2*v4*w3 - 172*v5*w3 - 84*CA2*v5*w3 +
     &        48*v6*w3 + 60*CA2*v6*w3 + 10*v7*w3 +
     &        30*CA2*v7*w3 - 32*CA2*v8*w3 + 2*v2*w4 -
     &        2*CA2*v2*w4 + 4*v3*w4 - 61*v4*w4 +
     &        17*CA2*v4*w4 + 84*v5*w4 - 10*CA2*v5*w4 -
     &        12*v6*w4 + 22*CA2*v6*w4 - 8*v7*w4 -
     &        84*CA2*v7*w4 - 8*v8*w4 + 56*CA2*v8*w4 +
     &        3*v4*w5 - 3*CA2*v4*w5 + 3*v5*w5 -
     &        3*CA2*v5*w5 - 11*v6*w5 - 13*CA2*v6*w5 -
     &        19*v7*w5 + 79*CA2*v7*w5 + 24*v8*w5 -
     &        60*CA2*v8*w5 - v4*w6 + CA2*v4*w6 -
     &        3*v5*w6 + CA2*v5*w6 + 2*v6*w6 +
     &        4*CA2*v6*w6 + 23*v7*w6 - 43*CA2*v7*w6 -
     &        26*v8*w6 + 42*CA2*v8*w6 - 6*v7*w7 +
     &        10*CA2*v7*w7 + 12*v8*w7 - 16*CA2*v8*w7 -
     &        2*v8*w8 + 2*CA2*v8*w8))/
     &    ((1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      struv9= part1 + part2 + part3 + part4 +
     &        part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV10(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (-2*CF*lmss*(2 - 4*v + 2*v2 + 2*v*w - 
     &     2*v2*w + v2*w2)*
     &      (1 - 2*v + 2*v2 + 2*v*w - 4*v2*w + 2*v2*w2)*
     &      (2*CA2 - 4*CA2*v + 2*CA2*v2 + 2*CA2*v*w -
     &        2*CA2*v2*w - v2*w2 + CA2*v2*w2))/
     &    ((1 - v)*v2*w2*(1 - v + v*w)**2) 
      part2 = (4*CF*l1vw*(1 - w)*(2 - CA2 - 8*v + 2*CA2*v + 10*v2 -
     &        CA2*v2 - 4*v3 + 7*v*w - 2*CA2*v*w - 19*v2*w +
     &        5*CA2*v2*w + 12*v3*w - 3*CA2*v3*w + 
     &     8*v2*w2 -
     &        3*CA2*v2*w2 - 10*v3*w2 + 4*CA2*v3*w2 +
     &        2*v3*w3 - CA2*v3*w3))/(v*w*(1 - v + v*w)**2) 
      part3 = (16*CF**2*lvw*(1 - w)*(CA - 2*CA*v + 2*CA*v2 + 
     &     v*w + 2*CA*v*w -
     &        2*v2*w - 7*CA*v2*w + 2*v3*w + 4*CA*v3*w +
     &        4*CA*v2*w2 - 2*v3*w2 - 5*CA*v3*w2 +
     &        2*CA*v4*w2 + v3*w3 + CA*v3*w3 - 
     &     2*CA*v4*w3 +
     &        CA*v4*w4))/((1 - v)*v2*w2) 
      part4 = -(4*CF*lw*(2*CA - 4*CA**3 - 6*CA*v + 18*CA**3*v + 
     &   8*CA*v2 -
     &        38*CA**3*v2 - 4*CA*v3 + 44*CA**3*v3 - 
     &     28*CA**3*v4 +
     &        8*CA**3*v5 + w - CA2*w - 3*v*w + 3*CA2*v*w -
     &        6*CA**3*v*w + 4*v2*w - 4*CA2*v2*w + 
     &     26*CA**3*v2*w -
     &        2*v3*w - 2*CA*v3*w + 2*CA2*v3*w - 
     &     38*CA**3*v3*w +
     &        4*CA*v4*w + 16*CA**3*v4*w + 8*CA**3*v5*w -
     &        8*CA**3*v6*w + v*w2 + CA2*v*w2 - 6*v2*w2 -
     &        5*CA*v2*w2 - 4*CA2*v2*w2 - 3*CA**3*v2*w2 +
     &        10*v3*w2 + 19*CA*v3*w2 + 10*CA2*v3*w2 +
     &        CA**3*v3*w2 - 8*v4*w2 - 26*CA*v4*w2 -
     &        12*CA2*v4*w2 + 26*CA**3*v4*w2 + 4*v5*w2 +
     &        8*CA*v5*w2 + 4*CA2*v5*w2 - 40*CA**3*v5*w2 +
     &        20*CA**3*v6*w2 + v2*w3 + CA2*v2*w3 -
     &        4*v3*w3 - 8*CA*v3*w3 - 6*CA2*v3*w3 +
     &        4*CA**3*v3*w3 + 6*v4*w3 + 23*CA*v4*w3 +
     &        14*CA2*v4*w3 - 27*CA**3*v4*w3 - 6*v5*w3 -
     &        12*CA*v5*w3 - 6*CA2*v5*w3 + 
     &     40*CA**3*v5*w3 -
     &        20*CA**3*v6*w3 + v3*w4 + CA2*v3*w4 -
     &        2*v4*w4 - 5*CA*v4*w4 - 8*CA2*v4*w4 +
     &        7*CA**3*v4*w4 + 4*v5*w4 + 4*CA*v5*w4 +
     &        4*CA2*v5*w4 - 16*CA**3*v5*w4 + 
     &     10*CA**3*v6*w4 +
     &        2*CA2*v4*w5 - v5*w5 - CA2*v5*w5 +
     &        2*CA**3*v5*w5 - 2*CA**3*v6*w5))/
     &    (CA*(1 - v)*v2*w2*(1 - v*w)*(1 - v + v*w)) 
      part5 = -(4*CF*l1v*(2*CA**3 - 8*CA**3*v + 14*CA**3*v2 - 
     &     12*CA**3*v3 +
     &        4*CA**3*v4 + w + CA2*w - 5*v*w +
     &      6*CA*v*w - 3*CA2*v*w +
     &        CA**3*v*w + 10*v2*w - 24*CA*v2*w + 4*CA2*v2*w -
     &        4*CA**3*v2*w - 10*v3*w + 46*CA*v3*w - 
     &     2*CA2*v3*w -
     &        CA**3*v3*w + 4*v4*w - 44*CA*v4*w + 
     &     12*CA**3*v4*w +
     &        16*CA*v5*w - 8*CA**3*v5*w + v*w2 + CA2*v*w2 -
     &        4*v2*w2 + 10*CA*v2*w2 - 4*CA2*v2*w2 +
     &        CA**3*v2*w2 + 6*v3*w2 - 46*CA*v3*w2 +
     &        4*CA2*v3*w2 + 5*CA**3*v3*w2 + 72*CA*v4*w2 -
     &        18*CA**3*v4*w2 - 4*v5*w2 - 36*CA*v5*w2 +
     &        8*CA**3*v5*w2 + 4*CA**3*v6*w2 + v2*w3 +
     &        CA2*v2*w3 - 4*v3*w3 + 14*CA*v3*w3 -
     &        4*CA2*v3*w3 - 3*CA**3*v3*w3 - 44*CA*v4*w3 +
     &        2*CA2*v4*w3 + 12*CA**3*v4*w3 + 6*v5*w3 +
     &        30*CA*v5*w3 - 2*CA2*v5*w3 - CA**3*v5*w3 -
     &        8*CA**3*v6*w3 + v3*w4 + CA2*v3*w4 +
     &        10*CA*v4*w4 - 3*CA**3*v4*w4 - 4*v5*w4 -
     &        10*CA*v5*w4 + 2*CA2*v5*w4 - 3*CA**3*v5*w4 +
     &        6*CA**3*v6*w4 + v5*w5 - CA2*v5*w5 +
     &        2*CA**3*v5*w5 - 2*CA**3*v6*w5))/
     &    (CA*(1 - v)*v2*w2*(1 - v*w)*(1 - v + v*w)) 
      part6 = (2*CF*lms*(4*CA - 4*CA**3 - 8*CA*v + 12*CA**3*v + 
     &     8*CA*v2 -
     &        20*CA**3*v2 + 16*CA**3*v3 - 8*CA**3*v4 - 2*CA*w +
     &        2*CA**3*w + 2*v*w - 2*CA*v*w - 2*CA2*v*w - 
     &     2*CA**3*v*w -
     &        4*v2*w + 4*CA*v2*w + 4*CA2*v2*w + 4*v3*w -
     &        12*CA*v3*w - 4*CA2*v3*w + 20*CA**3*v3*w -
     &        24*CA**3*v4*w + 16*CA**3*v5*w + CA*w2 - 
     &     CA**3*w2 +
     &        2*CA*v*w2 - 4*v2*w2 + 4*CA2*v2*w2 +
     &        2*CA**3*v2*w2 + 4*v3*w2 + 4*CA*v3*w2 -
     &        4*CA2*v3*w2 - 16*CA**3*v3*w2 - 8*v4*w2 +
     &        8*CA*v4*w2 + 8*CA2*v4*w2 + 8*CA**3*v4*w2 -
     &        8*CA**3*v6*w2 - 2*CA*v*w3 + 2*CA**3*v*w3 +
     &        2*CA*v2*w3 - 6*CA**3*v2*w3 + 4*v3*w3 -
     &        4*CA*v3*w3 - 4*CA2*v3*w3 + 12*CA**3*v3*w3 +
     &        4*v4*w3 - 4*CA*v4*w3 - 4*CA2*v4*w3 +
     &        4*v5*w3 - 4*CA*v5*w3 - 4*CA2*v5*w3 -
     &        4*CA**3*v5*w3 + 8*CA**3*v6*w3 + CA*v2*w4 -
     &        CA**3*v2*w4 - 2*CA*v3*w4 + 4*CA**3*v3*w4 -
     &        4*v4*w4 + 6*CA*v4*w4 + 4*CA2*v4*w4 -
     &        12*CA**3*v4*w4 - 4*v5*w4 + 4*CA2*v5*w4 +
     &        8*CA**3*v5*w4 + 4*CA*v6*w4 - 8*CA**3*v6*w4 +
     &        2*v5*w5 - 2*CA*v5*w5 - 2*CA2*v5*w5 +
     &        2*CA**3*v5*w5 - 4*CA*v6*w5 + 4*CA**3*v6*w5 +
     &        2*CA*v6*w6 - 2*CA**3*v6*w6))/
     &    (CA*(1 - v)*v2*w2*(1 - v*w)**2) 
      part7 = -(2*CF*l1w*(4*CA - 12*CA**3 - 16*CA*v + 64*CA**3*v + 
     &     28*CA*v2 -
     &        156*CA**3*v2 - 24*CA*v3 + 216*CA**3*v3 + 
     &     8*CA*v4 -
     &        176*CA**3*v4 + 80*CA**3*v5 - 16*CA**3*v6 - 
     &     2*w - 2*CA*w +
     &        2*CA2*w + 2*CA**3*w + 10*v*w - 2*CA*v*w - 
     &     14*CA2*v*w -
     &        24*CA**3*v*w - 22*v2*w + 42*CA*v2*w + 
     &     38*CA2*v2*w +
     &        84*CA**3*v2*w + 26*v3*w - 130*CA*v3*w - 
     &     54*CA2*v3*w -
     &        108*CA**3*v3*w - 16*v4*w + 192*CA*v4*w +
     &        40*CA2*v4*w - 6*CA**3*v4*w + 4*v5*w - 
     &     132*CA*v5*w -
     &        12*CA2*v5*w + 132*CA**3*v5*w + 32*CA*v6*w -
     &        112*CA**3*v6*w + 32*CA**3*v7*w + CA*w2 - 
     &     CA**3*w2 -
     &        2*v*w2 - 4*CA*v*w2 - 2*CA2*v*w2 + 
     &     6*CA**3*v*w2 +
     &        8*v2*w2 - 17*CA*v2*w2 + 12*CA2*v2*w2 -
     &        13*CA**3*v2*w2 - 10*v3*w2 + 102*CA*v3*w2 -
     &        30*CA2*v3*w2 - 38*CA**3*v3*w2 - 4*v4*w2 -
     &        198*CA*v4*w2 + 52*CA2*v4*w2 + 
     &     232*CA**3*v4*w2 +
     &        16*v5*w2 + 124*CA*v5*w2 - 56*CA2*v5*w2 -
     &        354*CA**3*v5*w2 - 8*v6*w2 + 28*CA*v6*w2 +
     &        24*CA2*v6*w2 + 196*CA**3*v6*w2 - 
     &     32*CA*v7*w2 -
     &        16*CA**3*v7*w2 - 16*CA**3*v8*w2 + 
     &     6*CA*v2*w3 -
     &        4*CA2*v2*w3 - 6*CA**3*v2*w3 - 
     &     30*CA*v3*w3 +
     &        20*CA2*v3*w3 + 48*CA**3*v3*w3 + 12*v4*w3 +
     &        60*CA*v4*w3 - 52*CA2*v4*w3 - 
     &     182*CA**3*v4*w3 -
     &        24*v5*w3 + 32*CA*v5*w3 + 64*CA2*v5*w3 +
     &        242*CA**3*v5*w3 + 8*v6*w3 - 172*CA*v6*w3 -
     &        16*CA2*v6*w3 - 46*CA**3*v6*w3 + 4*v7*w3 +
     &        92*CA*v7*w3 - 12*CA2*v7*w3 - 
     &     92*CA**3*v7*w3 +
     &        48*CA**3*v8*w3 - 2*CA*v2*w4 + 2*CA**3*v2*w4 +
     &        6*CA*v3*w4 - 10*CA**3*v3*w4 - 8*v4*w4 -
     &        2*CA*v4*w4 + 8*CA2*v4*w4 + 40*CA**3*v4*w4 +
     &        20*v5*w4 - 76*CA*v5*w4 - 12*CA2*v5*w4 -
     &        34*CA**3*v5*w4 - 4*v6*w4 + 202*CA*v6*w4 -
     &        28*CA2*v6*w4 - 108*CA**3*v6*w4 - 8*v7*w4 -
     &        120*CA*v7*w4 + 32*CA2*v7*w4 + 
     &     158*CA**3*v7*w4 +
     &        8*CA*v8*w4 - 64*CA**3*v8*w4 + 2*v4*w5 -
     &        4*CA*v4*w5 - 2*CA2*v4*w5 + 4*CA**3*v4*w5 -
     &        10*v5*w5 + 32*CA*v5*w5 - 2*CA2*v5*w5 -
     &        24*CA**3*v5*w5 + 2*v6*w5 - 110*CA*v6*w5 +
     &        38*CA2*v6*w5 + 98*CA**3*v6*w5 + 6*v7*w5 +
     &        90*CA*v7*w5 - 34*CA2*v7*w5 - 
     &     118*CA**3*v7*w5 -
     &        24*CA*v8*w5 + 56*CA**3*v8*w5 + CA*v4*w6 -
     &        CA**3*v4*w6 + 2*v5*w6 - 2*CA*v5*w6 +
     &        2*CA2*v5*w6 + 4*CA**3*v5*w6 + 23*CA*v6*w6 -
     &        20*CA2*v6*w6 - 23*CA**3*v6*w6 - 2*v7*w6 -
     &        34*CA*v7*w6 + 18*CA2*v7*w6 + 
     &     40*CA**3*v7*w6 +
     &        28*CA*v8*w6 - 36*CA**3*v8*w6 + 
     &     4*CA2*v6*w7 +
     &        4*CA*v7*w7 - 4*CA2*v7*w7 - 4*CA**3*v7*w7 -
     &        16*CA*v8*w7 + 16*CA**3*v8*w7 + 4*CA*v8*w8 -
     &        4*CA**3*v8*w8))/
     &    (CA*(1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part8 = -(2*CF*lv*(4*CA - 16*CA**3 - 16*CA*v + 84*CA**3*v + 
     &     28*CA*v2 -
     &        200*CA**3*v2 - 24*CA*v3 + 268*CA**3*v3 + 
     &     8*CA*v4 -
     &        208*CA**3*v4 + 88*CA**3*v5 - 16*CA**3*v6 - 2*CA*w +
     &        2*CA**3*w + 2*v*w - 6*CA*v*w - 6*CA2*v*w - 
     &     30*CA**3*v*w -
     &        8*v2*w + 62*CA*v2*w + 24*CA2*v2*w + 
     &     110*CA**3*v2*w +
     &        14*v3*w - 166*CA*v3*w - 42*CA2*v3*w -
     &        150*CA**3*v3*w - 12*v4*w + 220*CA*v4*w +
     &        36*CA2*v4*w + 16*CA**3*v4*w + 4*v5*w - 
     &     140*CA*v5*w -
     &        12*CA2*v5*w + 140*CA**3*v5*w + 32*CA*v6*w -
     &        120*CA**3*v6*w + 32*CA**3*v7*w + CA*w2 - 
     &     CA**3*w2 -
     &        4*CA*v*w2 + 6*CA**3*v*w2 - 4*v2*w2 - 
     &     29*CA*v2*w2 -
     &        11*CA**3*v2*w2 + 16*v3*w2 + 146*CA*v3*w2 +
     &        4*CA2*v3*w2 - 56*CA**3*v3*w2 - 32*v4*w2 -
     &        250*CA*v4*w2 + 298*CA**3*v4*w2 + 36*v5*w2 +
     &        144*CA*v5*w2 - 20*CA2*v5*w2 - 
     &     452*CA**3*v5*w2 -
     &        16*v6*w2 + 28*CA*v6*w2 + 16*CA2*v6*w2 +
     &        252*CA**3*v6*w2 - 32*CA*v7*w2 - 
     &     24*CA**3*v7*w2 -
     &        16*CA**3*v8*w2 + 6*CA*v2*w3 - 6*CA**3*v2*w3 -
     &        2*v3*w3 - 38*CA*v3*w3 + 2*CA2*v3*w3 +
     &        56*CA**3*v3*w3 + 12*v4*w3 + 68*CA*v4*w3 -
     &        12*CA2*v4*w3 - 230*CA**3*v4*w3 - 24*v5*w3 +
     &        52*CA*v5*w3 + 40*CA2*v5*w3 + 
     &     332*CA**3*v5*w3 +
     &        4*v6*w3 - 200*CA*v6*w3 - 28*CA2*v6*w3 -
     &        96*CA**3*v6*w3 + 12*v7*w3 + 100*CA*v7*w3 -
     &        4*CA2*v7*w3 - 100*CA**3*v7*w3 + 
     &     56*CA**3*v8*w3 -
     &        2*CA*v2*w4 + 2*CA**3*v2*w4 + 6*CA*v3*w4 -
     &        10*CA**3*v3*w4 + 6*CA*v4*w4 + 
     &     48*CA**3*v4*w4 +
     &        4*v5*w4 - 116*CA*v5*w4 - 20*CA2*v5*w4 -
     &        58*CA**3*v5*w4 + 16*v6*w4 + 254*CA*v6*w4 +
     &        16*CA2*v6*w4 - 114*CA**3*v6*w4 - 28*v7*w4 -
     &        140*CA*v7*w4 + 12*CA2*v7*w4 + 
     &     204*CA**3*v7*w4 +
     &        8*CA*v8*w4 - 88*CA**3*v8*w4 - 4*CA*v4*w5 +
     &        4*CA**3*v4*w5 - 2*v5*w5 + 44*CA*v5*w5 +
     &        6*CA2*v5*w5 - 22*CA**3*v5*w5 - 12*v6*w5 -
     &        138*CA*v6*w5 - 4*CA2*v6*w5 + 
     &     116*CA**3*v6*w5 +
     &        26*v7*w5 + 106*CA*v7*w5 - 14*CA2*v7*w5 -
     &        166*CA**3*v7*w5 - 24*CA*v8*w5 + 
     &     84*CA**3*v8*w5 +
     &        CA*v4*w6 - CA**3*v4*w6 - 2*CA*v5*w6 +
     &        4*CA**3*v5*w6 + 4*v6*w6 + 27*CA*v6*w6 -
     &        29*CA**3*v6*w6 - 12*v7*w6 - 38*CA*v7*w6 +
     &        8*CA2*v7*w6 + 62*CA**3*v7*w6 + 
     &     28*CA*v8*w6 -
     &        52*CA**3*v8*w6 + 2*v7*w7 + 4*CA*v7*w7 -
     &        2*CA2*v7*w7 - 8*CA**3*v7*w7 - 16*CA*v8*w7 +
     &        20*CA**3*v8*w7 + 4*CA*v8*w8 - 
     &     4*CA**3*v8*w8))/
     &    (CA*(1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2) 
      part9 = -(2*CF*(4*CA*v - 12*CA**3*v - 12*CA*v2 + 52*CA**3*v2 +
     &        12*CA*v3 - 92*CA**3*v3 - 4*CA*v4 + 84*CA**3*v4 -
     &        40*CA**3*v5 + 8*CA**3*v6 + CA*w - 
     &     CA**3*w - 9*CA*v*w +
     &        17*CA**3*v*w + 4*v2*w + 27*CA*v2*w - 4*CA2*v2*w -
     &        75*CA**3*v2*w - 12*v3*w - 29*CA*v3*w + 
     &     12*CA2*v3*w +
     &        125*CA**3*v3*w + 12*v4*w + 8*CA*v4*w - 
     &     12*CA2*v4*w -
     &        72*CA**3*v4*w - 4*v5*w + 2*CA*v5*w + 
     &     4*CA2*v5*w -
     &        26*CA**3*v5*w + 48*CA**3*v6*w - 
     &     16*CA**3*v7*w - CA*w2 +
     &        CA**3*w2 - CA*v*w2 - CA**3*v*w2 - 4*v2*w2 -
     &        CA*v2*w2 + 4*CA2*v2*w2 + 7*CA**3*v2*w2 +
     &        12*v3*w2 + 41*CA*v3*w2 - 12*CA2*v3*w2 -
     &        3*CA**3*v3*w2 - 4*v4*w2 - 94*CA*v4*w2 +
     &        4*CA2*v4*w2 - 80*CA**3*v4*w2 - 12*v5*w2 +
     &        84*CA*v5*w2 + 12*CA2*v5*w2 + 
     &     176*CA**3*v5*w2 +
     &        8*v6*w2 - 28*CA*v6*w2 - 8*CA2*v6*w2 -
     &        132*CA**3*v6*w2 + 24*CA**3*v7*w2 + 
     &     8*CA**3*v8*w2 -
     &        4*CA*v2*w3 + 4*CA**3*v2*w3 - 38*CA*v3*w3 -
     &        6*CA**3*v3*w3 - 16*v4*w3 + 156*CA*v4*w3 +
     &        16*CA2*v4*w3 + 28*CA**3*v4*w3 + 32*v5*w3 -
     &        172*CA*v5*w3 - 32*CA2*v5*w3 - 
     &     84*CA**3*v5*w3 -
     &        12*v6*w3 + 48*CA*v6*w3 + 12*CA2*v6*w3 +
     &        60*CA**3*v6*w3 - 4*v7*w3 + 10*CA*v7*w3 +
     &        4*CA2*v7*w3 + 30*CA**3*v7*w3 - 
     &     32*CA**3*v8*w3 +
     &        2*CA*v2*w4 - 2*CA**3*v2*w4 + 4*CA*v3*w4 +
     &        8*v4*w4 - 61*CA*v4*w4 - 8*CA2*v4*w4 +
     &        17*CA**3*v4*w4 - 16*v5*w4 + 84*CA*v5*w4 +
     &        16*CA2*v5*w4 - 10*CA**3*v5*w4 - 4*v6*w4 -
     &        12*CA*v6*w4 + 4*CA2*v6*w4 + 
     &     22*CA**3*v6*w4 +
     &        12*v7*w4 - 8*CA*v7*w4 - 12*CA2*v7*w4 -
     &        84*CA**3*v7*w4 - 8*CA*v8*w4 +
     &      56*CA**3*v8*w4 +
     &        3*CA*v4*w5 - 3*CA**3*v4*w5 + 3*CA*v5*w5 -
     &        3*CA**3*v5*w5 + 12*v6*w5 - 11*CA*v6*w5 -
     &        12*CA2*v6*w5 - 13*CA**3*v6*w5 - 12*v7*w5 -
     &        19*CA*v7*w5 + 12*CA2*v7*w5 + 
     &     79*CA**3*v7*w5 +
     &        24*CA*v8*w5 - 60*CA**3*v8*w5 - CA*v4*w6 +
     &        CA**3*v4*w6 - 3*CA*v5*w6 + CA**3*v5*w6 -
     &        4*v6*w6 + 2*CA*v6*w6 + 4*CA2*v6*w6 +
     &        4*CA**3*v6*w6 + 4*v7*w6 + 23*CA*v7*w6 -
     &        4*CA2*v7*w6 - 43*CA**3*v7*w6 - 
     &     26*CA*v8*w6 +
     &        42*CA**3*v8*w6 - 6*CA*v7*w7 +
     &      10*CA**3*v7*w7 +
     &        12*CA*v8*w7 - 16*CA**3*v8*w7 - 2*CA*v8*w8 +
     &        2*CA**3*v8*w8))/
     &    (CA*(1 - v)*v2*w2*(1 - v*w)**2*(1 - v + v*w)**2)
      struv10= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV11(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (-4*CF*lvw*(CA - CA**3 + CA*v2 - CA**3*v2 + 
     &     2*w + 7*CA*w -
     &        7*CA**3*w - 4*v*w - 16*CA*v*w + 13*CA**3*v*w + 3*v2*w +
     &        8*CA*v2*w - 3*CA2*v2*w - 9*CA**3*v2*w - v3*w -
     &        3*CA*v3*w + 3*CA2*v3*w + 3*CA**3*v3*w - 
     &     CA**3*v*w2 +
     &        7*CA*v2*w2 + CA2*v2*w2 - 4*CA**3*v2*w2 -
     &        5*CA*v3*w2 - CA2*v3*w2 + 3*CA**3*v3*w2 +
     &        2*CA*v4*w2 - 2*CA**3*v4*w2 + 4*CA*v2*w3 -
     &        4*CA2*v2*w3 - 8*CA**3*v2*w3 - 8*CA*v3*w3 +
     &        4*CA2*v3*w3 + 8*CA**3*v3*w3 + 2*CA*v4*w4 -
     &        2*CA**3*v4*w4))/(CA*(1 - v)**2*v*w) 
      part2 = -(4*CF*l1v*(v - CA2*v - 2*v2 + 2*CA2*v2 + v3 - 
     &     CA2*v3 +
     &        4*CA*w - 4*CA**3*w - 6*v*w - 30*CA*v*w - 2*CA2*v*w +
     &        10*CA**3*v*w + 14*v2*w + 60*CA*v2*w + 
     &     4*CA2*v2*w -
     &        16*CA**3*v2*w - 9*v3*w - 48*CA*v3*w - 
     &     3*CA2*v3*w +
     &        14*CA**3*v3*w + v4*w + 14*CA*v4*w + CA2*v4*w -
     &        4*CA**3*v4*w - 2*v*w2 - 28*CA*v*w2 - 
     &     6*CA2*v*w2 +
     &        4*v2*w2 + 68*CA*v2*w2 + 12*CA2*v2*w2 -
     &        10*CA**3*v2*w2 - 6*v3*w2 - 76*CA*v3*w2 -
     &        8*CA2*v3*w2 + 15*CA**3*v3*w2 + 4*v4*w2 +
     &        48*CA*v4*w2 + 2*CA2*v4*w2 - 
     &     13*CA**3*v4*w2 -
     &        14*CA*v5*w2 + 4*CA**3*v5*w2 + 2*v2*w3 +
     &        8*CA*v2*w3 - 2*CA**3*v2*w3 - 3*v3*w3 -
     &        18*CA*v3*w3 - 2*CA2*v3*w3 + 3*CA**3*v3*w3 +
     &        v4*w3 + 10*CA*v4*w3 + 2*CA2*v4*w3 -
     &        CA**3*v4*w3 + 2*v3*w4 + 2*CA2*v3*w4 -
     &        2*v4*w4 + 4*CA*v4*w4 - 2*CA2*v4*w4 -
     &        2*CA*v5*w4))/(CA*(1 - v)**2*w*(1 - v*w)*
     &     (1 - v + v*w)) 
      part3 = (4*CF*l1vw*(2 + 8*CA - 8*CA**3 + 8*CA2*CF - 
     &    10*v - 41*CA*v +
     &        42*CA**3*v - 40*CA2*CF*v + 21*v2 + 89*CA*v2 -
     &        3*CA2*v2 - 96*CA**3*v2 + 84*CA2*CF*v2 - 
     &     24*v3 -
     &        106*CA*v3 + 12*CA2*v3 + 124*CA**3*v3 -
     &        96*CA2*CF*v3 + 16*v4 + 74*CA*v4 - 18*CA2*v4 -
     &        96*CA**3*v4 + 64*CA2*CF*v4 - 6*v5 - 29*CA*v5 +
     &        12*CA2*v5 + 42*CA**3*v5 - 24*CA2*CF*v5 + v6 +
     &        5*CA*v6 - 3*CA2*v6 - 8*CA**3*v6 + 
     &     4*CA2*CF*v6 +
     &        6*v*w + 25*CA*v*w - 22*CA**3*v*w + 24*CA2*CF*v*w -
     &        24*v2*w - 102*CA*v2*w + CA2*v2*w + 
     &     92*CA**3*v2*w -
     &        104*CA2*CF*v2*w + 39*v3*w + 170*CA*v3*w -
     &        13*CA2*v3*w - 160*CA**3*v3*w + 
     &     184*CA2*CF*v3*w -
     &        33*v4*w - 148*CA*v4*w + 33*CA2*v4*w +
     &        148*CA**3*v4*w - 168*CA2*CF*v4*w + 15*v5*w +
     &        69*CA*v5*w - 31*CA2*v5*w - 74*CA**3*v5*w +
     &        80*CA2*CF*v5*w - 3*v6*w - 14*CA*v6*w +
     &        10*CA2*v6*w + 16*CA**3*v6*w - 16*CA2*CF*v6*w +
     &        6*v2*w2 + 33*CA*v2*w2 - 4*CA2*v2*w2 -
     &        30*CA**3*v2*w2 + 8*CA*CF*v2*w2 +
     &        52*CA2*CF*v2*w2 - 18*v3*w2 - 101*CA*v3*w2 +
     &        19*CA2*v3*w2 + 94*CA**3*v3*w2 - 
     &     24*CA*CF*v3*w2 -
     &        176*CA2*CF*v3*w2 + 21*v4*w2 + 
     &     119*CA*v4*w2 -
     &        42*CA2*v4*w2 - 118*CA**3*v4*w2 + 
     &     28*CA*CF*v4*w2 +
     &        232*CA2*CF*v4*w2 - 12*v5*w2 - 67*CA*v5*w2 +
     &        43*CA2*v5*w2 + 74*CA**3*v5*w2 - 
     &     16*CA*CF*v5*w2 -
     &        144*CA2*CF*v5*w2 + 3*v6*w2 + 16*CA*v6*w2 -
     &        16*CA2*v6*w2 - 20*CA**3*v6*w2 + 
     &     4*CA*CF*v6*w2 +
     &        36*CA2*CF*v6*w2 + 2*v3*w3 + 29*CA*v3*w3 -
     &        12*CA2*v3*w3 - 32*CA**3*v3*w3 + 
     &     16*CA*CF*v3*w3 +
     &        64*CA2*CF*v3*w3 - 4*v4*w3 - 66*CA*v4*w3 +
     &        39*CA2*v4*w3 + 86*CA**3*v4*w3 - 
     &     40*CA*CF*v4*w3 -
     &        168*CA2*CF*v4*w3 + 3*v5*w3 + 53*CA*v5*w3 -
     &        45*CA2*v5*w3 - 84*CA**3*v5*w3 + 
     &     36*CA*CF*v5*w3 +
     &        156*CA2*CF*v5*w3 - v6*w3 - 16*CA*v6*w3 +
     &        18*CA2*v6*w3 + 30*CA**3*v6*w3 - 
     &     12*CA*CF*v6*w3 -
     &        52*CA2*CF*v6*w3 + 21*CA*v4*w4 - 
     &     12*CA2*v4*w4 -
     &        28*CA**3*v4*w4 + 12*CA*CF*v4*w4 +
     &        48*CA2*CF*v4*w4 - 36*CA*v5*w4 + 
     &     25*CA2*v5*w4 +
     &        52*CA**3*v5*w4 - 24*CA*CF*v5*w4 -
     &        88*CA2*CF*v5*w4 + 17*CA*v6*w4 - 
     &     13*CA2*v6*w4 -
     &        26*CA**3*v6*w4 + 12*CA*CF*v6*w4 +
     &        44*CA2*CF*v6*w4 + 10*CA*v5*w5 - 
     &     4*CA2*v5*w5 -
     &        10*CA**3*v5*w5 + 4*CA*CF*v5*w5 +
     &        20*CA2*CF*v5*w5 - 10*CA*v6*w5 + 
     &     4*CA2*v6*w5 +
     &        10*CA**3*v6*w5 - 4*CA*CF*v6*w5 -
     &        20*CA2*CF*v6*w5 + 2*CA*v6*w6 - 
     &     2*CA**3*v6*w6 +
     &        4*CA2*CF*v6*w6))/(CA*(1 - v)**2*v*
     &     (1 - v + v*w)**3) 
      part4 = -(2*CF*lmss*(2*CA*CF + 2*v - 2*CA2*v - 10*CA*CF*v - 
     &      11*v2 +
     &        11*CA2*v2 + 22*CA*CF*v2 + 27*v3 - 27*CA2*v3 -
     &        26*CA*CF*v3 - 37*v4 + 37*CA2*v4 + 16*CA*CF*v4 +
     &        29*v5 - 29*CA2*v5 - 4*CA*CF*v5 - 12*v6 +
     &        12*CA2*v6 + 2*v7 - 2*CA2*v7 - 2*v*w + 
     &     4*CA*CF*v*w +
     &        16*v2*w - 4*CA2*v2*w - 8*CF*v2*w - 
     &     40*CA*CF*v2*w -
     &        54*v3*w + 22*CA2*v3*w + 32*CF*v3*w + 
     &     116*CA*CF*v3*w +
     &        94*v4*w - 46*CA2*v4*w - 48*CF*v4*w - 
     &     156*CA*CF*v4*w -
     &        88*v5*w + 46*CA2*v5*w + 32*CF*v5*w + 
     &     112*CA*CF*v5*w +
     &        42*v6*w - 22*CA2*v6*w - 8*CF*v6*w - 
     &     44*CA*CF*v6*w -
     &        8*v7*w + 4*CA2*v7*w + 8*CA*CF*v7*w - 
     &     5*v2*w2 +
     &        CA2*v2*w2 - 8*CF*v2*w2 - 22*CA*CF*v2*w2 +
     &        35*v3*w2 - 11*CA2*v3*w2 + 8*CF*v3*w2 +
     &        30*CA*CF*v3*w2 - 88*v4*w2 + 32*CA2*v4*w2 +
     &        24*CF*v4*w2 + 36*CA*CF*v4*w2 + 106*v5*w2 -
     &        42*CA2*v5*w2 - 40*CF*v5*w2 - 
     &     88*CA*CF*v5*w2 -
     &        62*v6*w2 + 26*CA2*v6*w2 + 16*CF*v6*w2 +
     &        60*CA*CF*v6*w2 + 14*v7*w2 - 6*CA2*v7*w2 -
     &        16*CA*CF*v7*w2 - 8*v3*w3 + 2*CA2*v3*w3 -
     &        16*CF*v3*w3 - 52*CA*CF*v3*w3 + 38*v4*w3 -
     &        14*CA2*v4*w3 + 28*CF*v4*w3 +
     &      92*CA*CF*v4*w3 -
     &        70*v5*w3 + 32*CA2*v5*w3 - 8*CF*v5*w3 -
     &        48*CA*CF*v5*w3 + 56*v6*w3 - 28*CA2*v6*w3 -
     &        4*CF*v6*w3 - 16*v7*w3 + 8*CA2*v7*w3 +
     &        8*CA*CF*v7*w3 - 7*v4*w4 + 3*CA2*v4*w4 -
     &        12*CF*v4*w4 - 44*CA*CF*v4*w4 + 29*v5*w4 -
     &        13*CA2*v5*w4 + 20*CF*v5*w4 + 
     &     64*CA*CF*v5*w4 -
     &        36*v6*w4 + 16*CA2*v6*w4 - 8*CF*v6*w4 -
     &        28*CA*CF*v6*w4 + 14*v7*w4 - 6*CA2*v7*w4 -
     &        6*v5*w5 + 2*CA2*v5*w5 - 4*CF*v5*w5 -
     &        20*CA*CF*v5*w5 + 14*v6*w5 - 6*CA2*v6*w5 +
     &        4*CF*v6*w5 + 16*CA*CF*v6*w5 - 8*v7*w5 +
     &        4*CA2*v7*w5 - 2*v6*w6 + 2*CA2*v6*w6 -
     &        4*CA*CF*v6*w6 + 2*v7*w6 - 2*CA2*v7*w6))/
     &    ((1 - v)**2*v*w*(1 - v + v*w)**3) 
      part5 = -(2*CF*lms*(6*CA*CF - 8*CA*CF*v + 4*CF*v2 + 18*CA*CF*v2 -
     &        4*CF*v3 - 12*CA*CF*v3 + 4*CA*CF*v4 - 
     &     4*CA*CF*w + 2*v*w -
     &        8*CA*CF*v*w - 12*CF*v2*w - 16*CA*CF*v2*w + 2*v3*w -
     &        16*CA*CF*v3*w + 12*CF*v4*w + 24*CA*CF*v4*w -
     &        12*CA*CF*v5*w + 12*CA*CF*v*w2 - 5*v2*w2 +
     &        CA2*v2*w2 - 8*CF*v2*w2 - 34*CA*CF*v2*w2 -
     &        4*v3*w2 + 44*CF*v3*w2 + 116*CA*CF*v3*w2 -
     &        v4*w2 + CA2*v4*w2 - 24*CF*v4*w2 -
     &        58*CA*CF*v4*w2 - 12*CF*v5*w2 + 
     &     12*CA*CF*v6*w2 -
     &        12*CA*CF*v2*w3 + 8*v3*w3 - 2*CA2*v3*w3 +
     &        16*CF*v3*w3 + 72*CA*CF*v3*w3 + 2*v4*w3 -
     &        2*CA2*v4*w3 - 48*CF*v4*w3 - 
     &     164*CA*CF*v4*w3 +
     &        2*v5*w3 + 28*CF*v5*w3 + 100*CA*CF*v5*w3 +
     &        4*CF*v6*w3 - 24*CA*CF*v6*w3 - 4*CA*CF*v7*w3 +
     &        4*CA*CF*v3*w4 - 7*v4*w4 + 3*CA2*v4*w4 -
     &        12*CF*v4*w4 - 52*CA*CF*v4*w4 - 4*v5*w4 +
     &        20*CF*v5*w4 + 92*CA*CF*v5*w4 - v6*w4 +
     &        CA2*v6*w4 - 8*CF*v6*w4 - 48*CA*CF*v6*w4 +
     &        12*CA*CF*v7*w4 + 6*v5*w5 - 2*CA2*v5*w5 +
     &        4*CF*v5*w5 + 20*CA*CF*v5*w5 + 2*v6*w5 -
     &        2*CA2*v6*w5 - 4*CF*v6*w5 - 24*CA*CF*v6*w5 +
     &        4*CA*CF*v7*w5 - 2*v6*w6 + 2*CA2*v6*w6 -
     &        4*CA*CF*v6*w6 + 4*CA*CF*v7*w6))/
     &    ((1 - v)**2*v*w*(1 - v*w)**3) 
      part6 = -(CF*(4 + 10*CA - 4*CA2 - 10*CA**3 + 8*CA*CF + 8*CA2*CF -
     &   20*v - 46*CA*v + 20*CA2*v + 46*CA**3*v - 40*CA*CF*v -
     &   40*CA2*CF*v + 48*v2 + 88*CA*v2 - 48*CA2*v2 - 88*CA**3*v2 +
     &   88*CA*CF*v2 + 80*CA2*CF*v2 - 72*v3 - 88*CA*v3 + 72*CA2*v3 +
     &   88*CA**3*v3 - 112*CA*CF*v3 - 80*CA2*CF*v3 + 68*v4 + 46*CA*v4 -
     &   68*CA2*v4 - 46*CA**3*v4 + 88*CA*CF*v4 + 40*CA2*CF*v4 - 
     &   36*v5 - 10*CA*v5 + 36*CA2*v5 + 10*CA**3*v5 - 40*CA*CF*v5 -
     &   8*CA2*CF*v5 + 8*v6 - 8*CA2*v6 + 8*CA*CF*v6 - 4*w - 24*CA*w +
     &   4*CA2*w + 8*CA**3*w - 8*CA*CF*w - 8*CA2*CF*w + 20*v*w +
     &   118*CA*v*w - 20*CA2*v*w - 42*CA**3*v*w + 40*CA*CF*v*w +
     &   40*CA2*CF*v*w - 44*v2*w - 199*CA*v2*w + 28*CA2*v2*w +
     &   43*CA**3*v2*w - 80*CA*CF*v2*w - 40*CA2*CF*v2*w + 56*v3*w +
     &   104*CA*v3*w + 8*CA2*v3*w + 72*CA**3*v3*w + 80*CA*CF*v3*w -
     &   80*CA2*CF*v3*w - 20*v4*w + 80*CA*v4*w - 76*CA2*v4*w - 
     &   184*CA**3*v4*w - 16*CA*CF*v4*w + 200*CA2*CF*v4*w - 52*v5*w -
     &   128*CA*v5*w + 116*CA2*v5*w + 140*CA**3*v5*w - 64*CA*CF*v5*w -
     &   152*CA2*CF*v5*w + 68*v6*w + 55*CA*v6*w - 84*CA2*v6*w -
     &   35*CA**3*v6*w + 72*CA*CF*v6*w + 40*CA2*CF*v6*w - 24*v7*w -
     &   6*CA*v7*w + 24*CA2*v7*w - 2*CA**3*v7*w - 24*CA*CF*v7*w -
     &   4*v2*w2 + 54*CA*v2*w2 + 60*CA2*v2*w2 - 2*CA**3*v2*w2 -
     &   96*CA*CF*v2*w2 - 64*CA2*CF*v2*w2 + 16*v3*w2 - 234*CA*v3*w2 -
     &   240*CA2*v3*w2 + 18*CA**3*v3*w2 + 384*CA*CF*v3*w2 +
     &   256*CA2*CF*v3*w2 - 72*v4*w2 + 473*CA*v4*w2 + 388*CA2*v4*w2 - 
     &   109*CA**3*v4*w2 - 656*CA*CF*v4*w2 - 320*CA2*CF*v4*w2 + 
     &   160*v5*w2 - 567*CA*v5*w2 - 324*CA2*v5*w2 + 227*CA**3*v5*w2 +
     &   624*CA*CF*v5*w2 + 64*CA2*CF*v5*w2 - 124*v6*w2 + 411*CA*v6*w2 +
     &   120*CA2*v6*w2 - 195*CA**3*v6*w2 - 312*CA*CF*v6*w2 +
     &   128*CA2*CF*v6*w2 - 155*CA*v7*w2 + 20*CA2*v7*w2 + 
     &   55*CA**3*v7*w2 + 32*CA*CF*v7*w2 - 64*CA2*CF*v7*w2 + 
     &   24*v8*w2 + 18*CA*v8*w2 - 24*CA2*v8*w2 + 6*CA**3*v8*w2 +
     &   24*CA*CF*v8*w2 + 12*v2*w3 + 72*CA*v2*w3 - 12*CA2*v2*w3 -
     &   24*CA**3*v2*w3 + 24*CA*CF*v2*w3 + 24*CA2*CF*v2*w3 - 48*v3*w3 -
     &   282*CA*v3*w3 + 48*CA2*v3*w3 + 90*CA**3*v3*w3 - 
     &   96*CA*CF*v3*w3 - 96*CA2*CF*v3*w3 + 96*v4*w3 + 397*CA*v4*w3 -
     &   153*CA**3*v4*w3 + 56*CA*CF*v4*w3 - 120*v5*w3 - 233*CA*v5*w3 -
     &   168*CA2*v5*w3 + 169*CA**3*v5*w3 + 168*CA*CF*v5*w3 +
     &   336*CA2*CF*v5*w3 + 14*v6*w3 + 114*CA*v6*w3 + 270*CA2*v6*w3 - 
     &   170*CA**3*v6*w3 - 360*CA*CF*v6*w3 - 368*CA2*CF*v6*w3 + 
     &   116*v7*w3 - 171*CA*v7*w3 - 204*CA2*v7*w3 + 147*CA**3*v7*w3 +
     &   328*CA*CF*v7*w3 + 64*CA2*CF*v7*w3 - 62*v8*w3 + 125*CA*v8*w3 +
     &   58*CA2*v8*w3 - 57*CA**3*v8*w3 - 112*CA*CF*v8*w3 +
     &   40*CA2*CF*v8*w3 - 8*v9*w3 - 18*CA*v9*w3 + 8*CA2*v9*w3 -
     &   6*CA**3*v9*w3 - 8*CA*CF*v9*w3 + 20*v4*w4 + 118*CA*v4*w4 -
     &   84*CA2*v4*w4 + 18*CA**3*v4*w4 + 152*CA*CF*v4*w4 +
     &   104*CA2*CF*v4*w4 - 60*v5*w4 - 338*CA*v5*w4 + 252*CA2*v5*w4 -
     &   38*CA**3*v5*w4 - 456*CA*CF*v5*w4 - 312*CA2*CF*v5*w4 + 
     &   114*v6*w4 + 190*CA*v6*w4 - 278*CA2*v6*w4 + 86*CA**3*v6*w4 +
     &   536*CA*CF*v6*w4 + 200*CA2*CF*v6*w4 - 128*v7*w4 +
     &   192*CA*v7*w4 + 136*CA2*v7*w4 - 136*CA**3*v7*w4 -
     &   312*CA*CF*v7*w4 + 120*CA2*CF*v7*w4 + 24*v8*w4 - 167*CA*v8*w4 +
     &   4*CA2*v8*w4 + 55*CA**3*v8*w4 + 32*CA*CF*v8*w4 - 
     &   104*CA2*CF*v8*w4 + 30*v9*w4 - 11*CA*v9*w4 - 30*CA2*v9*w4 + 
     &   23*CA**3*v9*w4 + 48*CA*CF*v9*w4 - 8*CA2*CF*v9*w4 + 
     &   6*CA*v10*w4 + 2*CA**3*v10*w4 - 12*v4*w5 - 72*CA*v4*w5 +
     &   12*CA2*v4*w5 + 24*CA**3*v4*w5 - 24*CA*CF*v4*w5 -
     &   24*CA2*CF*v4*w5 + 36*v5*w5 + 214*CA*v5*w5 - 36*CA2*v5*w5 -
     &   82*CA**3*v5*w5 + 72*CA*CF*v5*w5 + 72*CA2*CF*v5*w5 - 12*v6*w5 -
     &   61*CA*v6*w5 - 20*CA2*v6*w5 + CA**3*v6*w5 + 16*CA*CF*v6*w5 +
     &   40*CA2*CF*v6*w5 - 36*v7*w5 - 250*CA*v7*w5 + 100*CA2*v7*w5 +
     &   158*CA**3*v7*w5 - 152*CA*CF*v7*w5 - 200*CA2*CF*v7*w5 + 
     &   56*v8*w5 + 130*CA*v8*w5 - 88*CA2*v8*w5 - 66*CA**3*v8*w5 +
     &   160*CA*CF*v8*w5 + 88*CA2*CF*v8*w5 - 32*v9*w5 +
     &   67*CA*v9*w5 + 32*CA2*v9*w5 - 47*CA**3*v9*w5 -
     &   72*CA*CF*v9*w5 + 24*CA2*CF*v9*w5 - 12*CA*v10*w5 -
     &   4*CA**3*v10*w5 - 36*v6*w6 - 118*CA*v6*w6 + 44*CA2*v6*w6 +
     &   82*CA**3*v6*w6 - 80*CA*CF*v6*w6 - 48*CA2*CF*v6*w6 + 72*v7*w6 +
     &   250*CA*v7*w6 - 88*CA2*v7*w6 - 154*CA**3*v7*w6 +
     &   160*CA*CF*v7*w6 + 96*CA2*CF*v7*w6 - 26*v8*w6 - 103*CA*v8*w6 +
     &   34*CA2*v8*w6 + 31*CA**3*v8*w6 - 96*CA*CF*v8*w6 - 
     &   24*CA2*CF*v8*w6 - 10*v9*w6 - 61*CA*v9*w6 + 10*CA2*v9*w6 + 
     &   57*CA**3*v9*w6 + 16*CA*CF*v9*w6 - 24*CA2*CF*v9*w6 + 
     &   6*CA*v10*w6 + 10*CA**3*v10*w6 + 4*v6*w7 + 24*CA*v6*w7 -
     &   4*CA2*v6*w7 - 8*CA**3*v6*w7 + 8*CA*CF*v6*w7 +
     &   8*CA2*CF*v6*w7 - 8*v7*w7 - 54*CA*v7*w7 + 8*CA2*v7*w7 +
     &   6*CA**3*v7*w7 - 16*CA*CF*v7*w7 - 16*CA2*CF*v7*w7 - 32*v8*w7 +
     &   23*CA*v8*w7 + 32*CA2*v8*w7 + 13*CA**3*v8*w7 - 24*CA*CF*v8*w7 +
     &        36*v9*w7 + 35*CA*v9*w7 - 36*CA2*v9*w7 -
     &        23*CA**3*v9*w7 + 32*CA*CF*v9*w7 +
     &        8*CA2*CF*v9*w7 - 16*CA**3*v10*w7 + 
     &     16*v8*w8 -
     &        16*CA2*v8*w8 + 8*CA**3*v8*w8 + 
     &     16*CA*CF*v8*w8 -
     &        16*v9*w8 - 16*CA*v9*w8 + 16*CA2*v9*w8 -
     &        16*CA*CF*v9*w8 + 8*CA**3*v10*w8 + 
     &     4*CA*v9*w9 -
     &        4*CA**3*v9*w9))/
     &    (CA*(1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w)**3) 
      part7 = -(2*CF*lv*(4*CA - 4*CA**3 - 20*CA*v + 20*CA**3*v - 2*v2 +
     &        52*CA*v2 + 2*CA2*v2 - 52*CA**3*v2 + 8*v3 -
     &        88*CA*v3 - 8*CA2*v3 + 88*CA**3*v3 - 12*v4 +
     &        96*CA*v4 + 12*CA2*v4 - 96*CA**3*v4 + 8*v5 -
     &        64*CA*v5 - 8*CA2*v5 + 64*CA**3*v5 - 2*v6 +
     &        24*CA*v6 + 2*CA2*v6 - 24*CA**3*v6 - 4*CA*v7 +
     &        4*CA**3*v7 - 2*CA*w + 2*CA**3*w + 10*CA*v*w - 
     &     10*CA**3*v*w +
     &        14*v2*w + 39*CA*v2*w + 18*CA2*v2*w + 
     &     35*CA**3*v2*w -
     &        56*v3*w - 216*CA*v3*w - 72*CA2*v3*w - 
     &     80*CA**3*v3*w +
     &        82*v4*w + 400*CA*v4*w + 114*CA2*v4*w +
     &        72*CA**3*v4*w - 50*v5*w - 402*CA*v5*w - 
     &     90*CA2*v5*w +
     &        22*CA**3*v5*w + 8*v6*w + 239*CA*v6*w + 
     &     36*CA2*v6*w -
     &        81*CA**3*v6*w + 2*v7*w - 80*CA*v7*w - 
     &     6*CA2*v7*w +
     &        52*CA**3*v7*w + 12*CA*v8*w - 12*CA**3*v8*w +
     &        70*CA*v2*w2 + 32*CA2*v2*w2 + 
     &     38*CA**3*v2*w2 -
     &        280*CA*v3*w2 - 128*CA2*v3*w2 - 
     &     152*CA**3*v3*w2 +
     &        36*v4*w2 + 571*CA*v4*w2 + 244*CA2*v4*w2 +
     &        329*CA**3*v4*w2 - 108*v5*w2 - 733*CA*v5*w2 -
     &        284*CA2*v5*w2 - 455*CA**3*v5*w2 + 
     &     110*v6*w2 +
     &        607*CA*v6*w2 + 194*CA2*v6*w2 + 
     &     341*CA**3*v6*w2 -
     &        40*v7*w2 - 319*CA*v7*w2 - 64*CA2*v7*w2 -
     &        101*CA**3*v7*w2 + 2*v8*w2 + 96*CA*v8*w2 +
     &        6*CA2*v8*w2 - 12*CA**3*v8*w2 - 
     &     12*CA*v9*w2 +
     &        12*CA**3*v9*w2 + 6*CA*v2*w3 - 6*CA**3*v2*w3 -
     &        24*CA*v3*w3 + 24*CA**3*v3*w3 - 18*v4*w3 +
     &        135*CA*v4*w3 + 26*CA2*v4*w3 - 
     &     23*CA**3*v4*w3 +
     &        54*v5*w3 - 321*CA*v5*w3 - 78*CA2*v5*w3 -
     &        15*CA**3*v5*w3 - 30*v6*w3 + 416*CA*v6*w3 +
     &        126*CA2*v6*w3 + 154*CA**3*v6*w3 - 
     &     30*v7*w3 -
     &        325*CA*v7*w3 - 122*CA2*v7*w3 -
     &      255*CA**3*v7*w3 +
     &        26*v8*w3 + 161*CA*v8*w3 + 50*CA2*v8*w3 +
     &        157*CA**3*v8*w3 - 2*v9*w3 - 48*CA*v9*w3 -
     &        2*CA2*v9*w3 - 36*CA**3*v9*w3 + 
     &     4*CA*v10*w3 -
     &        4*CA**3*v10*w3 - 12*v4*w4 - 172*CA*v4*w4 -
     &        68*CA2*v4*w4 - 60*CA**3*v4*w4 + 36*v5*w4 +
     &        516*CA*v5*w4 + 204*CA2*v5*w4 + 
     &     180*CA**3*v5*w4 -
     &        64*v6*w4 - 602*CA*v6*w4 - 240*CA2*v6*w4 -
     &        306*CA**3*v6*w4 + 68*v7*w4 + 344*CA*v7*w4 +
     &        140*CA2*v7*w4 + 312*CA**3*v7*w4 - 
     &     26*v8*w4 -
     &        91*CA*v8*w4 - 22*CA2*v8*w4 - 
     &     125*CA**3*v8*w4 -
     &        2*v9*w4 + 5*CA*v9*w4 - 14*CA2*v9*w4 -
     &        CA**3*v9*w4 + 8*CA*v10*w4 + 20*CA**3*v10*w4 -
     &        6*CA*v4*w5 + 6*CA**3*v4*w5 + 18*CA*v5*w5 -
     &        18*CA**3*v5*w5 - 14*v6*w5 - 193*CA*v6*w5 -
     &        66*CA2*v6*w5 + 15*CA**3*v6*w5 + 28*v7*w5 +
     &        356*CA*v7*w5 + 132*CA2*v7*w5 - 22*v8*w5 -
     &        278*CA*v8*w5 - 94*CA2*v8*w5 - 
     &     72*CA**3*v8*w5 +
     &        8*v9*w5 + 103*CA*v9*w5 + 28*CA2*v9*w5 +
     &        69*CA**3*v9*w5 - 20*CA*v10*w5 - 
     &     36*CA**3*v10*w5 +
     &        16*v6*w6 + 74*CA*v6*w6 + 48*CA2*v6*w6 +
     &        26*CA**3*v6*w6 - 32*v7*w6 - 148*CA*v7*w6 -
     &        96*CA2*v7*w6 - 52*CA**3*v7*w6 + 14*v8*w6 +
     &        105*CA*v8*w6 + 42*CA2*v8*w6 + 
     &     63*CA**3*v8*w6 +
     &        2*v9*w6 - 31*CA*v9*w6 + 6*CA2*v9*w6 -
     &        37*CA**3*v9*w6 - 4*CA*v10*w6 + 
     &     28*CA**3*v10*w6 +
     &        2*CA*v6*w7 - 2*CA**3*v6*w7 - 4*CA*v7*w7 +
     &        4*CA**3*v7*w7 + 10*v8*w7 + 11*CA*v8*w7 +
     &        30*CA2*v8*w7 + 21*CA**3*v8*w7 - 10*v9*w7 -
     &        9*CA*v9*w7 - 30*CA2*v9*w7 - 
     &     23*CA**3*v9*w7 +
     &        20*CA*v10*w7 - 12*CA**3*v10*w7 - 4*v8*w8 +
     &        8*CA*v8*w8 - 12*CA2*v8*w8 - 
     &     16*CA**3*v8*w8 +
     &        4*v9*w8 - 8*CA*v9*w8 + 12*CA2*v9*w8 +
     &        16*CA**3*v9*w8 - 12*CA*v10*w8 + 
     &     8*CA**3*v10*w8 +
     &        4*CA*v10*w9 - 4*CA**3*v10*w9))/
     &    (CA*(1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w)**3) 
      part8 = -(2*CF*l1w*(4*CA - 4*CA**3 - 20*CA*v + 20*CA**3*v + 
     &     52*CA*v2 -
     &        52*CA**3*v2 - 88*CA*v3 + 88*CA**3*v3 + 96*CA*v4 -
     &        96*CA**3*v4 - 64*CA*v5 + 64*CA**3*v5 + 24*CA*v6 -
     &        24*CA**3*v6 - 4*CA*v7 + 4*CA**3*v7 - 
     &     2*CA*w + 2*CA**3*w +
     &        10*CA*v*w - 10*CA**3*v*w - 6*v2*w - 45*CA*v2*w +
     &        14*CA2*v2*w + 51*CA**3*v2*w + 24*v3*w + 
     &     120*CA*v3*w -
     &        56*CA2*v3*w - 144*CA**3*v3*w - 32*v4*w -
     &        140*CA*v4*w + 84*CA2*v4*w + 176*CA**3*v4*w +
     &        12*v5*w + 42*CA*v5*w - 56*CA2*v5*w - 
     &     66*CA**3*v5*w +
     &        6*v6*w + 47*CA*v6*w + 14*CA2*v6*w - 
     &     41*CA**3*v6*w -
     &        4*v7*w - 44*CA*v7*w + 44*CA**3*v7*w + 
     &     12*CA*v8*w -
     &        12*CA**3*v8*w - 12*v2*w2 - 18*CA*v2*w2 +
     &        20*CA2*v2*w2 + 50*CA**3*v2*w2 + 48*v3*w2 +
     &        72*CA*v3*w2 - 80*CA2*v3*w2 - 
     &     200*CA**3*v3*w2 -
     &        100*v4*w2 - 253*CA*v4*w2 + 164*CA2*v4*w2 +
     &        459*CA**3*v4*w2 + 132*v5*w2 + 507*CA*v5*w2 -
     &        212*CA2*v5*w2 - 677*CA**3*v5*w2 - 
     &     88*v6*w2 -
     &        477*CA*v6*w2 + 152*CA2*v6*w2 + 
     &     551*CA**3*v6*w2 +
     &        12*v7*w2 + 193*CA*v7*w2 - 44*CA2*v7*w2 -
     &        207*CA**3*v7*w2 + 8*v8*w2 - 12*CA*v8*w2 +
     &        12*CA**3*v8*w2 - 12*CA*v9*w2 + 
     &     12*CA**3*v9*w2 +
     &        6*CA*v2*w3 - 6*CA**3*v2*w3 - 24*CA*v3*w3 +
     &        24*CA**3*v3*w3 + 115*CA*v4*w3 + 
     &     6*CA2*v4*w3 -
     &        29*CA**3*v4*w3 - 261*CA*v5*w3 - 
     &     18*CA2*v5*w3 +
     &        3*CA**3*v5*w3 - 44*v6*w3 + 88*CA*v6*w3 +
     &        62*CA2*v6*w3 + 188*CA**3*v6*w3 + 88*v7*w3 +
     &        231*CA*v7*w3 - 94*CA2*v7*w3 - 
     &     353*CA**3*v7*w3 -
     &        40*v8*w3 - 215*CA*v8*w3 + 44*CA2*v8*w3 +
     &        233*CA**3*v8*w3 - 4*v9*w3 + 60*CA*v9*w3 -
     &        60*CA**3*v9*w3 + 4*CA*v10*w3 - 
     &     4*CA**3*v10*w3 +
     &        20*v4*w4 + 12*CA*v4*w4 - 40*CA2*v4*w4 -
     &        84*CA**3*v4*w4 - 60*v5*w4 - 36*CA*v5*w4 +
     &        120*CA2*v5*w4 + 252*CA**3*v5*w4 + 
     &     106*v6*w4 +
     &        266*CA*v6*w4 - 162*CA2*v6*w4 - 
     &     446*CA**3*v6*w4 -
     &        112*v7*w4 - 472*CA*v7*w4 + 124*CA2*v7*w4 +
     &        472*CA**3*v7*w4 + 24*v8*w4 + 241*CA*v8*w4 -
     &        28*CA2*v8*w4 - 199*CA**3*v8*w4 + 22*v9*w4 -
     &        11*CA*v9*w4 - 14*CA2*v9*w4 + 
     &     5*CA**3*v9*w4 -
     &        28*CA*v10*w4 + 28*CA**3*v10*w4 - 6*CA*v4*w5 +
     &        6*CA**3*v4*w5 + 18*CA*v5*w5 - 
     &     18*CA**3*v5*w5 +
     &        2*v6*w5 - 61*CA*v6*w5 - 30*CA2*v6*w5 +
     &        3*CA**3*v6*w5 - 4*v7*w5 + 92*CA*v7*w5 +
     &        60*CA2*v7*w5 + 24*CA**3*v7*w5 + 38*v8*w5 +
     &        70*CA*v8*w5 - 56*CA2*v8*w5 - 
     &     126*CA**3*v8*w5 -
     &        36*v9*w5 - 113*CA*v9*w5 + 26*CA2*v9*w5 +
     &        111*CA**3*v9*w5 + 52*CA*v10*w5 - 
     &     52*CA**3*v10*w5 -
     &        12*v6*w6 - 30*CA*v6*w6 + 28*CA2*v6*w6 +
     &        38*CA**3*v6*w6 + 24*v7*w6 + 60*CA*v7*w6 -
     &        56*CA2*v7*w6 - 76*CA**3*v7*w6 - 30*v8*w6 -
     &        103*CA*v8*w6 + 30*CA2*v8*w6 + 
     &     97*CA**3*v8*w6 +
     &        18*v9*w6 + 73*CA*v9*w6 - 2*CA2*v9*w6 -
     &        59*CA**3*v9*w6 - 36*CA*v10*w6 + 
     &     36*CA**3*v10*w6 +
     &        2*CA*v6*w7 - 2*CA**3*v6*w7 - 4*CA*v7*w7 +
     &        4*CA**3*v7*w7 - 4*v8*w7 - 17*CA*v8*w7 +
     &        18*CA2*v8*w7 + 23*CA**3*v8*w7 + 4*v9*w7 +
     &        19*CA*v9*w7 - 18*CA2*v9*w7 - 
     &     25*CA**3*v9*w7 +
     &        12*CA*v10*w7 - 12*CA**3*v10*w7 + 4*v8*w8 +
     &        16*CA*v8*w8 - 8*CA2*v8*w8 - 
     &     16*CA**3*v8*w8 -
     &        4*v9*w8 - 16*CA*v9*w8 + 8*CA2*v9*w8 +
     &        16*CA**3*v9*w8 - 8*CA*v10*w8 + 
     &     8*CA**3*v10*w8 +
     &        4*CA*v10*w9 - 4*CA**3*v10*w9))/
     &    (CA*(1 - v)**2*v*w*(1 - v*w)**3*(1 - v + v*w)**3) 
      part9 = -(4*CF*lw*(3*CA - 3*CA**3 + 2*CA2*CF - 15*CA*v + 
     &   15*CA**3*v - 10*CA2*CF*v + 2*v2 + 39*CA*v2 - 2*CA2*v2 -
     &   39*CA**3*v2 + 2*CA*CF*v2 + 26*CA2*CF*v2 - 8*v3 -
     &   66*CA*v3 + 8*CA2*v3 + 66*CA**3*v3 - 8*CA*CF*v3 -
     &   44*CA2*CF*v3 + 12*v4 + 72*CA*v4 - 12*CA2*v4 -
     &   72*CA**3*v4 + 12*CA*CF*v4 + 48*CA2*CF*v4 - 8*v5 -
     &   48*CA*v5 + 8*CA2*v5 + 48*CA**3*v5 - 8*CA*CF*v5 -
     &   32*CA2*CF*v5 + 2*v6 + 18*CA*v6 - 2*CA2*v6 - 18*CA**3*v6 +
     &   2*CA*CF*v6 + 12*CA2*CF*v6 - 3*CA*v7 + 3*CA**3*v7 -
     &   2*CA2*CF*v7 - 10*CA*w - 2*CA**3*w + 8*CA2*CF*w + 50*CA*v*w +
     &   10*CA**3*v*w - 40*CA2*CF*v*w - v2*w - 105*CA*v2*w -
     &   3*CA2*v2*w - 28*CA**3*v2*w + 90*CA2*CF*v2*w + 4*v3*w +
     &   120*CA*v3*w + 12*CA2*v3*w + 52*CA**3*v3*w - 120*CA2*CF*v3*w -
     &   55*CA*v4*w - 24*CA2*v4*w - 85*CA**3*v4*w +
     &   6*CA*CF*v4*w + 118*CA2*CF*v4*w - 14*v5*w - 45*CA*v5*w +
     &   30*CA2*v5*w + 115*CA**3*v5*w - 18*CA*CF*v5*w -
     &   102*CA2*CF*v5*w + 17*v6*w + 79*CA*v6*w - 21*CA2*v6*w -
     &   98*CA**3*v6*w + 18*CA*CF*v6*w + 70*CA2*CF*v6*w - 6*v7*w -
     &   43*CA*v7*w + 6*CA2*v7*w + 45*CA**3*v7*w - 6*CA*CF*v7*w -
     &   30*CA2*CF*v7*w + 9*CA*v8*w - 9*CA**3*v8*w + 6*CA2*CF*v8*w +
     &   CA*w2 - CA**3*w2 + 2*CA2*CF*w2 - 5*CA*v*w2 + 5*CA**3*v*w2 -
     &   10*CA2*CF*v*w2 - v2*w2 - 20*CA*v2*w2 + CA2*v2*w2 - 
     &   9*CA**3*v2*w2 + 2*CA*CF*v2*w2 + 36*CA2*CF*v2*w2 + 4*v3*w2 +
     &   110*CA*v3*w2 - 4*CA2*v3*w2 + 6*CA**3*v3*w2 - 8*CA*CF*v3*w2 -
     &   84*CA2*CF*v3*w2 - 16*v4*w2 - 229*CA*v4*w2 + 4*CA2*v4*w2 +
     &   19*CA**3*v4*w2 + 6*CA*CF*v4*w2 + 106*CA2*CF*v4*w2 + 34*v5*w2 +
     &   281*CA*v5*w2 + 2*CA2*v5*w2 - 57*CA**3*v5*w2 + 10*CA*CF*v5*w2 -
     &   66*CA2*CF*v5*w2 - 25*v6*w2 - 180*CA*v6*w2 - 11*CA2*v6*w2 +
     &   37*CA**3*v6*w2 - 10*CA*CF*v6*w2 + 28*CA2*CF*v6*w2 - 2*v7*w2 +
     &   30*CA*v7*w2 + 14*CA2*v7*w2 + 18*CA**3*v7*w2 - 6*CA*CF*v7*w2 -
     &   24*CA2*CF*v7*w2 + 6*v8*w2 + 21*CA*v8*w2 - 6*CA2*v8*w2 -
     &   27*CA**3*v8*w2 + 6*CA*CF*v8*w2 + 18*CA2*CF*v8*w2 -
     &   9*CA*v9*w2 + 9*CA**3*v9*w2 - 6*CA2*CF*v9*w2 - 2*v2*w3 +
     &   3*CA*v2*w3 - 10*CA2*v2*w3 - 5*CA**3*v2*w3 + 8*CA*CF*v2*w3 +
     &   10*CA2*CF*v2*w3 + 8*v3*w3 - 12*CA*v3*w3 + 40*CA2*v3*w3 +
     &   20*CA**3*v3*w3 - 32*CA*CF*v3*w3 - 40*CA2*CF*v3*w3 - 4*v4*w3 +
     &   37*CA*v4*w3 - 55*CA2*v4*w3 - 41*CA**3*v4*w3 + 62*CA*CF*v4*w3 +
     &   86*CA2*CF*v4*w3 - 16*v5*w3 - 69*CA*v5*w3 + 25*CA2*v5*w3 +
     &   53*CA**3*v5*w3 - 74*CA*CF*v5*w3 - 118*CA2*CF*v5*w3 + 
     &   6*v6*w3 - 9*CA*v6*w3 + 9*CA2*v6*w3 + 4*CA**3*v6*w3 + 
     &   40*CA*CF*v6*w3 + 74*CA2*CF*v6*w3 + 24*v7*w3 + 119*CA*v7*w3 -
     &   13*CA2*v7*w3 - 73*CA**3*v7*w3 + 6*CA*CF*v7*w3 +
     &   2*CA2*CF*v7*w3 - 14*v8*w3 - 84*CA*v8*w3 + 2*CA2*v8*w3 +
     &   51*CA**3*v8*w3 - 8*CA*CF*v8*w3 - 20*CA2*CF*v8*w3 - 2*v9*w3 +
     &   15*CA*v9*w3 + 2*CA2*v9*w3 - 9*CA**3*v9*w3 - 2*CA*CF*v9*w3 +
     &   6*CA2*CF*v9*w3 + 3*CA*v10*w3 - 3*CA**3*v10*w3 + 
     &   2*CA2*CF*v10*w3 - 3*CA*v2*w4 + 3*CA**3*v2*w4 -
     &   6*CA2*CF*v2*w4 + 12*CA*v3*w4 - 12*CA**3*v3*w4 +
     &   24*CA2*CF*v3*w4 - 14*v4*w4 - 62*CA*v4*w4 - 17*CA2*v4*w4 +
     &   32*CA**3*v4*w4 - 2*CA*CF*v4*w4 - 52*CA2*CF*v4*w4 + 42*v5*w4 +
     &   144*CA*v5*w4 + 51*CA2*v5*w4 - 54*CA**3*v5*w4 + 
     &   6*CA*CF*v5*w4 + 72*CA2*CF*v5*w4 - 21*v6*w4 - 55*CA*v6*w4 -
     &   46*CA2*v6*w4 + 12*CA**3*v6*w4 + 10*CA*CF*v6*w4 -
     &   38*CA2*CF*v6*w4 - 28*v7*w4 - 116*CA*v7*w4 + 7*CA2*v7*w4 +
     &   52*CA**3*v7*w4 - 30*CA*CF*v7*w4 - 16*CA2*CF*v7*w4 + 14*v8*w4 +
     &   75*CA*v8*w4 + 8*CA2*v8*w4 - 26*CA**3*v8*w4 + 12*CA*CF*v8*w4 +
     &   18*CA2*CF*v8*w4 + 7*v9*w4 + 5*CA*v9*w4 - 3*CA2*v9*w4 -
     &   7*CA**3*v9*w4 + 4*CA*CF*v9*w4 - 2*CA2*CF*v9*w4 -
     &   11*CA*v10*w4 + 9*CA**3*v10*w4 - 6*CA2*CF*v10*w4 + 10*v4*w5 +
     &   38*CA*v4*w5 + 20*CA2*v4*w5 + 6*CA**3*v4*w5 - 12*CA*CF*v4*w5 -
     &   24*CA2*CF*v4*w5 - 30*v5*w5 - 114*CA*v5*w5 - 60*CA2*v5*w5 -
     &   18*CA**3*v5*w5 + 36*CA*CF*v5*w5 + 72*CA2*CF*v5*w5 + 2*v6*w5 +
     &   25*CA*v6*w5 + 48*CA2*v6*w5 + 48*CA**3*v6*w5 - 54*CA*CF*v6*w5 -
     &   86*CA2*CF*v6*w5 + 46*v7*w5 + 140*CA*v7*w5 + 4*CA2*v7*w5 -
     &   66*CA**3*v7*w5 + 48*CA*CF*v7*w5 + 52*CA2*CF*v7*w5 - 17*v8*w5 -
     &   62*CA*v8*w5 - 10*CA2*v8*w5 + 17*CA**3*v8*w5 - 16*CA*CF*v8*w5 -
     &   16*CA2*CF*v8*w5 - 11*v9*w5 - 27*CA*v9*w5 - 2*CA2*v9*w5 +
     &   13*CA**3*v9*w5 - 2*CA*CF*v9*w5 + 2*CA2*CF*v9*w5 +
     &   19*CA*v10*w5 - 13*CA**3*v10*w5 + 10*CA2*CF*v10*w5 +
     &   3*CA*v4*w6 - 3*CA**3*v4*w6 + 6*CA2*CF*v4*w6 - 9*CA*v5*w6 + 
     &   9*CA**3*v5*w6 - 18*CA2*CF*v5*w6 + 27*v6*w6 + 90*CA*v6*w6 +
     &   23*CA2*v6*w6 - 11*CA**3*v6*w6 - 2*CA*CF*v6*w6 -
     &   16*CA2*CF*v6*w6 - 54*v7*w6 - 165*CA*v7*w6 - 46*CA2*v7*w6 +
     &   7*CA**3*v7*w6 + 4*CA*CF*v7*w6 + 62*CA2*CF*v7*w6 + 14*v8*w6 +
     &   42*CA*v8*w6 + 15*CA2*v8*w6 + 15*CA**3*v8*w6 -
     &   36*CA2*CF*v8*w6 + 13*v9*w6 + 39*CA*v9*w6 + 8*CA2*v9*w6 -
     &   17*CA**3*v9*w6 - 2*CA*CF*v9*w6 + 2*CA2*CF*v9*w6 -
     &   23*CA*v10*w6 + 15*CA**3*v10*w6 - 14*CA2*CF*v10*w6 - 10*v6*w7 -
     &   23*CA*v6*w7 - 14*CA2*v6*w7 - 11*CA**3*v6*w7 + 8*CA*CF*v6*w7 +
     &   30*CA2*CF*v6*w7 + 20*v7*w7 + 46*CA*v7*w7 + 28*CA2*v7*w7 +
     &   22*CA**3*v7*w7 - 16*CA*CF*v7*w7 - 60*CA2*CF*v7*w7 + 3*v8*w7 +
     &   15*CA*v8*w7 - 2*CA2*v8*w7 - 21*CA**3*v8*w7 +
     &   18*CA2*CF*v8*w7 - 13*v9*w7 - 38*CA*v9*w7 - 12*CA2*v9*w7 +
     &   10*CA**3*v9*w7 + 8*CA*CF*v9*w7 + 12*CA2*CF*v9*w7 +
     &   21*CA*v10*w7 - 13*CA**3*v10*w7 + 14*CA2*CF*v10*w7 - 
     &   CA*v6*w8 + CA**3*v6*w8 - 2*CA2*CF*v6*w8 + 2*CA*v7*w8 -
     &   2*CA**3*v7*w8 + 4*CA2*CF*v7*w8 - 8*v8*w8 - 16*CA*v8*w8 -
     &   11*CA2*v8*w8 - 4*CA**3*v8*w8 + 10*CA*CF*v8*w8 +
     &   20*CA2*CF*v8*w8 + 8*v9*w8 + 15*CA*v9*w8 + 11*CA2*v9*w8 + 
     &   5*CA**3*v9*w8 - 10*CA*CF*v9*w8 - 22*CA2*CF*v9*w8 - 
     &   13*CA*v10*w8 + 7*CA**3*v10*w8 - 10*CA2*CF*v10*w8 + 
     &   2*v8*w9 + 4*CA2*v8*w9 + 4*CA**3*v8*w9 - 4*CA*CF*v8*w9 -
     &   8*CA2*CF*v8*w9 - 2*v9*w9 - 4*CA2*v9*w9 - 4*CA**3*v9*w9 +
     &   4*CA*CF*v9*w9 + 8*CA2*CF*v9*w9 + 5*CA*v10*w9 -
     &   3*CA**3*v10*w9 + 6*CA2*CF*v10*w9 - CA*v10*w10 +
     &   CA**3*v10*w10 - 2*CA2*CF*v10*w10))/
     &   (CA*(1 - v)**2*v*(1 - w)*w*(1 - v*w)**3*(1 - v + v*w)**3)
      struv11= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV12(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (8*CA*CF**2*l1w*Nf*(1 - 2*v + v2 + v2*w2)*
     &      (1 + v2 - 2*v2*w + v2*w2))/(1 - v + v*w)**4 
      part2 = -(8*CA*CF**2*lmss*Nf*(1 - 2*v + v2 + v2*w2)*
     &      (1 + v2 - 2*v2*w + v2*w2))/(1 - v + v*w)**4 
      part3 = (8*CA*CF**2*lv*Nf*(1 - 2*v + v2 + v2*w2)*
     &      (1 + v2 - 2*v2*w + v2*w2))/(1 - v + v*w)**4 
      part4 = (16*CA*CF**2*l1vw*Nf*(1 - w)*(1 - 2*v + v2 + v2*w2)*
     &      (1 + v2 - 2*v2*w + v2*w2))/(1 - v + v*w)**4 
      part5 = (4*CF*lvw*(4*CA**3 + 4*CA**3*v2 + 4*w - 2*CA2*w - 
     &     4*CA**3*w +
     &        2*CA4*w - 6*v*w - 2*CA*v*w + 5*CA2*v*w - 
     &     8*CA**3*v*w -
     &        4*CA4*v*w + 4*v2*w - 4*CA2*v2*w - 
     &     10*CA**3*v2*w +
     &        3*CA4*v2*w - 2*v3*w + 2*CA*v3*w + CA2*v3*w -
     &        2*CA**3*v3*w - CA4*v3*w - 4*v*w2 + 2*CA*v*w2 +
     &        5*CA2*v*w2 + 8*CA**3*v*w2 - 8*CA4*v*w2 +
     &        3*v2*w2 - CA*v2*w2 - 9*CA2*v2*w2 +
     &        17*CA**3*v2*w2 + 15*CA4*v2*w2 - 
     &     4*CA*v3*w2 +
     &        2*CA2*v3*w2 + 6*CA**3*v3*w2 - 
     &     11*CA4*v3*w2 +
     &        v4*w2 - CA*v4*w2 + 2*CA2*v4*w2 +
     &        CA**3*v4*w2 + 4*CA4*v4*w2 + 4*v2*w3 +
     &        CA*v2*w3 - 12*CA2*v2*w3 - 11*CA**3*v2*w3 +
     &        13*CA4*v2*w3 - 2*v3*w3 + 4*CA*v3*w3 +
     &        23*CA2*v3*w3 - 10*CA**3*v3*w3 - 
     &     20*CA4*v3*w3 -
     &        2*v4*w3 + 3*CA*v4*w3 - 15*CA2*v4*w3 -
     &        3*CA**3*v4*w3 + 11*CA4*v4*w3 + 
     &     4*CA2*v5*w3 -
     &        4*CA4*v5*w3 - v3*w4 - 2*CA*v3*w4 +
     &        6*CA2*v3*w4 + 6*CA**3*v3*w4 - 
     &     9*CA4*v3*w4 +
     &        v4*w4 - 4*CA*v4*w4 - 6*CA2*v4*w4 +
     &        4*CA**3*v4*w4 + 9*CA4*v4*w4 + 2*CA*v4*w5 -
     &        4*CA2*v4*w5 - 2*CA**3*v4*w5 + 
     &     4*CA4*v4*w5 +
     &        4*CA2*v5*w5 - 4*CA4*v5*w5))/
     &    (CA*(1 - v)**2*v2*w2) 
      part6 = -(4*CF*Nf*(-2*CA*CF - 4*v + 4*CA2*v + 7*v2 - 7*CA2*v2 -
     &        8*v3 + 8*CA2*v3 + 7*v4 - 7*CA2*v4 - 4*v5 +
     &        4*CA2*v5 + v6 - CA2*v6 + 4*v*w - 4*CA2*v*w -
     &        13*v2*w + 24*CA2*v2*w + 14*v3*w - 
     &     58*CA2*v3*w -
     &        4*v4*w + 84*CA2*v4*w - 2*v5*w - 88*CA2*v5*w +
     &        v6*w + 64*CA2*v6*w - 26*CA2*v7*w + 
     &     4*CA2*v8*w +
     &        6*v2*w2 - 6*CA2*v2*w2 - 10*v3*w2 +
     &        23*CA2*v3*w2 + 30*v4*w2 - 81*CA2*v4*w2 -
     &        54*v5*w2 + 153*CA2*v5*w2 + 28*v6*w2 -
     &        137*CA2*v6*w2 + 60*CA2*v7*w2 - 
     &     12*CA2*v8*w2 +
     &        4*v3*w3 - 4*CA2*v3*w3 - 34*v4*w3 +
     &        41*CA2*v4*w3 + 58*v5*w3 - 82*CA2*v5*w3 -
     &        28*v6*w3 + 77*CA2*v6*w3 - 44*CA2*v7*w3 +
     &        12*CA2*v8*w3 + v4*w4 - CA2*v4*w4 +
     &        2*v5*w4 - 3*CA2*v5*w4 - v6*w4 -
     &        6*CA2*v6*w4 + 12*CA2*v7*w4 - 
     &     4*CA2*v8*w4 -
     &        v6*w5 + 3*CA2*v6*w5 - 2*CA2*v7*w5))/
     &    (3.*(1 - v)*v2*w*(1 - v + v*w)**4) 
      part7 = -(4*CF*l1v*(2*CA**3 - 4*CA**3*v + 4*CA**3*v2 -
     &   4*CA**3*v3 +
     &        2*CA**3*v4 + w - CA2*w + 2*CA4*w -
     &      2*v*w + 2*CA*v*w +
     &        5*CA2*v*w - 4*CA**3*v*w - 7*CA4*v*w + 2*v2*w -
     &        2*CA*v2*w - 8*CA2*v2*w + 4*CA**3*v2*w +
     &        11*CA4*v2*w - 4*v3*w - 2*CA*v3*w + 
     &     6*CA2*v3*w +
     &        2*CA**3*v3*w - 10*CA4*v3*w + 5*v4*w + 
     &     2*CA*v4*w -
     &        3*CA2*v4*w + 5*CA4*v4*w - 2*v5*w + 
     &     CA2*v5*w -
     &        2*CA**3*v5*w - CA4*v5*w - v*w2 - CA2*v*w2 +
     &        4*v2*w2 - 6*CA*v2*w2 + 7*CA2*v2*w2 +
     &        4*CA**3*v2*w2 - 7*CA4*v2*w2 - 2*v3*w2 +
     &        6*CA*v3*w2 - 22*CA2*v3*w2 - 6*CA**3*v3*w2 +
     &        19*CA4*v3*w2 + 2*v4*w2 + 2*CA*v4*w2 +
     &        34*CA2*v4*w2 - 6*CA**3*v4*w2 - 
     &     17*CA4*v4*w2 -
     &        5*v5*w2 - 2*CA*v5*w2 - 23*CA2*v5*w2 +
     &        4*CA**3*v5*w2 + 4*CA4*v5*w2 + 2*v6*w2 +
     &        5*CA2*v6*w2 + CA4*v6*w2 - 3*v2*w3 -
     &        CA2*v2*w3 - 2*CA4*v2*w3 + 2*v3*w3 +
     &        6*CA*v3*w3 + 2*CA2*v3*w3 - 2*CA**3*v3*w3 +
     &        16*CA4*v3*w3 - 6*v4*w3 - 6*CA*v4*w3 +
     &        2*CA2*v4*w3 + 10*CA**3*v4*w3 - 
     &     33*CA4*v4*w3 +
     &        8*v5*w3 - 14*CA2*v5*w3 + 28*CA4*v5*w3 -
     &        v6*w3 + 15*CA2*v6*w3 - 9*CA4*v6*w3 -
     &        4*CA2*v7*w3 - v3*w4 -
     &      CA2*v3*w4 + 8*v4*w4 -
     &        2*CA*v4*w4 + CA2*v4*w4 - 2*CA**3*v4*w4 -
     &        CA4*v4*w4 - 6*v5*w4 + 2*CA*v5*w4 +
     &        5*CA2*v5*w4 - 4*CA**3*v5*w4 + 
     &     2*CA4*v5*w4 -
     &        v6*w4 - 5*CA2*v6*w4 -
     &      CA4*v6*w4 - 2*v5*w5 +
     &        3*CA2*v5*w5 + 2*CA**3*v5*w5 - CA4*v5*w5 +
     &        2*v6*w5 - 7*CA2*v6*w5 + CA4*v6*w5 +
     &        4*CA2*v7*w5))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)*(1 - v + v*w)) 
      part8 = -(4*CF*lw*(2*CA - 4*CA**3 - 6*CA*v + 6*CA**3*v + 
     &   8*CA*v2 -
     &        8*CA**3*v2 - 8*CA*v3 + 8*CA**3*v3 + 6*CA*v4 -
     &        4*CA**3*v4 - 2*CA*v5 + 2*CA**3*v5 - 2*w - 2*CA*w +
     &        5*CA2*w + 4*CA**3*w - 3*CA4*w + 7*v*w + 6*CA*v*w -
     &        20*CA2*v*w + 13*CA4*v*w - 10*v2*w - 8*CA*v2*w +
     &        37*CA2*v2*w + 4*CA**3*v2*w - 27*CA4*v2*w + 
     &     8*v3*w +
     &        10*CA*v3*w - 43*CA2*v3*w - 14*CA**3*v3*w +
     &        35*CA4*v3*w - 4*v4*w - 8*CA*v4*w +
     &      34*CA2*v4*w +
     &        2*CA**3*v4*w - 30*CA4*v4*w + v5*w - 
     &     17*CA2*v5*w -
     &        2*CA**3*v5*w + 16*CA4*v5*w + 2*CA*v6*w +
     &        4*CA2*v6*w - 2*CA**3*v6*w - 4*CA4*v6*w + w2 -
     &        CA2*w2 - 6*v*w2 + 2*CA2*v*w2 - 6*CA**3*v*w2 +
     &        5*CA4*v*w2 + 8*v2*w2 + 3*CA*v2*w2 -
     &        3*CA2*v2*w2 - CA**3*v2*w2 - 
     &     14*CA4*v2*w2 -
     &        5*v3*w2 - CA*v3*w2 + 8*CA2*v3*w2 +
     &        17*CA**3*v3*w2 + 19*CA4*v3*w2 - 
     &     5*CA*v4*w2 -
     &        7*CA2*v4*w2 + 15*CA**3*v4*w2 - 
     &     20*CA4*v4*w2 +
     &        3*v5*w2 + 9*CA*v5*w2 - 4*CA2*v5*w2 -
     &        3*CA**3*v5*w2 + 18*CA4*v5*w2 - v6*w2 -
     &        6*CA*v6*w2 + 9*CA2*v6*w2 + 6*CA**3*v6*w2 -
     &        12*CA4*v6*w2 - 4*CA2*v7*w2 + 
     &     4*CA4*v7*w2 +
     &        v*w3 + CA2*v*w3 - v2*w3 - 3*CA*v2*w3 -
     &        5*CA2*v2*w3 + 5*CA**3*v2*w3 - 
     &     7*CA4*v2*w3 +
     &        4*v3*w3 - 9*CA*v3*w3 + 7*CA2*v3*w3 -
     &        11*CA**3*v3*w3 + 14*CA4*v3*w3 + v4*w3 +
     &        14*CA*v4*w3 - 15*CA2*v4*w3 - 
     &     30*CA**3*v4*w3 -
     &        7*CA4*v4*w3 - 5*v5*w3 - 9*CA*v5*w3 +
     &        27*CA2*v5*w3 + 3*CA**3*v5*w3 - 
     &     4*CA4*v5*w3 +
     &        7*CA*v6*w3 - 23*CA2*v6*w3 - 7*CA**3*v6*w3 +
     &        8*CA4*v6*w3 + 8*CA2*v7*w3 -
     &      4*CA4*v7*w3 +
     &        v2*w4 + CA2*v2*w4 - 6*v3*w4 + 
     &     8*CA*v3*w4 +
     &        CA2*v3*w4 + 6*CA4*v3*w4 - v4*w4 -
     &        4*CA*v4*w4 + 6*CA2*v4*w4 + 22*CA**3*v4*w4 -
     &        12*CA4*v4*w4 + 3*v5*w4 - 21*CA2*v5*w4 +
     &        4*CA**3*v5*w4 + 9*CA4*v5*w4 + 3*v6*w4 -
     &        4*CA*v6*w4 + 21*CA2*v6*w4 + 4*CA**3*v6*w4 -
     &        7*CA4*v6*w4 - 8*CA2*v7*w4 + 
     &     4*CA4*v7*w4 +
     &        v3*w5 + CA2*v3*w5 + 3*v4*w5 -
     &      3*CA*v4*w5 -
     &        10*CA2*v4*w5 - 5*CA**3*v4*w5 + 2*CA*v5*w5 +
     &        21*CA2*v5*w5 - 6*CA**3*v5*w5 - 
     &     2*CA4*v5*w5 -
     &        4*v6*w5 + CA*v6*w5 - 20*CA2*v6*w5 -
     &        CA**3*v6*w5 + 6*CA4*v6*w5 + 8*CA2*v7*w5 -
     &        4*CA4*v7*w5 + 2*CA2*v4*w6 - 2*v5*w6 -
     &        7*CA2*v5*w6 + 2*CA**3*v5*w6 - CA4*v5*w6 +
     &        2*v6*w6 + 9*CA2*v6*w6 + CA4*v6*w6 -
     &        4*CA2*v7*w6))/
     &    (CA*(1 - v)**2*v2*(1 - w)*w2*(1 - v*w)*(1 - v + v*w)) 
      part9 = (2*CF*lms*(2 - 6*CA2 - 4*v + 4*v2 - 
     &     8*CA2*v2 - 4*v3 +
     &        2*v4 - 2*CA2*v4 - 2*w + 2*CA2*w + 6*CF*w -
     &        6*CA2*CF*w + 20*CA2*v*w - 
     &     14*CF*v*w + 18*CA2*CF*v*w +
     &        4*v2*w + 12*CA2*v2*w + 14*CF*v2*w -
     &        30*CA2*CF*v2*w - 4*v3*w + 
     &     16*CA2*v3*w - 6*CF*v3*w +
     &        34*CA2*CF*v3*w + 6*v4*w + 2*CA2*v4*w -
     &        24*CA2*CF*v4*w - 4*v5*w + 4*CA2*v5*w +
     &        8*CA2*CF*v5*w + w2 - CA2*w2 - 2*CF*w2 +
     &        2*CA2*CF*w2 + 2*v*w2 - 4*CA2*v*w2 -
     &      6*CF*v*w2 +
     &        2*CA2*CF*v*w2 - 2*v2*w2 - 36*CA2*v2*w2 +
     &        14*CF*v2*w2 + 10*CA2*CF*v2*w2 + 2*v3*w2 -
     &        24*CA2*v3*w2 - 18*CF*v3*w2 - 
     &     22*CA2*CF*v3*w2 -
     &        v4*w2 - 13*CA2*v4*w2 + 12*CF*v4*w2 -
     &        4*CA2*v5*w2 + 24*CA2*CF*v5*w2 + 2*v6*w2 -
     &        2*CA2*v6*w2 - 16*CA2*CF*v6*w2 - 2*v*w3 +
     &        2*CA2*v*w3 + 4*CF*v*w3 - 4*CA2*CF*v*w3 +
     &        2*v2*w3 + 2*CA2*v2*w3 - 6*CF*v2*w3 +
     &        14*CA2*CF*v2*w3 - 6*v3*w3 + 
     &     40*CA2*v3*w3 +
     &        10*CF*v3*w3 - 90*CA2*CF*v3*w3 - 4*v4*w3 +
     &        20*CA2*v4*w3 + 160*CA2*CF*v4*w3 +
     &        6*CA2*v5*w3 - 8*CF*v5*w3 - 
     &     112*CA2*CF*v5*w3 -
     &        2*v6*w3 + 2*CA2*v6*w3 + 24*CA2*CF*v6*w3 +
     &        8*CA2*CF*v7*w3 + v2*w4 - CA2*v2*w4 -
     &        2*CF*v2*w4 + 2*CA2*CF*v2*w4 - 2*v3*w4 +
     &        6*CF*v3*w4 - 10*CA2*CF*v3*w4 + 11*v4*w4 -
     &        27*CA2*v4*w4 - 12*CF*v4*w4 + 
     &     88*CA2*CF*v4*w4 +
     &        2*v5*w4 - 8*CA2*v5*w4 + 8*CF*v5*w4 -
     &        144*CA2*CF*v5*w4 + 2*v6*w4 - 
     &     2*CA2*v6*w4 +
     &        88*CA2*CF*v6*w4 - 24*CA2*CF*v7*w4 - 
     &     6*v5*w5 +
     &        10*CA2*v5*w5 - 32*CA2*CF*v5*w5 - 
     &     2*v6*w5 +
     &        2*CA2*v6*w5 + 40*CA2*CF*v6*w5 -
     &        8*CA2*CF*v7*w5 + 2*v6*w6 - 2*CA2*v6*w6 +
     &        8*CA2*CF*v6*w6 - 8*CA2*CF*v7*w6))/
     &    ((1 - v)**2*v2*w2*(1 - v*w)**2) 
      part10 = -(4*CF*l1vw*(3 - 2*CA2 + 3*CA4 - 8*CA**3*CF - 11*v -
     &      6*CA*v +
     &        8*CA2*v - 12*CA4*v + 32*CA**3*CF*v + 15*v2 + 
     &     18*CA*v2 -
     &        12*CA2*v2 + 18*CA4*v2 - 48*CA**3*CF*v2 - 
     &     7*v3 -
     &        12*CA*v3 + 10*CA2*v3 - 13*CA4*v3 + 
     &     32*CA**3*CF*v3 -
     &        7*v4 - 12*CA*v4 - 10*CA2*v4 + 7*CA4*v4 -
     &        8*CA**3*CF*v4 + 15*v5 + 18*CA*v5 + 12*CA2*v5 -
     &        6*CA4*v5 - 11*v6 - 6*CA*v6 - 8*CA2*v6 +
     &        4*CA4*v6 + 3*v7 + 2*CA2*v7 - 
     &     CA4*v7 + 9*v*w +
     &        6*CA*v*w - 3*CA2*v*w + 6*CA4*v*w - 16*CA**3*CF*v*w -
     &        26*v2*w - 36*CA*v2*w + 16*CA2*v2*w - 
     &     8*CA4*v2*w +
     &        48*CA**3*CF*v2*w + 19*v3*w + 40*CA*v3*w -
     &        37*CA2*v3*w - 4*CA**3*v3*w - 19*CA4*v3*w -
     &        48*CA**3*CF*v3*w + 24*v4*w + 36*CA*v4*w +
     &        60*CA2*v4*w + 12*CA**3*v4*w + 52*CA4*v4*w +
     &        16*CA**3*CF*v4*w - 61*v5*w - 78*CA*v5*w -
     &        73*CA2*v5*w - 12*CA**3*v5*w - 52*CA4*v5*w +
     &        50*v6*w + 32*CA*v6*w + 52*CA2*v6*w + 
     &     4*CA**3*v6*w +
     &        28*CA4*v6*w - 15*v7*w - 15*CA2*v7*w -
     &        7*CA4*v7*w + 13*v2*w2 + 18*CA*v2*w2 -
     &        9*CA2*v2*w2 - 4*CA4*v2*w2 - 
     &     16*CA**3*CF*v2*w2 -
     &        19*v3*w2 - 44*CA*v3*w2 + 47*CA2*v3*w2 +
     &        8*CA**3*v3*w2 + 23*CA4*v3*w2 +
     &        32*CA**3*CF*v3*w2 - 30*v4*w2 - 40*CA*v4*w2 -
     &        110*CA2*v4*w2 - 32*CA**3*v4*w2 - 
     &     50*CA4*v4*w2 -
     &        16*CA**3*CF*v4*w2 + 98*v5*w2 + 140*CA*v5*w2 +
     &        162*CA2*v5*w2 + 40*CA**3*v5*w2 + 
     &     64*CA4*v5*w2 -
     &        95*v6*w2 - 74*CA*v6*w2 - 137*CA2*v6*w2 -
     &        16*CA**3*v6*w2 - 50*CA4*v6*w2 + 33*v7*w2 +
     &        47*CA2*v7*w2 + 17*CA4*v7*w2 + 7*v3*w3 +
     &        16*CA*v3*w3 - 20*CA2*v3*w3 -
     &      4*CA**3*v3*w3 +
     &        9*CA4*v3*w3 - 16*CA**3*CF*v3*w3 + 
     &     16*v4*w3 +
     &        20*CA*v4*w3 + 80*CA2*v4*w3 + 
     &     28*CA**3*v4*w3 -
     &        20*CA4*v4*w3 + 16*CA**3*CF*v4*w3 - 
     &     78*v5*w3 -
     &        132*CA*v5*w3 - 168*CA2*v5*w3 - 
     &     48*CA**3*v5*w3 +
     &        20*CA4*v5*w3 + 100*v6*w3 + 96*CA*v6*w3 +
     &        188*CA2*v6*w3 + 24*CA**3*v6*w3 - 45*v7*w3 -
     &        80*CA2*v7*w3 - 9*CA4*v7*w3 - 3*v4*w4 -
     &        4*CA*v4*w4 - 20*CA2*v4*w4 - 8*CA**3*v4*w4 +
     &        11*CA4*v4*w4 - 8*CA**3*CF*v4*w4 +
     &      31*v5*w4 +
     &        66*CA*v5*w4 + 82*CA2*v5*w4 + 
     &     24*CA**3*v5*w4 -
     &        26*CA4*v5*w4 - 65*v6*w4 - 74*CA*v6*w4 -
     &        142*CA2*v6*w4 - 16*CA**3*v6*w4 + 
     &     8*CA4*v6*w4 +
     &        45*v7*w4 + 80*CA2*v7*w4 + 9*CA4*v7*w4 -
     &        5*v5*w5 - 14*CA*v5*w5 - 15*CA2*v5*w5 -
     &        4*CA**3*v5*w5 + 26*v6*w5 + 32*CA*v6*w5 +
     &        56*CA2*v6*w5 + 4*CA**3*v6*w5 + 
     &     20*CA4*v6*w5 -
     &        33*v7*w5 - 47*CA2*v7*w5 - 17*CA4*v7*w5 -
     &        5*v6*w6 - 6*CA*v6*w6 - 9*CA2*v6*w6 -
     &        10*CA4*v6*w6 + 15*v7*w6 + 15*CA2*v7*w6 +
     &        7*CA4*v7*w6 - 3*v7*w7 - 2*CA2*v7*w7 +
     &        CA4*v7*w7))/(CA*(1 - v)*v2*w*(1 - v + v*w)**4) 
      part11 = (4*CF*lmss*(CA - CA**3 - 6*CA*v + 6*CA**3*v + 16*CA*v2 -
     &   16*CA**3*v2 - 26*CA*v3 + 26*CA**3*v3 + 30*CA*v4 -
     &   30*CA**3*v4 - 26*CA*v5 + 26*CA**3*v5 + 16*CA*v6 -
     &   16*CA**3*v6 - 6*CA*v7 + 6*CA**3*v7 + CA*v8 - CA**3*v8 -
     &   4*CA**3*CF*w + v*w + 3*CA*v*w -  5*CA2*v*w -
     &   3*CA**3*v*w + 4*CA4*v*w + 20*CA**3*CF*v*w - 5*v2*w -
     &   17*CA*v2*w + 31*CA2*v2*w + 17*CA**3*v2*w - 26*CA4*v2*w -
     &   40*CA**3*CF*v2*w + 11*v3*w + 45*CA*v3*w - 85*CA2*v3*w -
     &   45*CA**3*v3*w + 74*CA4*v3*w + 40*CA**3*CF*v3*w - 15*v4*w -
     &   75*CA*v4*w + 137*CA2*v4*w + 75*CA**3*v4*w - 122*CA4*v4*w -
     &   20*CA**3*CF*v4*w + 15*v5*w + 85*CA*v5*w - 145*CA2*v5*w -
     &   85*CA**3*v5*w + 130*CA4*v5*w + 4*CA**3*CF*v5*w - 11*v6*w -
     &   63*CA*v6*w + 105*CA2*v6*w + 63*CA**3*v6*w - 94*CA4*v6*w +
     &   5*v7*w + 27*CA*v7*w - 51*CA2*v7*w - 27*CA**3*v7*w +
     &   46*CA4*v7*w - v8*w - 5*CA*v8*w + 15*CA2*v8*w + 5*CA**3*v8*w -
     &   14*CA4*v8*w - 2*CA2*v9*w + 2*CA4*v9*w - 16*CA**3*CF*v*w2 +
     &   v2*w2 + 5*CA*v2*w2 - 15*CA2*v2*w2 - 5*CA**3*v2*w2 +
     &   22*CA4*v2*w2 + 80*CA**3*CF*v2*w2 - 6*v3*w2 - 28*CA*v3*w2 +
     &   88*CA2*v3*w2 + 28*CA**3*v3*w2 - 130*CA4*v3*w2 -
     &   168*CA**3*CF*v3*w2 + 17*v4*w2 + 73*CA*v4*w2 - 227*CA2*v4*w2 -
     &   73*CA**3*v4*w2 + 338*CA4*v4*w2 + 192*CA**3*CF*v4*w2 - 
     &   28*v5*w2 - 112*CA*v5*w2 + 338*CA2*v5*w2 + 112*CA**3*v5*w2 -
     &   506*CA4*v5*w2 - 128*CA**3*CF*v5*w2 + 27*v6*w2 +
     &   103*CA*v6*w2 - 317*CA2*v6*w2 - 103*CA**3*v6*w2 +
     &   474*CA4*v6*w2 + 48*CA**3*CF*v6*w2 - 14*v7*w2 - 52*CA*v7*w2 +
     &   188*CA2*v7*w2 + 52*CA**3*v7*w2 - 278*CA4*v7*w2 -
     &   8*CA**3*CF*v7*w2 + 3*v8*w2 + 11*CA*v8*w2 - 65*CA2*v8*w2 - 
     &   11*CA**3*v8*w2 + 94*CA4*v8*w2 + 10*CA2*v9*w2 - 14*CA4*v9*w2 -
     &   24*CA**3*CF*v2*w3 + v3*w3 + 7*CA*v3*w3 - 25*CA2*v3*w3 -
     &   7*CA**3*v3*w3 + 32*CA4*v3*w3 + 112*CA**3*CF*v3*w3 - 7*v4*w3 -
     &   35*CA*v4*w3 + 133*CA2*v4*w3 + 35*CA**3*v4*w3 - 182*CA4*v4*w3 -
     &   216*CA**3*CF*v4*w3 + 18*v5*w3 + 78*CA*v5*w3 - 304*CA2*v5*w3 - 
     &   78*CA**3*v5*w3 + 442*CA4*v5*w3 + 216*CA**3*CF*v5*w3 - 
     &   22*v6*w3 - 94*CA*v6*w3 + 386*CA2*v6*w3 + 94*CA**3*v6*w3 -
     &   588*CA4*v6*w3 - 112*CA**3*CF*v6*w3 + 13*v7*w3 + 59*CA*v7*w3 -
     &   289*CA2*v7*w3 - 59*CA**3*v7*w3 + 452*CA4*v7*w3 +
     &   24*CA**3*CF*v7*w3 - 3*v8*w3 - 15*CA*v8*w3 + 121*CA2*v8*w3 + 
     &   15*CA**3*v8*w3 - 190*CA4*v8*w3 - 22*CA2*v9*w3 + 
     &   34*CA4*v9*w3 - 16*CA**3*CF*v3*w4 + 2*v4*w4 + 8*CA*v4*w4 -
     &   32*CA2*v4*w4 - 8*CA**3*v4*w4 + 30*CA4*v4*w4 +
     &   64*CA**3*CF*v4*w4 - 6*v5*w4 - 32*CA*v5*w4 + 140*CA2*v5*w4 +
     &   32*CA**3*v5*w4 - 138*CA4*v5*w4 - 104*CA**3*CF*v5*w4 + 
     &   8*v6*w4 + 56*CA*v6*w4 - 260*CA2*v6*w4 - 56*CA**3*v6*w4 +
     &   276*CA4*v6*w4 + 80*CA**3*CF*v6*w4 - 6*v7*w4 - 48*CA*v7*w4 +
     &   258*CA2*v7*w4 + 48*CA**3*v7*w4 - 296*CA4*v7*w4 -
     &   24*CA**3*CF*v7*w4 + 2*v8*w4 + 16*CA*v8*w4 - 136*CA2*v8*w4 - 
     &   16*CA**3*v8*w4 + 166*CA4*v8*w4 + 30*CA2*v9*w4 - 
     &   38*CA4*v9*w4 - 4*CA**3*CF*v4*w5 + v5*w5 + 7*CA*v5*w5 -
     &        29*CA2*v5*w5 - 7*CA**3*v5*w5 + 
     &     16*CA4*v5*w5 +
     &        12*CA**3*CF*v5*w5 - 3*v6*w5 - 23*CA*v6*w5 +
     &        105*CA2*v6*w5 + 23*CA**3*v6*w5 -
     &      62*CA4*v6*w5 -
     &        16*CA**3*CF*v6*w5 + 5*v7*w5 + 31*CA*v7*w5 -
     &        155*CA2*v7*w5 - 31*CA**3*v7*w5 + 
     &     98*CA4*v7*w5 +
     &        8*CA**3*CF*v7*w5 - 3*v8*w5 - 15*CA*v8*w5 +
     &        109*CA2*v8*w5 + 15*CA**3*v8*w5 - 
     &     74*CA4*v8*w5 -
     &        30*CA2*v9*w5 + 22*CA4*v9*w5 + v6*w6 +
     &        5*CA*v6*w6 - 19*CA2*v6*w6 - 5*CA**3*v6*w6 +
     &        10*CA4*v6*w6 - 4*v7*w6 - 14*CA*v7*w6 +
     &        58*CA2*v7*w6 + 14*CA**3*v7*w6 -
     &      26*CA4*v7*w6 +
     &        3*v8*w6 + 11*CA*v8*w6 - 61*CA2*v8*w6 -
     &        11*CA**3*v8*w6 + 26*CA4*v8*w6 + 
     &     22*CA2*v9*w6 -
     &        10*CA4*v9*w6 + v7*w7 + 3*CA*v7*w7 -
     &        9*CA2*v7*w7 - 3*CA**3*v7*w7 + 
     &     4*CA4*v7*w7 -
     &        v8*w7 - 5*CA*v8*w7 + 19*CA2*v8*w7 +
     &        5*CA**3*v8*w7 - 10*CA4*v8*w7 - 
     &     10*CA2*v9*w7 +
     &        6*CA4*v9*w7 + CA*v8*w8 - 2*CA2*v8*w8 -
     &        CA**3*v8*w8 + 2*CA4*v8*w8 + 2*CA2*v9*w8 -
     &        2*CA4*v9*w8))/(CA*(1 - v)**2*v2*w2*
     &     (1 - v + v*w)**4)
      part12 = -(2*CF*(6*CA - 18*CA**3 - 48*CA*v + 120*CA**3*v +
     &   168*CA*v2 - 348*CA**3*v2 - 336*CA*v3 + 576*CA**3*v3 + 
     &   420*CA*v4 - 600*CA**3*v4 - 336*CA*v5 + 408*CA**3*v5 + 
     &   168*CA*v6 - 180*CA**3*v6 - 48*CA*v7 + 48*CA**3*v7 + 6*CA*v8 -
     &   6*CA**3*v8 - 9*w + 3*CA*w + 7*CA2*w - 3*CA**3*w + 2*CA4*w +
     &   57*v*w + 3*CA*v*w - 59*CA2*v*w - 27*CA**3*v*w + 2*CA4*v*w -
     &   177*v2*w - 81*CA*v2*w + 239*CA2*v2*w +
     &   177*CA**3*v2*w - 62*CA4*v2*w + 351*v3*w + 249*CA*v3*w -
     &   585*CA2*v3*w - 369*CA**3*v3*w + 234*CA4*v3*w - 465*v4*w -
     &   333*CA*v4*w + 927*CA2*v4*w + 333*CA**3*v4*w - 462*CA4*v4*w +
     &   399*v5*w + 177*CA*v5*w - 977*CA2*v5*w - 57*CA**3*v5*w +
     &   578*CA4*v5*w - 207*v6*w + 57*CA*v6*w + 689*CA2*v6*w -
     &   153*CA**3*v6*w - 482*CA4*v6*w + 57*v7*w - 129*CA*v7*w -
     &   319*CA2*v7*w + 153*CA**3*v7*w + 262*CA4*v7*w - 6*v8*w +
     &   66*CA*v8*w + 90*CA2*v8*w - 66*CA**3*v8*w - 84*CA4*v8*w -
     &   12*CA*v9*w - 12*CA2*v9*w + 12*CA**3*v9*w + 12*CA4*v9*w -
     &   3*CA*w2 + 6*CA2*w2 + 3*CA**3*w2 - 6*CA4*w2 - 27*v*w2 +
     &   15*CA*v*w2 - 37*CA2*v*w2 - 27*CA**3*v*w2 + 52*CA4*v*w2 +
     &   135*v2*w2 + 9*CA*v2*w2 + 62*CA2*v2*w2 + 63*CA**3*v2*w2 -
     &   237*CA4*v2*w2 - 273*v3*w2 - 81*CA*v3*w2 + 80*CA2*v3*w2 - 
     &   135*CA**3*v3*w2 + 639*CA4*v3*w2 + 225*v4*w2 - 39*CA*v4*w2 -
     &   411*CA2*v4*w2 + 399*CA**3*v4*w2 - 1136*CA4*v4*w2 + 87*v5*w2 +
     &   465*CA*v5*w2 + 539*CA2*v5*w2 - 765*CA**3*v5*w2 +
     &   1466*CA4*v5*w2 - 351*v6*w2 - 747*CA*v6*w2 - 224*CA2*v6*w2 +
     &   819*CA**3*v6*w2 - 1421*CA4*v6*w2 + 297*v7*w2 + 549*CA*v7*w2 -
     &   158*CA2*v7*w2 - 501*CA**3*v7*w2 + 1003*CA4*v7*w2 - 105*v8*w2 -
     &   186*CA*v8*w2 + 239*CA2*v8*w2 + 162*CA**3*v8*w2 -
     &   488*CA4*v8*w2 + 12*v9*w2 + 12*CA*v9*w2 - 120*CA2*v9*w2 -
     &   12*CA**3*v9*w2 + 152*CA4*v9*w2 + 6*CA*v10*w2 + 
     &   24*CA2*v10*w2 - 6*CA**3*v10*w2 - 24*CA4*v10*w2 - 6*CA*v*w3 +
     &   12*CA2*v*w3 + 6*CA**3*v*w3 - 12*CA4*v*w3 - 9*v2*w3 +
     &   3*CA*v2*w3 - 85*CA2*v2*w3 - 27*CA**3*v2*w3 + 70*CA4*v2*w3 -
     &   63*v3*w3 + 21*CA*v3*w3 + 238*CA2*v3*w3 + 51*CA**3*v3*w3 -
     &   61*CA4*v3*w3 + 450*v4*w3 + 174*CA*v4*w3 - 435*CA2*v4*w3 -
     &   198*CA**3*v4*w3 - 179*CA4*v4*w3 - 1044*v5*w3 - 720*CA*v5*w3 +
     &   804*CA2*v5*w3 + 540*CA**3*v5*w3 + 472*CA4*v5*w3 + 1185*v6*w3 +
     &   1005*CA*v6*w3 - 1341*CA2*v6*w3 - 693*CA**3*v6*w3 -
     &   636*CA4*v6*w3 - 651*v7*w3 - 615*CA*v7*w3 + 1494*CA2*v7*w3 +
     &   423*CA**3*v7*w3 + 663*CA4*v7*w3 + 108*v8*w3 + 96*CA*v8*w3 -
     &   989*CA2*v8*w3 - 72*CA**3*v8*w3 - 467*CA4*v8*w3 + 30*v9*w3 +
     &   72*CA*v9*w3 + 344*CA2*v9*w3 - 60*CA**3*v9*w3 +
     &   190*CA4*v9*w3 - 6*v10*w3 - 30*CA*v10*w3 - 30*CA2*v10*w3 +
     &   30*CA**3*v10*w3 - 52*CA4*v10*w3 - 12*CA2*v11*w3 +
     &   12*CA4*v11*w3 + 3*CA*v2*w4 - 6*CA2*v2*w4 - 3*CA**3*v2*w4 +
     &   6*CA4*v2*w4 + 45*v3*w4 - 33*CA*v3*w4 + 47*CA2*v3*w4 + 
     &   45*CA**3*v3*w4 - 80*CA4*v3*w4 - 303*v4*w4 - 42*CA*v4*w4 -
     &   83*CA2*v4*w4 - 102*CA**3*v4*w4 + 240*CA4*v4*w4 + 702*v5*w4 +
     &   396*CA*v5*w4 - 198*CA2*v5*w4 - 48*CA**3*v5*w4 - 
     &   588*CA4*v5*w4 - 660*v6*w4 - 597*CA*v6*w4 + 922*CA2*v6*w4 +
     &   285*CA**3*v6*w4 + 1126*CA4*v6*w4 + 57*v7*w4 + 243*CA*v7*w4 -
     &   1337*CA2*v7*w4 - 195*CA**3*v7*w4 - 1222*CA4*v7*w4 + 
     &   321*v8*w4 + 180*CA*v8*w4 + 869*CA2*v8*w4 - 84*CA**3*v8*w4 +
     &   654*CA4*v8*w4 - 180*v9*w4 - 222*CA*v9*w4 - 148*CA2*v9*w4 +
     &   174*CA**3*v9*w4 - 150*CA4*v9*w4 + 18*v10*w4 + 72*CA*v10*w4 -
     &   126*CA2*v10*w4 - 72*CA**3*v10*w4 + 30*CA4*v10*w4 +
     &   60*CA2*v11*w4 - 16*CA4*v11*w4 + 12*CA*v3*w5 - 24*CA2*v3*w5 -
     &   12*CA**3*v3*w5 + 24*CA4*v3*w5 + 45*v4*w5 - 15*CA*v4*w5 +
     &   149*CA2*v4*w5 + 63*CA**3*v4*w5 - 146*CA4*v4*w5 - 111*v5*w5 -
     &   93*CA*v5*w5 - 245*CA2*v5*w5 - 27*CA**3*v5*w5 +
     &   408*CA4*v5*w5 - 75*v6*w5 + 165*CA*v6*w5 - 45*CA2*v6*w5 -
     &   165*CA**3*v6*w5 - 248*CA4*v6*w5 + 519*v7*w5 + 57*CA*v7*w5 +
     &   419*CA2*v7*w5 + 147*CA**3*v7*w5 - 600*CA4*v7*w5 - 630*v8*w5 -
     &   342*CA*v8*w5 - 208*CA2*v8*w5 + 126*CA**3*v8*w5 +
     &   1106*CA4*v8*w5 + 270*v9*w5 + 330*CA*v9*w5 - 280*CA2*v9*w5 -
     &   246*CA**3*v9*w5 - 748*CA4*v9*w5 - 18*v10*w5 - 114*CA*v10*w5 +
     &   366*CA2*v10*w5 + 114*CA**3*v10*w5 + 204*CA4*v10*w5 -
     &   132*CA2*v11*w5 + 3*CA*v4*w6 - 6*CA2*v4*w6 - 3*CA**3*v4*w6 +
     &   6*CA4*v4*w6 - 9*v5*w6 + 21*CA*v5*w6 + 17*CA2*v5*w6 -
     &   9*CA**3*v5*w6 + 4*CA4*v5*w6 + 129*v6*w6 - 33*CA*v6*w6 +
     &   76*CA2*v6*w6 + 105*CA**3*v6*w6 - 329*CA4*v6*w6 - 381*v7*w6 -
     &   81*CA*v7*w6 - 216*CA2*v7*w6 - 111*CA**3*v7*w6 + 
     &   835*CA4*v7*w6 + 465*v8*w6 + 282*CA*v8*w6 + 35*CA2*v8*w6 -
     &   102*CA**3*v8*w6 - 940*CA4*v8*w6 - 216*v9*w6 - 324*CA*v9*w6 +
     &   370*CA2*v9*w6 + 240*CA**3*v9*w6 + 580*CA4*v9*w6 + 12*v10*w6 +
     &   132*CA*v10*w6 - 456*CA2*v10*w6 - 132*CA**3*v10*w6 -
     &   108*CA4*v10*w6 + 180*CA2*v11*w6 - 48*CA4*v11*w6 - 6*CA*v5*w7 +
     &   12*CA2*v5*w7 + 6*CA**3*v5*w7 - 12*CA4*v5*w7 - 27*v6*w7 +
     &   9*CA*v6*w7 - 71*CA2*v6*w7 - 33*CA**3*v6*w7 + 74*CA4*v6*w7 +
     &   111*v7*w7 + 27*CA*v7*w7 + 144*CA2*v7*w7 + 45*CA**3*v7*w7 -
     &   101*CA4*v7*w7 - 192*v8*w7 - 132*CA*v8*w7 - 41*CA2*v8*w7 +
     &   60*CA**3*v8*w7 + 109*CA4*v8*w7 + 126*v9*w7 + 216*CA*v9*w7 -
     &   254*CA2*v9*w7 - 168*CA**3*v9*w7 - 34*CA4*v9*w7 - 18*v10*w7 -
     &   114*CA*v10*w7 + 390*CA2*v10*w7 + 114*CA**3*v10*w7 - 
     &   172*CA4*v10*w7 - 180*CA2*v11*w7 + 136*CA4*v11*w7 - 
     &   3*CA*v6*w8 + 6*CA2*v6*w8 + 3*CA**3*v6*w8 - 6*CA4*v6*w8 -
     &   9*v7*w8 - 3*CA*v7*w8 - 27*CA2*v7*w8 - 9*CA**3*v7*w8 +
     &   24*CA4*v7*w8 + 39*v8*w8 + 30*CA*v8*w8 + 5*CA2*v8*w8 - 
     &   18*CA**3*v8*w8 - 54*CA4*v8*w8 - 48*v9*w8 - 90*CA*v9*w8 +
     &   130*CA2*v9*w8 + 78*CA**3*v9*w8 - 26*CA4*v9*w8 + 18*v10*w8 +
     &   72*CA*v10*w8 - 246*CA2*v10*w8 - 72*CA**3*v10*w8 +
     &   194*CA4*v10*w8 + 132*CA2*v11*w8 - 132*CA4*v11*w8 + 6*v9*w9 +
     &   18*CA*v9*w9 - 30*CA2*v9*w9 - 18*CA**3*v9*w9 + 24*CA4*v9*w9 -
     &   6*v10*w9 - 30*CA*v10*w9 + 90*CA2*v10*w9 + 30*CA**3*v10*w9 -
     &   84*CA4*v10*w9 - 60*CA2*v11*w9 + 60*CA4*v11*w9 + 6*CA*v10*w10 -
     &   12*CA2*v10*w10 - 6*CA**3*v10*w10 + 12*CA4*v10*w10 +
     &   12*CA2*v11*w10 - 12*CA4*v11*w10))/
     &    (3.*CA*(1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**4) 
      part13 = -(2*CF*lv*(4*CA - 16*CA**3 - 24*CA*v + 76*CA**3*v +
     &   64*CA*v2 - 164*CA**3*v2 - 104*CA*v3 + 228*CA**3*v3 + 
     &   120*CA*v4 - 236*CA**3*v4 - 104*CA*v5 + 180*CA**3*v5 + 
     &   64*CA*v6 - 92*CA**3*v6 - 24*CA*v7 + 28*CA**3*v7 + 4*CA*v8 -
     &   4*CA**3*v8 - 3*w - 2*CA*w + 18*CA2*w + 2*CA**3*w -
     &   15*CA4*w + 21*v*w + 2*CA*v*w - 126*CA2*v*w - 14*CA**3*v*w +
     &   101*CA4*v*w - 63*v2*w + 6*CA*v2*w + 396*CA2*v2*w +
     &   18*CA**3*v2*w - 313*CA4*v2*w + 107*v3*w - 10*CA*v3*w -
     &   746*CA2*v3*w - 42*CA**3*v3*w + 595*CA4*v3*w - 115*v4*w +
     &   6*CA*v4*w + 948*CA2*v4*w + 130*CA**3*v4*w - 773*CA4*v4*w +
     &   83*v5*w - 10*CA*v5*w - 858*CA2*v5*w - 138*CA**3*v5*w +
     &   715*CA4*v5*w - 41*v6*w + 26*CA*v6*w + 556*CA2*v6*w +
     &   14*CA**3*v6*w - 471*CA4*v6*w + 13*v7*w - 38*CA*v7*w -
     &   246*CA2*v7*w + 58*CA**3*v7*w + 213*CA4*v7*w - 2*v8*w +
     &   28*CA*v8*w + 66*CA2*v8*w - 36*CA**3*v8*w - 60*CA4*v8*w -
     &   8*CA*v9*w - 8*CA2*v9*w + 8*CA**3*v9*w + 8*CA4*v9*w + w2 +
     &   CA*w2 - 2*CA2*w2 - CA**3*w2 + CA4*w2 - 13*v*w2 -
     &   10*CA*v*w2 + 52*CA2*v*w2 + 8*CA**3*v*w2 - 39*CA4*v*w2 +
     &   45*v2*w2 + 32*CA*v2*w2 - 268*CA2*v2*w2 + 6*CA**3*v2*w2 + 
     &   259*CA4*v2*w2 - 73*v3*w2 - 54*CA*v3*w2 + 678*CA2*v3*w2 -
     &   16*CA**3*v3*w2 - 841*CA4*v3*w2 + 59*v4*w2 - 2*CA*v4*w2 - 
     &   1046*CA2*v4*w2 - 64*CA**3*v4*w2 + 1627*CA4*v4*w2 - 3*v5*w2 +
     &   138*CA*v5*w2 + 1040*CA2*v5*w2 - 20*CA**3*v5*w2 - 
     &   1985*CA4*v5*w2 - 49*v6*w2 - 188*CA*v6*w2 - 608*CA2*v6*w2 +
     &   302*CA**3*v6*w2 + 1509*CA4*v6*w2 + 53*v7*w2 + 134*CA*v7*w2 +
     &   102*CA2*v7*w2 - 312*CA**3*v7*w2 - 639*CA4*v7*w2 - 24*v8*w2 -
     &   71*CA*v8*w2 + 116*CA2*v8*w2 + 113*CA**3*v8*w2 + 76*CA4*v8*w2 +
     &   4*v9*w2 + 16*CA*v9*w2 - 80*CA2*v9*w2 - 12*CA**3*v9*w2 +
     &   48*CA4*v9*w2 + 4*CA*v10*w2 + 16*CA2*v10*w2 - 4*CA**3*v10*w2 -
     &   16*CA4*v10*w2 + 2*v*w3 + 2*CA*v*w3 - 4*CA2*v*w3 -
     &   2*CA**3*v*w3 + 2*CA4*v*w3 - 7*v2*w3 - 6*CA*v2*w3 +
     &   6*CA2*v2*w3 + 2*CA**3*v2*w3 + CA4*v2*w3 + v3*w3 +
     &   10*CA*v3*w3 + 62*CA2*v3*w3 - 4*CA**3*v3*w3 - 71*CA4*v3*w3 +
     &   34*v4*w3 + 72*CA*v4*w3 - 304*CA2*v4*w3 + 56*CA**3*v4*w3 +
     &   386*CA4*v4*w3 - 110*v5*w3 - 224*CA*v5*w3 + 758*CA2*v5*w3 + 
     &   122*CA**3*v5*w3 - 1148*CA4*v5*w3 + 175*v6*w3 + 150*CA*v6*w3 -
     &   1276*CA2*v6*w3 - 514*CA**3*v6*w3 + 2089*CA4*v6*w3 -
     &   127*v7*w3 + 18*CA*v7*w3 + 1434*CA2*v7*w3 + 400*CA**3*v7*w3 -
     &   2359*CA4*v7*w3 + 24*v8*w3 - 4*CA*v8*w3 - 984*CA2*v8*w3 -
     &   44*CA**3*v8*w3 + 1616*CA4*v8*w3 + 10*v9*w3 + 2*CA*v9*w3 +
     &   350*CA2*v9*w3 - 36*CA**3*v9*w3 - 608*CA4*v9*w3 - 2*v10*w3 -
     &   20*CA*v10*w3 - 34*CA2*v10*w3 + 20*CA**3*v10*w3 +
     &   84*CA4*v10*w3 - 8*CA2*v11*w3 + 8*CA4*v11*w3 - v2*w4 -
     &   CA*v2*w4 + 2*CA2*v2*w4 + CA**3*v2*w4 - CA4*v2*w4 + 23*v3*w4 +
     &   18*CA*v3*w4 - 96*CA2*v3*w4 - 16*CA**3*v3*w4 + 73*CA4*v3*w4 -
     &   78*v4*w4 - 93*CA*v4*w4 + 448*CA2*v4*w4 + 21*CA**3*v4*w4 -
     &   482*CA4*v4*w4 + 158*v5*w4 + 166*CA*v5*w4 - 1076*CA2*v5*w4 -
     &   148*CA**3*v5*w4 + 1550*CA4*v5*w4 - 195*v6*w4 + 81*CA*v6*w4 +
     &   1720*CA2*v6*w4 + 391*CA**3*v6*w4 - 2833*CA4*v6*w4 + 63*v7*w4 -
     &   318*CA*v7*w4 - 1816*CA2*v7*w4 - 140*CA**3*v7*w4 +
     &   3013*CA4*v7*w4 + 88*v8*w4 + 105*CA*v8*w4 + 1084*CA2*v8*w4 -
     &   193*CA**3*v8*w4 - 1756*CA4*v8*w4 - 60*v9*w4 - 2*CA*v9*w4 -
     &   204*CA2*v9*w4 + 112*CA**3*v9*w4 + 372*CA4*v9*w4 + 2*v10*w4 +
     &   44*CA*v10*w4 - 110*CA2*v10*w4 - 44*CA**3*v10*w4 +
     &   136*CA4*v10*w4 + 48*CA2*v11*w4 - 72*CA4*v11*w4 - 4*v3*w5 -
     &   4*CA*v3*w5 + 8*CA2*v3*w5 + 4*CA**3*v3*w5 - 4*CA4*v3*w5 +
     &   23*v4*w5 + 18*CA*v4*w5 - 66*CA2*v4*w5 - 10*CA**3*v4*w5 +
     &   43*CA4*v4*w5 - 59*v5*w5 - 18*CA*v5*w5 + 188*CA2*v5*w5 + 
     &   38*CA**3*v5*w5 - 189*CA4*v5*w5 + 59*v6*w5 - 198*CA*v6*w5 -
     &   338*CA2*v6*w5 - 98*CA**3*v6*w5 + 275*CA4*v6*w5 + 95*v7*w5 +
     &   350*CA*v7*w5 + 306*CA2*v7*w5 - 134*CA**3*v7*w5 +
     &   135*CA4*v7*w5 - 214*v8*w5 - 20*CA*v8*w5 + 122*CA2*v8*w5 +
     &   308*CA**3*v8*w5 - 856*CA4*v8*w5 + 86*v9*w5 - 68*CA*v9*w5 -
     &   522*CA2*v9*w5 - 112*CA**3*v9*w5 + 1100*CA4*v9*w5 + 14*v10*w5 -
     &   60*CA*v10*w5 + 414*CA2*v10*w5 + 60*CA**3*v10*w5 -
     &   688*CA4*v10*w5 - 112*CA2*v11*w5 + 184*CA4*v11*w5 - v4*w6 -
     &   CA*v4*w6 +2*CA2*v4*w6 + CA**3*v4*w6 - CA4*v4*w6 -
     &   7*v5*w6 - 6*CA*v5*w6 + 36*CA2*v5*w6 + 8*CA**3*v5*w6 -
     &   29*CA4*v5*w6 + 37*v6*w6 + 86*CA*v6*w6 - 156*CA2*v6*w6 - 
     &   16*CA**3*v6*w6 + 251*CA4*v6*w6 - 141*v7*w6 - 130*CA*v7*w6 +
     &   350*CA2*v7*w6 + 128*CA**3*v7*w6 - 797*CA4*v7*w6 + 192*v8*w6 -
     &   125*CA*v8*w6 - 632*CA2*v8*w6 - 181*CA**3*v8*w6 +
     &   1316*CA4*v8*w6 - 40*v9*w6 + 120*CA*v9*w6 + 772*CA2*v9*w6 + 
     &   36*CA**3*v9*w6 - 1268*CA4*v9*w6 - 40*v10*w6 + 64*CA*v10*w6 -
     &   508*CA2*v10*w6 - 64*CA**3*v10*w6 + 728*CA4*v10*w6 +
     &   136*CA2*v11*w6 - 200*CA4*v11*w6 + 2*v5*w7 + 2*CA*v5*w7 -
     &   4*CA2*v5*w7 - 2*CA**3*v5*w7 + 2*CA4*v5*w7 - 13*v6*w7 -
     &   10*CA*v6*w7 + 42*CA2*v6*w7 + 6*CA**3*v6*w7 - 29*CA4*v6*w7 +
     &   51*v7*w7 + 10*CA*v7*w7 - 122*CA2*v7*w7 - 28*CA**3*v7*w7 +
     &   135*CA4*v7*w7 - 68*v8*w7 + 108*CA*v8*w7 + 260*CA2*v8*w7 + 
     &   36*CA**3*v8*w7 - 292*CA4*v8*w7 - 18*v9*w7 - 74*CA*v9*w7 -
     &   394*CA2*v9*w7 + 12*CA**3*v9*w7 + 348*CA4*v9*w7 + 46*v10*w7 -
     &   60*CA*v10*w7 + 322*CA2*v10*w7 + 60*CA**3*v10*w7 -
     &   268*CA4*v10*w7 - 104*CA2*v11*w7 + 104*CA4*v11*w7 + v6*w8 +
     &   CA*v6*w8 - 2*CA2*v6*w8 - CA**3*v6*w8 + CA4*v6*w8 - 3*v7*w8 -
     &   2*CA*v7*w8 + 8*CA2*v7*w8 - 5*CA4*v7*w8 + 4*v8*w8 - 
     &   25*CA*v8*w8 - 32*CA2*v8*w8 + CA**3*v8*w8 + 4*CA4*v8*w8 +
     &   24*v9*w8 + 10*CA*v9*w8 + 104*CA2*v9*w8 - 8*CA**3*v9*w8 -
     &   26*v10*w8 + 44*CA*v10*w8 - 142*CA2*v10*w8 - 44*CA**3*v10*w8 +
     &   40*CA4*v10*w8 + 64*CA2*v11*w8 - 40*CA4*v11*w8 - 6*v9*w9 +
     &   4*CA*v9*w9 - 18*CA2*v9*w9 + 6*v10*w9 - 20*CA*v10*w9 +
     &   50*CA2*v10*w9 + 20*CA**3*v10*w9 - 24*CA4*v10*w9 - 
     &   32*CA2*v11*w9 + 24*CA4*v11*w9 + 4*CA*v10*w10 - 8*CA2*v10*w10 -
     &   4*CA**3*v10*w10 + 8*CA4*v10*w10 + 8*CA2*v11*w10 - 
     &   8*CA4*v11*w10))/(CA*(1 - v)**2*v2*w2*(1 - v*w)**2*
     &   (1 - v + v*w)**4)
      part14 = -(2*CF*l1w*(4*CA - 12*CA**3 - 24*CA*v + 56*CA**3*v +
     &   64*CA*v2 - 120*CA**3*v2 - 104*CA*v3 + 168*CA**3*v3 + 
     &   120*CA*v4 - 176*CA**3*v4 - 104*CA*v5 + 136*CA**3*v5 + 
     &   64*CA*v6 - 72*CA**3*v6 - 24*CA*v7 + 24*CA**3*v7 + 4*CA*v8 -
     &   4*CA**3*v8 - 5*w - 2*CA*w + 16*CA2*w + 2*CA**3*w - 11*CA4*w +
     &   31*v*w + 6*CA*v*w - 114*CA2*v*w - 14*CA**3*v*w +
     &   75*CA4*v*w - 85*v2*w - 10*CA*v2*w + 364*CA2*v2*w +
     &   22*CA**3*v2*w - 237*CA4*v2*w + 133*v3*w + 10*CA*v3*w -
     &   698*CA2*v3*w - 46*CA**3*v3*w + 463*CA4*v3*w - 125*v4*w +
     &   6*CA*v4*w + 908*CA2*v4*w + 114*CA**3*v4*w - 623*CA4*v4*w +
     &   65*v5*w - 30*CA*v5*w - 846*CA2*v5*w - 114*CA**3*v5*w + 
     &   601*CA4*v5*w - 11*v6*w + 42*CA*v6*w + 564*CA2*v6*w +
     &   18*CA**3*v6*w - 415*CA4*v6*w - 5*v7*w - 42*CA*v7*w - 
     &   254*CA2*v7*w + 38*CA**3*v7*w + 197*CA4*v7*w + 2*v8*w +
     &   28*CA*v8*w + 68*CA2*v8*w - 28*CA**3*v8*w - 58*CA4*v8*w -
     &   8*CA*v9*w - 8*CA2*v9*w + 8*CA**3*v9*w + 8*CA4*v9*w + w2 +
     &   CA*w2 - 2*CA2*w2 -  CA**3*w2 + CA4*w2 - 15*v*w2 - 10*CA*v*w2 +
     &   50*CA2*v*w2 + 8*CA**3*v*w2 - 31*CA4*v*w2 + 49*v2*w2 +
     &   28*CA*v2*w2 - 270*CA2*v2*w2 - 2*CA**3*v2*w2 + 205*CA4*v2*w2 -
     &   67*v3*w2 - 26*CA*v3*w2 + 732*CA2*v3*w2 - 675*CA4*v3*w2 + 
     &   23*v4*w2 - 66*CA*v4*w2 - 1236*CA2*v4*w2 - 68*CA**3*v4*w2 +
     &   1335*CA4*v4*w2 + 63*v5*w2 + 194*CA*v5*w2 + 1374*CA2*v5*w2 +
     &   8*CA**3*v5*w2 - 1681*CA4*v5*w2 - 101*v6*w2 - 192*CA*v6*w2 -
     &   942*CA2*v6*w2 + 206*CA**3*v6*w2 + 1335*CA4*v6*w2 + 55*v7*w2 +
     &   114*CA*v7*w2 + 288*CA2*v7*w2 - 224*CA**3*v7*w2 - 
     &   601*CA4*v7*w2 - 4*v8*w2 - 63*CA*v8*w2 + 66*CA2*v8*w2 +
     &   93*CA**3*v8*w2 + 84*CA4*v8*w2 - 4*v9*w2 + 16*CA*v9*w2 -
     &   76*CA2*v9*w2 - 16*CA**3*v9*w2 + 44*CA4*v9*w2 + 4*CA*v10*w2 +
     &   16*CA2*v10*w2 - 4*CA**3*v10*w2 - 16*CA4*v10*w2 + 2*v*w3 +
     &   2*CA*v*w3 - 4*CA2*v*w3 - 2*CA**3*v*w3 + 2*CA4*v*w3 + 3*v2*w3 -
     &   6*CA*v2*w3 + 12*CA2*v2*w3 + 2*CA**3*v2*w3 - 3*CA4*v2*w3 -
     &   57*v3*w3 - 2*CA*v3*w3 + 20*CA2*v3*w3 - 43*CA4*v3*w3 +
     &   174*v4*w3 + 96*CA*v4*w3 - 188*CA2*v4*w3 + 44*CA**3*v4*w3 +
     &   276*CA4*v4*w3 - 286*v5*w3 - 204*CA*v5*w3 + 612*CA2*v5*w3 +
     &   70*CA**3*v5*w3 - 868*CA4*v5*w3 + 269*v6*w3 + 70*CA*v6*w3 -
     &   1260*CA2*v6*w3 - 362*CA**3*v6*w3 + 1651*CA4*v6*w3 - 97*v7*w3 +
     &   78*CA*v7*w3 + 1620*CA2*v7*w3 + 300*CA**3*v7*w3 - 
     &   1963*CA4*v7*w3 - 32*v8*w3 - 12*CA*v8*w3 - 1204*CA2*v8*w3 -
     &   56*CA**3*v8*w3 + 1430*CA4*v8*w3 + 22*v9*w3 - 2*CA*v9*w3 +
     &   448*CA2*v9*w3 - 16*CA**3*v9*w3 - 576*CA4*v9*w3 + 2*v10*w3 -
     &   20*CA*v10*w3 - 48*CA2*v10*w3 + 20*CA**3*v10*w3 +
     &   86*CA4*v10*w3 - 8*CA2*v11*w3 + 8*CA4*v11*w3 - v2*w4 -
     &   CA*v2*w4 + 2*CA2*v2*w4 + CA**3*v2*w4 - CA4*v2*w4 + 41*v3*w4 +
     &   18*CA*v3*w4 - 86*CA2*v3*w4 - 16*CA**3*v3*w4 + 57*CA4*v3*w4 -
     &   146*v4*w4 - 81*CA*v4*w4 + 430*CA2*v4*w4 + 21*CA**3*v4*w4 -
     &   368*CA4*v4*w4 + 242*v5*w4 + 106*CA*v5*w4 - 1144*CA2*v5*w4 -
     &   100*CA**3*v5*w4 + 1204*CA4*v5*w4 - 193*v6*w4 + 161*CA*v6*w4 +
     &   2028*CA2*v6*w4 + 291*CA**3*v6*w4 - 2285*CA4*v6*w4 - 31*v7*w4 -
     &   334*CA*v7*w4 - 2326*CA2*v7*w4 - 144*CA**3*v7*w4 +
     &   2565*CA4*v7*w4 + 152*v8*w4 + 77*CA*v8*w4 + 1482*CA2*v8*w4 -
     &   101*CA**3*v8*w4 - 1618*CA4*v8*w4 - 52*v9*w4 + 10*CA*v9*w4 -
     &   316*CA2*v9*w4 + 76*CA**3*v9*w4 + 406*CA4*v9*w4 - 12*v10*w4 +
     &   44*CA*v10*w4 - 126*CA2*v10*w4 - 44*CA**3*v10*w4 +
     &   112*CA4*v10*w4 + 56*CA2*v11*w4 - 72*CA4*v11*w4 - 4*v3*w5 -
     &   4*CA*v3*w5 + 8*CA2*v3*w5 + 4*CA**3*v3*w5 - 4*CA4*v3*w5 +
     &   25*v4*w5 + 18*CA*v4*w5 - 60*CA2*v4*w5 - 10*CA**3*v4*w5 +
     &   39*CA4*v4*w5 - 33*v5*w5 - 6*CA*v5*w5 + 186*CA2*v5*w5 + 
     &   26*CA**3*v5*w5 - 165*CA4*v5*w5 - 47*v6*w5 - 198*CA*v6*w5 -
     &   380*CA2*v6*w5 - 86*CA**3*v6*w5 + 265*CA4*v6*w5 + 229*v7*w5 +
     &   298*CA*v7*w5 + 374*CA2*v7*w5 - 42*CA**3*v7*w5 + 15*CA4*v7*w5 -
     &   262*v8*w5 + 28*CA*v8*w5 + 160*CA2*v8*w5 + 188*CA**3*v8*w5 -
     &   606*CA4*v8*w5 + 62*v9*w5 - 76*CA*v9*w5 - 688*CA2*v9*w5 -
     &   92*CA**3*v9*w5 + 902*CA4*v9*w5 + 30*v10*w5 - 60*CA*v10*w5 +
     &   536*CA2*v10*w5 + 60*CA**3*v10*w5 - 630*CA4*v10*w5 - 
     &   136*CA2*v11*w5 + 184*CA4*v11*w5 - v4*w6 - CA*v4*w6 +
     &   2*CA2*v4*w6 + CA**3*v4*w6 - CA4*v4*w6 - 21*v5*w6 - 
     &   6*CA*v5*w6 + 38*CA2*v5*w6 + 8*CA**3*v5*w6 - 21*CA4*v5*w6 +
     &   97*v6*w6 + 74*CA*v6*w6 - 178*CA2*v6*w6 - 8*CA**3*v6*w6 +
     &   185*CA4*v6*w6 - 215*v7*w6 - 94*CA*v7*w6 + 436*CA2*v7*w6 +
     &   64*CA**3*v7*w6 - 603*CA4*v7*w6 + 212*v8*w6 - 141*CA*v8*w6 -
     &   814*CA2*v8*w6 - 121*CA**3*v8*w6 + 1052*CA4*v8*w6 - 32*v9*w6 +
     &   112*CA*v9*w6 + 984*CA2*v9*w6 + 56*CA**3*v9*w6 - 
     &   1092*CA4*v9*w6 - 40*v10*w6 + 64*CA*v10*w6 - 620*CA2*v10*w6 -
     &   64*CA**3*v10*w6 + 680*CA4*v10*w6 + 152*CA2*v11*w6 -
     &   200*CA4*v11*w6 + 2*v5*w7 + 2*CA*v5*w7 - 4*CA2*v5*w7 -
     &   2*CA**3*v5*w7 + 2*CA4*v5*w7 - 23*v6*w7 - 10*CA*v6*w7 +
     &   36*CA2*v6*w7 + 6*CA**3*v6*w7 - 25*CA4*v6*w7 + 69*v7*w7 +
     &   6*CA*v7*w7 - 104*CA2*v7*w7 - 16*CA**3*v7*w7 + 107*CA4*v7*w7 -
     &   76*v8*w7 + 100*CA*v8*w7 + 236*CA2*v8*w7 + 32*CA**3*v8*w7 -
     &   246*CA4*v8*w7 - 2*v9*w7 - 62*CA*v9*w7 - 372*CA2*v9*w7 -
     &   24*CA**3*v9*w7 + 320*CA4*v9*w7 + 30*v10*w7 - 60*CA*v10*w7 +
     &   296*CA2*v10*w7 + 60*CA**3*v10*w7 - 262*CA4*v10*w7 -
     &   88*CA2*v11*w7 + 104*CA4*v11*w7 + v6*w8 + CA*v6*w8 -
     &   2*CA2*v6*w8 - CA**3*v6*w8 + CA4*v6*w8 - 5*v7*w8 - 2*CA*v7*w8 -
     &   2*CA2*v7*w8 - 5*CA4*v7*w8 + 8*v8*w8 - 21*CA*v8*w8 +
     &   10*CA2*v8*w8 - 3*CA**3*v8*w8 + 10*CA4*v8*w8 + 8*v9*w8 + 
     &   6*CA*v9*w8 + 32*CA2*v9*w8 + 12*CA**3*v9*w8 - 14*CA4*v9*w8 -
     &   12*v10*w8 + 44*CA*v10*w8 - 78*CA2*v10*w8 - 44*CA**3*v10*w8 +
     &        48*CA4*v10*w8 + 40*CA2*v11*w8 -
     &        40*CA4*v11*w8 - 4*CA2*v8*w9 - 2*v9*w9 +
     &        4*CA*v9*w9 - 4*CA2*v9*w9 - 4*CA**3*v9*w9 +
     &        2*CA4*v9*w9 + 2*v10*w9 - 20*CA*v10*w9 +
     &        32*CA2*v10*w9 + 20*CA**3*v10*w9 -
     &        26*CA4*v10*w9 - 24*CA2*v11*w9 +
     &        24*CA4*v11*w9 + 4*CA*v10*w10 - 8*CA2*v10*w10 -
     &        4*CA**3*v10*w10 + 8*CA4*v10*w10 +
     &        8*CA2*v11*w10 - 8*CA4*v11*w10))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**2*(1 - v + v*w)**4)
      struv12= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + 
     &         part9 + part10 + part11 + part12 +
     &         part13 + part14
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION STRUV13(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (8*CA*CF**2*l1w*Nf*(1 + v2*w2)*
     &      (1 - 2*v + 2*v2 - 2*v2*w + v2*w2))/
     &     ((1 - v)*(1 - v*w)**4)
      part2 = -(8*CA*CF**2*lms*Nf*(1 + v2*w2)*
     &      (1 - 2*v + 2*v2 - 2*v2*w + v2*w2))/
     &     ((1 - v)*(1 - v*w)**4)
      part3 = (8*CA*CF**2*lv*Nf*(1 + v2*w2)*
     &      (1 - 2*v + 2*v2 - 2*v2*w + v2*w2))/
     &     ((1 - v)*(1 - v*w)**4)
      part4 = -(4*CF*lvw*(2*CA - 2*CA**3 - 6*CA*v + 6*CA**3*v + 
     &   8*CA*v2 - 8*CA**3*v2 - 4*CA*v3 + 4*CA**3*v3 - 3*w - 2*CA*w +
     &        2*CA2*w + 2*CA**3*w - 3*CA4*w + 11*v*w + 2*CA*v*w +
     &        2*CA2*v*w - 8*CA**3*v*w + 4*CA4*v*w - 16*v2*w +
     &        6*CA*v2*w - 6*CA2*v2*w + 18*CA**3*v2*w -
     &        6*CA4*v2*w + 12*v3*w - 14*CA*v3*w + 10*CA2*v3*w -
     &        16*CA**3*v3*w + 3*CA4*v3*w - 4*v4*w + 8*CA*v4*w -
     &        4*CA2*v4*w + 4*CA**3*v4*w - 2*CA4*v4*w - v*w2 +
     &        4*CA*v*w2 - 5*CA2*v*w2 + 2*CA**3*v*w2 -
     &        6*CA4*v*w2 + 3*v2*w2 - 14*CA*v2*w2 +
     &        CA2*v2*w2 - 16*CA**3*v2*w2 + 14*CA4*v2*w2 -
     &        6*v3*w2 + 16*CA*v3*w2 - 10*CA2*v3*w2 +
     &        26*CA**3*v3*w2 - 13*CA4*v3*w2 + 4*v4*w2 -
     &        2*CA*v4*w2 + 6*CA2*v4*w2 - 16*CA**3*v4*w2 +
     &        5*CA4*v4*w2 - 4*CA*v5*w2 + 4*CA**3*v5*w2 -
     &        v2*w3 + 7*CA2*v2*w3 + 6*CA**3*v2*w3 -
     &        14*CA4*v2*w3 + 5*v3*w3 + 4*CA*v3*w3 +
     &        CA2*v3*w3 - 16*CA**3*v3*w3 + 15*CA4*v3*w3 -
     &        4*v4*w3 - 12*CA*v4*w3 + 18*CA**3*v4*w3 -
     &        9*CA4*v4*w3 + 8*CA*v5*w3 - 8*CA**3*v5*w3 -
     &        v3*w4 - 2*CA*v3*w4 - 4*CA2*v3*w4 +
     &        2*CA**3*v3*w4 - 7*CA4*v3*w4 + v4*w4 +
     &        8*CA*v4*w4 - 4*CA2*v4*w4 - 8*CA**3*v4*w4 +
     &        7*CA4*v4*w4 - 6*CA*v5*w4 + 6*CA**3*v5*w4 -
     &        2*CA*v4*w5 + 4*CA2*v4*w5 + 2*CA**3*v4*w5 -
     &        4*CA4*v4*w5 + 2*CA*v5*w5 - 2*CA**3*v5*w5))/
     &    (CA*(1 - v)**2*v2*w2) 
      part5 = -(4*CF*l1v*(2*CA**3 - 10*CA**3*v + 22*CA**3*v2 - 
     &   26*CA**3*v3 + 16*CA**3*v4 - 4*CA**3*v5 - w - 3*CA2*w + 
     &   2*CA4*w + 5*v*w - 6*CA*v*w + 10*CA2*v*w + 4*CA**3*v*w - 
     &     7*CA4*v*w - 11*v2*w + 30*CA*v2*w - 15*CA2*v2*w - 
     &     20*CA**3*v2*w + 11*CA4*v2*w + 11*v3*w - 70*CA*v3*w + 
     &     12*CA2*v3*w + 38*CA**3*v3*w - 10*CA4*v3*w - 4*v4*w + 
     &     90*CA*v4*w - 4*CA2*v4*w - 32*CA**3*v4*w + 5*CA4*v4*w -
     &        60*CA*v5*w + 10*CA**3*v5*w - CA4*v5*w + 16*CA*v6*w -
     &        v*w2 - CA2*v*w2 + 4*v2*w2 - 10*CA*v2*w2 -
     &        CA2*v2*w2 + 6*CA**3*v2*w2 - 9*CA4*v2*w2 -
     &        5*v3*w2 + 56*CA*v3*w2 + 2*CA2*v3*w2 -
     &        24*CA**3*v3*w2 + 16*CA4*v3*w2 - 118*CA*v4*w2 -
     &        4*CA2*v4*w2 + 36*CA**3*v4*w2 - 11*CA4*v4*w2 +
     &        2*v5*w2 + 108*CA*v5*w2 + 4*CA2*v5*w2 -
     &        28*CA**3*v5*w2 + 3*CA4*v5*w2 - 36*CA*v6*w2 +
     &        14*CA**3*v6*w2 + CA4*v6*w2 - 4*CA**3*v7*w2 +
     &        3*v2*w3 + 5*CA2*v2*w3 - 2*CA4*v2*w3 -
     &        7*v3*w3 - 14*CA*v3*w3 - 19*CA2*v3*w3 +
     &        4*CA**3*v3*w3 - 8*CA4*v3*w3 + 6*v4*w3 +
     &        58*CA*v4*w3 + 22*CA2*v4*w3 - 14*CA**3*v4*w3 +
     &        7*CA4*v4*w3 - 74*CA*v5*w3 - 12*CA2*v5*w3 +
     &        24*CA**3*v5*w3 - 3*CA4*v5*w3 - 2*v6*w3 +
     &        30*CA*v6*w3 - 22*CA**3*v6*w3 - 2*CA4*v6*w3 +
     &        8*CA**3*v7*w3 - v3*w4 - CA2*v3*w4 + 4*v4*w4 -
     &        10*CA*v4*w4 + 3*CA2*v4*w4 + 2*CA**3*v4*w4 -
     &        3*CA4*v4*w4 - 8*v5*w4 + 20*CA*v5*w4 -
     &        6*CA2*v5*w4 - 10*CA**3*v5*w4 + 3*CA4*v5*w4 +
     &        5*v6*w4 - 10*CA*v6*w4 + 4*CA2*v6*w4 +
     &        14*CA**3*v6*w4 - 6*CA**3*v7*w4 + 2*v5*w5 +
     &        7*CA2*v5*w5 + 2*CA**3*v5*w5 - CA4*v5*w5 -
     &        2*v6*w5 - 3*CA2*v6*w5 - 4*CA**3*v6*w5 +
     &        CA4*v6*w5 + 2*CA**3*v7*w5))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)*(1 - v + v*w))
      part6 = - (4*CF*Nf*(-2*CA*CF - 4*v + 4*CA2*v + 7*v2 - 7*CA2*v2 -
     &        6*v3 + 6*CA2*v3 + 2*v4 - 2*CA2*v4 - 4*v*w +
     &        4*CA2*v*w + 12*v2*w - CA2*v2*w - 11*v3*w -
     &        11*CA2*v3*w - v4*w + 26*CA2*v4*w + 4*v5*w -
     &        14*CA2*v5*w + 6*v2*w2 - 6*CA2*v2*w2 -
     &        26*v3*w2 + 13*CA2*v3*w2 + 76*v4*w2 -
     &        62*CA2*v4*w2 - 92*v5*w2 + 67*CA2*v5*w2 +
     &        36*v6*w2 - 24*CA2*v6*w2 - 4*v3*w3 +
     &        4*CA2*v3*w3 - 20*v4*w3 + 27*CA2*v4*w3 +
     &        50*v5*w3 - 54*CA2*v5*w3 - 30*v6*w3 +
     &        49*CA2*v6*w3 + 4*v7*w3 - 14*CA2*v7*w3 +
     &        v4*w4 - CA2*v4*w4 - 6*v5*w4 + 7*CA2*v5*w4 +
     &        17*v6*w4 - 27*CA2*v6*w4 - 14*v7*w4 +
     &        19*CA2*v7*w4 + 2*v8*w4 - 2*CA2*v8*w4 -
     &        4*v6*w5 + 6*CA2*v6*w5 + 5*v7*w5 -
     &        7*CA2*v7*w5 - v8*w5 + CA2*v8*w5))/
     &    (3.*(1 - v)**2*v2*w*(1 - v*w)**4) 
      part7 = -(4*CF*l1vw*(4 - 2*CA2 + 2*CA4 - 8*CA2*CF**2 - 
     &   22*v + 2*CA*v + 9*CA2*v - 10*CA4*v - 8*CA*CF*v +
     &   8*CA**3*CF*v + 24*CA2*CF**2*v + 52*v2 - 12*CA*v2 - 16*CA2*v2 +
     &        2*CA**3*v2 + 21*CA4*v2 + 28*CA*CF*v2 -
     &        28*CA**3*CF*v2 - 28*CA2*CF**2*v2 - 68*v3 + 28*CA*v3 +
     &        14*CA2*v3 - 8*CA**3*v3 - 24*CA4*v3 - 40*CA*CF*v3 +
     &        40*CA**3*CF*v3 + 16*CA2*CF**2*v3 + 52*v4 - 32*CA*v4 -
     &        6*CA2*v4 + 12*CA**3*v4 + 16*CA4*v4 + 30*CA*CF*v4 -
     &        30*CA**3*CF*v4 - 4*CA2*CF**2*v4 - 22*v5 + 18*CA*v5 +
     &        CA2*v5 - 8*CA**3*v5 - 6*CA4*v5 - 12*CA*CF*v5 +
     &        12*CA**3*CF*v5 + 4*v6 - 4*CA*v6 + 2*CA**3*v6 +
     &        CA4*v6 + 2*CA*CF*v6 - 2*CA**3*CF*v6 + 12*v*w -
     &        2*CA*v*w - 9*CA2*v*w + 12*CA4*v*w - 8*CA**3*CF*v*w -
     &        24*CA2*CF**2*v*w - 57*v2*w + 17*CA*v2*w +
     &        35*CA2*v2*w - 3*CA**3*v2*w - 49*CA4*v2*w -
     &        24*CA*CF*v2*w + 56*CA**3*CF*v2*w + 56*CA2*CF**2*v2*w +
     &        113*v3*w - 51*CA*v3*w - 50*CA2*v3*w + 13*CA**3*v3*w +
     &        82*CA4*v3*w + 68*CA*CF*v3*w - 120*CA**3*CF*v3*w -
     &        48*CA2*CF**2*v3*w - 117*v4*w + 71*CA*v4*w +
     &        30*CA2*v4*w - 21*CA**3*v4*w - 72*CA4*v4*w -
     &        76*CA*CF*v4*w + 120*CA**3*CF*v4*w + 16*CA2*CF**2*v4*w +
     &        63*v5*w - 47*CA*v5*w - 5*CA2*v5*w + 15*CA**3*v5*w +
     &        34*CA4*v5*w + 40*CA*CF*v5*w - 60*CA**3*CF*v5*w -
     &        14*v6*w + 12*CA*v6*w - CA2*v6*w - 4*CA**3*v6*w -
     &        7*CA4*v6*w - 8*CA*CF*v6*w + 12*CA**3*CF*v6*w +
     &        16*v2*w2 - 5*CA*v2*w2 - 24*CA2*v2*w2 +
     &        CA**3*v2*w2 + 31*CA4*v2*w2 - 32*CA**3*CF*v2*w2 -
     &        36*CA2*CF**2*v2*w2 - 66*v3*w2 + 27*CA*v3*w2 +
     &        74*CA2*v3*w2 - 5*CA**3*v3*w2 - 101*CA4*v3*w2 -
     &        36*CA*CF*v3*w2 + 140*CA**3*CF*v3*w2 +
     &        56*CA2*CF**2*v3*w2 + 106*v4*w2 - 53*CA*v4*w2 -
     &        79*CA2*v4*w2 + 9*CA**3*v4*w2 + 130*CA4*v4*w2 +
     &        74*CA*CF*v4*w2 - 206*CA**3*CF*v4*w2 -
     &        28*CA2*CF**2*v4*w2 - 78*v5*w2 + 45*CA*v5*w2 +
     &        32*CA2*v5*w2 - 7*CA**3*v5*w2 - 81*CA4*v5*w2 -
     &        56*CA*CF*v5*w2 + 136*CA**3*CF*v5*w2 + 22*v6*w2 -
     &        14*CA*v6*w2 - 3*CA2*v6*w2 + 2*CA**3*v6*w2 +
     &        21*CA4*v6*w2 + 14*CA*CF*v6*w2 -
     &        34*CA**3*CF*v6*w2 + 13*v3*w3 - 4*CA*v3*w3 -
     &        35*CA2*v3*w3 + 43*CA4*v3*w3 -
     &        60*CA**3*CF*v3*w3 - 24*CA2*CF**2*v3*w3 -
     &        44*v4*w3 + 16*CA*v4*w3 + 79*CA2*v4*w3 -
     &        108*CA4*v4*w3 - 24*CA*CF*v4*w3 +
     &        172*CA**3*CF*v4*w3 + 24*CA2*CF**2*v4*w3 +
     &        50*v5*w3 - 20*CA*v5*w3 - 59*CA2*v5*w3 +
     &        99*CA4*v5*w3 + 36*CA*CF*v5*w3 -
     &        168*CA**3*CF*v5*w3 - 19*v6*w3 + 8*CA*v6*w3 +
     &        15*CA2*v6*w3 - 34*CA4*v6*w3 - 12*CA*CF*v6*w3 +
     &        56*CA**3*CF*v6*w3 + 6*v4*w4 - 2*CA*v4*w4 -
     &        28*CA2*v4*w4 + 35*CA4*v4*w4 -
     &        60*CA**3*CF*v4*w4 - 8*CA2*CF**2*v4*w4 -
     &        14*v5*w4 + 4*CA*v5*w4 + 45*CA2*v5*w4 -
     &        63*CA4*v5*w4 - 8*CA*CF*v5*w4 +
     &        112*CA**3*CF*v5*w4 + 8*v6*w4 - 2*CA*v6*w4 -
     &        21*CA2*v6*w4 + 32*CA4*v6*w4 + 4*CA*CF*v6*w4 -
     &        56*CA**3*CF*v6*w4 + v5*w5 - 14*CA2*v5*w5 +
     &        17*CA4*v5*w5 - 32*CA**3*CF*v5*w5 - v6*w5 +
     &        14*CA2*v6*w5 - 17*CA4*v6*w5 +
     &        32*CA**3*CF*v6*w5 - 4*CA2*v6*w6 + 4*CA4*v6*w6 -
     &        8*CA**3*CF*v6*w6))/(CA*(1 - v)**2*v2*w*
     &     (1 - v + v*w)**2)
      part8 = -(2*CF*lmss*(4*CA2 - 28*CA2*v + 88*CA2*v2 -
     &        160*CA2*v3 + 180*CA2*v4 - 124*CA2*v5 +
     &        48*CA2*v6 - 8*CA2*v7 + 8*CA*CF**2*w + 16*CA2*v*w +
     &        8*CF*v*w - 8*CA2*CF*v*w - 24*CA*CF**2*v*w -
     &        104*CA2*v2*w - 28*CF*v2*w + 28*CA2*CF*v2*w +
     &        28*CA*CF**2*v2*w + 288*CA2*v3*w + 40*CF*v3*w -
     &        40*CA2*CF*v3*w - 16*CA*CF**2*v3*w - 432*CA2*v4*w -
     &        30*CF*v4*w + 30*CA2*CF*v4*w + 4*CA*CF**2*v4*w +
     &        368*CA2*v5*w + 12*CF*v5*w - 12*CA2*CF*v5*w -
     &        168*CA2*v6*w - 2*CF*v6*w + 2*CA2*CF*v6*w +
     &        32*CA2*v7*w + 16*CA*CF**2*v*w2 - 2*v2*w2 +
     &        32*CA2*v2*w2 + 16*CF*v2*w2 - 32*CA2*CF*v2*w2 -
     &        44*CA*CF**2*v2*w2 + 10*v3*w2 - 176*CA2*v3*w2 -
     &        52*CF*v3*w2 + 88*CA2*CF*v3*w2 +
     &        48*CA*CF**2*v3*w2 - 22*v4*w2 + 392*CA2*v4*w2 +
     &        70*CF*v4*w2 - 98*CA2*CF*v4*w2 -
     &        20*CA*CF**2*v4*w2 + 26*v5*w2 - 440*CA2*v5*w2 -
     &        44*CF*v5*w2 + 56*CA2*CF*v5*w2 - 16*v6*w2 +
     &        248*CA2*v6*w2 + 10*CF*v6*w2 - 14*CA2*CF*v6*w2 +
     &        4*v7*w2 - 56*CA2*v7*w2 + 8*CA*CF**2*v2*w3 -
     &        6*v3*w3 + 36*CA2*v3*w3 + 8*CF*v3*w3 -
     &        60*CA2*CF*v3*w3 - 24*CA*CF**2*v3*w3 + 28*v4*w3 -
     &        160*CA2*v4*w3 - 28*CF*v4*w3 +
     &        116*CA2*CF*v4*w3 + 24*CA*CF**2*v4*w3 - 50*v5*w3 +
     &        268*CA2*v5*w3 + 36*CF*v5*w3 - 88*CA2*CF*v5*w3 +
     &        40*v6*w3 - 200*CA2*v6*w3 - 12*CF*v6*w3 +
     &        28*CA2*CF*v6*w3 - 12*v7*w3 + 56*CA2*v7*w3 -
     &        9*v4*w4 + 25*CA2*v4*w4 - 60*CA2*CF*v4*w4 -
     &        8*CA*CF**2*v4*w4 + 31*v5*w4 - 83*CA2*v5*w4 -
     &        8*CF*v5*w4 + 80*CA2*CF*v5*w4 - 36*v6*w4 +
     &        92*CA2*v6*w4 + 4*CF*v6*w4 - 32*CA2*CF*v6*w4 +
     &        14*v7*w4 - 34*CA2*v7*w4 - 6*v5*w5 +
     &        10*CA2*v5*w5 - 32*CA2*CF*v5*w5 + 14*v6*w5 -
     &        22*CA2*v6*w5 + 24*CA2*CF*v6*w5 - 8*v7*w5 +
     &        12*CA2*v7*w5 - 2*v6*w6 + 2*CA2*v6*w6 -
     &        8*CA2*CF*v6*w6 + 2*v7*w6 - 2*CA2*v7*w6))/
     &    ((1 - v)**2*v2*w2*(1 - v + v*w)**2) 
      part9 = - (4*CF*lw*(2*CA - 4*CA**3 - 8*CA*v + 22*CA**3*v + 
     &   14*CA*v2 - 56*CA**3*v2 - 12*CA*v3 + 82*CA**3*v3 + 4*CA*v4 -
     &   72*CA**3*v4 + 36*CA**3*v5 - 8*CA**3*v6 - 2*w - 2*CA*w +
     &        5*CA2*w + 4*CA**3*w - 3*CA4*w + 7*v*w + 8*CA*v*w -
     &        15*CA2*v*w - 28*CA**3*v*w + 8*CA4*v*w - 10*v2*w -
     &        14*CA*v2*w + 22*CA2*v2*w + 88*CA**3*v2*w -
     &        12*CA4*v2*w + 7*v3*w + 10*CA*v3*w - 17*CA2*v3*w -
     &        146*CA**3*v3*w + 10*CA4*v3*w - 2*v4*w + 2*CA*v4*w +
     &        7*CA2*v4*w + 126*CA**3*v4*w - 5*CA4*v4*w -
     &        4*CA*v5*w - 2*CA2*v5*w - 44*CA**3*v5*w +
     &        2*CA4*v5*w - 8*CA**3*v6*w + 8*CA**3*v7*w + w2 -
     &        CA2*w2 - v*w2 + 5*CA2*v*w2 + 6*CA**3*v*w2 -
     &        5*CA4*v*w2 - 7*v2*w2 + 3*CA*v2*w2 -
     &        12*CA2*v2*w2 - 37*CA**3*v2*w2 + 16*CA4*v2*w2 +
     &        20*v3*w2 - 14*CA*v3*w2 + 12*CA2*v3*w2 +
     &        78*CA**3*v3*w2 - 24*CA4*v3*w2 - 25*v4*w2 +
     &        21*CA*v4*w2 - 47*CA**3*v4*w2 + 16*CA4*v4*w2 +
     &        16*v5*w2 - 18*CA*v5*w2 - 2*CA2*v5*w2 -
     &        44*CA**3*v5*w2 - 7*CA4*v5*w2 - 4*v6*w2 +
     &        8*CA*v6*w2 + 2*CA2*v6*w2 + 72*CA**3*v6*w2 -
     &        28*CA**3*v7*w2 - v*w3 - CA2*v*w3 + 5*v2*w3 -
     &        3*CA*v2*w3 + CA2*v2*w3 + 5*CA**3*v2*w3 -
     &        7*CA4*v2*w3 - 14*v3*w3 + 24*CA*v3*w3 +
     &        3*CA2*v3*w3 - 14*CA**3*v3*w3 + 21*CA4*v3*w3 +
     &        27*v4*w3 - 52*CA*v4*w3 - 17*CA2*v4*w3 -
     &        24*CA**3*v4*w3 - 21*CA4*v4*w3 - 27*v5*w3 +
     &        51*CA*v5*w3 + 11*CA2*v5*w3 + 103*CA**3*v5*w3 +
     &        11*CA4*v5*w3 + 10*v6*w3 - 20*CA*v6*w3 -
     &        5*CA2*v6*w3 - 110*CA**3*v6*w3 + 40*CA**3*v7*w3 +
     &        v2*w4 + CA2*v2*w4 + v3*w4 - 8*CA*v3*w4 -
     &        6*CA2*v3*w4 - 6*CA4*v3*w4 - 15*v4*w4 +
     &        28*CA*v4*w4 + 20*CA2*v4*w4 + 22*CA**3*v4*w4 +
     &        12*CA4*v4*w4 + 26*v5*w4 - 36*CA*v5*w4 -
     &        13*CA2*v5*w4 - 70*CA**3*v5*w4 - 9*CA4*v5*w4 -
     &        13*v6*w4 + 16*CA*v6*w4 + 6*CA2*v6*w4 +
     &        78*CA**3*v6*w4 - CA4*v6*w4 - 30*CA**3*v7*w4 -
     &        v3*w5 - CA2*v3*w5 + 7*v4*w5 - 3*CA*v4*w5 -
     &        6*CA2*v4*w5 - 5*CA**3*v4*w5 - 15*v5*w5 +
     &        7*CA*v5*w5 + 3*CA2*v5*w5 + 21*CA**3*v5*w5 +
     &        2*CA4*v5*w5 + 9*v6*w5 - 4*CA*v6*w5 -
     &        4*CA2*v6*w5 - 28*CA**3*v6*w5 + 2*CA4*v6*w5 +
     &        12*CA**3*v7*w5 + 2*CA2*v4*w6 + 2*v5*w6 +
     &        CA2*v5*w6 - 2*CA**3*v5*w6 + CA4*v5*w6 -
     &        2*v6*w6 + CA2*v6*w6 + 4*CA**3*v6*w6 -
     &        CA4*v6*w6 - 2*CA**3*v7*w6))/
     &    (CA*(1 - v)**2*v2*(1 - w)*w2*(1 - v*w)*(1 - v + v*w))
      part10 = (2*CF*lms*(4*CA - 4*CA**3 - 12*CA*v + 16*CA**3*v + 
     &     16*CA*v2 - 32*CA**3*v2 - 8*CA*v3 + 36*CA**3*v3 - 
     &     24*CA**3*v4 + 8*CA**3*v5 - w - 2*CA*w + 6*CA2*w + 
     &     2*CA**3*w - 5*CA4*w - 8*CA*v*w - 8*CA2*v*w + 4*CA**3*v*w + 
     &     8*CA4*v*w + 4*v2*w +
     &        30*CA*v2*w + 10*CA2*v2*w - 30*CA**3*v2*w -
     &        14*CA4*v2*w - 6*v3*w - 48*CA*v3*w - 2*CA2*v3*w +
     &        84*CA**3*v3*w + 8*CA4*v3*w + 3*v4*w + 28*CA*v4*w +
     &        2*CA2*v4*w - 116*CA**3*v4*w - 5*CA4*v4*w +
     &        88*CA**3*v5*w - 32*CA**3*v6*w + w2 + CA*w2 -
     &        2*CA2*w2 - CA**3*w2 + CA4*w2 + 2*v*w2 +
     &        5*CA*v*w2 - 22*CA2*v*w2 - 3*CA**3*v*w2 +
     &        20*CA4*v*w2 - 4*v2*w2 + 2*CA*v2*w2 +
     &        22*CA2*v2*w2 + 6*CA**3*v2*w2 - 2*CA4*v2*w2 +
     &        4*v3*w2 - 20*CA*v3*w2 - 36*CA2*v3*w2 -
     &        6*CA**3*v3*w2 + 16*CA4*v3*w2 - 3*v4*w2 +
     &        52*CA*v4*w2 + 10*CA2*v4*w2 - 48*CA**3*v4*w2 +
     &        9*CA4*v4*w2 - 40*CA*v5*w2 - 12*CA2*v5*w2 +
     &        116*CA**3*v5*w2 + 4*CA4*v5*w2 - 112*CA**3*v6*w2 +
     &        48*CA**3*v7*w2 - 4*v*w3 - 4*CA*v*w3 + 8*CA2*v*w3 +
     &        4*CA**3*v*w3 - 4*CA4*v*w3 + 2*v2*w3 +
     &        28*CA2*v2*w3 - 8*CA**3*v2*w3 - 30*CA4*v2*w3 +
     &        2*v3*w3 - 2*CA*v3*w3 - 2*CA2*v3*w3 +
     &        10*CA**3*v3*w3 - 16*CA4*v3*w3 - 6*v4*w3 -
     &        2*CA*v4*w3 + 30*CA2*v4*w3 + 26*CA**3*v4*w3 -
     &        40*CA4*v4*w3 + 12*v5*w3 - 24*CA*v5*w3 -
     &        32*CA**3*v5*w3 - 4*CA4*v5*w3 - 6*v6*w3 +
     &        32*CA*v6*w3 + 20*CA2*v6*w3 - 16*CA**3*v6*w3 -
     &        14*CA4*v6*w3 + 48*CA**3*v7*w3 - 32*CA**3*v8*w3 +
     &        6*v2*w4 + 6*CA*v2*w4 - 12*CA2*v2*w4 -
     &        6*CA**3*v2*w4 + 6*CA4*v2*w4 - 8*v3*w4 -
     &        10*CA*v3*w4 - 12*CA2*v3*w4 + 22*CA**3*v3*w4 +
     &        20*CA4*v3*w4 + 8*v4*w4 + 18*CA*v4*w4 -
     &        48*CA2*v4*w4 - 50*CA**3*v4*w4 + 40*CA4*v4*w4 -
     &        12*v5*w4 - 2*CA*v5*w4 + 12*CA2*v5*w4 +
     &        26*CA**3*v5*w4 + 8*CA4*v5*w4 + 6*v6*w4 +
     &        8*CA*v6*w4 - 28*CA2*v6*w4 + 16*CA**3*v6*w4 +
     &        38*CA4*v6*w4 - 20*CA*v7*w4 - 12*CA2*v7*w4 -
     &        24*CA**3*v7*w4 + 4*CA4*v7*w4 + 8*CA**3*v8*w4 +
     &        8*CA**3*v9*w4 - 4*v3*w5 - 4*CA*v3*w5 +
     &        8*CA2*v3*w5 + 4*CA**3*v3*w5 - 4*CA4*v3*w5 +
     &        7*v4*w5 + 10*CA*v4*w5 - 2*CA2*v4*w5 -
     &        18*CA**3*v4*w5 - 5*CA4*v4*w5 - 8*v5*w5 -
     &        24*CA*v5*w5 + 64*CA2*v5*w5 + 52*CA**3*v5*w5 -
     &        32*CA4*v5*w5 + 8*v6*w5 + 10*CA*v6*w5 -
     &        18*CA2*v6*w5 - 50*CA**3*v6*w5 - 6*CA4*v6*w5 -
     &        6*v7*w5 - 4*CA*v7*w5 + 26*CA2*v7*w5 +
     &        24*CA**3*v7*w5 - 12*CA4*v7*w5 + 3*v8*w5 +
     &        12*CA*v8*w5 + 2*CA2*v8*w5 - 4*CA**3*v8*w5 -
     &        5*CA4*v8*w5 - 8*CA**3*v9*w5 + v4*w6 +
     &        CA*v4*w6 - 2*CA2*v4*w6 - CA**3*v4*w6 +
     &        CA4*v4*w6 - 2*v5*w6 - 3*CA*v5*w6 +
     &        2*CA2*v5*w6 + 5*CA**3*v5*w6 + 4*v6*w6 +
     &        14*CA*v6*w6 - 42*CA2*v6*w6 - 22*CA**3*v6*w6 +
     &        22*CA4*v6*w6 - 4*CA*v7*w6 + 18*CA**3*v7*w6 -
     &        8*CA4*v7*w6 - 3*v8*w6 - 4*CA*v8*w6 -
     &        6*CA2*v8*w6 - 8*CA**3*v8*w6 + 9*CA4*v8*w6 -
     &        4*CA*v9*w6 + 8*CA**3*v9*w6 - 2*v7*w7 -
     &        6*CA*v7*w7 + 18*CA2*v7*w7 + 6*CA**3*v7*w7 -
     &        8*CA4*v7*w7 + 2*v8*w7 + 2*CA*v8*w7 +
     &        2*CA2*v8*w7 - 2*CA**3*v8*w7 - 4*CA4*v8*w7 +
     &        4*CA*v9*w7 - 4*CA**3*v9*w7 + 2*CA*v8*w8 -
     &        4*CA2*v8*w8 - 2*CA**3*v8*w8 + 4*CA4*v8*w8 -
     &        2*CA*v9*w8 + 2*CA**3*v9*w8))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**4)
      part11 = -(2*CF*(12*CA*v - 36*CA**3*v - 48*CA*v2 + 192*CA**3*v2 +
     &   72*CA*v3 - 432*CA**3*v3 - 48*CA*v4 + 528*CA**3*v4 +
     &   12*CA*v5 - 372*CA**3*v5 + 144*CA**3*v6 - 24*CA**3*v7 -
     &   9*w + 3*CA*w + 7*CA2*w - 3*CA**3*w + 2*CA4*w + 48*v*w -
     &   30*CA*v*w - 24*CA2*v*w + 54*CA**3*v*w - 24*CA4*v*w - 
     &   132*v2*w + 84*CA*v2*w + 64*CA2*v2*w - 204*CA**3*v2*w +
     &   68*CA4*v2*w + 228*v3*w - 72*CA*v3*w - 132*CA2*v3*w +
     &   216*CA**3*v3*w - 96*CA4*v3*w - 237*v4*w - 33*CA*v4*w +
     &   159*CA2*v4*w + 273*CA**3*v4*w + 78*CA4*v4*w + 132*v5*w +
     &   78*CA*v5*w - 100*CA2*v5*w - 918*CA**3*v5*w - 32*CA4*v5*w -
     &   30*v6*w - 30*CA*v6*w + 26*CA2*v6*w + 966*CA**3*v6*w + 
     &   4*CA4*v6*w - 480*CA**3*v7*w + 96*CA**3*v8*w - 3*CA*w2 +
     &   6*CA2*w2 + 3*CA**3*w2 - 6*CA4*w2 + 27*v*w2 + 18*CA*v*w2 -
     &   29*CA2*v*w2 - 6*CA**3*v*w2 + 14*CA4*v*w2 - 129*v2*w2 - 
     &   12*CA*v2*w2 + 16*CA2*v2*w2 - 66*CA**3*v2*w2 - 47*CA4*v2*w2 +
     &   243*v3*w2 - 150*CA*v3*w2 + 67*CA2*v3*w2 + 504*CA**3*v3*w2 +
     &   144*CA4*v3*w2 - 285*v4*w2 + 387*CA*v4*w2 - 53*CA2*v4*w2 -
     &   1365*CA**3*v4*w2 - 296*CA4*v4*w2 + 282*v5*w2 - 336*CA*v5*w2 -
     &   82*CA2*v5*w2 + 1698*CA**3*v5*w2 + 354*CA4*v5*w2 - 198*v6*w2 +
     &   84*CA*v6*w2 + 131*CA2*v6*w2 - 780*CA**3*v6*w2 -
     &   227*CA4*v6*w2 + 60*v7*w2 + 12*CA*v7*w2 - 56*CA2*v7*w2 -
     &   324*CA**3*v7*w2 + 64*CA4*v7*w2 + 480*CA**3*v8*w2 - 
     &   144*CA**3*v9*w2 + 6*CA*v*w3 - 12*CA2*v*w3 - 6*CA**3*v*w3 +
     &   12*CA4*v*w3 - 9*v2*w3 - 57*CA*v2*w3 + 35*CA2*v2*w3 +
     &   33*CA**3*v2*w3 - 50*CA4*v2*w3 + 144*v3*w3 + 228*CA*v3*w3 -
     &   13*CA2*v3*w3 - 114*CA**3*v3*w3 - 29*CA4*v3*w3 - 378*v4*w3 - 
     &   330*CA*v4*w3 - 151*CA2*v4*w3 + 78*CA**3*v4*w3 +
     &   413*CA4*v4*w3 + 384*v5*w3 + 66*CA*v5*w3 + 227*CA2*v5*w3 +
     &   624*CA**3*v5*w3 - 871*CA4*v5*w3 - 225*v6*w3 + 183*CA*v6*w3 -
     &   76*CA2*v6*w3 - 1815*CA**3*v6*w3 + 817*CA4*v6*w3 + 144*v7*w3 -
     &   84*CA*v7*w3 - 58*CA2*v7*w3 + 2052*CA**3*v7*w3 -
     &   328*CA4*v7*w3 - 60*v8*w3 - 12*CA*v8*w3 + 60*CA2*v8*w3 -
     &   948*CA**3*v8*w3 + 24*CA4*v8*w3 + 96*CA**3*v10*w3 +
     &   3*CA*v2*w4 - 6*CA2*v2*w4 - 3*CA**3*v2*w4 + 6*CA4*v2*w4 -
     &   45*v3*w4 + 6*CA*v3*w4 + 7*CA2*v3*w4 - 18*CA**3*v3*w4 +
     &   26*CA4*v3*w4 + 51*v4*w4 - 183*CA*v4*w4 + 83*CA2*v4*w4 +
     &   183*CA**3*v4*w4 - 184*CA4*v4*w4 + 213*v5*w4 + 525*CA*v5*w4 -
     &   87*CA2*v5*w4 - 597*CA**3*v5*w4 + 644*CA4*v5*w4 - 423*v6*w4 -
     &   549*CA*v6*w4 - CA2*v6*w4 + 993*CA**3*v6*w4 - 1086*CA4*v6*w4 +
     &   288*v7*w4 + 282*CA*v7*w4 - 52*CA2*v7*w4 -
     &   678*CA**3*v7*w4 + 856*CA4*v7*w4 - 144*v8*w4 - 156*CA*v8*w4 +
     &   52*CA2*v8*w4 - 312*CA**3*v8*w4 - 310*CA4*v8*w4 + 60*v9*w4 +
     &   72*CA*v9*w4 - 56*CA2*v9*w4 + 696*CA**3*v9*w4 + 64*CA4*v9*w4 -
     &        240*CA**3*v10*w4 - 24*CA**3*v11*w4 - 12*CA*v3*w5 +
     &        24*CA2*v3*w5 + 12*CA**3*v3*w5 - 24*CA4*v3*w5 +
     &        45*v4*w5 + 81*CA*v4*w5 - 43*CA2*v4*w5 -
     &        33*CA**3*v4*w5 + 46*CA4*v4*w5 - 216*v5*w5 -
     &        144*CA*v5*w5 - 114*CA2*v5*w5 - 24*CA**3*v5*w5 -
     &        58*CA4*v5*w5 + 264*v6*w5 + 54*CA*v6*w5 +
     &        210*CA2*v6*w5 + 138*CA**3*v6*w5 +
     &        478*CA4*v6*w5 - 132*v7*w5 - 48*CA*v7*w5 +
     &        24*CA2*v7*w5 - 348*CA**3*v7*w5 -
     &        850*CA4*v7*w5 + 81*v8*w5 + 93*CA*v8*w5 -
     &        41*CA2*v8*w5 + 819*CA**3*v8*w5 +
     &        620*CA4*v8*w5 - 12*v9*w5 + 54*CA*v9*w5 +
     &        46*CA2*v9*w5 - 810*CA**3*v9*w5 -
     &        216*CA4*v9*w5 - 30*v10*w5 - 78*CA*v10*w5 +
     &        26*CA2*v10*w5 + 150*CA**3*v10*w5 +
     &        4*CA4*v10*w5 + 96*CA**3*v11*w5 + 3*CA*v4*w6 -
     &        6*CA2*v4*w6 - 3*CA**3*v4*w6 + 6*CA4*v4*w6 +
     &        9*v5*w6 - 42*CA*v5*w6 + 25*CA2*v5*w6 +
     &        30*CA**3*v5*w6 - 46*CA4*v5*w6 + 69*v6*w6 +
     &        132*CA*v6*w6 + 58*CA2*v6*w6 - 18*CA**3*v6*w6 -
     &        179*CA4*v6*w6 - 123*v7*w6 - 114*CA*v7*w6 -
     &        215*CA2*v7*w6 + 12*CA**3*v7*w6 +
     &        540*CA4*v7*w6 + 81*v8*w6 + 165*CA*v8*w6 +
     &        31*CA2*v8*w6 - 339*CA**3*v8*w6 -
     &        600*CA4*v8*w6 - 114*v9*w6 - 312*CA*v9*w6 -
     &        2*CA2*v9*w6 + 426*CA**3*v9*w6 + 334*CA4*v9*w6 +
     &        78*v10*w6 + 144*CA*v10*w6 - 71*CA2*v10*w6 +
     &        60*CA**3*v10*w6 - 7*CA4*v10*w6 + 24*CA*v11*w6 -
     &        168*CA**3*v11*w6 + 6*CA*v5*w7 - 12*CA2*v5*w7 -
     &        6*CA**3*v5*w7 + 12*CA4*v5*w7 - 27*v6*w7 -
     &        27*CA*v6*w7 + CA2*v6*w7 + 3*CA**3*v6*w7 +
     &        2*CA4*v6*w7 + 24*v7*w7 + 24*CA*v7*w7 +
     &        31*CA2*v7*w7 - 6*CA**3*v7*w7 - 89*CA4*v7*w7 -
     &        42*v8*w7 - 114*CA*v8*w7 + 89*CA2*v8*w7 +
     &        174*CA**3*v8*w7 + 205*CA4*v8*w7 + 120*v9*w7 +
     &        210*CA*v9*w7 - 23*CA2*v9*w7 - 180*CA**3*v9*w7 -
     &        247*CA4*v9*w7 - 75*v10*w7 - 27*CA*v10*w7 +
     &        94*CA2*v10*w7 - 165*CA**3*v10*w7 -
     &        19*CA4*v10*w7 - 72*CA*v11*w7 + 180*CA**3*v11*w7 -
     &        3*CA*v6*w8 + 6*CA2*v6*w8 + 3*CA**3*v6*w8 -
     &        6*CA4*v6*w8 + 9*v7*w8 + 18*CA*v7*w8 -
     &        3*CA2*v7*w8 - 6*CA**3*v7*w8 + 6*CA4*v7*w8 +
     &        9*v8*w8 - 3*CA*v8*w8 - 49*CA2*v8*w8 -
     &        33*CA**3*v8*w8 - 18*CA4*v8*w8 - 45*v9*w8 -
     &        15*CA*v9*w8 - 13*CA2*v9*w8 + 3*CA**3*v9*w8 +
     &        104*CA4*v9*w8 + 27*v10*w8 - 75*CA*v10*w8 -
     &        73*CA2*v10*w8 + 159*CA**3*v10*w8 +
     &        46*CA4*v10*w8 + 78*CA*v11*w8 - 126*CA**3*v11*w8 -
     &       6*CA*v9*w9 + 24*CA2*v9*w9 + 18*CA**3*v9*w9 -
     &        24*CA4*v9*w9 + 42*CA*v10*w9 + 36*CA2*v10*w9 -
     &        66*CA**3*v10*w9 - 36*CA4*v10*w9 - 36*CA*v11*w9 +
     &        48*CA**3*v11*w9 - 6*CA*v10*w10 - 12*CA2*v10*w10 +
     &        6*CA**3*v10*w10 + 12*CA4*v10*w10 + 6*CA*v11*w10 -
     &        6*CA**3*v11*w10))/
     &    (3.*CA*(1 - v)**2*v2*w2*(1 - v*w)**4*(1 - v + v*w)**2)
      part12 = -(2*CF*lv*(4*CA - 16*CA**3 - 20*CA*v + 100*CA**3*v +
     &   44*CA*v2 - 284*CA**3*v2 - 52*CA*v3 + 468*CA**3*v3 + 
     &   32*CA*v4 - 476*CA**3*v4 - 8*CA*v5 + 296*CA**3*v5 - 
     &   104*CA**3*v6 + 16*CA**3*v7 - 3*w - 2*CA*w + 18*CA2*w + 
     &   2*CA**3*w - 15*CA4*w + 12*v*w + 20*CA*v*w - 72*CA2*v*w -
     &   8*CA**3*v*w + 64*CA4*v*w - 18*v2*w - 84*CA*v2*w + 
     &   126*CA2*v2*w - 12*CA**3*v2*w - 128*CA4*v2*w + 10*v3*w +
     &   196*CA*v3*w - 118*CA2*v3*w + 180*CA**3*v3*w + 152*CA4*v3*w +
     &   3*v4*w - 278*CA*v4*w + 56*CA2*v4*w - 578*CA**3*v4*w -
     &   111*CA4*v4*w - 6*v5*w + 248*CA*v5*w - 10*CA2*v5*w +
     &   908*CA**3*v5*w + 48*CA4*v5*w + 2*v6*w - 132*CA*v6*w -
     &   772*CA**3*v6*w - 10*CA4*v6*w + 32*CA*v7*w + 344*CA**3*v7*w -
     &   64*CA**3*v8*w + w2 + CA*w2 - 2*CA2*w2 - CA**3*w2 + CA4*w2 +
     &   2*v*w2 - CA*v*w2 - 30*CA2*v*w2 + 3*CA**3*v*w2 + 28*CA4*v*w2 -
     &   30*v2*w2 - 13*CA*v2*w2 + 142*CA2*v2*w2 + 31*CA**3*v2*w2 -
     &   76*CA4*v2*w2 + 88*v3*w2 + 51*CA*v3*w2 - 276*CA2*v3*w2 -
     &   233*CA**3*v3*w2 + 100*CA4*v3*w2 - 135*v4*w2 - 152*CA*v4*w2 +
     &   310*CA2*v4*w2 + 654*CA**3*v4*w2 - 127*CA4*v4*w2 + 122*v5*w2 +
     &   338*CA*v5*w2 - 186*CA2*v5*w2 - 806*CA**3*v5*w2 +
     &   116*CA4*v5*w2 - 60*v6*w2 - 452*CA*v6*w2 + 46*CA2*v6*w2 +
     &   252*CA**3*v6*w2 - 62*CA4*v6*w2 + 12*v7*w2 + 324*CA*v7*w2 -
     &   4*CA2*v7*w2 + 356*CA**3*v7*w2 + 20*CA4*v7*w2 - 96*CA*v8*w2 -
     &   352*CA**3*v8*w2 + 96*CA**3*v9*w2 - 2*v*w3 - 2*CA*v*w3 +
     &   4*CA2*v*w3 + 2*CA**3*v*w3 - 2*CA4*v*w3 + 13*v2*w3 +
     &   14*CA*v2*w3 - 34*CA2*v2*w3 - 18*CA**3*v2*w3 + 21*CA4*v2*w3 -
     &   28*v3*w3 - 46*CA*v3*w3 + 64*CA2*v3*w3 + 76*CA**3*v3*w3 - 
     &   28*CA4*v3*w3 + 30*v4*w3 + 176*CA*v4*w3 - 72*CA2*v4*w3 - 
     &   144*CA**3*v4*w3 + 94*CA4*v4*w3 + 12*v5*w3 - 476*CA*v5*w3 -
     &   30*CA2*v5*w3 - 150*CA**3*v5*w3 - 70*CA4*v5*w3 - 93*v6*w3 +
     &   626*CA*v6*w3 + 108*CA2*v6*w3 + 918*CA**3*v6*w3 - 
     &   39*CA4*v6*w3 + 104*v7*w3 - 292*CA*v7*w3 - 40*CA2*v7*w3 -
     &   1172*CA**3*v7*w3 + 48*CA4*v7*w3 - 36*v8*w3 - 96*CA*v8*w3 +
     &   8*CA2*v8*w3 + 520*CA**3*v8*w3 - 32*CA4*v8*w3 + 96*CA*v9*w3 +
     &   32*CA**3*v9*w3 -64*CA**3*v10*w3 - v2*w4 - CA*v2*w4 +
     &   2*CA2*v2*w4 + CA**3*v2*w4 - CA4*v2*w4 - 14*v3*w4 -
     &   9*CA*v3*w4 + 78*CA2*v3*w4 + 7*CA**3*v3*w4 - 64*CA4*v3*w4 +
     &   70*v4*w4 + 15*CA*v4*w4 - 248*CA2*v4*w4 - 71*CA**3*v4*w4 +
     &   66*CA4*v4*w4 - 172*v5*w4 + 65*CA*v5*w4 + 460*CA2*v5*w4 +
     &   365*CA**3*v5*w4 - 136*CA4*v5*w4 + 277*v6*w4 + 6*CA*v6*w4 -
     &   452*CA2*v6*w4 - 826*CA**3*v6*w4 + 307*CA4*v6*w4 - 212*v7*w4 -
     &   456*CA*v7*w4 + 144*CA2*v7*w4 + 664*CA**3*v7*w4 -
     &   212*CA4*v7*w4 + 24*v8*w4 + 632*CA*v8*w4 - 28*CA2*v8*w4 +
     &   120*CA**3*v8*w4 + 100*CA4*v8*w4 + 28*v9*w4 - 220*CA*v9*w4 -
     &   4*CA2*v9*w4 - 412*CA**3*v9*w4 + 12*CA4*v9*w4 -
     &   32*CA*v10*w4 + 136*CA**3*v10*w4 + 16*CA**3*v11*w4 + 4*v3*w5 +
     &   4*CA*v3*w5 - 8*CA2*v3*w5 - 4*CA**3*v3*w5 + 4*CA4*v3*w5 -
     &   9*v4*w5 - 14*CA*v4*w5 - 2*CA2*v4*w5 + 22*CA**3*v4*w5 +
     &   11*CA4*v4*w5 + 10*v5*w5 + 4*CA*v5*w5 + 50*CA2*v5*w5 -
     &   80*CA**3*v5*w5 - 36*v6*w5 - 152*CA*v6*w5 - 148*CA2*v6*w5 +
     &   144*CA**3*v6*w5 - 180*CA4*v6*w5 - 30*v7*w5 + 560*CA*v7*w5 +
     &   314*CA2*v7*w5 + 124*CA**3*v7*w5 + 100*CA4*v7*w5 + 157*v8*w5 -
     &   554*CA*v8*w5 - 136*CA2*v8*w5 - 574*CA**3*v8*w5 -
     &   65*CA4*v8*w5 - 90*v9*w5 + 12*CA*v9*w5 + 42*CA2*v9*w5 +
     &   500*CA**3*v9*w5 - 48*CA4*v9*w5 - 6*v10*w5 + 140*CA*v10*w5 -
     &   76*CA**3*v10*w5 - 6*CA4*v10*w5 - 56*CA**3*v11*w5 - v4*w6 -
     &        CA*v4*w6 + 2*CA2*v4*w6 + CA**3*v4*w6 -
     &        CA4*v4*w6 + 14*v5*w6 + 13*CA*v5*w6 -
     &        50*CA2*v5*w6 - 15*CA**3*v5*w6 + 36*CA4*v5*w6 -
     &        26*v6*w6 + 29*CA*v6*w6 + 102*CA2*v6*w6 +
     &        53*CA**3*v6*w6 + 56*CA4*v6*w6 + 96*v7*w6 -
     &        175*CA*v7*w6 - 180*CA2*v7*w6 - 203*CA**3*v7*w6 +
     &        12*CA4*v7*w6 - 177*v8*w6 + 60*CA*v8*w6 -
     &        2*CA2*v8*w6 + 366*CA**3*v8*w6 + 23*CA4*v8*w6 +
     &        66*v9*w6 + 286*CA*v9*w6 + 2*CA2*v9*w6 -
     &        242*CA**3*v9*w6 + 48*CA4*v9*w6 + 28*v10*w6 -
     &        204*CA*v10*w6 - 10*CA2*v10*w6 - 48*CA**3*v10*w6 +
     &        26*CA4*v10*w6 - 8*CA*v11*w6 + 88*CA**3*v11*w6 -
     &        2*v5*w7 - 2*CA*v5*w7 + 4*CA2*v5*w7 +
     &        2*CA**3*v5*w7 - 2*CA4*v5*w7 - v6*w7 +
     &        2*CA*v6*w7 + 18*CA2*v6*w7 - 6*CA**3*v6*w7 -
     &        17*CA4*v6*w7 - 16*v7*w7 + 10*CA*v7*w7 -
     &        28*CA2*v7*w7 + 28*CA**3*v7*w7 - 20*CA4*v7*w7 +
     &        46*v8*w7 + 88*CA*v8*w7 + 112*CA2*v8*w7 -
     &        56*CA**3*v8*w7 - 2*CA4*v8*w7 + 16*v9*w7 -
     &        240*CA*v9*w7 - 14*CA2*v9*w7 + 18*CA**3*v9*w7 -
     &        22*CA4*v9*w7 - 43*v10*w7 + 118*CA*v10*w7 +
     &        12*CA2*v10*w7 + 98*CA**3*v10*w7 -
     &        41*CA4*v10*w7 + 24*CA*v11*w7 - 84*CA**3*v11*w7 +
     &        v6*w8 + CA*v6*w8 - 2*CA2*v6*w8 -
     &        CA**3*v6*w8 + CA4*v6*w8 - 2*v7*w8 -
     &        3*CA*v7*w8 + 2*CA2*v7*w8 + 5*CA**3*v7*w8 +
     &        2*v8*w8 - 23*CA*v8*w8 - 20*CA2*v8*w8 -
     &        9*CA**3*v8*w8 - 6*CA4*v8*w8 - 28*v9*w8 +
     &        67*CA*v9*w8 - 36*CA2*v9*w8 + 15*CA**3*v9*w8 +
     &        8*CA4*v9*w8 + 27*v10*w8 - 14*CA*v10*w8 -
     &        8*CA2*v10*w8 - 62*CA**3*v10*w8 +
     &        37*CA4*v10*w8 - 28*CA*v11*w8 + 52*CA**3*v11*w8 +
     &        6*v9*w9 - 4*CA*v9*w9 + 18*CA2*v9*w9 -
     &        6*v10*w9 - 12*CA*v10*w9 + 14*CA2*v10*w9 +
     &        20*CA**3*v10*w9 - 24*CA4*v10*w9 + 16*CA*v11*w9 -
     &        20*CA**3*v11*w9 + 4*CA*v10*w10 - 8*CA2*v10*w10 -
     &        4*CA**3*v10*w10 + 8*CA4*v10*w10 - 4*CA*v11*w10 +
     &        4*CA**3*v11*w10))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**4*(1 - v + v*w)**2)
      part13 = -(2*CF*l1w*(4*CA - 12*CA**3 - 20*CA*v + 76*CA**3*v + 
     &   44*CA*v2 - 220*CA**3*v2 - 52*CA*v3 + 372*CA**3*v3 + 
     &   32*CA*v4 - 392*CA**3*v4 - 8*CA*v5 + 256*CA**3*v5 - 
     &   96*CA**3*v6 + 16*CA**3*v7 - 5*w - 2*CA*w + 16*CA2*w +
     &   2*CA**3*w - 11*CA4*w + 24*v*w + 16*CA*v*w - 62*CA2*v*w - 
     &   8*CA**3*v*w + 46*CA4*v*w - 50*v2*w - 60*CA*v2*w + 
     &   104*CA2*v2*w - 8*CA**3*v2*w - 92*CA4*v2*w + 62*v3*w + 
     &   140*CA*v3*w - 88*CA2*v3*w + 148*CA**3*v3*w + 110*CA4*v3*w -
     &   51*v4*w - 214*CA*v4*w + 28*CA2*v4*w - 482*CA**3*v4*w - 
     &   81*CA4*v4*w + 26*v5*w + 212*CA*v5*w + 6*CA2*v5*w + 
     &   772*CA**3*v5*w + 36*CA4*v5*w - 6*v6*w - 124*CA*v6*w -
     &   4*CA2*v6*w - 680*CA**3*v6*w - 8*CA4*v6*w + 32*CA*v7*w +
     &   320*CA**3*v7*w - 64*CA**3*v8*w + w2 + CA*w2 - 2*CA2*w2 - 
     &   CA**3*w2 + CA4*w2 + 4*v*w2 - CA*v*w2 - 28*CA2*v*w2 + 
     &   3*CA**3*v*w2 + 20*CA4*v*w2 - 46*v2*w2 - 17*CA*v2*w2 +
     &   120*CA2*v2*w2 + 23*CA**3*v2*w2 - 50*CA4*v2*w2 + 136*v3*w2 +
     &   59*CA*v3*w2 - 222*CA2*v3*w2 - 177*CA**3*v3*w2 + 60*CA4*v3*w2 -
     &   219*v4*w2 - 136*CA*v4*w2 + 240*CA2*v4*w2 + 490*CA**3*v4*w2 -
     &   75*CA4*v4*w2 + 224*v5*w2 + 282*CA*v5*w2 - 114*CA2*v5*w2 -
     &   582*CA**3*v5*w2 + 64*CA4*v5*w2 - 136*v6*w2 - 400*CA*v6*w2 -
     &   6*CA2*v6*w2 + 128*CA**3*v6*w2 - 36*CA4*v6*w2 + 36*v7*w2 +
     &   308*CA*v7*w2 + 12*CA2*v7*w2 + 356*CA**3*v7*w2 + 16*CA4*v7*w2 -
     &   96*CA*v8*w2 - 336*CA**3*v8*w2 + 96*CA**3*v9*w2 - 2*v*w3 -
     &   2*CA*v*w3 + 4*CA2*v*w3 + 2*CA**3*v*w3 - 2*CA4*v*w3 + 
     &   23*v2*w3 + 14*CA*v2*w3 - 28*CA2*v2*w3 - 18*CA**3*v2*w3 +
     &   17*CA4*v2*w3 - 60*v3*w3 - 34*CA*v3*w3 + 52*CA2*v3*w3 +
     &   72*CA**3*v3*w3 - 20*CA4*v3*w3 + 66*v4*w3 + 104*CA*v4*w3 -
     &   76*CA2*v4*w3 - 124*CA**3*v4*w3 + 64*CA4*v4*w3 - 8*v5*w3 -
     &   328*CA*v5*w3 - 24*CA2*v5*w3 - 126*CA**3*v5*w3 - 28*CA4*v5*w3 -
     &   103*v6*w3 + 498*CA*v6*w3 + 88*CA2*v6*w3 + 730*CA**3*v6*w3 - 
     &   43*CA4*v6*w3 + 144*v7*w3 - 252*CA*v7*w3 + 8*CA2*v7*w3 - 
     &   912*CA**3*v7*w3 + 36*CA4*v7*w3 - 60*v8*w3 - 96*CA*v8*w3 -
     &   16*CA2*v8*w3 + 392*CA**3*v8*w3 - 32*CA4*v8*w3 + 96*CA*v9*w3 +
     &   48*CA**3*v9*w3 - 64*CA**3*v10*w3 - v2*w4 - CA*v2*w4 +
     &   2*CA2*v2*w4 + CA**3*v2*w4 - CA4*v2*w4 - 32*v3*w4 -
     &   9*CA*v3*w4 + 68*CA2*v3*w4 + 7*CA**3*v3*w4 - 48*CA4*v3*w4 +
     &   146*v4*w4 + 27*CA*v4*w4 - 186*CA2*v4*w4 - 71*CA**3*v4*w4 +
     &   52*CA4*v4*w4 - 284*v5*w4 + 41*CA*v5*w4 + 374*CA2*v5*w4 +
     &   317*CA**3*v5*w4 - 140*CA4*v5*w4 + 363*v6*w4 - 22*CA*v6*w4 -
     &   370*CA2*v6*w4 - 638*CA**3*v6*w4 + 277*CA4*v6*w4 - 268*v7*w4 -
     &   360*CA*v7*w4 + 64*CA2*v7*w4 + 448*CA**3*v7*w4 - 
     &   184*CA4*v7*w4 + 40*v8*w4 + 560*CA*v8*w4 - 20*CA2*v8*w4 +
     &   156*CA**3*v8*w4 + 100*CA4*v8*w4 + 36*v9*w4 - 204*CA*v9*w4 +
     &   12*CA2*v9*w4 - 348*CA**3*v9*w4 + 16*CA4*v9*w4 -
     &   32*CA*v10*w4 + 112*CA**3*v10*w4 + 16*CA**3*v11*w4 +
     &   4*v3*w5 + 4*CA*v3*w5 - 8*CA2*v3*w5 - 4*CA**3*v3*w5 + 
     &   4*CA4*v3*w5 - 7*v4*w5 - 14*CA*v4*w5 + 4*CA2*v4*w5 + 
     &   22*CA**3*v4*w5 + 7*CA4*v4*w5 - 30*v5*w5 - 8*CA*v5*w5 +
     &   10*CA2*v5*w5 - 68*CA**3*v5*w5 + 4*CA4*v5*w5 + 56*v6*w5 -
     &   80*CA*v6*w5 - 76*CA2*v6*w5 + 84*CA**3*v6*w5 - 130*CA4*v6*w5 -
     &   94*v7*w5 + 432*CA*v7*w5 + 276*CA2*v7*w5 + 152*CA**3*v7*w5 +
     &        50*CA4*v7*w5 + 175*v8*w5 - 474*CA*v8*w5 -
     &        76*CA2*v8*w5 - 446*CA**3*v8*w5 - 55*CA4*v8*w5 -
     &        98*v9*w5 + 8*CA*v9*w5 + 10*CA2*v9*w5 +
     &        348*CA**3*v9*w5 - 56*CA4*v9*w5 - 6*v10*w5 +
     &        132*CA*v10*w5 - 4*CA2*v10*w5 - 40*CA**3*v10*w5 -
     &        8*CA4*v10*w5 - 48*CA**3*v11*w5 - v4*w6 -
     &        CA*v4*w6 + 2*CA2*v4*w6 + CA**3*v4*w6 -
     &        CA4*v4*w6 + 28*v5*w6 + 13*CA*v5*w6 -
     &        52*CA2*v5*w6 - 15*CA**3*v5*w6 + 28*CA4*v5*w6 -
     &        50*v6*w6 + 17*CA*v6*w6 + 92*CA2*v6*w6 +
     &        61*CA**3*v6*w6 + 38*CA4*v6*w6 + 80*v7*w6 -
     &        151*CA*v7*w6 - 186*CA2*v7*w6 - 179*CA**3*v7*w6 +
     &        28*CA4*v7*w6 - 133*v8*w6 + 68*CA*v8*w6 -
     &        20*CA2*v8*w6 + 250*CA**3*v8*w6 + 35*CA4*v8*w6 +
     &        52*v9*w6 + 246*CA*v9*w6 + 10*CA2*v9*w6 -
     &        138*CA**3*v9*w6 + 40*CA4*v9*w6 + 24*v10*w6 -
     &        184*CA*v10*w6 + 2*CA2*v10*w6 - 44*CA**3*v10*w6 +
     &        32*CA4*v10*w6 - 8*CA*v11*w6 + 64*CA**3*v11*w6 -
     &        2*v5*w7 - 2*CA*v5*w7 + 4*CA2*v5*w7 +
     &        2*CA**3*v5*w7 - 2*CA4*v5*w7 - 11*v6*w7 +
     &        2*CA*v6*w7 + 12*CA2*v6*w7 - 6*CA**3*v6*w7 -
     &        13*CA4*v6*w7 + 16*v7*w7 + 14*CA*v7*w7 -
     &        16*CA2*v7*w7 + 16*CA**3*v7*w7 - 12*CA4*v7*w7 +
     &        10*v8*w7 + 64*CA*v8*w7 + 100*CA2*v8*w7 -
     &        12*CA**3*v8*w7 - 28*CA4*v8*w7 + 16*v9*w7 -
     &        204*CA*v9*w7 - 12*CA2*v9*w7 - 6*CA**3*v9*w7 -
     &        4*CA4*v9*w7 - 29*v10*w7 + 102*CA*v10*w7 +
     &        62*CA**3*v10*w7 - 45*CA4*v10*w7 + 24*CA*v11*w7 -
     &        56*CA**3*v11*w7 + v6*w8 + CA*v6*w8 -
     &        2*CA2*v6*w8 - CA**3*v6*w8 + CA4*v6*w8 -
     &        3*CA*v7*w8 + 12*CA2*v7*w8 + 5*CA**3*v7*w8 -
     &        2*v8*w8 - 19*CA*v8*w8 - 18*CA2*v8*w8 -
     &        13*CA**3*v8*w8 - 12*v9*w8 + 59*CA*v9*w8 -
     &        30*CA2*v9*w8 + 7*CA**3*v9*w8 + 4*CA4*v9*w8 +
     &        13*v10*w8 - 10*CA*v10*w8 - 2*CA2*v10*w8 -
     &        34*CA**3*v10*w8 + 35*CA4*v10*w8 - 28*CA*v11*w8 +
     &        36*CA**3*v11*w8 - 4*CA2*v8*w9 + 2*v9*w9 -
     &        4*CA*v9*w9 + 16*CA2*v9*w9 + 4*CA**3*v9*w9 -
     &        2*CA4*v9*w9 - 2*v10*w9 - 12*CA*v10*w9 +
     &        12*CA2*v10*w9 + 12*CA**3*v10*w9 -
     &        22*CA4*v10*w9 + 16*CA*v11*w9 - 16*CA**3*v11*w9 +
     &        4*CA*v10*w10 - 8*CA2*v10*w10 - 4*CA**3*v10*w10 +
     &        8*CA4*v10*w10 - 4*CA*v11*w10 + 4*CA**3*v11*w10))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**4*(1 - v + v*w)**2)
      struv13= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + 
     &         part9 + part10 + part11 + part12 +
     &         part13 
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION STRUV14(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 = (4*CF*lvw*(16*CA4 - 32*CA4*v + 48*CA4*v2 - 
     &     32*CA4*v3 +
     &        16*CA4*v4 - CA2*v*w - 20*CA4*v*w + v2*w +
     &        2*CA2*v2*w - 15*CA4*v2*w - 3*v3*w + 
     &     21*CA4*v3*w +
     &        4*v4*w - 3*CA2*v4*w - 40*CA4*v4*w - 2*v5*w +
     &        2*CA2*v5*w + 6*CA4*v5*w - 6*CA2*v2*w2 +
     &        40*CA4*v2*w2 + 9*CA2*v3*w2 - 
     &     29*CA4*v3*w2 -
     &        3*v4*w2 - 18*CA2*v4*w2 + 72*CA4*v4*w2 +
     &        3*v5*w2 + 15*CA2*v5*w2 - 27*CA4*v5*w2 -
     &        8*CA2*v6*w2 + 8*CA4*v6*w2 - 
     &     15*CA4*v3*w3 +
     &        5*v4*w3 + 25*CA2*v4*w3 - 42*CA4*v4*w3 -
     &        5*v5*w3 - 25*CA2*v5*w3 + 25*CA4*v5*w3 +
     &        16*CA2*v6*w3 - 16*CA4*v6*w3 - 2*v4*w4 -
     &        10*CA2*v4*w4 + 16*CA4*v4*w4 + 2*v5*w4 +
     &        10*CA2*v5*w4 - 8*CA4*v5*w4 -
     &      12*CA2*v6*w4 +
     &        12*CA4*v6*w4 + 4*CA2*v6*w5 - 
     &     4*CA4*v6*w5))/
     &    (CA*(1 - v)**2*v3*w2)
      part2 = -(4*CF*l1v*(8*CA4 - 32*CA4*v + 64*CA4*v2 - 
     &     80*CA4*v3 +
     &        64*CA4*v4 - 32*CA4*v5 + 8*CA4*v6 + 
     &     8*CA4*w -
     &        34*CA4*v*w + 51*CA4*v2*w - 39*CA4*v3*w
     &      + v4*w -
     &        CA2*v4*w + 18*CA4*v4*w - 2*v5*w + 
     &     2*CA2*v5*w -
     &        9*CA4*v5*w + v6*w - CA2*v6*w + 
     &     13*CA4*v6*w -
     &        8*CA4*v7*w + 2*CA4*v3*w2 - 2*v4*w2 +
     &        2*CA2*v4*w2 + 3*CA4*v4*w2 + 4*v5*w2 -
     &        2*CA2*v5*w2 - 22*CA4*v5*w2 - v6*w2 -
     &        CA2*v6*w2 + 14*CA4*v6*w2 - v7*w2 +
     &        CA2*v7*w2 - 13*CA4*v7*w2 +
     &      8*CA4*v8*w2 -
     &        8*CA4*v2*w3 + 26*CA4*v3*w3 + 2*v4*w3 -
     &        2*CA2*v4*w3 - 41*CA4*v4*w3 - 4*v5*w3 +
     &        56*CA4*v5*w3 - v6*w3 + 3*CA2*v6*w3 -
     &        38*CA4*v6*w3 + 3*v7*w3 - CA2*v7*w3 +
     &        37*CA4*v7*w3 - 16*CA4*v8*w3 + 
     &     4*CA2*v5*w4 -
     &        8*CA4*v5*w4 + 4*v6*w4 - 4*CA2*v6*w4 +
     &        10*CA4*v6*w4 - 4*v7*w4 - 26*CA4*v7*w4 +
     &        12*CA4*v8*w4 - 2*v6*w5 - 2*CA2*v6*w5 +
     &        2*v7*w5 + 2*CA2*v7*w5 + 8*CA4*v7*w5 -
     &        4*CA4*v8*w5))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v*w)*(1 - v + v*w))
      part3 = -(4*CF*l1vw*(CA2 + 4*CA4 - 3*CA2*v - 12*CA4*v +
     &        3*CA2*v2 + 12*CA4*v2 - 3*CA2*v4 - 
     &     12*CA4*v4 +
     &        3*CA2*v5 + 12*CA4*v5 - CA2*v6 - 
     &     4*CA4*v6 +
     &        v*w + 3*CA2*v*w + 5*CA4*v*w - 2*v2*w -
     &      5*CA2*v2*w -
     &        15*CA4*v2*w + v3*w - 8*CA2*v3*w + 
     &     3*CA4*v3*w +
     &        v4*w + 30*CA2*v4*w + 29*CA4*v4*w - 2*v5*w -
     &        31*CA2*v5*w - 32*CA4*v5*w + v6*w + 
     &     11*CA2*v6*w +
     &        10*CA4*v6*w + 2*v2*w2 + 3*CA2*v2*w2 +
     &        7*CA4*v2*w2 - v3*w2 + 3*CA2*v3*w2 -
     &        2*CA4*v3*w2 - 8*v4*w2 - 47*CA2*v4*w2 -
     &        33*CA4*v4*w2 + 13*v5*w2 + 71*CA2*v5*w2 +
     &        34*CA4*v5*w2 - 6*v6*w2 - 30*CA2*v6*w2 -
     &        6*CA4*v6*w2 + v3*w3 + 5*CA2*v3*w3 -
     &        2*CA4*v3*w3 + 11*v4*w3 + 22*CA2*v4*w3 +
     &        23*CA4*v4*w3 - 24*v5*w3 - 69*CA2*v5*w3 -
     &        21*CA4*v5*w3 + 14*v6*w3 + 36*CA2*v6*w3 -
     &        2*CA4*v6*w3 - 4*v4*w4 - 2*CA2*v4*w4 -
     &        7*CA4*v4*w4 + 17*v5*w4 + 34*CA2*v5*w4 +
     &        8*CA4*v5*w4 - 16*v6*w4 - 23*CA2*v6*w4 +
     &        2*CA4*v6*w4 - 4*v5*w5 - 8*CA2*v5*w5 -
     &        CA4*v5*w5 + 9*v6*w5 + 9*CA2*v6*w5 -
     &        2*v6*w6 - 2*CA2*v6*w6))/
     &    (CA*(1 - v)*v2*w*(1 - v + v*w)**3) 
      part4 = (4*CF*lw*(24*CA4 - 96*CA4*v + 208*CA4*v2 -
     &      288*CA4*v3 +
     &        272*CA4*v4 - 176*CA4*v5 + 72*CA4*v6 -
     &        16*CA4*v7 - 8*CA4*w + 32*CA4*v*w + 2*v2*w -
     &        CA2*v2*w - 129*CA4*v2*w - 6*v3*w + 
     &     3*CA2*v3*w +
     &        275*CA4*v3*w + 8*v4*w - 5*CA2*v4*w -
     &        323*CA4*v4*w - 6*v5*w + 5*CA2*v5*w +
     &        225*CA4*v5*w + 2*v6*w - 2*CA2*v6*w -
     &        56*CA4*v6*w - 16*CA4*v7*w + 16*CA4*v8*w +
     &        4*CA2*v2*w2 + 29*CA4*v2*w2 - 
     &     12*CA2*v3*w2 -
     &        87*CA4*v3*w2 - v4*w2 + 19*CA2*v4*w2 +
     &        62*CA4*v4*w2 + 2*v5*w2 - 18*CA2*v5*w2 +
     &        21*CA4*v5*w2 + v6*w2 + 7*CA2*v6*w2 -
     &        169*CA4*v6*w2 - 2*v7*w2 + 144*CA4*v7*w2 -
     &        56*CA4*v8*w2 + 2*CA4*v2*w3 - 
     &     6*CA4*v3*w3 +
     &        4*v4*w3 - 7*CA2*v4*w3 + 85*CA4*v4*w3 -
     &        8*v5*w3 + 14*CA2*v5*w3 - 160*CA4*v5*w3 -
     &        v6*w3 - 7*CA2*v6*w3 + 299*CA4*v6*w3 +
     &        5*v7*w3 - 220*CA4*v7*w3 + 80*CA4*v8*w3 -
     &        2*v4*w4 + 4*CA2*v4*w4 - 41*CA4*v4*w4 +
     &        4*v5*w4 - 8*CA2*v5*w4 + 82*CA4*v5*w4 +
     &        7*v6*w4 + 4*CA2*v6*w4 - 194*CA4*v6*w4 -
     &        9*v7*w4 + 153*CA4*v7*w4 - 60*CA4*v8*w4 +
     &        6*CA4*v4*w5 - 12*CA4*v5*w5 - 8*v6*w5 +
     &        57*CA4*v6*w5 + 8*v7*w5 - 51*CA4*v7*w5 +
     &        24*CA4*v8*w5 + 2*v6*w6 - 6*CA4*v6*w6 -
     &        2*v7*w6 + 6*CA4*v7*w6 - 4*CA4*v8*w6))/
     &    (CA*(1 - v)**2*v3*(1 - w)*w2*(1 - v*w)*(1 - v + v*w))
      part5 = -(2*CF*lms*(32*CA4 - 80*CA4*v + 144*CA4*v2 - 
     &     144*CA4*v3 +
     &        112*CA4*v4 - 48*CA4*v5 + 16*CA4*v6 - 
     &     4*CA2*w +
     &        4*CA4*w + 12*CA2*v*w - 108*CA4*v*w + 3*v2*w -
     &        20*CA2*v2*w + 177*CA4*v2*w - 7*v3*w +
     &        22*CA2*v3*w - 319*CA4*v3*w + 7*v4*w -
     &        16*CA2*v4*w + 297*CA4*v4*w - 3*v5*w +
     &        6*CA2*v5*w - 275*CA4*v5*w + 128*CA4*v6*w -
     &        48*CA4*v7*w + 4*CA2*w2 + 4*CA4*w2 -
     &        24*CA4*v*w2 - 2*v2*w2 - 18*CA2*v2*w2 +
     &        164*CA4*v2*w2 - 4*v3*w2 + 44*CA2*v3*w2 -
     &        96*CA4*v3*w2 + 13*v4*w2 - 66*CA2*v4*w2 +
     &        213*CA4*v4*w2 - 14*v5*w2 + 54*CA2*v5*w2 -
     &        128*CA4*v5*w2 + 7*v6*w2 - 26*CA2*v6*w2 +
     &        195*CA4*v6*w2 - 96*CA4*v7*w2 +
     &      48*CA4*v8*w2 -
     &        12*CA2*v*w3 - 12*CA4*v*w3 + 24*CA2*v2*w3 +
     &        48*CA4*v2*w3 + 6*v3*w3 - 18*CA2*v3*w3 -
     &        164*CA4*v3*w3 - 3*v4*w3 + 2*CA2*v4*w3 -
     &        23*CA4*v4*w3 - 6*v5*w3 + 44*CA2*v5*w3 -
     &        54*CA4*v5*w3 + 10*v6*w3 - 42*CA2*v6*w3 -
     &        80*CA4*v6*w3 - 7*v7*w3 + 34*CA2*v7*w3 -
     &        11*CA4*v7*w3 - 16*CA4*v9*w3 + 
     &     12*CA2*v2*w4 +
     &        12*CA4*v2*w4 - 32*CA2*v3*w4 - 
     &     40*CA4*v3*w4 -
     &        6*v4*w4 + 42*CA2*v4*w4 + 108*CA4*v4*w4 +
     &        8*v5*w4 - 56*CA2*v5*w4 + 24*CA4*v5*w4 -
     &        7*v6*w4 + 2*CA2*v6*w4 + 77*CA4*v6*w4 +
     &        4*v7*w4 + 2*CA2*v7*w4 + 26*CA4*v7*w4 +
     &        v8*w4 - 22*CA2*v8*w4 - 3*CA4*v8*w4 +
     &        16*CA4*v9*w4 - 4*CA2*v3*w5 - 
     &     4*CA4*v3*w5 +
     &        12*CA2*v4*w5 + 12*CA4*v4*w5 + 2*v5*w5 -
     &        18*CA2*v5*w5 - 32*CA4*v5*w5 - 2*v6*w5 +
     &        38*CA2*v6*w5 - 20*CA4*v6*w5 + v7*w5 -
     &        2*CA2*v7*w5 - 55*CA4*v7*w5 - v8*w5 +
     &        10*CA2*v8*w5 + 15*CA4*v8*w5 +
     &      8*CA2*v9*w5 -
     &        16*CA4*v9*w5 - 8*CA2*v7*w6 + 
     &     16*CA4*v7*w6 -
     &        4*CA2*v8*w6 + 4*CA4*v8*w6 -
     &      8*CA2*v9*w6 +
     &        8*CA4*v9*w6 + 4*CA2*v9*w7 - 
     &     4*CA4*v9*w7))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v*w)**3)
      part6 = -(2*CF*lmss*(16*CA4 - 112*CA4*v + 368*CA4*v2 -
     &        752*CA4*v3 + 1056*CA4*v4 - 1056*CA4*v5 +
     &        752*CA4*v6 - 368*CA4*v7 + 112*CA4*v8 -
     &        16*CA4*v9 + 16*CA4*w - 48*CA4*v*w + v2*w -
     &        2*CA2*v2*w - 79*CA4*v2*w - 4*v3*w + 
     &     8*CA2*v3*w +
     &        684*CA4*v3*w + 7*v4*w - 14*CA2*v4*w -
     &        1721*CA4*v4*w - 8*v5*w + 16*CA2*v5*w +
     &        2488*CA4*v5*w + 7*v6*w - 14*CA2*v6*w -
     &        2313*CA4*v6*w - 4*v7*w + 8*CA2*v7*w +
     &        1388*CA4*v7*w + v8*w - 2*CA2*v8*w -
     &        495*CA4*v8*w + 80*CA4*v9*w + 48*CA4*v*w2 -
     &        192*CA4*v2*w2 + 2*v3*w2 - 6*CA2*v3*w2 +
     &        148*CA4*v3*w2 - 8*v4*w2 + 10*CA2*v4*w2 +
     &        638*CA4*v4*w2 + 16*v5*w2 + 8*CA2*v5*w2 -
     &        2008*CA4*v5*w2 - 20*v6*w2 - 
     &     32*CA2*v6*w2 +
     &        2764*CA4*v6*w2 + 14*v7*w2 + 
     &     38*CA2*v7*w2 -
     &        2156*CA4*v7*w2 - 4*v8*w2 - 26*CA2*v8*w2 +
     &        934*CA4*v8*w2 + 8*CA2*v9*w2 -
     &      176*CA4*v9*w2 +
     &        48*CA4*v2*w3 - 176*CA4*v3*w3 + 3*v4*w3 +
     &        8*CA2*v4*w3 + 85*CA4*v4*w3 - 12*v5*w3 -
     &        60*CA2*v5*w3 + 664*CA4*v5*w3 + 22*v6*w3 +
     &        136*CA2*v6*w3 - 1654*CA4*v6*w3 
     &     - 20*v7*w3 -
     &        156*CA2*v7*w3 + 1792*CA4*v7*w3 +
     &      7*v8*w3 +
     &        104*CA2*v8*w3 - 983*CA4*v8*w3 -
     &        32*CA2*v9*w3 + 224*CA4*v9*w3 +
     &      16*CA4*v3*w4 -
     &        48*CA4*v4*w4 + 2*v5*w4 + 26*CA2*v5*w4 -
     &        76*CA4*v5*w4 - 10*v6*w4 - 110*CA2*v6*w4 +
     &        512*CA4*v6*w4 + 14*v7*w4 + 
     &     174*CA2*v7*w4 -
     &        844*CA4*v7*w4 - 6*v8*w4 - 142*CA2*v8*w4 +
     &        620*CA4*v8*w4 + 52*CA2*v9*w4 -
     &        180*CA4*v9*w4 + 2*v6*w5 + 22*CA2*v6*w5 -
     &        64*CA4*v6*w5 - 4*v7*w5 - 72*CA2*v7*w5 +
     &        204*CA4*v7*w5 + 2*v8*w5 + 86*CA2*v8*w5 -
     &        224*CA4*v8*w5 - 44*CA2*v9*w5 +
     &      92*CA4*v9*w5 +
     &        8*CA2*v7*w6 - 16*CA4*v7*w6 -
     &      20*CA2*v8*w6 +
     &        36*CA4*v8*w6 + 20*CA2*v9*w6 -
     &      28*CA4*v9*w6 -
     &        4*CA2*v9*w7 + 4*CA4*v9*w7))/
     &    (CA*(1 - v)**2*v3*w2*(1 - v + v*w)**3)
      part7 = -(CF*(12*CA2 - 28*CA4 - 72*CA2*v + 168*CA4*v - 6*v2 +
     &        208*CA2*v2 - 474*CA4*v2 + 8*CA**3*CF*v2 + 28*v3 -
     &        376*CA2*v3 + 828*CA4*v3 - 40*CA**3*CF*v3 - 54*v4 +
     &        452*CA2*v4 - 974*CA4*v4 - 4*CA*CF*v4 +
     &        84*CA**3*CF*v4 + 56*v5 - 360*CA2*v5 + 784*CA4*v5 +
     &        16*CA*CF*v5 - 96*CA**3*CF*v5 - 34*v6 + 184*CA2*v6 -
     &        422*CA4*v6 - 24*CA*CF*v6 + 64*CA**3*CF*v6 + 12*v7 -
     &        56*CA2*v7 + 140*CA4*v7 + 16*CA*CF*v7 -
     &        24*CA**3*CF*v7 - 2*v8 + 8*CA2*v8 - 22*CA4*v8 -
     &        4*CA*CF*v8 + 4*CA**3*CF*v8 + 16*CA2*w + 16*CA4*w -
     &        96*CA2*v*w - 96*CA4*v*w - 4*v2*w + 254*CA2*v2*w +
     &        146*CA4*v2*w - 8*CA**3*CF*v2*w + 22*v3*w -
     &        390*CA2*v3*w + 148*CA4*v3*w + 40*CA**3*CF*v3*w -
     &        75*v4*w + 420*CA2*v4*w - 781*CA4*v4*w +
     &        4*CA*CF*v4*w - 92*CA**3*CF*v4*w + 158*v5*w -
     &        384*CA2*v5*w + 1178*CA4*v5*w - 16*CA*CF*v5*w +
     &        128*CA**3*CF*v5*w - 196*v6*w + 302*CA2*v6*w -
     &        990*CA4*v6*w + 12*CA*CF*v6*w - 100*CA**3*CF*v6*w +
     &        144*v7*w - 174*CA2*v7*w + 538*CA4*v7*w +
     &        20*CA*CF*v7*w + 20*CA**3*CF*v7*w - 61*v8*w +
     &        64*CA2*v8*w - 199*CA4*v8*w - 32*CA*CF*v8*w +
     &        24*CA**3*CF*v8*w + 12*v9*w - 12*CA2*v9*w +
     &        40*CA4*v9*w + 12*CA*CF*v9*w - 12*CA**3*CF*v9*w +
     &        12*CA2*v2*w2 + 132*CA4*v2*w2 - 60*CA2*v3*w2 -
     &        660*CA4*v3*w2 + 4*v4*w2 + 70*CA2*v4*w2 +
     &        1250*CA4*v4*w2 - 80*CA**3*CF*v4*w2 - 16*v5*w2 +
     &   84*CA2*v5*w2 - 1044*CA4*v5*w2 + 320*CA**3*CF*v5*w2 +
     &   49*v6*w2 - 238*CA2*v6*w2 + 101*CA4*v6*w2 + 24*CA*CF*v6*w2 -
     &   576*CA**3*CF*v6*w2 - 93*v7*w2 + 186*CA2*v7*w2 +
     &   563*CA4*v7*w2 - 72*CA*CF*v7*w2 + 608*CA**3*CF*v7*w2 +
     &   81*v8*w2 - 64*CA2*v8*w2 - 525*CA4*v8*w2 + 60*CA*CF*v8*w2 -
     &   356*CA**3*CF*v8*w2 - 23*v9*w2 + 6*CA2*v9*w2 +
     &   241*CA4*v9*w2 + 72*CA**3*CF*v9*w2 - 2*v10*w2 + 4*CA2*v10*w2 -
     &   58*CA4*v10*w2 - 12*CA*CF*v10*w2 + 12*CA**3*CF*v10*w2 -
     &   48*CA2*v2*w3 - 48*CA4*v2*w3 + 240*CA2*v3*w3 +
     &   240*CA4*v3*w3 + 12*v4*w3 - 438*CA2*v4*w3 - 234*CA4*v4*w3 +
     &   24*CA**3*CF*v4*w3 - 44*v5*w3 + 300*CA2*v5*w3 - 496*CA4*v5*w3 -
     &   96*CA**3*CF*v5*w3 - 47*v6*w3 - 28*CA2*v6*w3 + 1423*CA4*v6*w3 -
     &   12*CA*CF*v6*w3 + 100*CA**3*CF*v6*w3 + 295*v7*w3 +
     &   32*CA2*v7*w3 - 1515*CA4*v7*w3 + 36*CA*CF*v7*w3 +
     &   36*CA**3*CF*v7*w3 - 368*v8*w3 - 114*CA2*v8*w3 +
     &   790*CA4*v8*w3 - 232*CA**3*CF*v8*w3 + 189*v9*w3 +
     &   104*CA2*v9*w3 - 169*CA4*v9*w3 - 60*CA*CF*v9*w3 + 
     &   292*CA**3*CF*v9*w3 - 45*v10*w3 - 56*CA2*v10*w3 - 
     &   23*CA4*v10*w3 + 32*CA*CF*v10*w3 - 120*CA**3*CF*v10*w3 +
     &   8*v11*w3 + 8*CA2*v11*w3 + 32*CA4*v11*w3 + 4*CA*CF*v11*w3 -
     &   4*CA**3*CF*v11*w3 - 60*CA2*v4*w4 - 180*CA4*v4*w4 +
     &   240*CA2*v5*w4 + 720*CA4*v5*w4 + 102*v6*w4 - 244*CA2*v6*w4 -
     &   1018*CA4*v6*w4 + 136*CA**3*CF*v6*w4 - 320*v7*w4 -
     &   88*CA2*v7*w4 + 528*CA4*v7*w4 - 408*CA**3*CF*v7*w4 +
     &   408*v8*w4 + 250*CA2*v8*w4 + 186*CA4*v8*w4 - 36*CA*CF*v8*w4 +
     &   540*CA**3*CF*v8*w4 - 272*v9*w4 - 160*CA2*v9*w4 -
     &   336*CA4*v9*w4 + 72*CA*CF*v9*w4 - 400*CA**3*CF*v9*w4 +
     &   109*v10*w4 + 28*CA2*v10*w4 + 155*CA4*v10*w4 - 
     &   20*CA*CF*v10*w4 + 84*CA**3*CF*v10*w4 - 27*v11*w4 + 
     &   50*CA2*v11*w4 - 71*CA4*v11*w4 - 16*CA*CF*v11*w4 +
     &   48*CA**3*CF*v11*w4 - 16*CA2*v**12*w4 + 48*CA2*v4*w5 +
     &   48*CA4*v4*w5 - 192*CA2*v5*w5 - 192*CA4*v5*w5 - 12*v6*w5 +
     &   210*CA2*v6*w5 + 126*CA4*v6*w5 - 24*CA**3*CF*v6*w5 + 42*v7*w5 +
     &   34*CA2*v7*w5 + 296*CA4*v7*w5 + 72*CA**3*CF*v7*w5 - 101*v8*w5 -
     &   92*CA2*v8*w5 - 515*CA4*v8*w5 + 12*CA*CF*v8*w5 -
     &   36*CA**3*CF*v8*w5 + 138*v9*w5 - 54*CA2*v9*w5 +
     &   264*CA4*v9*w5 - 24*CA*CF*v9*w5 - 48*CA**3*CF*v9*w5 -
     &   106*v10*w5 + 228*CA2*v10*w5 - 74*CA4*v10*w5 - 
     &   12*CA*CF*v10*w5 + 140*CA**3*CF*v10*w5 + 39*v11*w5 -
     &   246*CA2*v11*w5 + 111*CA4*v11*w5 + 24*CA*CF*v11*w5 -
     &   104*CA**3*CF*v11*w5 + 64*CA2*v**12*w5 - 16*CA4*v**12*w5 +
     &   36*CA2*v6*w6 + 76*CA4*v6*w6 - 108*CA2*v7*w6 -
     &   228*CA4*v7*w6 + 12*v8*w6 + 22*CA2*v8*w6 + 202*CA4*v8*w6 -
     &   64*CA**3*CF*v8*w6 - 36*v9*w6 + 140*CA2*v9*w6 - 16*CA4*v9*w6 +
     &   128*CA**3*CF*v9*w6 + 53*v10*w6 - 344*CA2*v10*w6 +
     &   47*CA4*v10*w6 + 16*CA*CF*v10*w6 - 144*CA**3*CF*v10*w6 -
     &   29*v11*w6 + 358*CA2*v11*w6 - 185*CA4*v11*w6 -
     &   16*CA*CF*v11*w6 + 80*CA**3*CF*v11*w6 - 96*CA2*v**12*w6 +
     &   48*CA4*v**12*w6 - 16*CA2*v6*w7 - 16*CA4*v6*w7 + 48*CA2*v7*w7 +
     &   48*CA4*v7*w7 + 4*v8*w7 - 26*CA2*v8*w7 - 38*CA4*v8*w7 +
     &   8*CA**3*CF*v8*w7 - 4*v9*w7 - 32*CA2*v9*w7 - 4*CA4*v9*w7 -
     &   16*CA**3*CF*v9*w7 - 9*v10*w7 + 164*CA2*v10*w7 -
     &   71*CA4*v10*w7 - 4*CA*CF*v10*w7 + 28*CA**3*CF*v10*w7 +
     &   9*v11*w7 - 226*CA2*v11*w7 + 169*CA4*v11*w7 + 4*CA*CF*v11*w7 -
     &   20*CA**3*CF*v11*w7 + 64*CA2*v**12*w7 - 48*CA4*v**12*w7 -
     &   24*CA2*v10*w8 + 24*CA4*v10*w8 + 64*CA2*v11*w8 -
     &        64*CA4*v11*w8 - 16*CA2*v**12*w8 +
     &        16*CA4*v**12*w8 - 8*CA2*v11*w9 + 8*CA4*v11*w9))/
     &    (CA*(1 - v)**2*v3*w*(1 - v*w)**3*(1 - v + v*w)**3)
      part8 = (2*CF*l1w*(64*CA4 - 384*CA4*v + 1120*CA4*v2 -
     &        2080*CA4*v3 + 2688*CA4*v4 - 2496*CA4*v5 +
     &        1664*CA4*v6 - 768*CA4*v7 + 224*CA4*v8 -
     &        32*CA4*v9 - 4*CA2*w + 52*CA4*w + 24*CA2*v*w -
     &        312*CA4*v*w + 4*v2*w - 70*CA2*v2*w + 864*CA4*v2*w -
     &        20*v3*w + 130*CA2*v3*w - 1460*CA4*v3*w + 46*v4*w -
     &        170*CA2*v4*w + 1906*CA4*v4*w - 64*v5*w +
     &        164*CA2*v5*w - 2296*CA4*v5*w + 56*v6*w -
     &        114*CA2*v6*w + 2524*CA4*v6*w - 28*v7*w +
     &        50*CA2*v7*w - 2236*CA4*v7*w + 6*v8*w -
     &        10*CA2*v8*w + 1390*CA4*v8*w - 528*CA4*v9*w +
     &        96*CA4*v10*w + 4*CA2*w2 + 4*CA4*w2 -
     &        24*CA2*v*w2 - 24*CA4*v*w2 - 2*v2*w2 +
     &        54*CA2*v2*w2 + 32*CA4*v2*w2 + 10*v3*w2 -
     &        50*CA2*v3*w2 + 60*CA4*v3*w2 - 20*v4*w2 -
     &        44*CA2*v4*w2 - 598*CA4*v4*w2 + 20*v5*w2 +
     &        212*CA2*v5*w2 + 1768*CA4*v5*w2 + 6*v6*w2 -
     &        336*CA2*v6*w2 - 2752*CA4*v6*w2 - 46*v7*w2 +
     &        320*CA2*v7*w2 + 2584*CA4*v7*w2 + 48*v8*w2 -
     &        182*CA2*v8*w2 - 1262*CA4*v8*w2 - 16*v9*w2 +
     &        46*CA2*v9*w2 + 28*CA4*v9*w2 +
     &        256*CA4*v10*w2 - 96*CA4*v11*w2 +
     &        24*CA2*v2*w3 - 144*CA4*v2*w3 -
     &        120*CA2*v3*w3 + 720*CA4*v3*w3 - 8*v4*w3 +
     &        312*CA2*v4*w3 - 1462*CA4*v4*w3 + 32*v5*w3 -
     &        528*CA2*v5*w3 + 1528*CA4*v5*w3 - 76*v6*w3 +
     &        550*CA2*v6*w3 - 1118*CA4*v6*w3 + 116*v7*w3 -
     &        306*CA2*v7*w3 + 1030*CA4*v7*w3 - 76*v8*w3 -
     &        4*CA2*v8*w3 - 1592*CA4*v8*w3 - 4*v9*w3 +
     &        142*CA2*v9*w3 + 1810*CA4*v9*w3 + 16*v10*w3 -
     &        70*CA2*v10*w3 - 948*CA4*v10*w3 +
     &        176*CA4*v11*w3 + 32*CA4*v**12*w3 -
     &        12*CA2*v2*w4 - 12*CA4*v2*w4 + 60*CA2*v3*w4 +
     &        60*CA4*v3*w4 + 6*v4*w4 - 102*CA2*v4*w4 -
     &        252*CA4*v4*w4 - 24*v5*w4 + 48*CA2*v5*w4 +
     &        648*CA4*v5*w4 + 43*v6*w4 + 202*CA2*v6*w4 -
     &        545*CA4*v6*w4 - 45*v7*w4 - 522*CA2*v7*w4 -
     &        381*CA4*v7*w4 - 29*v8*w4 + 584*CA2*v8*w4 +
     &        1583*CA4*v8*w4 + 105*v9*w4 - 362*CA2*v9*w4 -
     &        1895*CA4*v9*w4 - 52*v10*w4 + 54*CA2*v10*w4 +
     &        722*CA4*v10*w4 - 4*v11*w4 + 50*CA2*v11*w4 +
     &        72*CA4*v11*w4 - 128*CA4*v**12*w4 -
     &        36*CA2*v4*w5 + 132*CA4*v4*w5 +
     &        144*CA2*v5*w5 - 528*CA4*v5*w5 + 
     &     12*v6*w5 -
     &        330*CA2*v6*w5 + 688*CA4*v6*w5 - 
     &     36*v7*w5 +
     &        486*CA2*v7*w5 - 216*CA4*v7*w5 + 
     &     116*v8*w5 -
     &        320*CA2*v8*w5 - 518*CA4*v8*w5 - 
     &     172*v9*w5 -
     &        2*CA2*v9*w5 + 780*CA4*v9*w5 + 
     &     62*v10*w5 +
     &        186*CA2*v10*w5 + 102*CA4*v10*w5 + 
     &     18*v11*w5 -
     &        128*CA2*v11*w5 - 440*CA4*v11*w5 -
     &        16*CA2*v**12*w5 + 224*CA4*v**12*w5 +
     &        12*CA2*v4*w6 + 12*CA4*v4*w6 -
     &      48*CA2*v5*w6 -
     &        48*CA4*v5*w6 - 6*v6*w6 + 66*CA2*v6*w6 +
     &        176*CA4*v6*w6 + 18*v7*w6 - 30*CA2*v7*w6 -
     &        360*CA4*v7*w6 - 72*v8*w6 - 
     &     172*CA2*v8*w6 +
     &        382*CA4*v8*w6 + 114*v9*w6 + 
     &     338*CA2*v9*w6 -
     &        220*CA4*v9*w6 - 19*v10*w6 -
     &      244*CA2*v10*w6 -
     &        465*CA4*v10*w6 - 35*v11*w6 + 
     &     78*CA2*v11*w6 +
     &        523*CA4*v11*w6 + 64*CA2*v**12*w6 -
     &        240*CA4*v**12*w6 + 16*CA2*v6*w7 -
     &        40*CA4*v6*w7 - 48*CA2*v7*w7 + 
     &     120*CA4*v7*w7 +
     &        12*v8*w7 + 128*CA2*v8*w7 - 
     &     110*CA4*v8*w7 -
     &        24*v9*w7 - 176*CA2*v9*w7 + 20*CA4*v9*w7 -
     &        24*v10*w7 + 28*CA2*v10*w7 +
     &      348*CA4*v10*w7 +
     &        36*v11*w7 + 52*CA2*v11*w7 -
     &      338*CA4*v11*w7 -
     &        104*CA2*v**12*w7 + 184*CA4*v**12*w7 -
     &        4*CA2*v6*w8 - 4*CA4*v6*w8 +
     &      12*CA2*v7*w8 +
     &        12*CA4*v7*w8 + 2*v8*w8 - 18*CA2*v8*w8 -
     &        20*CA4*v8*w8 - 4*v9*w8 + 16*CA2*v9*w8 +
     &        20*CA4*v9*w8 + 21*v10*w8 + 
     &     70*CA2*v10*w8 -
     &        131*CA4*v10*w8 - 19*v11*w8 -
     &      76*CA2*v11*w8 +
     &        123*CA4*v11*w8 + 88*CA2*v**12*w8 -
     &        104*CA4*v**12*w8 - 4*v10*w9 - 
     &     24*CA2*v10*w9 +
     &        20*CA4*v10*w9 + 4*v11*w9 + 
     &     24*CA2*v11*w9 -
     &        20*CA4*v11*w9 - 40*CA2*v**12*w9 +
     &        40*CA4*v**12*w9 + 8*CA2*v**12*w10 
     &     - 8*CA4*v**12*w10))/
     &   (CA*(1 - v)**2*v3*w2*(1 - v*w)**3*(1 - v + v*w)**3)
      part9 = (2*CF*lv*(80*CA4 - 480*CA4*v + 1392*CA4*v2 -
     &        2560*CA4*v3 + 3264*CA4*v4 - 2976*CA4*v5 +
     &        1936*CA4*v6 - 864*CA4*v7 + 240*CA4*v8 -
     &        32*CA4*v9 - 4*CA2*w + 52*CA4*w + 24*CA2*v*w -
     &        312*CA4*v*w + 4*v2*w - 70*CA2*v2*w +
     &      866*CA4*v2*w - 20*v3*w + 130*CA2*v3*w - 
     &     1470*CA4*v3*w + 48*v4*w -
     &        172*CA2*v4*w + 1948*CA4*v4*w - 72*v5*w +
     &        172*CA2*v5*w - 2404*CA4*v5*w + 68*v6*w -
     &        126*CA2*v6*w + 2698*CA4*v6*w - 36*v7*w +
     &        58*CA2*v7*w - 2422*CA4*v7*w + 8*v8*w -
     &        12*CA2*v8*w + 1508*CA4*v8*w - 560*CA4*v9*w +
     &        96*CA4*v10*w + 4*CA2*w2 + 4*CA4*w2 -
     &        24*CA2*v*w2 - 24*CA4*v*w2 - 2*v2*w2 +
     &        54*CA2*v2*w2 - 4*CA4*v2*w2 + 10*v3*w2 -
     &        50*CA2*v3*w2 + 240*CA4*v3*w2 - 24*v4*w2 -
     &        36*CA2*v4*w2 - 1060*CA4*v4*w2 + 36*v5*w2 +
     &        180*CA2*v5*w2 + 2536*CA4*v5*w2 - 12*v6*w2 -
     &        290*CA2*v6*w2 - 3686*CA4*v6*w2 - 48*v7*w2 +
     &        294*CA2*v7*w2 + 3454*CA4*v7*w2 + 62*v8*w2 -
     &        180*CA2*v8*w2 - 1790*CA4*v8*w2 - 22*v9*w2 +
     &        48*CA2*v9*w2 + 170*CA4*v9*w2 +
     &        256*CA4*v10*w2 - 96*CA4*v11*w2 +
     &        24*CA2*v2*w3 - 144*CA4*v2*w3 -
     &        120*CA2*v3*w3 + 720*CA4*v3*w3 - 4*v4*w3 +
     &        296*CA2*v4*w3 - 1436*CA4*v4*w3 + 16*v5*w3 -
     &        464*CA2*v5*w3 + 1424*CA4*v5*w3 - 70*v6*w3 +
     &        468*CA2*v6*w3 - 874*CA4*v6*w3 + 154*v7*w3 -
     &        284*CA2*v7*w3 + 662*CA4*v7*w3 - 120*v8*w3 +
     &        24*CA2*v8*w3 - 1460*CA4*v8*w3 + 2*v9*w3 +
     &        124*CA2*v9*w3 + 2038*CA4*v9*w3 + 22*v10*w3 -
     &        68*CA2*v10*w3 - 1138*CA4*v10*w3 +
     &        208*CA4*v11*w3 + 32*CA4*v**12*w3 -
     &        12*CA2*v2*w4 - 12*CA4*v2*w4 + 60*CA2*v3*w4 +
     &        60*CA4*v3*w4 + 6*v4*w4 - 102*CA2*v4*w4 -
     &        240*CA4*v4*w4 - 24*v5*w4 + 48*CA2*v5*w4 +
     &        600*CA4*v5*w4 + 67*v6*w4 + 154*CA2*v6*w4 -
     &        437*CA4*v6*w4 - 117*v7*w4 - 378*CA2*v7*w4 -
     &        537*CA4*v7*w4 + 19*v8*w4 + 436*CA2*v8*w4 +
     &        2073*CA4*v8*w4 + 129*v9*w4 - 306*CA2*v9*w4 -
     &        2671*CA4*v9*w4 - 74*v10*w4 + 52*CA2*v10*w4 +
     &        1130*CA4*v10*w4 - 6*v11*w4 + 48*CA2*v11*w4 +
     &        34*CA4*v11*w4 - 144*CA4*v**12*w4 -
     &        36*CA2*v4*w5 + 132*CA4*v4*w5 +
     &        144*CA2*v5*w5 - 528*CA4*v5*w5 -
     &        298*CA2*v6*w5 + 642*CA4*v6*w5 +
     &        390*CA2*v7*w5 - 78*CA4*v7*w5 + 122*v8*w5 -
     &        246*CA2*v8*w5 - 972*CA4*v8*w5 - 244*v9*w5 +
     &        10*CA2*v9*w5 + 1458*CA4*v9*w5 + 94*v10*w5 +
     &        154*CA2*v10*w5 - 124*CA4*v10*w5 + 28*v11*w5 -
     &        118*CA2*v11*w5 - 530*CA4*v11*w5 -
     &        16*CA2*v**12*w5 + 288*CA4*v**12*w5 +
     &        12*CA2*v4*w6 + 12*CA4*v4*w6 - 48*CA2*v5*w6 -
     &        48*CA4*v5*w6 - 6*v6*w6 + 66*CA2*v6*w6 +
     &        196*CA4*v6*w6 + 18*v7*w6 - 30*CA2*v7*w6 -
     &        420*CA4*v7*w6 - 108*v8*w6 - 132*CA2*v8*w6 +
     &        552*CA4*v8*w6 + 186*v9*w6 + 258*CA2*v9*w6 -
     &        460*CA4*v9*w6 - 33*v10*w6 - 190*CA2*v10*w6 -
     &        573*CA4*v10*w6 - 57*v11*w6 + 64*CA2*v11*w6 +
     &        741*CA4*v11*w6 + 64*CA2*v**12*w6 -
     &        344*CA4*v**12*w6 + 16*CA2*v6*w7 -
     &        40*CA4*v6*w7 - 48*CA2*v7*w7 + 120*CA4*v7*w7 +
     &        24*v8*w7 + 112*CA2*v8*w7 - 104*CA4*v8*w7 -
     &        48*v9*w7 - 144*CA2*v9*w7 + 8*CA4*v9*w7 -
     &        38*v10*w7 + 6*CA2*v10*w7 + 540*CA4*v10*w7 +
     &        62*v11*w7 + 58*CA2*v11*w7 - 524*CA4*v11*w7 -
     &        104*CA2*v**12*w7 + 272*CA4*v**12*w7 -
     &        4*CA2*v6*w8 - 4*CA4*v6*w8 + 12*CA2*v7*w8 +
     &        12*CA4*v7*w8 + 2*v8*w8 - 18*CA2*v8*w8 -
     &        32*CA4*v8*w8 - 4*v9*w8 + 16*CA2*v9*w8 +
     &        44*CA4*v9*w8 + 37*v10*w8 + 70*CA2*v10*w8 -
     &        219*CA4*v10*w8 - 35*v11*w8 - 76*CA2*v11*w8 +
     &        199*CA4*v11*w8 + 88*CA2*v**12*w8 -
     &        144*CA4*v**12*w8 - 8*v10*w9 - 24*CA2*v10*w9 +
     &        32*CA4*v10*w9 + 8*v11*w9 + 24*CA2*v11*w9 -
     &        32*CA4*v11*w9 - 40*CA2*v**12*w9 +
     &        48*CA4*v**12*w9 + 8*CA2*v**12*w10 -
     &      8*CA4*v**12*w10))/
     &  (CA*(1 - v)**2*v3*w2*(1 - v*w)**3*(1 - v + v*w)**3)
      struv14= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + 
     &         part9
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV15(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...

      part1 =   (-256*CA**3*lmss*Nf*(1 - 2*v + v2 + v2*w2)*
     &      (1 + v2 - 2*v2*w + v2*w2)*
     &      (4 - 8*v + 4*v2 - v*w + v2*w + 4*v2*w2))/
     &    (243.*(1 - v)*v*w*(1 - v + v*w)**4) 
      part2 = (32*CA**3*lvw*Nf*(81 - 179*v + 115*v2 - 34*v3 - 145*v*w +
     &        130*v2*w - 13*v3*w + 79*v2*w2 - 51*v3*w2 +
     &        17*v3*w3))/(243.*(1 - v)*v2*w) -
     &   (32*CA**3*l1v*Nf*(18 - 45*v + 54*v2 - 44*v3 + 17*v4 +
     &        16*v3*w - 17*v5*w + 20*v2*w2 - 46*v3*w2 +
     &      14*v4*w2 + 24*v5*w2 + 12*v4*w3 - 24*v5*w3 +
     &        2*v4*w4 - v5*w4))/
     &    (243.*(1 - v)*v2*w*(1 - v*w)*(1 - v + v*w))
      part3 = -(32*CA**3*lw*Nf*(88 - 176*v + 104*v2 - 16*v3 + 
     &   13*w -26*v*w + 82*v2*w - 69*v3*w + 18*v4*w - 20*w2 + 40*v*w2 -
     &      164*v2*w2 + 144*v3*w2 - 45*v4*w2 -2*v2*w3 +
     &        2*v3*w3 - 10*v4*w3 - 20*v2*w4 +20*v3*w4 +
     &        20*v4*w4 + 17*v4*w5))/
     &    (243.*(1 - v)*v*(1 - w)*w*(1 - v*w)*(1 - v + v*w))
      part4 = (128*CA**3*lvw*(4 - 8*v + 12*v2 - 8*v3 + 4*v4 - 5*v*w -
     &      4*v2*w + 6*v3*w - 11*v4*w + 2*v5*w + 18*v2*w2 -
     &        20*v3*w2 + 35*v4*w2 - 13*v5*w2 +4*v6*w2 -
     &        8*v3*w3 - 21*v4*w3 + 9*v5*w3 - 8*v6*w3 +
     &        19*v4*w4 - 7*v5*w4 + 12*v6*w4 - 4*v5*w5 -
     &        8*v6*w5 + 4*v6*w6))/((1 - v)**2*v3*w2) 
      part5 = -(128*CA**3*l1v*(2 - 8*v + 16*v2 - 20*v3 + 16*v4 - 8*v5 +
     &        2*v6 + 2*w - 8*v*w + 11*v2*w - 7*v3*w + 2*v4*w -
     &      v5*w + 3*v6*w - 2*v7*w + 2*v3*w2 - 5*v4*w2 +
     &        4*v5*w2 - 4*v6*w2 - v7*w2 + 2*v8*w2 -
     &        2*v2*w3 + 8*v3*w3 - 19*v4*w3 + 25*v5*w3 -
     &        16*v6*w3 + 14*v7*w3 - 6*v8*w3 + v5*w4 -
     &        3*v6*w4 - 10*v7*w4 + 6*v8*w4 - 2*v5*w5 +
     &        3*v6*w5 + 7*v7*w5 - 4*v8*w5 - 4*v7*w6 +
     &  2*v8*w6))/((1 - v)**2*v3*w2*(1 - v*w)*(1 - v + v*w))
      part6 = -(256*CA**3*lms*Nf*(18 - 36*v + 39*v2 - 25*v3 + 
     &   12*v4 -36*w + 74*v2*w - 105*v3*w + 50*v4*w - 14*v5*w +36*w2 +
     &        72*v*w2 - 110*v2*w2 + 30*v3*w2 + 94*v4*w2 -
     &        88*v5*w2 + 40*v6*w2 - 144*v*w3 + 72*v2*w3 +
     &       80*v3*w3 - 140*v4*w3 + 34*v5*w3 + 2*v6*w3-
     &   14*v7*w3 + 216*v2*w4 - 288*v3*w4 + 150*v4*w4 +
     &       40*v5*w4 - 41*v6*w4 + 5*v7*w4 + 12*v8*w4 -
     &   144*v3*w5 + 252*v4*w5 - 208*v5*w5 + 66*v6*w5 +
     &       3*v7*w5 - 16*v8*w5 + 36*v4*w6 - 72*v5*w6 +
     &        70*v6*w6 - 34*v7*w6 + 12*v8*w6))/
     &    (243.*(1 - v)*v3*w*(1 - v*w)**4) 
      part7 = (32*CA**3*l1vw*Nf*(81 - 324*v + 486*v2 - 243*v3 - 
     &   243*v4 + 486*v5 - 324*v6 + 81*v7 + 405*v*w - 1728*v2*w +
     &        2750*v3*w - 1716*v4*w - 147*v5*w + 652*v6*w -
     &  216*v7*w + 1161*v2*w2 - 3808*v3*w2 + 4684*v4*w2 -
     &        2354*v5*w2 + 83*v6*w2 + 234*v7*w2 +
     &        1301*v3*w3 - 3824*v4*w3 + 3660*v5*w3 -
     &        948*v6*w3 - 189*v7*w3 + 1099*v4*w4 -
     &        1812*v5*w4 + 686*v6*w4 + 189*v7*w4 +
     &   167*v5*w5 - 176*v6*w5 - 234*v7*w5 + 27*v6*w6 +
     &        216*v7*w6 - 81*v7*w7))/
     &    (243.*(1 - v)*v2*w*(1 - v + v*w)**4)
      part8 = -(128*CA**3*l1vw*(1 - 4*v + 6*v2 - 3*v3 - 3*v4 + 6*v5 -
     &     4*v6 + v7 + 2*v*w - 10*v2*w + 17*v3*w - 12*v4*w +
     &      4*v5*w - 2*v6*w + v7*w + 6*v2*w2 - 17*v3*w2 +
     &     18*v4*w2 - 20*v5*w2 + 24*v6*w2 - 11*v7*w2 +
     &        3*v3*w3 - 4*v4*w3 + 16*v5*w3 - 36*v6*w3 +
     &        21*v7*w3 + v4*w4 - 10*v5*w4 + 32*v6*w4 -
     &     21*v7*w4 + 4*v5*w5 - 18*v6*w5 + 11*v7*w5 +
     &        4*v6*w6 - v7*w6 - v7*w7))/
     &    ((1 - v)*v2*w*(1 - v + v*w)**4)
      part9 = (128*CA**3*lw*(6 - 24*v + 52*v2 - 72*v3 + 68*v4 - 44*v5 +
     &        18*v6 - 4*v7 - 2*w + 8*v*w - 32*v2*w + 68*v3*w -
     &        80*v4*w + 56*v5*w - 14*v6*w - 4*v7*w + 4*v8*w +
     &     13*v2*w2 - 39*v3*w2 + 48*v4*w2 - 31*v5*w2 -
     &        21*v6*w2 + 30*v7*w2 - 14*v8*w2 - v2*w3 +
     &        3*v3*w3 + v4*w3 - 7*v5*w3 + 65*v6*w3 -
     &        61*v7*w3 + 26*v8*w3 - 3*v4*w4 + 6*v5*w4 -
     &        69*v6*w4 + 66*v7*w4 - 30*v8*w4 + v4*w5 -
     &       2*v5*w5 + 43*v6*w5 - 42*v7*w5 + 22*v8*w5 -
     &       15*v6*w6 + 15*v7*w6 - 10*v8*w6 + 2*v6*w7 -
     &        2*v7*w7 + 2*v8*w7))/
     &    ((1 - v)**2*v3*(1 - w)*w2*(1 - v*w)*(1 - v + v*w)) 
      part10 = -(256*CA**3*lms*(2 - 5*v + 9*v2 - 9*v3 + 7*v4 -3*v5 +
     &   v6 - 8*v*w + 15*v2*w - 28*v3*w + 27*v4*w - 24*v5*w +
     &        11*v6*w - 4*v7*w + w2 - 3*v*w2 + 18*v2*w2 -
     &      17*v3*w2 + 34*v4*w2 - 27*v5*w2 + 30*v6*w2 -
     &       14*v7*w2 + 6*v8*w2 - w3 - v*w3 + 6*v2*w3 -
     &        25*v3*w3 + 2*v4*w3 - 15*v5*w3 - 2*v6*w3 -
     &        13*v7*w3 + 6*v8*w3 - 4*v9*w3 + 4*v*w4 -
     &        6*v2*w4 + 6*v3*w4 + 10*v4*w4 + 27*v5*w4 -
     &        3*v6*w4 + 20*v7*w4 + v9*w4 + v10*w4 -
     &      6*v2*w5 + 14*v3*w5 - 24*v4*w5 + 18*v5*w5 -
     &        53*v6*w5 + 12*v7*w5 - 18*v8*w5 - v10*w5 +
     &      4*v3*w6 - 11*v4*w6 + 21*v5*w6 - 22*v6*w6 +
     &        47*v7*w6 - 9*v8*w6 + 9*v9*w6 + v10*w6 -
     &        v4*w7 + 3*v5*w7 - 6*v6*w7 + 7*v7*w7 -
     &        21*v8*w7 + 2*v9*w7 - 3*v10*w7 + 5*v9*w8 +
     &  v10*w8 - v10*w9))/((1 - v)**2*v3*w2*(1 - v*w)**4)
      part11 = -(256*CA**3*lmss*(1 - 8*v + 30*v2 - 70*v3 + 113*v4 - 
     &   132*v5 +
     &     113*v6 - 70*v7 + 30*v8 - 8*v9 + v10 + w - 3*v*w -
     &   9*v2*w + 71*v3*w - 198*v4*w + 330*v5*w - 367*v6*w +
     &    279*v7*w - 141*v8*w + 43*v9*w - 6*v10*w + 4*v*w2 -
     &    18*v2*w2 + 16*v3*w2 + 76*v4*w2 - 282*v5*w2 +
     &  471*v6*w2 - 472*v7*w2 + 294*v8*w2 - 106*v9*w2 +
     &      17*v10*w2 + 6*v2*w3 - 26*v3*w3 + 18*v4*w3 +
     &   114*v5*w3 - 355*v6*w3 + 497*v7*w3 -393*v8*w3 +
     &    171*v9*w3 - 32*v10*w3 + 4*v3*w4 - 15*v4*w4 -
     &    18*v5*w4 + 170*v6*w4 - 357*v7*w4 +372*v8*w4 -
     &        201*v9*w4 + 45*v10*w4 + v4*w5 - 3*v5*w5 -
     &    39*v6*w5 + 155*v7*w5 - 234*v8*w5 +167*v9*w5 -
     &    47*v10*w5 - 30*v7*w6 + 87*v8*w6 - 92*v9*w6 +
     &   35*v10*w6 - 15*v8*w7 + 31*v9*w7 - 18*v10*w7 -
     &        5*v9*w8 + 6*v10*w8 - v10*w9))/
     &    ((1 - v)**2*v3*w2*(1 - v + v*w)**4)
      part12 = (64*CA**3*lv*Nf*(72 - 432*v + 1180*v2 - 1940*v3 + 
     &   2136*v4 -
     &    1656*v5 + 908*v6 - 324*v7 + 56*v8 - 144*w + 864*v*w -
     &        2008*v2*w + 2120*v3*w - 524*v4*w - 1120*v5*w +
     &   1399*v6*w - 877*v7*w + 377*v8*w - 87*v9*w + 144*w2 -
     &   864*v*w2 + 1432*v2*w2 + 760*v3*w2 - 5168*v4*w2 +
     &        6608*v5*w2 - 3166*v6*w2 - 934*v7*w2 +
     &        2011*v8*w2 - 1028*v9*w2 + 205*v10*w2 +
     &        1152*v2*w3 - 5760*v3*w3 + 10880*v4*w3 -
     &        8960*v5*w3 - 162*v6*w3 + 7654*v7*w3 -
     &        7149*v8*w3 + 2608*v9*w3 - 178*v10*w3 -
     &        85*v11*w3 - 576*v2*w4 + 2880*v3*w4 -
     &        3280*v4*w4 - 4160*v5*w4 + 14136*v6*w4 -
     &        15752*v7*w4 + 8900*v8*w4 - 2160*v9*w4 +
     &        60*v10*w4 - 48*v11*w4 + 55*v**12*w4 -
     &        2592*v4*w5 + 10368*v5*w5 - 14832*v6*w5 +
     &        8208*v7*w5 + 286*v8*w5 - 2156*v9*w5 +
     &        57*v10*w5 + 661*v11*w5 - 239*v**12*w5 +
     &        864*v4*w6 - 3456*v5*w6 + 2256*v6*w6 +
     &        5328*v7*w6 - 9360*v8*w6 + 5808*v9*w6 +
     &        82*v10*w6 - 1522*v11*w6 + 533*v**12*w6 +
     &        2304*v6*w7 - 6912*v7*w7 + 7360*v8*w7 -
     &        3200*v9*w7 - 1638*v10*w7 + 2086*v11*w7 -
     &        707*v**12*w7 - 576*v6*w8 + 1728*v7*w8 -
     &        760*v8*w8 - 1360*v9*w8 + 2604*v10*w8 -
     &        1636*v11*w8 + 560*v**12*w8 - 720*v8*w9 +
     &        1440*v9*w9 - 1400*v10*w9 + 680*v11*w9 -
     &        266*v**12*w9 + 144*v8*w10 - 288*v9*w10 +
     &        280*v10*w10 - 136*v11*w10 + 64*v**12*w10))/
     &    (243.*(1 - v)*v3*w*(1 - v*w)**4*(1 - v + v*w)**4)
      part13 = (32*CA**3*l1w*Nf*(144 - 864*v + 2369*v2 - 3925*v3 + 
     &     4361*v4 -
     &        3398*v5 + 1855*v6 - 653*v7 + 111*v8 - 288*w +
     &        1728*v*w - 3996*v2*w + 4140*v3*w - 781*v4*w -
     &        2708*v5*w + 3298*v6*w - 2036*v7*w + 815*v8*w -
     &        172*v9*w + 288*w2 - 1728*v*w2 + 2864*v2*w2 +
     &        1520*v3*w2 - 10317*v4*w2 + 13140*v5*w2 -
     &        6044*v6*w2 - 2466*v7*w2 + 4563*v8*w2 -
     &        2230*v9*w2 + 410*v10*w2 + 2304*v2*w3 -
     &        11520*v3*w3 + 21720*v4*w3 - 17760*v5*w3 -
     &        744*v6*w3 + 16008*v7*w3 - 14674*v8*w3 +
     &        4988*v9*w3 - 150*v10*w3 - 172*v11*w3 -
     &        1152*v2*w4 + 5760*v3*w4 - 6560*v4*w4 -
     &        8320*v5*w4 + 28298*v6*w4 - 31582*v7*w4 +
     &        17442*v8*w4 - 3474*v9*w4 - 207*v10*w4 -
     &        205*v11*w4 + 111*v**12*w4 - 5184*v4*w5 +
     &        20736*v5*w5 - 29664*v6*w5 + 16416*v7*w5 +
     &        926*v8*w5 - 5020*v9*w5 + 118*v10*w5 +
     &        1672*v11*w5 - 457*v**12*w5 + 1728*v4*w6 -
     &        6912*v5*w6 + 4512*v6*w6 + 10656*v7*w6 -
     &        18882*v8*w6 + 11940*v9*w6 + 452*v10*w6 -
     &        3494*v11*w6 + 973*v**12*w6 + 4608*v6*w7 -
     &        13824*v7*w7 + 14760*v8*w7 - 6480*v9*w7 -
     &        3552*v10*w7 + 4488*v11*w7 - 1254*v**12*w7 -
     &        1152*v6*w8 + 3456*v7*w8 - 1520*v8*w8 -
     &        2720*v9*w8 + 5333*v10*w8 - 3397*v11*w8 +
     &        973*v**12*w8 - 1440*v8*w9 + 2880*v9*w9 -
     &        2820*v10*w9 + 1380*v11*w9 - 457*v**12*w9 +
     &        288*v8*w10 - 576*v9*w10 + 560*v10*w10 -
     &        272*v11*w10 + 111*v**12*w10))/
     &    (243.*(1 - v)*v3*w*(1 - v*w)**4*(1 - v + v*w)**4)
      part14 = -(64*CA**3*Nf*(126 - 882*v + 2882*v2 - 5842*v3 + 
     &   8172*v4 - 8236*v5 + 6002*v6 - 3042*v7 + 962*v8 - 142*v9 -
     &   630*w + 4410*v*w - 13105*v2*w + 21368*v3*w - 20644*v4*w +
     &        12127*v5*w - 4721*v6*w + 1654*v7*w - 444*v8*w -
     &        87*v9*w + 72*v10*w + 684*w2 - 4788*v*w2 +
     &        11890*v2*w2 - 9096*v3*w2 - 12922*v4*w2 +
     &        33660*v5*w2 - 26657*v6*w2 + 1132*v7*w2 +
     &        14112*v8*w2 - 11608*v9*w2 + 4245*v10*w2 -
     &        652*v11*w2 + 5256*v2*w3 - 31536*v3*w3 +
     &        77264*v4*w3 - 96952*v5*w3 + 55397*v6*w3 +
     &        12308*v7*w3 - 44489*v8*w3 + 31249*v9*w3 -
     &        9058*v10*w3 + 219*v11*w3 + 342*v**12*w3 -
     &        2736*v2*w4 + 16416*v3*w4 - 30820*v4*w4 +
     &        3620*v5*w4 + 68334*v6*w4 - 115112*v7*w4 +
     &        95658*v8*w4 - 48090*v9*w4 + 16335*v10*w4 -
     &        4656*v11*w4 + 1409*v**12*w4 - 358*v**13*w4 -
     &        11988*v4*w5 + 59940*v5*w5 - 116946*v6*w5 +
     &        108584*v7*w5 - 45969*v8*w5 + 9015*v9*w5 -
     &        9443*v10*w5 + 10929*v11*w5 - 5060*v**12*w5 +
     &        938*v**13*w5 + 54*v**14*w5 + 4104*v4*w6 -
     &        20520*v5*w6 + 28068*v6*w6 + 10848*v7*w6 -
     &        59922*v8*w6 + 54966*v9*w6 - 14329*v10*w6 -
     &        8688*v11*w6 + 6651*v**12*w6 - 1178*v**13*w6 -
     &        162*v**14*w6 + 10728*v6*w7 - 42912*v7*w7 +
     &        65464*v8*w7 - 45912*v9*w7 + 7619*v10*w7 +
     &        11162*v11*w7 - 7507*v**12*w7 + 1358*v**13*w7 +
     &        162*v**14*w7 - 2736*v6*w8 + 10944*v7*w8 -
     &        11866*v8*w8 - 2706*v9*w8 + 16200*v10*w8 -
     &        15338*v11*w8 + 6788*v**12*w8 - 1286*v**13*w8 -
     &        54*v**14*w8 - 3366*v8*w9 + 10098*v9*w9 -
     &        12677*v10*w9 + 8592*v11*w9 - 3315*v**12*w9 +
     &        668*v**13*w9 + 684*v8*w10 - 2052*v9*w10 +
     &        2602*v10*w10 - 1784*v11*w10 + 692*v**12*w10 -
     &        142*v**13*w10))/
     &    (243.*(1 - v)**2*v3*w*(1 - v*w)**4*(1 - v + v*w)**4)
      part15 = (64*CA**3*(17 - 119*v + 391*v2 - 799*v3 + 1122*v4 -
     &   1122*v5 + 799*v6 - 391*v7 + 119*v8 - 17*v9 - 38*w +
     &   266*v*w - 758*v2*w + 1090*v3*w - 683*v4*w - 237*v5*w +
     &        833*v6*w - 788*v7*w + 449*v8*w - 163*v9*w +
     &        29*v10*w + 40*w2 - 280*v*w2 + 670*v2*w2 -
     &        380*v3*w2 - 1053*v4*w2 + 2075*v5*w2 -
     &        1015*v6*w2 - 1130*v7*w2 + 2111*v8*w2 -
     &        1495*v9*w2 + 543*v10*w2 - 86*v11*w2 +
     &        312*v2*w3 - 1872*v3*w3 + 4560*v4*w3 -
     &        5640*v5*w3 + 3131*v6*w3 + 724*v7*w3 -
     &        2216*v8*w3 + 1018*v9*w3 + 233*v10*w3 -
     &        334*v11*w3 + 84*v**12*w3 - 160*v2*w4 +
     &        960*v3*w4 - 1802*v4*w4 + 210*v5*w4 +
     &        3970*v6*w4 - 6580*v7*w4 + 5059*v8*w4 -
     &        1825*v9*w4 + 181*v10*w4 - 137*v11*w4 +
     &        185*v**12*w4 - 61*v**13*w4 - 708*v4*w5 +
     &        3540*v5*w5 - 7028*v6*w5 + 6872*v7*w5 -
     &        3309*v8*w5 + 743*v9*w5 - 719*v10*w5 +
     &        1137*v11*w5 - 651*v**12*w5 + 123*v**13*w5 +
     &        11*v**14*w5 + 240*v4*w6 - 1200*v5*w6 +
     &        1696*v6*w6 + 416*v7*w6 - 3350*v8*w6 +
     &        3554*v9*w6 - 1121*v10*w6 - 796*v11*w6 +
     &        642*v**12*w6 - 81*v**13*w6 - 33*v**14*w6 +
     &        632*v6*w7 - 2528*v7*w7 + 4048*v8*w7 -
     &        3296*v9*w7 + 869*v10*w7 + 806*v11*w7 -
     &        602*v**12*w7 + 71*v**13*w7 + 33*v**14*w7 -
     &160*v6*w8 + 640*v7*w8 - 751*v8*w8 + 13*v9*w8 +
     &        935*v10*w8 - 1145*v11*w8 + 571*v**12*w8 -
     &        103*v**13*w8 - 11*v**14*w8 - 198*v8*w9 +
     &        594*v9*w9 - 822*v10*w9 + 654*v11*w9 -
     &        296*v**12*w9 + 68*v**13*w9 + 40*v8*w10 -
     &        120*v9*w10 + 170*v10*w10 - 140*v11*w10 +
     &        67*v**12*w10 - 17*v**13*w10))/
     &    (3.*(1 - v)**2*v3*w*(1 - v*w)**4*(1 - v + v*w)**4)
      part16 = (128*CA**3*l1w*(8 - 56*v + 188*v2 - 400*v3 + 596*v4 -
     & 648*v5 + 520*v6 - 304*v7 + 124*v8 - 32*v9 + 4*v10 +
     &        6*w - 42*v*w + 144*v2*w - 318*v3*w + 535*v4*w -
     &        761*v5*w + 920*v6*w - 896*v7*w + 657*v8*w -
     &335*v9*w + 106*v10*w - 16*v11*w + 2*w2 - 14*v*w2 +
     & 37*v2*w2 - 40*v3*w2 - 62*v4*w2 + 343*v5*w2 -
     &657*v6*w2 + 702*v7*w2 - 388*v8*w2 - 13*v9*w2 +
     & 172*v10*w2 - 106*v11*w2 + 24*v**12*w2 - 2*w3 +
     & 14*v*w3 - 64*v2*w3 + 202*v3*w3 - 451*v4*w3 +
     &        737*v5*w3 - 1005*v6*w3 + 1248*v7*w3 -
     &        1375*v8*w3 + 1231*v9*w3 - 741*v10*w3 +
     & 220*v11*w3 + 2*v**12*w3 - 16*v**13*w3 - 16*v2*w4 +
     & 96*v3*w4 - 290*v4*w4 + 570*v5*w4 - 739*v6*w4 +
     & 592*v7*w4 - 221*v8*w4 - 71*v9*w4 - 22*v10*w4 +
     &        241*v11*w4 - 190*v**12*w4 + 50*v**13*w4 +
     &  4*v**14*w4 + 8*v2*w5 - 48*v3*w5 + 144*v4*w5 -
     &  280*v5*w5 + 370*v6*w5 - 328*v7*w5 + 166*v8*w5 +
     &        2*v9*w5 + 127*v10*w5 - 344*v11*w5 +
     &  214*v**12*w5 - 31*v**13*w5 - 20*v**14*w5 + 36*v4*w6 -
     &  180*v5*w6 + 466*v6*w6 - 784*v7*w6 + 975*v8*w6 -
     &     937*v9*w6 + 481*v10*w6 + 45*v11*w6 -
     &  88*v**12*w6 - 14*v**13*w6 + 42*v**14*w6 - 12*v4*w7 +
     &  60*v5*w7 - 132*v6*w7 + 168*v7*w7 - 122*v8*w7 +
     &        30*v9*w7 + 183*v10*w7 - 340*v11*w7 +
     &  166*v**12*w7 - v**13*w7 - 62*v**14*w7 - 32*v6*w8 +
     &        128*v7*w8 - 284*v8*w8 + 404*v9*w8 -
     &        437*v10*w8 + 350*v11*w8 - 115*v**12*w8 -
     &  14*v**13*w8 + 90*v**14*w8 + 8*v6*w9 - 32*v7*w9 +
     &  58*v8*w9 - 62*v9*w9 + 22*v10*w9 + 22*v11*w9 -
     &        75*v**12*w9 + 59*v**13*w9 - 110*v**14*w9 +
     &        10*v8*w10 - 30*v9*w10 + 65*v10*w10 -
     &        80*v11*w10 + 71*v**12*w10 - 36*v**13*w10 +
     & 92*v**14*w10 - 2*v8*w11 + 6*v9*w11 - 12*v10*w11 +
     & 14*v11*w11 - 7*v**12*w11 + v**13*w11 - 52*v**14*w11 -
     & 2*v**12*w**12 + 2*v**13*w**12 + 20*v**14*w**12 - 4*v**14*w**13))
     &     /((1 - v)**2*v3*w2*(1 - v*w)**4*(1 - v + v*w)**4)
      part17 = (128*CA**3*lv*(10 - 70*v + 234*v2 - 494*v3 + 728*v4 -
     & 780*v5 + 614*v6 - 350*v7 + 138*v8 - 34*v9 + 4*v10 +
     &  6*w - 42*v*w + 146*v2*w - 330*v3*w + 572*v4*w -
     &        836*v5*w + 1026*v6*w - 1002*v7*w + 730*v8*w -
     &  366*v9*w + 112*v10*w - 16*v11*w + 2*w2 - 14*v*w2 +
     &  32*v2*w2 - 10*v3*w2 - 152*v4*w2 + 518*v5*w2 -
     &  896*v6*w2 + 938*v7*w2 - 548*v8*w2 + 46*v9*w2 +
     &        170*v10*w2 - 110*v11*w2 + 24*v**12*w2 -2*w3 +
     &        14*v*w3 - 64*v2*w3 + 202*v3*w3 -448*v4*w3 +
     &        722*v5*w3 - 973*v6*w3 + 1210*v7*w3 -
     &        1380*v8*w3 + 1316*v9*w3 - 839*v10*w3 +
     &  260*v11*w3 - 2*v**12*w3 - 16*v**13*w3 - 16*v2*w4 +
     &  96*v3*w4 - 288*v4*w4 + 560*v5*w4 - 703*v6*w4 +
     &  508*v7*w4 - 41*v8*w4 - 359*v9*w4 + 214*v10*w4 +
     &        171*v11*w4 - 198*v**12*w4 + 56*v**13*w4 +
     &        4*v**14*w4 + 8*v2*w5 - 48*v3*w5+144*v4*w5 -
     &  280*v5*w5 + 358*v6*w5 - 280*v7*w5 + 30*v8*w5 +
     &        242*v9*w5 - 23*v10*w5 - 388*v11*w5 +
     &  293*v**12*w5 - 56*v**13*w5 - 22*v**14*w5 + 36*v4*w6 -
     &        180*v5*w6 + 468*v6*w6 - 792*v7*w6 +
     &        1005*v8*w6 - 999*v9*w6 + 430*v10*w6 +
     &        241*v11*w6 - 228*v**12*w6 + 19*v**13*w6 +
     &  54*v**14*w6 - 12*v4*w7 + 60*v5*w7 - 132*v6*w7 +
     &  168*v7*w7 - 124*v8*w7 + 36*v9*w7 + 265*v10*w7 -
     &        514*v11*w7 + 261*v**12*w7 - 8*v**13*w7 -
     &  92*v**14*w7 - 32*v6*w8 + 128*v7*w8 - 282*v8*w8 +
     &  398*v9*w8 - 473*v10*w8 + 432*v11*w8 -
     &  137*v**12*w8 - 34*v**13*w8 + 132*v**14*w8 + 8*v6*w9 -
     &  32*v7*w9 + 58*v8*w9 - 62*v9*w9 + 40*v10*w9 -
     &        14*v11*w9 - 80*v**12*w9 + 82*v**13*w9 -
     &        148*v**14*w9 + 10*v8*w10 - 30*v9*w10 +
     &        60*v10*w10 - 70*v11*w10 + 83*v**12*w10 -
     &        53*v**13*w10 + 116*v**14*w10 - 2*v8*w11 +
     &        6*v9*w11 - 12*v10*w11 + 14*v11*w11 -
     &        16*v**12*w11 + 10*v**13*w11 - 62*v**14*w11 +
     &        22*v**14*w**12 - 4*v**14*w**13))/
     &    ((1 - v)**2*v3*w2*(1 - v*w)**4*(1 - v + v*w)**4)
      struv15= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + 
     &         part9 + part10 + part11 + part12 +
     &         part13 + part14 + part15 + part16 + part17 
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRUV16(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,CA,GTR,CF,VC
      COMMON / PRECOLOR / CACF,CA2,CA4
      COMMON / PREV / V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12
      COMMON / PREW / W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12
      COMMON / PRELOG / L1V,LV,L1W,LW,LVW,L1VW,LMU,LMS,LMSS
*
cc      CALL PRECALC(V,W,S)
      Nf = 2.*GTR 
C...
      part1 =  (-4*CF*lvw*(4*CA2 - 4*CA4 - 12*CA2*v + 
     &        12*CA4*v +
     &        16*CA2*v2 - 16*CA4*v2 - 8*CA2*v3 + 
     &        8*CA4*v3 -
     &    2*w - 2*CA2*w + 7*v*w + 11*CA2*v*w - 9*CA4*v*w -
     &    8*v2*w - 16*CA2*v2*w + 28*CA4*v2*w + 5*v3*w +
     &  9*CA2*v3*w - 29*CA4*v3*w - 2*v4*w + 2*CA2*v4*w +
     &  6*CA4*v4*w - 3*v*w2 - 5*CA2*v*w2 + 19*CA4*v*w2 +
     &        5*v2*w2 + 4*CA2*v2*w2 - 54*CA4*v2*w2 -
     &        6*v3*w2 - 6*CA2*v3*w2 + 87*CA4*v3*w2 +
     &        4*v4*w2 - 5*CA2*v4*w2 - 48*CA4*v4*w2 +
     &  16*CA4*v5*w2 + 5*CA2*v2*w3 - 13*CA4*v2*w3 +
     &        3*v3*w3 + 8*CA2*v3*w3 - 36*CA4*v3*w3 -
     &        3*v4*w3 + 3*CA2*v4*w3 + 33*CA4*v4*w3 -
     &        32*CA4*v5*w3 - v3*w4 - 7*CA2*v3*w4 +
     &        43*CA4*v3*w4 + v4*w4 - CA2*v4*w4 -
     &  27*CA4*v4*w4 + 48*CA4*v5*w4 - 16*CA4*v4*w5 -
     &        32*CA4*v5*w5 + 16*CA4*v5*w6))/
     &    (CA*(1 - v)**2*v2*w2)
      part2 = -(4*CF*l1vw*(3 - 2*CA2 + 3*CA4 - 8*CA2*CF**2 - 14*v +
     &        11*CA2*v - 15*CA4*v - 8*CA*CF*v + 8*CA**3*CF*v +
     &  24*CA2*CF**2*v + 27*v2 - 25*CA2*v2 + 35*CA4*v2 +
     &        28*CA*CF*v2 - 28*CA**3*CF*v2 - 28*CA2*CF**2*v2 -
     &        28*v3 + 30*CA2*v3 - 50*CA4*v3 - 40*CA*CF*v3 +
     &        40*CA**3*CF*v3 + 16*CA2*CF**2*v3 + 17*v4 -
     &        20*CA2*v4 + 45*CA4*v4 + 30*CA*CF*v4 -
     & 30*CA**3*CF*v4 - 4*CA2*CF**2*v4 - 6*v5 + 7*CA2*v5 -
     &        23*CA4*v5 - 12*CA*CF*v5 + 12*CA**3*CF*v5 + v6 -
     & CA2*v6 + 5*CA4*v6 + 2*CA*CF*v6 - 2*CA**3*CF*v6 +
     &        12*v*w - 6*CA2*v*w + CA4*v*w + 8*CA*CF*v*w -
     &        8*CA2*CF**2*v*w - 44*v2*w + 29*CA2*v2*w -
     &        8*CA4*v2*w - 40*CA*CF*v2*w + 8*CA**3*CF*v2*w +
     &        24*CA2*CF**2*v2*w + 64*v3*w - 56*CA2*v3*w +
     &        18*CA4*v3*w + 80*CA*CF*v3*w - 28*CA**3*CF*v3*w -
     &        24*CA2*CF**2*v3*w - 48*v4*w + 54*CA2*v4*w -
     &        16*CA4*v4*w - 80*CA*CF*v4*w + 36*CA**3*CF*v4*w +
     &        8*CA2*CF**2*v4*w + 20*v5*w - 26*CA2*v5*w +
     &        5*CA4*v5*w + 40*CA*CF*v5*w - 20*CA**3*CF*v5*w -
     &        4*v6*w + 5*CA2*v6*w - 8*CA*CF*v6*w +
     &        4*CA**3*CF*v6*w + 19*v2*w2 - 9*CA2*v2*w2 +
     &  3*CA4*v2*w2 + 12*CA*CF*v2*w2 - 4*CA**3*CF*v2*w2 -
     & 12*CA2*CF**2*v2*w2 - 52*v3*w2 + 38*CA2*v3*w2 -
     &        6*CA4*v3*w2 - 56*CA*CF*v3*w2 +
     &        24*CA**3*CF*v3*w2 + 16*CA2*CF**2*v3*w2 +
     &        54*v4*w2 - 60*CA2*v4*w2 - 10*CA4*v4*w2 +
     &        84*CA*CF*v4*w2 - 36*CA**3*CF*v4*w2 -
     &  8*CA2*CF**2*v4*w2 - 28*v5*w2 + 42*CA2*v5*w2 +
     &        26*CA4*v5*w2 - 56*CA*CF*v5*w2 +
     &        24*CA**3*CF*v5*w2 + 7*v6*w2 -11*CA2*v6*w2 -
     &        13*CA4*v6*w2 + 14*CA*CF*v6*w2 -
     &        6*CA**3*CF*v6*w2 + 16*v3*w3 -12*CA2*v3*w3 +
     &  6*CA4*v3*w3 + 16*CA*CF*v3*w3 - 4*CA**3*CF*v3*w3 -
     &  8*CA2*CF**2*v3*w3 - 32*v4*w3 + 38*CA2*v4*w3 -
     &        48*CA*CF*v4*w3 + 20*CA**3*CF*v4*w3 +
     &  8*CA2*CF**2*v4*w3 + 24*v5*w3 - 40*CA2*v5*w3 -
     &        22*CA4*v5*w3 + 48*CA*CF*v5*w3 -
     &  24*CA**3*CF*v5*w3 - 8*v6*w3 + 14*CA2*v6*w3 +
     &        16*CA4*v6*w3 - 16*CA*CF*v6*w3 +
     &        8*CA**3*CF*v6*w3 + 9*v4*w4 - 12*CA2*v4*w4 -
     &  3*CA4*v4*w4 + 14*CA*CF*v4*w4 - 6*CA**3*CF*v4*w4 -
     &  4*CA2*CF**2*v4*w4 - 14*v5*w4 + 23*CA2*v5*w4 +
     &        21*CA4*v5*w4 - 28*CA*CF*v5*w4 +
     &        12*CA**3*CF*v5*w4 + 7*v6*w4 -11*CA2*v6*w4 -
     &        13*CA4*v6*w4 + 14*CA*CF*v6*w4 -
     &        6*CA**3*CF*v6*w4 + 4*v5*w5 - 6*CA2*v5*w5 -
     &  7*CA4*v5*w5 + 8*CA*CF*v5*w5 - 4*CA**3*CF*v5*w5 -
     &        4*v6*w5 + 5*CA2*v6*w5 - 8*CA*CF*v6*w5 +
     &        4*CA**3*CF*v6*w5 + v6*w6 - CA2*v6*w6 +
     & 5*CA4*v6*w6 + 2*CA*CF*v6*w6 - 2*CA**3*CF*v6*w6))/
     &    (CA*(1 - v)*v2*w*(1 - v + v*w)**3)
      part3 = -(4*CF*l1v*(4*CA4 - 20*CA4*v + 44*CA4*v2 - 52*CA4*v3 +
     & 32*CA4*v4 - 8*CA4*v5 - 2*CA2*w + 2*CA4*w - v*w +
     &        5*CA2*v*w - 3*CA4*v*w + 2*v2*w - 6*CA2*v2*w -
     & 10*CA4*v2*w - 2*v3*w + 6*CA2*v3*w + 24*CA4*v3*w +
     &        2*v4*w - 4*CA2*v4*w - 8*CA4*v4*w - v5*w +
     &        CA2*v5*w - 13*CA4*v5*w + 8*CA4*v6*w +
     &        2*CA2*v2*w2 - CA4*v2*w2 + v3*w2 -
     &        7*CA2*v3*w2 + 6*CA4*v3*w2 - 3*v4*w2 +
     &        5*CA2*v4*w2 - 24*CA4*v4*w2 + v5*w2 +
     &        CA2*v5*w2 + 22*CA4*v5*w2 + v6*w2 -
     &        CA2*v6*w2 + 5*CA4*v6*w2 - 8*CA4*v7*w2 -
     &10*CA4*v2*w3 + 2*CA2*v3*w3 + 42*CA4*v3*w3 +
     &        v4*w3 + 3*CA2*v4*w3 - 56*CA4*v4*w3 +
     &        2*v5*w3 - 8*CA2*v5*w3 + 52*CA4*v5*w3 -
     &        3*v6*w3 + 3*CA2*v6*w3 - 56*CA4*v6*w3 +
     &24*CA4*v7*w3 - 6*CA2*v4*w4 - 17*CA4*v4*w4 -
     &        3*v5*w4 + 9*CA2*v5*w4 + 13*CA4*v5*w4 +
     &        3*v6*w4 - 3*CA2*v6*w4 + 40*CA4*v6*w4 -
     &24*CA4*v7*w4 + 2*CA2*v4*w5 + 8*CA4*v4*w5 +
     &        v5*w5 - 3*CA2*v5*w5 - 13*CA4*v5*w5 -
     &        v6*w5 + CA2*v6*w5 - 27*CA4*v6*w5 +
     & 16*CA4*v7*w5 + 16*CA4*v6*w6 - 8*CA4*v7*w6))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)*(1 - v + v*w))
      part4 = -(4*CF*lw*(4*CA2 - 8*CA4 - 16*CA2*v + 44*CA4*v +
     &        28*CA2*v2 - 112*CA4*v2 - 24*CA2*v3 +
     &        164*CA4*v3 + 8*CA2*v4 - 144*CA4*v4 +
     &        72*CA4*v5 - 16*CA4*v6 - w + 5*CA4*w + 3*v*w +
     &        5*CA2*v*w - 48*CA4*v*w - 3*v2*w - 13*CA2*v2*w +
     & 164*CA4*v2*w + v3*w + 9*CA2*v3*w - 282*CA4*v3*w +
     &        9*CA2*v4*w + 247*CA4*v4*w - 10*CA2*v5*w -
     & 86*CA4*v5*w - 16*CA4*v6*w + 16*CA4*v7*w - v*w2 -
     &        6*CA2*v*w2 + 28*CA4*v*w2 + v2*w2 +
     &        9*CA2*v2*w2 - 129*CA4*v2*w2 + v3*w2 +
     &        3*CA2*v3*w2 + 266*CA4*v3*w2 - v4*w2 -
     &29*CA2*v4*w2 - 227*CA4*v4*w2 + 27*CA2*v5*w2 -
     & 4*CA4*v5*w2 + 118*CA4*v6*w2 - 56*CA4*v7*w2 -
     &        6*CA4*v*w3 + 6*CA2*v2*w3 + 29*CA4*v2*w3 -
     &        17*CA2*v3*w3 - 97*CA4*v3*w3 - v4*w3 +
     &        35*CA2*v4*w3 + 67*CA4*v4*w3 + v5*w3 -
     &        40*CA2*v5*w3 + 162*CA4*v5*w3 -
     & 239*CA4*v6*w3 + 104*CA4*v7*w3 + 3*CA2*v3*w4 +
     &  25*CA4*v3*w4 - 11*CA2*v4*w4 - 4*CA4*v4*w4 +
     &        v5*w4 + 34*CA2*v5*w4 - 202*CA4*v5*w4 -
     &        v6*w4 + 2*CA2*v6*w4 + 257*CA4*v6*w4 -
     &        120*CA4*v7*w4 - 2*CA4*v3*w5 + v4*w5 -
     &        2*CA2*v4*w5 - 8*CA4*v4*w5 - 3*v5*w5 -
     &        19*CA2*v5*w5 + 137*CA4*v5*w5 + 2*v6*w5 -
     &  3*CA2*v6*w5 - 163*CA4*v6*w5 + 88*CA4*v7*w5 +
     &        v5*w6 + 7*CA2*v5*w6 - 51*CA4*v5*w6 -
     &        v6*w6 + CA2*v6*w6 + 59*CA4*v6*w6 -
     &    40*CA4*v7*w6 + 8*CA4*v5*w7 - 8*CA4*v6*w7 +
     &        8*CA4*v7*w7))/
     &    (CA*(1 - v)**2*v2*(1 - w)*w2*(1 - v*w)*(1 - v + v*w))
      part5 = -(4*CF*lmss*(4*CA**3 - 32*CA**3*v + 116*CA**3*v2 - 
     &   248*CA**3*v3 +
     &        340*CA**3*v4 - 304*CA**3*v5 + 172*CA**3*v6 -
     & 56*CA**3*v7 + 8*CA**3*v8 + 4*CA*CF**2*w + 20*CA**3*v*w +
     &        4*CF*v*w - 4*CA2*CF*v*w - 16*CA*CF**2*v*w -
     &        148*CA**3*v2*w - 18*CF*v2*w + 18*CA2*CF*v2*w +
     &        26*CA*CF**2*v2*w + 480*CA**3*v3*w + 34*CF*v3*w -
     & 34*CA2*CF*v3*w - 22*CA*CF**2*v3*w - 880*CA**3*v4*w -
     &        35*CF*v4*w + 35*CA2*CF*v4*w + 10*CA*CF**2*v4*w +
     &        980*CA**3*v5*w + 21*CF*v5*w - 21*CA2*CF*v5*w -
     &        2*CA*CF**2*v5*w - 660*CA**3*v6*w - 7*CF*v6*w +
     &        7*CA2*CF*v6*w + 248*CA**3*v7*w + CF*v7*w -
     &        CA2*CF*v7*w - 40*CA**3*v8*w + 12*CA*CF**2*v*w2 +
     & 56*CA**3*v2*w2 + 6*CF*v2*w2 - 14*CA2*CF*v2*w2 -
     & 54*CA*CF**2*v2*w2 - 360*CA**3*v3*w2 - 26*CF*v3*w2 +
     &        72*CA2*CF*v3*w2 + 98*CA*CF**2*v3*w2 +
     &        984*CA**3*v4*w2 + 45*CF*v4*w2 -
     &        157*CA2*CF*v4*w2 - 90*CA*CF**2*v4*w2 -
     &        1456*CA**3*v5*w2 - 39*CF*v5*w2 +
     &        187*CA2*CF*v5*w2 + 42*CA*CF**2*v5*w2 +
     &        1224*CA**3*v6*w2 + 17*CF*v6*w2 -
     &        129*CA2*CF*v6*w2 - 8*CA*CF**2*v6*w2 -
     & 552*CA**3*v7*w2 - 3*CF*v7*w2 + 49*CA2*CF*v7*w2 +
     &        104*CA**3*v8*w2 - 8*CA2*CF*v8*w2 +
     & 12*CA*CF**2*v2*w3 + 100*CA**3*v3*w3 + 8*CF*v3*w3 -
     &        14*CA2*CF*v3*w3 - 48*CA*CF**2*v3*w3 -
     & 548*CA**3*v4*w3 - 25*CF*v4*w3 + 65*CA2*CF*v4*w3 +
     &        78*CA*CF**2*v4*w3 + 1220*CA**3*v5*w3 +
     &        31*CF*v5*w3 - 131*CA2*CF*v5*w3 -
     &        58*CA*CF**2*v5*w3 - 1372*CA**3*v6*w3 -
     &        18*CF*v6*w3 + 138*CA2*CF*v6*w3 +
     & 16*CA*CF**2*v6*w3 + 776*CA**3*v7*w3 + 4*CF*v7*w3 -
     &        74*CA2*CF*v7*w3 - 176*CA**3*v8*w3 +
     &        16*CA2*CF*v8*w3 + 4*CA*CF**2*v3*w4 +
     & 128*CA**3*v4*w4 + 7*CF*v4*w4 - 7*CA2*CF*v4*w4 -
     & 14*CA*CF**2*v4*w4 - 568*CA**3*v5*w4 - 17*CF*v5*w4 +
     &        23*CA2*CF*v5*w4 + 18*CA*CF**2*v5*w4 +
     & 960*CA**3*v6*w4 + 14*CF*v6*w4 - 34*CA2*CF*v6*w4 -
     & 8*CA*CF**2*v6*w4 - 728*CA**3*v7*w4 - 4*CF*v7*w4 +
     &        26*CA2*CF*v7*w4 + 208*CA**3*v8*w4 -
     & 8*CA2*CF*v8*w4 + 116*CA**3*v5*w5 + 4*CF*v5*w5 -
     & 2*CA2*CF*v5*w5 - 396*CA**3*v6*w5 - 7*CF*v6*w5 +
     & 3*CA2*CF*v6*w5 + 456*CA**3*v7*w5 + 3*CF*v7*w5 -
     & CA2*CF*v7*w5 - 176*CA**3*v8*w5 + 76*CA**3*v6*w6 +
     &        CF*v6*w6 - CA2*CF*v6*w6 - 176*CA**3*v7*w6 -
     &        CF*v7*w6 + CA2*CF*v7*w6 + 104*CA**3*v8*w6 +
     &  32*CA**3*v7*w7 - 40*CA**3*v8*w7 + 8*CA**3*v8*w8))/
     &    ((1 - v)**2*v2*w2*(1 - v + v*w)**3)
      part6 = (2*CF*lms*(8*CA2 - 8*CA4 - 24*CA2*v + 32*CA4*v +
     &   32*CA2*v2 - 64*CA4*v2 - 16*CA2*v3 +72*CA4*v3 -
     &   48*CA4*v4 + 16*CA4*v5 - w - 2*CA2*w + 3*CA4*w -
     &        8*CA**3*CF*w + 2*v*w - 14*CA2*v*w + 4*CA4*v*w +
     &8*CA**3*CF*v*w - 2*v2*w + 52*CA2*v2*w - 42*CA4*v2*w +
     &        2*v3*w - 74*CA2*v3*w + 112*CA4*v3*w - v4*w +
     &        46*CA2*v4*w - 165*CA4*v4*w + 128*CA4*v5*w -
     &        48*CA4*v6*w + 2*w2 - 2*CA4*w2 - v*w2 -
     &        2*CA2*v*w2 + 11*CA4*v*w2 + 24*CA**3*CF*v*w2 -
     &        4*v2*w2 + 8*CA2*v2*w2 + 4*CA4*v2*w2 -
     &   16*CA**3*CF*v2*w2 + 8*v3*w2 - 46*CA2*v3*w2 -
     &        10*CA4*v3*w2 - 12*v4*w2 + 62*CA2*v4*w2 +
     &        6*CA4*v4*w2 + 7*v5*w2 - 58*CA2*v5*w2 +
     & 83*CA4*v5*w2 - 96*CA4*v6*w2 + 48*CA4*v7*w2 -
     &        2*w3 + 2*CA4*w3 - 2*v*w3 + 8*CA2*v*w3 -
     &        14*CA4*v*w3 + 5*v2*w3 + 10*CA2*v2*w3 -
     &7*CA4*v2*w3 - 24*CA**3*CF*v2*w3 + 4*CA2*v3*w3 -
     &        68*CA4*v3*w3 + 8*CA**3*CF*v3*w3 + 2*v4*w3 +
     &        22*CA2*v4*w3 + 72*CA4*v4*w3 + 4*v5*w3 +
     &        2*CA2*v5*w3 - 134*CA4*v5*w3 - 7*v6*w3 +
     & 34*CA2*v6*w3 + 37*CA4*v6*w3 - 16*CA4*v8*w3 +
     &  6*v*w4 - 6*CA4*v*w4 - 6*v2*w4 - 24*CA2*v2*w4 +
     &        54*CA4*v2*w4 + v3*w4 + 2*CA2*v3*w4 -
     &        75*CA4*v3*w4 + 8*CA**3*CF*v3*w4 - 5*v4*w4 -
     &        18*CA2*v4*w4 + 263*CA4*v4*w4 - 3*v5*w4 -
     &        24*CA2*v5*w4 - 149*CA4*v5*w4 + 4*v6*w4 -
     &  30*CA2*v6*w4 + 138*CA4*v6*w4 + 3*v7*w4 -
     &   6*CA2*v7*w4 - 13*CA4*v7*w4 + 16*CA4*v8*w4 -
     &        6*v2*w5 + 6*CA4*v2*w5 + 10*v3*w5 +
     &        24*CA2*v3*w5 - 58*CA4*v3*w5 - 8*v4*w5 -
     &        16*CA2*v4*w5 + 112*CA4*v4*w5 + 10*v5*w5 +
     &        22*CA2*v5*w5 - 352*CA4*v5*w5 - 2*v6*w5 +
     &        26*CA2*v6*w5 + 152*CA4*v6*w5 - 4*v7*w5 +
     & 8*CA2*v7*w5 - 100*CA4*v7*w5 - 16*CA4*v8*w5 +
     &        2*v3*w6 - 2*CA4*v3*w6 - 4*v4*w6 -
     &        8*CA2*v4*w6 + 20*CA4*v4*w6 + 4*v5*w6 +
     &        8*CA2*v5*w6 - 44*CA4*v5*w6 - 5*v6*w6 -
     &        10*CA2*v6*w6 + 207*CA4*v6*w6 + 3*v7*w6 -
     & 6*CA2*v7*w6 - 45*CA4*v7*w6 + 48*CA4*v8*w6 -
     & 64*CA4*v7*w7 - 16*CA4*v8*w7 + 16*CA4*v8*w8))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**3)
      part7 = -(CF*(16*CA2*v - 48*CA4*v - 80*CA2*v2 + 304*CA4*v2 +
     &        160*CA2*v3 - 832*CA4*v3 - 160*CA2*v4 +
     &        1280*CA4*v4 + 80*CA2*v5 - 1200*CA4*v5 -
     &        16*CA2*v6 + 688*CA4*v6 - 224*CA4*v7 +
     &        32*CA4*v8 - w + 6*CA2*w - 5*CA4*w + 7*v*w -
     &58*CA2*v*w + 83*CA4*v*w - 24*v2*w + 232*CA2*v2*w -
     &        464*CA4*v2*w + 50*v3*w - 452*CA2*v3*w +
     &        1170*CA4*v3*w - 65*v4*w + 438*CA2*v4*w -
     &        1397*CA4*v4*w + 51*v5*w - 178*CA2*v5*w +
     &        447*CA4*v5*w - 22*v6*w - 4*CA2*v6*w +
     &        794*CA4*v6*w + 4*v7*w + 16*CA2*v7*w -
     &        1044*CA4*v7*w + 512*CA4*v8*w - 96*CA4*v9*w -
     &        14*w2 + 8*CA2*w2 + 6*CA4*w2 + 72*v*w2 -
     &        22*CA2*v*w2 - 54*CA4*v*w2 - 161*v2*w2 +
     &        12*CA2*v2*w2 + 201*CA4*v2*w2 +186*v3*w2 -
     &        38*CA2*v3*w2 - 264*CA4*v3*w2 -74*v4*w2 +
     &        248*CA2*v4*w2 - 510*CA4*v4*w2 -80*v5*w2 -
     & 466*CA2*v5*w2 + 2342*CA4*v5*w2 + 119*v6*w2 +
     &        348*CA2*v6*w2 - 3471*CA4*v6*w2-58*v7*w2 -
     &        82*CA2*v7*w2 + 2400*CA4*v7*w2 +10*v8*w2 -
     &  8*CA2*v8*w2 - 554*CA4*v8*w2 - 192*CA4*v9*w2 +
     &96*CA4*v10*w2 + 16*w3 - 8*CA2*w3 - 8*CA4*w3 -
     &80*v*w3 + 12*CA2*v*w3 + 76*CA4*v*w3 + 127*v2*w3 +
     &        34*CA2*v2*w3 - 217*CA4*v2*w3 - 2*v3*w3 -
     &  74*CA2*v3*w3 + 176*CA4*v3*w3 - 315*v4*w3 +
     &      118*CA2*v4*w3 + 429*CA4*v4*w3 + 609*v5*w3 -
     &  414*CA2*v5*w3 - 1571*CA4*v5*w3 - 590*v6*w3 +
     &  824*CA2*v6*w3 + 2054*CA4*v6*w3 + 293*v7*w3 -
     &  756*CA2*v7*w3 - 661*CA4*v7*w3 - 58*v8*w3 +
     &        288*CA2*v8*w3 - 1094*CA4*v8*w3 -
     &        24*CA2*v9*w3 + 1104*CA4*v9*w3 -
     &  256*CA4*v10*w3 - 32*CA4*v11*w3 + 90*v2*w4 -
     &  48*CA2*v2*w4 - 42*CA4*v2*w4 - 366*v3*w4 +
     &  54*CA2*v3*w4 + 348*CA4*v3*w4 + 624*v4*w4 +
     &  124*CA2*v4*w4 - 908*CA4*v4*w4 - 658*v5*w4 -
     &  126*CA2*v5*w4 + 1404*CA4*v5*w4 + 501*v6*w4 -
     &  96*CA2*v6*w4 - 1441*CA4*v6*w4 - 222*v7*w4 -
     &  42*CA2*v7*w4 + 4*CA4*v7*w4 + 21*v8*w4 +
     &  352*CA2*v8*w4 + 1919*CA4*v8*w4 + 8*v9*w4 -
     &  254*CA2*v9*w4 - 1646*CA4*v9*w4 + 2*v10*w4 +
     &        36*CA2*v10*w4 + 202*CA4*v10*w4 +
     &  160*CA4*v11*w4 - 48*v2*w5 + 24*CA2*v2*w5 +
     &        24*CA4*v2*w5 + 192*v3*w5 - 12*CA2*v3*w5 -
     &      204*CA4*v3*w5 - 177*v4*w5 - 150*CA2*v4*w5 +
     &      471*CA4*v4*w5 - 113*v5*w5 + 30*CA2*v5*w5 -
     &      397*CA4*v5*w5 + 262*v6*w5 + 244*CA2*v6*w5 +
     &      126*CA4*v6*w5 - 171*v7*w5 + 76*CA2*v7*w5 +
     &      771*CA4*v7*w5 + 57*v8*w5 - 466*CA2*v8*w5 -
     &      2183*CA4*v8*w5 + 12*v9*w5 + 262*CA2*v9*w5 +
     &   1642*CA4*v9*w5 - 14*v10*w5 + 20*CA2*v10*w5 +
     &        138*CA4*v10*w5 - 12*CA2*v11*w5 -
     &      404*CA4*v11*w5 - 138*v4*w6 + 72*CA2*v4*w6 +
     &        66*CA4*v4*w6 + 420*v5*w6 + 6*CA2*v5*w6 -
     &      486*CA4*v5*w6 - 429*v6*w6 - 172*CA2*v6*w6 +
     &      941*CA4*v6*w6 + 156*v7*w6 - 66*CA2*v7*w6 -
     &      1102*CA4*v7*w6 + 35*v8*w6 + 220*CA2*v8*w6 +
     &      1581*CA4*v8*w6 - 76*v9*w6 - 28*CA2*v9*w6 -
     &   1016*CA4*v9*w6 + 32*v10*w6 - 120*CA2*v10*w6 -
     &        600*CA4*v10*w6 + 24*CA2*v11*w6 +
     &      680*CA4*v11*w6 + 48*v4*w7 - 24*CA2*v4*w7 -
     &        24*CA4*v4*w7 - 144*v5*w7 - 12*CA2*v5*w7 +
     &      180*CA4*v5*w7 + 73*v6*w7 + 126*CA2*v6*w7 -
     &      319*CA4*v6*w7 + 100*v7*w7 + 26*CA2*v7*w7 +
     &      102*CA4*v7*w7 - 137*v8*w7 - 82*CA2*v8*w7 -
     &        277*CA4*v8*w7 + 92*v9*w7 - 10*CA2*v9*w7 +
     &    202*CA4*v9*w7 - 32*v10*w7 + 84*CA2*v10*w7 +
     &        828*CA4*v10*w7 - 12*CA2*v11*w7 -
     &      820*CA4*v11*w7 + 62*v6*w8 - 32*CA2*v6*w8 -
     &        30*CA4*v6*w8 - 126*v7*w8 - 38*CA2*v7*w8 +
     &        192*CA4*v7*w8 + 94*v8*w8 + 36*CA2*v8*w8 -
     &        138*CA4*v8*w8 - 44*v9*w8 - 10*CA2*v9*w8 +
     &      82*CA4*v9*w8 + 14*v10*w8 - 20*CA2*v10*w8 -
     &   618*CA4*v10*w8 + 704*CA4*v11*w8 - 16*v6*w9 +
     &        8*CA2*v6*w9 + 8*CA4*v6*w9 + 32*v7*w9 +
     &        12*CA2*v7*w9 - 52*CA4*v7*w9 - 22*v8*w9 -
     &        16*CA2*v8*w9 + 70*CA4*v8*w9 + 8*v9*w9 +
     &        12*CA2*v9*w9 - 60*CA4*v9*w9 - 2*v10*w9 +
     &        242*CA4*v10*w9 - 416*CA4*v11*w9 -
     &        32*CA4*v10*w10 + 160*CA4*v11*w10 -
     &        32*CA4*v11*w11))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**3*(1 - v + v*w)**3)
      part8 = -(2*CF*lv*(8*CA2 - 32*CA4 - 48*CA2*v + 232*CA4*v +
     &        128*CA2*v2 - 768*CA4*v2 - 192*CA2*v3 +
     &        1504*CA4*v3 + 168*CA2*v4 - 1888*CA4*v4 -
     &        80*CA2*v5 + 1544*CA4*v5 + 16*CA2*v6 -
     &        800*CA4*v6 + 240*CA4*v7 - 32*CA4*v8 - 3*w +
     &  14*CA2*w - 11*CA4*w + 17*v*w - 68*CA2*v*w + 3*CA4*v*w -
     &   42*v2*w + 162*CA2*v2*w + 176*CA4*v2*w + 60*v3*w -
     &        268*CA2*v3*w - 496*CA4*v3*w - 55*v4*w +
     &        354*CA2*v4*w + 341*CA4*v4*w + 33*v5*w -
     &        348*CA2*v5*w + 651*CA4*v5*w - 12*v6*w +
     &        206*CA2*v6*w - 1530*CA4*v6*w + 2*v7*w -
     &        52*CA2*v7*w + 1330*CA4*v7*w - 560*CA4*v8*w +
     &        96*CA4*v9*w + 2*w2 - 2*CA4*w2 - 10*v*w2 -
     &        8*CA2*v*w2 + 26*CA4*v*w2 + 14*v2*w2 +
     &        22*CA2*v2*w2 - 76*CA4*v2*w2 + 12*v3*w2 +
     &        12*CA2*v3*w2 - 112*CA4*v3*w2 - 68*v4*w2 -
     &      96*CA2*v4*w2 + 1192*CA4*v4*w2 + 110*v5*w2 +
     &      100*CA2*v5*w2 - 3066*CA4*v5*w2 - 98*v6*w2 +
     &      42*CA2*v6*w2 + 3800*CA4*v6*w2 + 48*v7*w2 -
     &      136*CA2*v7*w2 - 2256*CA4*v7*w2 - 10*v8*w2 +
     &        64*CA2*v8*w2 + 334*CA4*v8*w2 +
     &        256*CA4*v9*w2 - 96*CA4*v10*w2 - 2*w3 +
     &    2*CA4*w3 + 10*v*w3 + 8*CA2*v*w3 - 26*CA4*v*w3 -
     &      7*v2*w3 - 74*CA2*v2*w3 + 137*CA4*v2*w3 -
     &      41*v3*w3 + 218*CA2*v3*w3 - 337*CA4*v3*w3 +
     &      120*v4*w3 - 444*CA2*v4*w3 + 300*CA4*v4*w3 -
     &      166*v5*w3 + 716*CA2*v5*w3 + 322*CA4*v5*w3 +
     &     129*v6*w3 - 818*CA2*v6*w3 - 471*CA4*v6*w3 -
     &       37*v7*w3 + 574*CA2*v7*w3 - 889*CA4*v7*w3 -
     &     16*v8*w3 - 160*CA2*v8*w3 + 1896*CA4*v8*w3 +
     &      10*v9*w3 - 20*CA2*v9*w3 - 1174*CA4*v9*w3 +
     &      208*CA4*v10*w3 + 32*CA4*v11*w3 - 12*v2*w4 +
     &      12*CA4*v2*w4 + 48*v3*w4 + 48*CA2*v3*w4 -
     &      144*CA4*v3*w4 - 70*v4*w4 - 144*CA2*v4*w4 +
     &      550*CA4*v4*w4 + 48*v5*w4 + 142*CA2*v5*w4 -
     &      1286*CA4*v5*w4 + 2*v6*w4 - 132*CA2*v6*w4 +
     &      1434*CA4*v6*w4 - 62*v7*w4 + 234*CA2*v7*w4 +
     &      196*CA4*v7*w4 + 68*v8*w4 - 340*CA2*v8*w4 -
     &      1652*CA4*v8*w4 - 18*v9*w4 + 216*CA2*v9*w4 +
     &      1010*CA4*v9*w4 - 4*v10*w4 - 16*CA2*v10*w4 +
     &       16*CA4*v10*w4 - 144*CA4*v11*w4 + 6*v2*w5 -
     &       6*CA4*v2*w5 - 24*v3*w5 - 24*CA2*v3*w5 +
     &        72*CA4*v3*w5 + 15*v4*w5 + 114*CA2*v4*w5 -
     &       273*CA4*v4*w5 + 37*v5*w5 - 148*CA2*v5*w5 +
     &       639*CA4*v5*w5 - 74*v6*w5 + 206*CA2*v6*w5 -
     &       772*CA4*v6*w5 + 88*v7*w5 - 346*CA2*v7*w5 -
     &       374*CA4*v7*w5 - 60*v8*w5 + 386*CA2*v8*w5 +
     &       1546*CA4*v8*w5 - 4*v9*w5 - 206*CA2*v9*w5 -
     &      862*CA4*v9*w5 + 16*v10*w5 - 30*CA2*v10*w5 -
     &        218*CA4*v10*w5 + 8*CA2*v11*w5 +
     &      296*CA4*v11*w5 + 18*v4*w6 - 18*CA4*v4*w6 -
     &        54*v5*w6 - 72*CA2*v5*w6 + 198*CA4*v5*w6 +
     &       58*v6*w6 + 122*CA2*v6*w6 - 628*CA4*v6*w6 -
     &      42*v7*w6 - 36*CA2*v7*w6 + 1486*CA4*v7*w6 +
     &      12*v8*w6 - 26*CA2*v8*w6 - 1934*CA4*v8*w6 +
     &        36*v9*w6 + 14*CA2*v9*w6 + 926*CA4*v9*w6 -
     &    28*v10*w6 + 102*CA2*v10*w6 + 290*CA4*v10*w6 -
     &     16*CA2*v11*w6 - 464*CA4*v11*w6 - 6*v4*w7 +
     &        6*CA4*v4*w7 + 18*v5*w7 + 24*CA2*v5*w7 -
     &        66*CA4*v5*w7 - 9*v6*w7 - 62*CA2*v6*w7 +
     &        191*CA4*v6*w7 - 7*v7*w7 + 22*CA2*v7*w7 -
     &        391*CA4*v7*w7 + 20*v8*w7 - 12*CA2*v8*w7 +
     &        472*CA4*v8*w7 - 44*v9*w7 - 6*CA2*v9*w7 -
     &      30*CA4*v9*w7 + 28*v10*w7 - 78*CA2*v10*w7 -
     &        526*CA4*v10*w7 + 8*CA2*v11*w7 +
     &        616*CA4*v11*w7 - 8*v6*w8 + 8*CA4*v6*w8 +
     &        16*v7*w8 + 32*CA2*v7*w8 - 80*CA4*v7*w8 -
     &        18*v8*w8 - 24*CA2*v8*w8 + 186*CA4*v8*w8 +
     &        26*v9*w8 + 26*CA2*v9*w8 - 316*CA4*v9*w8 -
     &     16*v10*w8 + 30*CA2*v10*w8 + 434*CA4*v10*w8 -
     &        576*CA4*v11*w8 + 2*v6*w9 - 2*CA4*v6*w9 -
     &        4*v7*w9 - 8*CA2*v7*w9 + 20*CA4*v7*w9 +
     &        4*v8*w9 + 8*CA2*v8*w9 - 44*CA4*v8*w9 -
     &        6*v9*w9 - 8*CA2*v9*w9 + 70*CA4*v9*w9 +
     &      4*v10*w9 - 8*CA2*v10*w9 - 108*CA4*v10*w9 +
     &        352*CA4*v11*w9 - 144*CA4*v11*w10 +
     &        32*CA4*v11*w11))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**3*(1 - v + v*w)**3)
      part9 = -(2*CF*l1w*(8*CA2 - 24*CA4 - 48*CA2*v + 176*CA4*v +
     &        128*CA2*v2 - 592*CA4*v2 - 192*CA2*v3 +
     &        1184*CA4*v3 + 168*CA2*v4 - 1528*CA4*v4 -
     &        80*CA2*v5 + 1296*CA4*v5 + 16*CA2*v6 -
     &        704*CA4*v6 + 224*CA4*v7 - 32*CA4*v8 - 3*w +
     &        14*CA2*w - 11*CA4*w + 19*v*w - 66*CA2*v*w +
     &  21*CA4*v*w - 50*v2*w + 154*CA2*v2*w + 72*CA4*v2*w +
     &    74*v3*w - 254*CA2*v3*w - 254*CA4*v3*w - 71*v4*w +
     &        338*CA2*v4*w + 85*CA4*v4*w + 47*v5*w -
     &        334*CA2*v5*w + 705*CA4*v5*w - 20*v6*w +
     &        198*CA2*v6*w - 1394*CA4*v6*w + 4*v7*w -
     &        50*CA2*v7*w + 1208*CA4*v7*w - 528*CA4*v8*w +
     &        96*CA4*v9*w + 2*w2 - 2*CA4*w2 - 10*v*w2 -
     &        8*CA2*v*w2 + 14*CA4*v*w2 + 14*v2*w2 +
     &        26*CA2*v2*w2 - 18*CA4*v2*w2 + 14*v3*w2 -
     &        2*CA2*v3*w2 - 188*CA4*v3*w2 - 70*v4*w2 -
     &      70*CA2*v4*w2 + 1058*CA4*v4*w2 + 110*v5*w2 +
     &      68*CA2*v5*w2 - 2494*CA4*v5*w2 - 106*v6*w2 +
     &      62*CA2*v6*w2 + 3030*CA4*v6*w2 + 62*v7*w2 -
     &      138*CA2*v7*w2 - 1780*CA4*v7*w2 - 16*v8*w2 +
     &        62*CA2*v8*w2 + 220*CA4*v8*w2 +
     &        256*CA4*v9*w2 - 96*CA4*v10*w2 - 2*w3 +
     &    2*CA4*w3 + 10*v*w3 + 8*CA2*v*w3 - 26*CA4*v*w3 -
     &        7*v2*w3 - 74*CA2*v2*w3 + 137*CA4*v2*w3 -
     &      45*v3*w3 + 218*CA2*v3*w3 - 385*CA4*v3*w3 +
     &      130*v4*w3 - 446*CA2*v4*w3 + 592*CA4*v4*w3 -
     &      184*v5*w3 + 722*CA2*v5*w3 - 412*CA4*v5*w3 +
     &      169*v6*w3 - 810*CA2*v6*w3 + 409*CA4*v6*w3 -
     &      77*v7*w3 + 546*CA2*v7*w3 - 1309*CA4*v7*w3 -
     &      10*v8*w3 - 142*CA2*v8*w3 + 1828*CA4*v8*w3 +
     &      16*v9*w3 - 22*CA2*v9*w3 - 1044*CA4*v9*w3 +
     &      176*CA4*v10*w3 + 32*CA4*v11*w3 - 12*v2*w4 +
     &        12*CA4*v2*w4 + 48*v3*w4 + 48*CA2*v3*w4 -
     &     124*CA4*v3*w4 - 70*v4*w4 - 148*CA2*v4*w4 +
     &     436*CA4*v4*w4 + 54*v5*w4 + 156*CA2*v5*w4 -
     &     1000*CA4*v5*w4 - 28*v6*w4 - 174*CA2*v6*w4 +
     &     1210*CA4*v6*w4 - 42*v7*w4 + 302*CA2*v7*w4 -
     &      90*CA4*v7*w4 + 96*v8*w4 - 380*CA2*v8*w4 -
     &     1004*CA4*v8*w4 - 40*v9*w4 + 218*CA2*v9*w4 +
     &     574*CA4*v9*w4 - 6*v10*w4 - 14*CA2*v10*w4 +
     &     106*CA4*v10*w4 - 128*CA4*v11*w4 + 6*v2*w5 -
     &        6*CA4*v2*w5 - 24*v3*w5 - 24*CA2*v3*w5 +
     &        72*CA4*v3*w5 + 15*v4*w5 + 114*CA2*v4*w5 -
     &       273*CA4*v4*w5 + 37*v5*w5 - 152*CA2*v5*w5 +
     &       623*CA4*v5*w5 - 68*v6*w5 + 232*CA2*v6*w5 -
     &      866*CA4*v6*w5 + 106*v7*w5 - 400*CA2*v7*w5 +
     &      136*CA4*v7*w5 - 120*v8*w5 + 414*CA2*v8*w5 +
     &      800*CA4*v8*w5 + 22*v9*w5 - 192*CA2*v9*w5 -
     &      424*CA4*v9*w5 + 26*v10*w5 - 40*CA2*v10*w5 -
     &      238*CA4*v10*w5 + 8*CA2*v11*w5 +
     &      216*CA4*v11*w5 + 18*v4*w6 - 18*CA4*v4*w6 -
     &        54*v5*w6 - 72*CA2*v5*w6 + 210*CA4*v5*w6 +
     &       58*v6*w6 + 118*CA2*v6*w6 - 622*CA4*v6*w6 -
     &       60*v7*w6 - 22*CA2*v7*w6 + 1262*CA4*v7*w6 +
     &       54*v8*w6 - 20*CA2*v8*w6 - 1588*CA4*v8*w6 +
     &        32*v9*w6 - 22*CA2*v9*w6 + 800*CA4*v9*w6 -
     & 48*v10*w6 + 122*CA2*v10*w6 + 156*CA4*v10*w6 -
     &      16*CA2*v11*w6 - 304*CA4*v11*w6 - 6*v4*w7 +
     &        6*CA4*v4*w7 + 18*v5*w7 + 24*CA2*v5*w7 -
     &        66*CA4*v5*w7 - 9*v6*w7 - 62*CA2*v6*w7 +
     &        191*CA4*v6*w7 - 3*v7*w7 + 22*CA2*v7*w7 -
     &        287*CA4*v7*w7 + 10*v8*w7 - 26*CA2*v8*w7 +
     &        320*CA4*v8*w7 - 58*v9*w7 + 28*CA2*v9*w7 -
     &      48*CA4*v9*w7 + 48*v10*w7 - 98*CA2*v10*w7 -
     &        372*CA4*v10*w7 + 8*CA2*v11*w7 +
     &        440*CA4*v11*w7 - 8*v6*w8 + 8*CA4*v6*w8 +
     &        16*v7*w8 + 32*CA2*v7*w8 - 116*CA4*v7*w8 -
     &        18*v8*w8 - 20*CA2*v8*w8 + 244*CA4*v8*w8 +
     &        36*v9*w8 + 12*CA2*v9*w8 - 262*CA4*v9*w8 -
     &   26*v10*w8 + 40*CA2*v10*w8 + 334*CA4*v10*w8 -
     &        448*CA4*v11*w8 + 2*v6*w9 - 2*CA4*v6*w9 -
     &        4*v7*w9 - 8*CA2*v7*w9 + 20*CA4*v7*w9 +
     &        4*v8*w9 + 8*CA2*v8*w9 - 44*CA4*v8*w9 -
     &        8*v9*w9 - 6*CA2*v9*w9 + 12*CA4*v9*w9 +
     &       6*v10*w9 - 10*CA2*v10*w9 - 50*CA4*v10*w9 +
     &        288*CA4*v11*w9 + 16*CA4*v9*w10 -
     &        16*CA4*v10*w10 - 128*CA4*v11*w10 +
     &        32*CA4*v11*w11))/
     &    (CA*(1 - v)**2*v2*w2*(1 - v*w)**3*(1 - v + v*w)**3)
      struv16= part1 + part2 + part3 + part4 +
     &         part5 + part6 + part7 + part8 + 
     &         part9
      RETURN
      END
C----------- LES COMBINAISONS DE FONCTION DE STRUCTURE
C*===================================================================
C
C                 STRUCTURE FUNCTION AT SCALE Q2
C
C==================================================================  */
      SUBROUTINE STRU(XUHA,XUBHA,XDHA,XDBHA,XSHA,XCHA,XGPROA,
     #                XUHB,XUBHB,XDHB,XDBHB,XSHB,XCHB,XGPROB,
     #                XDUP,XDUBP,XDDP,XDDBP,XDSP,XDSBP,XDCP,XDCBP,XDGP,
     #                GPPV,GPPC)
      IMPLICIT REAL*8(A-H,M-Z)
      DIMENSION GPPV(16),GPPC(16)
C
      GPPV(1)=XUHA*(XDHB+XSHB+XCHB)*XDUP+XDHA*(XUHB+XSHB+XCHB)*XDDP+
     &        XSHA*(XUHB+XDHB+XCHB)*XDSP+XCHA*(XUHB+XDHB+XSHB)*XDCP+
     & XUBHA*(XDBHB+XSHB+XCHB)*XDUBP+XDBHA*(XUBHB+XSHB+XCHB)*XDDBP+
     & XSHA*(XUBHB+XDBHB+XCHB)*XDSBP+XCHA*(XUBHB+XDBHB+XSHB)*XDCBP
      GPPC(1)=XUHB*(XDHA+XSHA+XCHA)*XDUP+XDHB*(XUHA+XSHA+XCHA)*XDDP+
     &        XSHB*(XUHA+XDHA+XCHA)*XDSP+XCHB*(XUHA+XDHA+XSHA)*XDCP+
     & XUBHB*(XDBHA+XSHA+XCHA)*XDUBP+XDBHB*(XUBHA+XSHA+XCHA)*XDDBP+
     & XSHB*(XUBHA+XDBHA+XCHA)*XDSBP+XCHB*(XUBHA+XDBHA+XSHA)*XDCBP
C*********************************************************************
      GPPV(2)=(XUHA*(XDHB+XSHB+XCHB)+XDHA*(XSHB+XCHB)+XSHA*XCHB+
     &      XUBHA*(XDBHB+XSHB+XCHB)+XDBHA*(XSHB+XCHB)+XSHA*XCHB)*XDGP
      GPPC(2)=(XUHB*(XDHA+XSHA+XCHA)+XDHB*(XSHA+XCHA)+XSHB*XCHA+
     &      XUBHB*(XDBHA+XSHA+XCHA)+XDBHB*(XSHA+XCHA)+XSHB*XCHA)*XDGP
C*********************************************************************
      GPPV(3)=XUHA*(XDBHB+XSHB+XCHB)*XDUP+XDHA*(XUBHB+XSHB+XCHB)*XDDP+
     &      XSHA*(XUBHB+XDBHB+XCHB)*XDSP+XCHA*(XUBHB+XDBHB+XSHB)*XDCP+
     &      XUBHA*(XDHB+XSHB+XCHB)*XDUBP+XDBHA*(XUHB+XSHB+XCHB)*XDDBP+
     &      XSHA*(XUHB+XDHB+XCHB)*XDSBP+XCHA*(XUHB+XDHB+XSHB)*XDCBP
      GPPC(3)=XUHB*(XDBHA+XSHA+XCHA)*XDUP+XDHB*(XUBHA+XSHA+XCHA)*XDDP+
     &      XSHB*(XUBHA+XDBHA+XCHA)*XDSP+XCHB*(XUBHA+XDBHA+XSHA)*XDCP+
     &      XUBHB*(XDHA+XSHA+XCHA)*XDUBP+XDBHB*(XUHA+XSHA+XCHA)*XDDBP+
     &      XSHB*(XUHA+XDHA+XCHA)*XDSBP+XCHB*(XUHA+XDHA+XSHA)*XDCBP
C*********************************************************************
      GPPV(4)=(XUHA*(XDBHB+XSHB+XCHB)+XDHA*(XUBHB+XSHB+XCHB)+
     &         XSHA*(XUBHB+XDBHB+XCHB)+XCHA*(XUBHB+XDBHB+XSHB))*XDGP
      GPPC(4)=(XUHB*(XDBHA+XSHA+XCHA)+XDHB*(XUBHA+XSHA+XCHA)+
     &         XSHB*(XUBHA+XDBHA+XCHA)+XCHB*(XUBHA+XDBHA+XSHA))*XDGP
C*********************************************************************
      GPPV(5)=(XDHA*XDBHB+XSHA*XSHB+XCHA*XCHB)*XDUP+
     &        (XUHA*XUBHB+XSHA*XSHB+XCHA*XCHB)*XDDP+
     &        (XUHA*XUBHB+XDHA*XDBHB+XCHA*XCHB)*XDSP+
     &        (XUHA*XUBHB+XDHA*XDBHB+XSHA*XSHB)*XDCP+
     &        (XDBHA*XDHB+XSHA*XSHB+XCHA*XCHB)*XDUBP+
     &        (XUBHA*XUHB+XSHA*XSHB+XCHA*XCHB)*XDDBP+
     &        (XUBHA*XUHB+XDBHA*XDHB+XCHA*XCHB)*XDSBP+
     &        (XUBHA*XUHB+XDBHA*XDHB+XSHA*XSHB)*XDCBP
      GPPC(5)=(XDHB*XDBHA+XSHB*XSHA+XCHB*XCHA)*XDUP+
     &        (XUHB*XUBHA+XSHB*XSHA+XCHB*XCHA)*XDDP+
     &        (XUHB*XUBHA+XDHB*XDBHA+XCHB*XCHA)*XDSP+
     &        (XUHB*XUBHA+XDHB*XDBHA+XSHB*XSHA)*XDCP+
     &        (XDBHB*XDHA+XSHB*XSHA+XCHB*XCHA)*XDUBP+
     &        (XUBHB*XUHA+XSHB*XSHA+XCHB*XCHA)*XDDBP+
     &        (XUBHB*XUHA+XDBHB*XDHA+XCHB*XCHA)*XDSBP+
     &        (XUBHB*XUHA+XDBHB*XDHA+XSHB*XSHA)*XDCBP
C*********************************************************************
      GPPV(6)=XUHA*XUHB*XDUP+XDHA*XDHB*XDDP+XSHA*XSHB*XDSP
     &       +XCHA*XCHB*XDCP
     &       +XUBHA*XUBHB*XDUBP+XDBHA*XDBHB*XDDBP+XSHA*XSHB*XDSBP
     &       +XCHA*XCHB*XDCBP
      GPPV(6)=GPPV(6)/2.D0
      GPPC(6)=XUHB*XUHA*XDUP+XDHB*XDHA*XDDP+XSHB*XSHA*XDSP
     &+XCHB*XCHA*XDCP
     &+XUBHB*XUBHA*XDUBP+XDBHB*XDBHA*XDDBP+XSHB*XSHA*XDSBP
     &+XCHB*XCHA*XDCBP
      GPPC(6)=GPPC(6)/2.D0
C*********************************************************************
      GPPV(7)=(XUHA*XUHB+XDHA*XDHB+XSHA*XSHB+XCHA*XCHB
     &      +XUBHA*XUBHB+XDBHA*XDBHB+XSHA*XSHB+XCHA*XCHB)*XDGP
      GPPV(7)=GPPV(7)/2.D0
      GPPC(7)=(XUHB*XUHA+XDHB*XDHA+XSHB*XSHA+XCHB*XCHA
     &      +XUBHB*XUBHA+XDBHB*XDBHA+XSHB*XSHA+XCHB*XCHA)*XDGP
      GPPC(7)=GPPC(7)/2.D0
C*********************************************************************
      GPPV(8)=((XDHA+XSHA+XCHA)*XDUP+(XUHA+XSHA+XCHA)*XDDP+
     &         (XUHA+XDHA+XCHA)*XDSP+(XUHA+XDHA+XSHA)*XDCP+
     &       (XDBHA+XSHA+XCHA)*XDUBP+(XUBHA+XSHA+XCHA)*XDDBP+
     &       (XUBHA+XDBHA+XCHA)*XDSBP+(XUBHA+XDBHA+XSHA)*XDCBP)*XGPROB
      GPPC(8)=((XDHB+XSHB+XCHB)*XDUP+(XUHB+XSHB+XCHB)*XDDP+
     &         (XUHB+XDHB+XCHB)*XDSP+(XUHB+XDHB+XSHB)*XDCP+
     &       (XDBHB+XSHB+XCHB)*XDUBP+(XUBHB+XSHB+XCHB)*XDDBP+
     &       (XUBHB+XDBHB+XCHB)*XDSBP+(XUBHB+XDBHB+XSHB)*XDCBP)*XGPROA
C*********************************************************************
      GPPV(9)=((XDHA+XSHA+XCHA)*XDUBP+(XUHA+XSHA+XCHA)*XDDBP+
     &         (XUHA+XDHA+XCHA)*XDSBP+(XUHA+XDHA+XSHA)*XDCBP+
     &         (XDBHA+XSHA+XCHA)*XDUP+(XUBHA+XSHA+XCHA)*XDDP+
     &         (XUBHA+XDBHA+XCHA)*XDSP+(XUBHA+XDBHA+XSHA)*XDCP)*XGPROB
      GPPC(9)=((XDHB+XSHB+XCHB)*XDUBP+(XUHB+XSHB+XCHB)*XDDBP+
     &         (XUHB+XDHB+XCHB)*XDSBP+(XUHB+XDHB+XSHB)*XDCBP+
     &         (XDBHB+XSHB+XCHB)*XDUP+(XUBHB+XSHB+XCHB)*XDDP+
     &         (XUBHB+XDBHB+XCHB)*XDSP+(XUBHB+XDBHB+XSHB)*XDCP)*XGPROA
C*********************************************************************
      GPPV(10)=(XUHA*XDUBP+XDHA*XDDBP+XSHA*XDSBP+XCHA*XDCBP+
     &          XUBHA*XDUP+XDBHA*XDDP+XSHA*XDSP+XCHA*XDCP)*XGPROB
      GPPC(10)=(XUHB*XDUBP+XDHB*XDDBP+XSHB*XDSBP+XCHB*XDCBP+
     &          XUBHB*XDUP+XDBHB*XDDP+XSHB*XDSP+XCHB*XDCP)*XGPROA
C*********************************************************************
      GPPV(11)=XUHA*XUBHB*XDUP+XDHA*XDBHB*XDDP+XSHA*XSHB*XDSP+
     &         XCHA*XCHB*XDCP+XUBHA*XUHB*XDUBP+XDBHA*XDHB*XDDBP+
     &         XSHA*XSHB*XDSBP+XCHA*XCHB*XDCBP
      GPPC(11)=XUHB*XUBHA*XDUP+XDHB*XDBHA*XDDP+XSHB*XSHA*XDSP+
     &         XCHB*XCHA*XDCP+XUBHB*XUHA*XDUBP+XDBHB*XDHA*XDDBP+
     &         XSHB*XSHA*XDSBP+XCHB*XCHA*XDCBP
C*********************************************************************
      GPPV(12)=(XUHA*XUBHB+XDHA*XDBHB+XSHA*XSHB+XCHA*XCHB)*XDGP
      GPPC(12)=(XUHB*XUBHA+XDHB*XDBHA+XSHB*XSHA+XCHB*XCHA)*XDGP
C*********************************************************************
      GPPV(13)=(XUHA*XDUP+XUBHA*XDUBP+XDHA*XDDP+XDBHA*XDDBP+
     &          XSHA*XDSP+XSHA*XDSBP+XCHA*XDCP+XCHA*XDCBP)*XGPROB
      GPPC(13)=(XUHB*XDUP+XUBHB*XDUBP+XDHB*XDDP+XDBHB*XDDBP+
     &          XSHB*XDSP+XSHB*XDSBP+XCHB*XDCP+XCHB*XDCBP)*XGPROA
C*********************************************************************
      GPPV(14)=(XUHA+XUBHA+XDHA+XDBHA+2.*XSHA+2.*XCHA)*XGPROB*XDGP
      GPPC(14)=(XUHB+XUBHB+XDHB+XDBHB+2.*XSHB+2.*XCHB)*XGPROA*XDGP
C*********************************************************************
      GPPV(15)=XGPROA*XGPROB*XDGP/2.
      GPPC(15)=XGPROB*XGPROA*XDGP/2.
C*********************************************************************
      GPPV(16)=XGPROA*XGPROB*(XDUP+XDUBP+XDDP+XDDBP+
     &         XDSP+XDSBP+XDCP+XDCBP)/2.
      GPPC(16)=XGPROB*XGPROA*(XDUP+XDUBP+XDDP+XDDBP+
     &         XDSP+XDSBP+XDCP+XDCBP)/2.
C
      RETURN
      END
