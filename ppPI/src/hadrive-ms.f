********************************************************
*     INCLUSIVE HADRON PRODUCTION AT O(ALPHA_S**3)     *
*       IN UNPOLARIZED HADRON-HADRON COLLISIONS        *
********************************************************
*
* STATUS: LAST CHANGES 24/OCT/07 (MS)
*         UPDATE by E.R. Nocera (NNPDF Collaboration)
*         for handling NNPDF parton sets SEP 2013
*
      PROGRAM HADRON
*
      IMPLICIT DOUBLE PRECISION (A-Z)
*
      INTEGER IIP, IIF, IPT, IY, IORD, IH1, IH2, ISIGM, JMAR, FINI
      INTEGER K1, J0
      INTEGER nmaxpt
*
      INTEGER IREP, FFREP, NPT, LREP, PDFREP
      INTEGER IDUM
      parameter(nmaxpt=300)
*
      double precision ptpt(nmaxpt)
*
      character*30 nameexp
      character*100 datafile
      character*100 theoryfile
      character*100 outfile
      character*10 cdum
*
      DIMENSION CC(16)
*
      COMMON / ORDER / IORD
      COMMON / IPART / IIP,IIF
      COMMON / IPTFL / IPT
      COMMON / IYFL / IY
      COMMON / HADR / IH1,IH2
      COMMON / XSECT / ISIGM
      COMMON / FRAGINI / FINI
      COMMON / KINVAR / SQS, PT, ETA
      COMMON / PTINT/ PTDO,PTUP
      COMMON / YINT / YDO,YUP
      COMMON / VERROR/ VERR
      COMMON / SCALESF / SCFAC,SCMU,SCFRAG
      COMMON / PARAMS / CF,CA,NF,PI,HC2
      COMMON / PREFAC / CC
      COMMON /CONTROL/ ZMIN
      COMMON /FRAGCUT/ ZCUT
*
************
* PARAMETERS
************
      PI=DACOS(-1.D0)
      CF=4.D0/3.D0
      CA=3.D0
      HC2=.38942957D0           !VALUE OF THE CONVERSION FACTOR: (HBAR*C)**2
                                !results in mbarn
*
      VC=CA**2-1.D0
      NC=CA
*
      DO 2 J0=1,16
         IF (J0.EQ.16.OR.J0.EQ.15) THEN
            CC(J0)=8.D0*VC**2
         ELSE IF (J0.EQ.14.OR.J0.EQ.13.OR.J0.EQ.10.OR.J0.EQ.9.OR.
     &            J0.EQ.8) THEN
            CC(J0)=8.D0*VC*NC
         ELSE
            CC(J0)=8.D0*NC**2
         ENDIF
 2    CONTINUE 
*      
      FINI = 0
      JMAR=2  ! FACT. SCHEME DEP. FOR FRAGM. PIECE (AVERSA), JMAR=2 FOR MSBAR
      
**********************
* REQUIRED USER INPUTS
**********************

      write(*,*) "Enter name of experiment"
      read(5,*) nameexp

      write(*,*) "Insert FF replica"
      read(5,*) FFREP

      write(*,*) "Insert PDF replica"
      read(5,*) PDFREP

      lrep=ffrep+1000

      write(datafile,701) "../data/",trim(nameexp),"/",
     1     trim(nameexp),".txt"
      write(theoryfile,701) "../theory/",trim(nameexp),"/",
     1     trim(nameexp),".txt"

 701  format(a,a,a,a,a)
      
      open(unit=10, status="old", file=datafile)
      read(10,*) cdum
      read(10,*) npt
      read(10,*) sqs
      read(10,*) ydo, yup
      read(10,*) cdum
      ydum = (yup - ydo)/2d0
      do ipt=1, npt
         read(10,*) idum, ptpt(ipt)
      enddo
      close(10)

      open(unit=20, status="old", file=theoryfile)
      read(20,*) ipt
      read(20,*) iy
      read(20,*) iord
      read(20,*) iip
      read(20,*) iif
      read(20,*) ih1, ih2
      read(20,*) isigm
      read(20,*) scfac
      close(20)
      
      nf = 5.D0     ! number of active flavours - do not change
      SCMU=SCFAC    ! do not change
      SCFRAG=SCFAC  ! do not change
      
***********************
* OUTPUT OF FLAGS, ETC.
***********************
*
      write(outfile,300) "../theory/",trim(nameexp),
     1     "/results/",iord, "_",iip,"_",iif,"_",scfac,
     1     "_",lrep,".res"

 300  format(a,a,a,i1,a,i1,a,i1,a,f3.1,a,i4,a)

      IF(IIP.EQ.0) THEN 
         WRITE(6,100) 
      ELSEIF(IIP.EQ.1) THEN
         WRITE(6,101)
      ELSEIF(IIP.EQ.2) THEN
         WRITE(6,102)
         IF(IORD.EQ.0) THEN
            CALL SetCtq6(4)     !LO
         ELSE
            CALL SetCtq6(1)     !NLO
         ENDIF
*
      ELSEIF(IIP.EQ.3) THEN
         WRITE(6,103)
         if(iord.eq.0)then
            CALL INITPDFSETBYNAMEM(4,"NNPDF31_lo_as_0118")
         elseif(iord.eq.1)then
            CALL INITPDFSETBYNAMEM(4,"NNPDF31_nlo_as_0118")
         endif

         IREP = PDFREP
         CALL INITPDFM(4,IREP)
      ELSEIF(IIP.EQ.4) THEN
         WRITE(6,105)
         CALL INITPDFSETBYNAMEM(4,"MSTW2008nlo68cl")
         !WRITE(6,104)
         !READ(*,*) IREP
         IREP = 0
         CALL INITPDFM(4,IREP)
      ENDIF

 100  FORMAT(4X,'CHOICE OF PDFs: GRV-98',/)
 101  FORMAT(4X,'CHOICE OF PDFs: CTEQ-5',/)
 102  FORMAT(4X,'CHOICE OF PDFs: CTEQ-6',/)
 103  FORMAT(4X,'CHOICE OF PDFs: NNPDF31',/)
 105  FORMAT(4X,'CHOICE OF PDFs: MSTW2008',/)
 104  FORMAT(4X,'INSERT NNPDF REPLICA NUMBER AND WAIT',/)
*
      IF(IIF.EQ.0) THEN
         WRITE(6,110)
      ELSEIF(IIF.EQ.1) THEN
         WRITE(6,111)
      ELSEIF(IIF.EQ.2) THEN
         WRITE(6,112)
         CALL INITPDFSETBYNAMEM(1,"DSS14_NLO_Pip")
         CALL INITPDFSETBYNAMEM(2,"DSS14_NLO_Pim")
         IREP = 0
         CALL INITPDFM(1,IREP)
         CALL INITPDFM(2,IREP)
      ELSEIF(IIF.EQ.3) THEN
         WRITE(6,113)
         CALL INITPDFSETBYNAMEM(1,"MAPFF10NLOPIp")
         CALL INITPDFSETBYNAMEM(2,"MAPFF10NLOPIm")
         IREP = FFREP
         CALL INITPDFM(1,IREP)
         CALL INITPDFM(2,IREP)
      ELSEIF(IIF.EQ.4) THEN
         WRITE(6,113)
         call initpdfsetbyname("MAPFF10NLOPIsum")          
         IREP=FFREP
         CALL INITPDF(IREP)
      ELSEIF(IIF.EQ.5) THEN
         WRITE(6,113)
         CALL INITPDFSETBYNAME("DSS07_NLO_HadronSum")
         IREP=0
         CALL INITPDF(IREP)
      ELSEIF(IIF.EQ.6) THEN
         WRITE(6,113)
         CALL INITPDFSETBYNAME("Ktz_NLO_HadSum")
         IREP=0
         CALL INITPDF(IREP)
      ELSEIF(IIF.EQ.7) THEN
         WRITE(6,114)
         CALL INITPDFSETBYNAMEM(1,"testff_p")
         CALL INITPDFSETBYNAMEM(2,"testff_m")
         IREP = FFREP
         CALL INITPDFM(1,IREP)
         CALL INITPDFM(2,IREP)
      ENDIF
*
 110  FORMAT(4X,'CHOICE OF FFs: DSS')
 111  FORMAT(4X,'CHOICE OF FFs: KRETZER')
 112  FORMAT(4X,'CHOICE OF FFs: DSS14 (LHAPDF)')
 113  FORMAT(4X,'CHOICE OF FFs: MAPFF1.0 (LHAPDF)')
 114  FORMAT(4X,'CHOICE OF FFs: testff')
*
***********************
* MAIN LOOP (PT-VALUES)
***********************
      OPEN(7,FILE=outfile)
*
CC      DO 1 PTDUM=1.25,17.25,0.5
c      DO 1 K1 = 1, 21, 1
      DO 1 K1 = 1, NPT, 1
cc      DO 1 K1 = 1, 12, 1
cc        DO 1 K1 = 1,4,1
*
c         ZMIN=1.D0
c         ZCUT=0.0D0
c         EPION=27.5D0+DBLE(K1)*2.5D0
c         PTDUM=EPION/DCOSH(YDUM)
c         PTDUM=PTB(K1)
c        PTDO = PTB(K1)
c        PTUP = PTB(K1+1)
c        PTDUM = (PTDO + PTUP)/2.D0
c        PTDIFF = PTUP - PTDO
c         PTDO = PTA(K1)
c         PTUP = PTA(K1+1)
c         PTDUM = (PTDO + PTUP)/2.D0
c         PTDIFF = PTUP - PTDO
c         PTDUM = PTB(K1)
*
c         PT = PTF(K1)
c         PT = PTALICE(K1)
         PT = ptpt(K1)
         ETA = YDUM
*************************************
*     CALCULATE CROSS SECTION PER BIN
*************************************     
         TH=FUNCDG(PTDUM,PTDO,PTUP,YDUM,YDO,YUP,IPT,IY) 
************
*     OUTPUT
************
         if(trim(nameexp).eq."UA2540")then
            WRITE(6,10) PT,4d0*TH !/1.D9/0.7D0
            WRITE(7,10) PT,4d0*TH !/1.D9/0.7D0
         else
            WRITE(6,10) PT,1d0*TH !/1.D9/0.7D0
            WRITE(7,10) PT,1d0*TH !/1.D9/0.7D0
         endif

c         WRITE(6,10) PT,TH/1.D9,ZMIN,ZCUT!/1.D9 !/30.D0
c         WRITE(7,10) PT,zmin,zcut,TH/1.D9/2.d0  !/0.7D0 
c         WRITE(6,10) YDO,YUP,PTDUM,TH,TH/PTDIFF
c         WRITE(7,10) YDO,YUP,PTDUM,TH,TH/PTDIFF
c         WRITE(6,10) PTDUM,TH,TH/PTDIFF,TH/PTDIFF/2.D0 
c         WRITE(7,10) PTDUM,TH,TH/PTDIFF,TH/PTDIFF/2.D0
 10      FORMAT(6(1PE14.6))   
*
 1    CONTINUE
*
      STOP
      END
*
***********************************************************
      FUNCTION FUNCDG (PT,PTDO,PTUP,ETA,ETADO,ETAUP,IPT,IY)
***********************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,L-Z)
*
      EXTERNAL DPLUS
*
      COMMON /RESULT/ERG1,ERG2,ERG3,ERG4
      COMMON /VERROR/VERR
*
      CALL VEGAS(DPLUS,1.D-9,5,50000,5,0,0)
      FUNCDG = ERG1
      VERR=ERG2
*     
      RETURN
      END
*
*************************
      FUNCTION ALPHAS(Q2)
*************************
*
*   RETURNS ALPHAS(Q**2) DIVIDED BY 4*PI !!
*      
      IMPLICIT DOUBLE PRECISION (A - Z)
*
      INTEGER NF, K, IIP, IIF, IORD
      DIMENSION LAMBDA (3:5), Q2THR (2)
*
      DOUBLE PRECISION ALPHANNPDF, ALPHAMZ, PI
      EXTERNAL ALPHAMZ
*
      COMMON /IPART/ IIP,IIF
      COMMON /ORDER/ IORD
*
********************************
*   SETUP HEAVY QUARK THRESHOLDS
********************************
      IF(IIP.EQ.0.OR.IIP.EQ.1.OR.IIP.EQ.2)THEN
*
         IF(IIP.EQ.0) THEN      !GRV-98
            Q2THR(1) = 1.96D0
            Q2THR(2) = 20.25D0
            LAMBDA(3)=0.2994D0
            LAMBDA(4)=0.246D0
            LAMBDA(5)=0.1677D0
         ELSE IF((IIP.EQ.1).OR.(IIP.EQ.2)) THEN !CTEQ-6M
            Q2THR(1) = 0.D0
            Q2THR(2) = 20.25D0
            LAMBDA(3) = 0.000D0
            IF(IORD.EQ.1) THEN
               LAMBDA(4)=0.326D0
               LAMBDA(5)=0.226D0
            ELSE
               LAMBDA(4)=0.215D0
               LAMBDA(5)=0.165D0
            ENDIF
*     ELSE IF(IIP.EQ.2) THEN    !MRST
*     Q2THR(1) = 2.6D0
*     Q2THR(1) = 2.24D0
*     Q2THR(1) = 1.D0
*     Q2THR(2) = 30.D0
*     LAMBDA(3)=0.000D0
*     LAMBDA(4)=0.231D0
*     LAMBDA(5)=0.151D0
         ENDIF
*     
         NF = 3
         DO 10 K = 1, 2
            IF (Q2 .GE. Q2THR (K)) THEN
               NF = NF + 1
            ELSE
               GO TO 20
            END IF
 10      CONTINUE
 20      B0 = 11.D0- 2.D0/3.D0* NF
         B0S = B0 * B0
         B1 = ( 102.D0- 38.D0/3.D0* NF ) * DBLE(IORD)
*     
         LAM2 = LAMBDA (NF) * LAMBDA (NF)
         LQ2 = DLOG (Q2 / LAM2)
         ALPHAS = 1.D0/ (B0 * LQ2) * (1.D0- B1 / B0S * DLOG (LQ2) / LQ2)
*
      ELSE IF(IIP.EQ.3.OR.IIP.EQ.4)THEN
c         IORD=1
c         CALL VFN(Q2,ALPHANNPDF,IORD,ALPHAMZ)
         PI=DACOS(-1.D0)
c         ALPHAS=ALPHANNPDF/4D0/PI 
         ALPHAS = alphasPDFM(4,dsqrt(Q2))/4D0/PI  ! Alphas from the PDF set
      ENDIF
      RETURN
      END
