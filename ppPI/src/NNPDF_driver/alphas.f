C     NNPDF subroutines for alpha computation

      SUBROUTINE VFN(qq2,alphas,ipt,FUNC)
*
      IMPLICIT none
*
      DOUBLE PRECISION EMC,ZETA2,ZETA3,ZETA4
      COMMON /CONSTS/ EMC,ZETA2,ZETA3,ZETA4
      DOUBLE PRECISION Q2TH(4:6),ASREF,Q2REF
      COMMON /ALPHASTHR/ Q2TH,ASREF,Q2REF
      DOUBLE PRECISION AS0,ASC,ASB,AST,ASQ
      COMMON /ASTHR/ AS0,ASC,ASB,AST,ASQ
      double precision pi
      parameter (pi=3.1415926535897932385d0)

      INTEGER nfi,nff,dnf,snf,ipt
      double precision q2,qq2,FUNC,alphasref,qq2ref
      double precision alphas,c2,ca,cf,tr,asi
      external FUNC

*     Color Factors

      CA= 3D0
      CF= 4D0/3D0
      TR= 0.5D0
      C2=14d0/3d0

*     Constants

      EMC   = 0.5772156649D0
      ZETA2 = 1.644934067D0  
      ZETA3 = 1.2020569031D0
      ZETA4 = 1.0823232337D0

*     Other parameters
      
      asref=.119d0
      q2ref=8.31744d3
      q2th(4)=2d0
      q2th(5)=22.5625d0
      q2th(6)=3.0625d4
      
      q2=qq2

      if(q2.gt.q2th(6))then
         nff=6
      elseif(q2.gt.q2th(5))then
         nff=5
      elseif(q2.gt.q2th(4))then
         nff=4
      else
         nff=3
      endif

      if(q2ref.gt.q2th(6))then
         nfi=6
      elseif(q2ref.gt.q2th(5))then
         nfi=5
      elseif(q2ref.gt.q2th(4))then
         nfi=4
      else
         nfi=3
      endif
*
      alphasref = asref
      qq2ref = q2ref

 10   if(nff.eq.nfi) then
         alphas = FUNC(nfi,qq2ref,alphasref,q2,ipt)
         return
      else
         if(nff.gt.nfi)then
            dnf=1
            snf=1
         else
            dnf=-1
            snf=0
         endif
         asi = FUNC(nfi,qq2ref,alphasref,q2th(nfi+snf),ipt)
         if(ipt.ge.2)then
            if(nff.gt.nfi) asi=asi+(c2/(4d0*pi)**2d0)*asi**3d0
            if(nff.lt.nfi) asi=asi-(c2/(4d0*pi)**2d0)*asi**3d0
         endif
         alphasref = asi
         qq2ref = q2th(nfi+snf)
         nfi = nfi+dnf
         goto 10
      endif
      end

      FUNCTION AlphaMZ(nfi,mz2,asz,q2,ipt)

      IMPLICIT none

      DOUBLE PRECISION EMC,ZETA2,ZETA3,ZETA4
      COMMON /CONSTS/ EMC,ZETA2,ZETA3,ZETA4
      DOUBLE PRECISION BETA0(3:6),BETA1(3:6),BETA2(3:6)
      DOUBLE PRECISION B1(3:6),B2(3:6)
      COMMON /BETA/ BETA0,BETA1,BETA2,B1,B2

      integer nfi,ipt,i
      double precision q2ref,q2,asi,asref
      double precision alo,t,as,den,mz2,asz,AlphaMZ
      double precision pi
      parameter (pi=3.1415926535897932385d0)

*     Beta function coefficients 
* 
      DO I=3,6         
         BETA0(I) = (33d0-2d0*I)/3d0  
         BETA1(I) = 102d0-38d0/3d0*I 
         BETA2(I) = 2857d0/2d0-5033d0/18d0*I+325d0/54d0*I**2d0 
         B1(I) = BETA1(I)/BETA0(I)
         B2(I) = BETA2(I)/BETA0(I)
      ENDDO

      q2ref=mz2
      asref=asz/4d0/pi

      asi=asref
      t=log(q2/q2ref)
      den=1+beta0(nfi)*asi*t
      alo=asi/den
*
*     LO
*
      as=alo
*
*     NLO
*
      if(ipt.ge.1)as=alo*(1-b1(nfi)*alo*log(den))
*
*     NNLO
*
      if(ipt.ge.2)then
         as=alo*(1d0+(alo*(alo-asi)*(b2(nfi)-b1(nfi)**2d0)
     #        +as*b1(nfi)*dlog(as/asi)))
      endif
*
      AlphaMZ=4d0*pi*as
*
      RETURN
      END
