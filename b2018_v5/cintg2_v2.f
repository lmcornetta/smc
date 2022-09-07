      module cintg_pp

c***
c Modified to remove common blocks CALLIN and FUNC (in consistency
c with the modifications done in MPLANE).
c***

      use mathlib, only: cderfm

      contains
 
c----------------------------------------------------------------
      double complex function cint(zval,c,ro,aion,sigma,iz,
     .                             xa,alpa,ntypa,xb,xc)
!    .                        xa,alpa,ntypa,ntypb,xb,alpb,xc,alpc,ntypc)

      use mplane1, only: nap, lap, map, iend, cterm, ll,
     .    nmax, lmax, mmax, ientr, jentr, cnorm, pi, twopi,
     .    pi3haf, pi5hf2, piquart, ifirst

cccc  IMPLICIT REAL*8 (A-H,O-Z)
      implicit none

      integer :: iz, ntypa!, ntypb, ntypc, ntypd
      real*8  :: zval, c(2), ro(2), aion(2,3,3), sigma(3,3)
      real*8  :: xa(3), alpa, xb(3)!, alpb
      real*8  :: xc(3)!, alpc, xd(3), alpd

cccc  COMPLEX*16 XINT,YINT,ZINT,PWCORE,PWION,UI,FATC,PFAB,PFCD
      COMPLEX*16 XINT,YINT,ZINT,UI,FATC

cccc  COMMON /FUNC/ NAP(60),LAP(60),MAP(60),IEND(60),CTERM(30),LL(60)
cccc &      ,NMAX(60),LMAX(60),MMAX(60),IENTR(60),JENTR(60),CNORM(60)
cccc &      ,PI,TWOPI,PI3HAF,PI5HF2,PIQUART,IFIRST

cccc  COMMON /CALLIN/XA,YA,ZA,ALPA,NTYPA,NTYPB,XB,YB,ZB,ALPB,
cccc & XC,YC,ZC,ALPC,NTYPC,NTYPD,XD,YD,ZD,ALPD,PFAB,PFCD

cccc  DIMENSION C(2),RO(2),CA(3),CK(3),LEXP(3),SIGMA(3,3),AION(2,3,3)
      integer :: lexp(3)
      real*8  :: ca(3), ck(3)

c***
c  Local variables
c
      real*8  :: alfa
c***

C--------------------------------------------------------------------
C
C     INITIALIZE THE INTEGRALS' VALUE
C
C--------------------------------------------------------------------
      CINT=(0.D0,0.D0)
      XINT=(0.D0,0.D0)
      YINT=(0.D0,0.D0)
      ZINT=(0.D0,0.D0)
C--------------------------------------------------------------------
C
C     DEFINE THE PARAMETERS
C
C--------------------------------------------------------------------
      UI=(0.D0,1.D0)
      ALFA=ALPA
      CA(:)=XA(:)-XC(:)
      CK(:)=XB(:)
      LEXP(1)=NAP(NTYPA)
      LEXP(2)=LAP(NTYPA)
      LEXP(3)=MAP(NTYPA)
      FATC=CDEXP(-UI*sum(CK(:)*XC(:)))
C--------------------------------------------------------------------
C     CALCULATION OF THE MATRIX ELEMENTS BEGINS
C--------------------------------------------------------------------
C
C     CORE PART INTEGRATION (XINT)
C
C--------------------------------------------------------------------
      XINT=PWCORE(ZVAL,ALFA,C,RO,CA,CK,LEXP)
C--------------------------------------------------------------------
C
C     ION PART INTEGRATION (YINT)
C
C--------------------------------------------------------------------
      YINT=PWION(AION,SIGMA,ALFA,CA,CK,LEXP,IZ)
C--------------------------------------------------------------------
C
C     INTERMEDIATE RESULT IS ZINT (CORE+ION)
C
C--------------------------------------------------------------------
      ZINT=XINT+YINT
      ZINT=FATC*ZINT
C--------------------------------------------------------------------
C
C     FINAL RESULT IS CINT (CORE+ION)
C
C--------------------------------------------------------------------
      CINT=DCONJG(ZINT)
C
ccc   RETURN
ccc   END
      end function cint
c----------------------------------------------------------------

c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION PWCORE(ZVAL,ALFA,BORE,RO,BA,BK,LEXP)
C!====================================================================!
C!                                                                    !
C!    CALCULA O ELEMENTO DE MATRIZ DA PARTE DO CORE DO PSEUDO-        !
C!    POTENCIAL CENTRADO NA ORIGEM ENTRE UMA ONDA PLANA E UMA         !
C!    GAUSSIANA CARTESIANA CENTRADA EM CA.                            !
C!                                                                    !
C!====================================================================!
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16(C)      
      DIMENSION RO(2),BORE(2),GAMA(2),GROOT(2),LEXP(3),BA(3),BK(3)
      PWCORE=(0.D0,0.D0)
      PI=4.D0*DATAN(1.D0)
      C1=(0.D0,1.D0)
      CVX=BK(1)/2.D0/ALFA/C1
      CVY=BK(2)/2.D0/ALFA/C1
      CVZ=BK(3)/2.D0/ALFA/C1
      AMOD=(BA(1)*BA(1)+BA(2)*BA(2)+BA(3)*BA(3))**0.5
      VSQ=CVX*CVX+CVY*CVY+CVZ*CVZ
      ASQ=AMOD*AMOD
      V1=ALFA*VSQ
      AK=BA(1)*BK(1)+BA(2)*BK(2)+BA(3)*BK(3)
      CAV=BA(1)*CVX+BA(2)*CVY+BA(3)*CVZ
      DO 2 I=1,2
    2 GAMA(I)=0.D0
      DO 3 I=1,2
    3 GAMA(I)=(RO(I)+ALFA)/(RO(I)*ALFA)/4.D0
      CFAT0=-ZVAL*(PI/ALFA)**1.5*CDEXP(-C1*AK)
      CMOD=CDSQRT(ASQ+VSQ+2.D0*CAV)
      CMODSQ=CMOD*CMOD
      LGAUS=LEXP(1)+LEXP(2)+LEXP(3)+1
      GO TO(101,102,103),LGAUS
C--------------------------------------------------------------------
C
C     101 - S SYMMETRY GAUSSIAN
C     102 - P SYMMETRY GAUSSIAN
C     103 - D SYMMETRY GAUSSIAN
C
C--------------------------------------------------------------------
  101 CSUM0=(0.D0,0.D0)
      DO 4 I=1,2
      GROOT(I)=1.D0/2.D0/GAMA(I)**0.5
    4 CSUM0=CSUM0+BORE(I)*GROOT(I)*CF0(CMOD*GROOT(I),V1)
C
      PWCORE=CFAT0*CSUM0
      RETURN
C
  102 CFAT1=2.D0*ALFA
      CFAT2=2.D0
      CSUM0=(0.D0,0.D0)
      CSUM1=(0.D0,0.D0)
      DO 5 I=1,2
      GROOT(I)=1.D0/2.D0/GAMA(I)**0.5
      CSUM0=CSUM0+BORE(I)*GROOT(I)*CF0(CMOD*GROOT(I),V1)
      CSUM1=CSUM1+BORE(I)*GROOT(I)*CF1(CMODSQ*GROOT(I)**2,V1)/
     &4.D0/GAMA(I)
    5 CONTINUE
      CPROD1=CFAT1*CSUM0
      CPROD2=CFAT2*CSUM1
C--------------------------------------------------------------------
C
C     X
C
C--------------------------------------------------------------------
      IF(LEXP(1).EQ.1)PWCORE=CPROD1*CVX+CPROD2*(BA(1)+CVX)
C--------------------------------------------------------------------
C
C     Y
C
C--------------------------------------------------------------------
      IF(LEXP(2).EQ.1)PWCORE=CPROD1*CVY+CPROD2*(BA(2)+CVY)
C--------------------------------------------------------------------
C
C     Z
C
C--------------------------------------------------------------------
      IF(LEXP(3).EQ.1)PWCORE=CPROD1*CVZ+CPROD2*(BA(3)+CVZ)
C
      PWCORE=CFAT0*PWCORE/2.D0/ALFA
      RETURN
C
  103 CFAT=2.D0
      CFAT1=4.D0*ALFA
      CFAT2=CFAT1*ALFA
      CFAT3=CFAT1/ALFA
      CSUM0=(0.D0,0.D0)
      CSUM1=(0.D0,0.D0)      
      CSUM2=(0.D0,0.D0)      
      DO 6 I=1,2
      GROOT(I)=1.D0/2.D0/GAMA(I)**0.5
      CSUM0=CSUM0+BORE(I)*GROOT(I)*CF0(CMOD*GROOT(I),V1)
      CSUM1=CSUM1+BORE(I)*GROOT(I)*CF1(CMODSQ*GROOT(I)**2,V1)/
     &4.D0/GAMA(I)
      CSUM2=CSUM2+BORE(I)*GROOT(I)*CF2(CMODSQ*GROOT(I)**2,V1)/
     &16.D0/GAMA(I)**2
    6 CONTINUE
      CPROD=CFAT*(ALFA*CSUM0+CSUM1)
      CPROD1=CFAT1*CSUM1
      CPROD2=CFAT2*CSUM0
      CPROD3=CFAT3*CSUM2
C--------------------------------------------------------------------
C
C     XY
C
C--------------------------------------------------------------------
      IF(LEXP(1).EQ.1.AND.LEXP(2).EQ.1)THEN
      PWCORE=CPROD1*(2.D0*CVX*CVY+CVX*BA(2)+CVY*BA(1))+
     &       CPROD2*CVX*CVY+
     &       CPROD3*(BA(1)+CVX)*(BA(2)+CVY)
      END IF
C--------------------------------------------------------------------
C
C     XZ
C
C--------------------------------------------------------------------
      IF(LEXP(1).EQ.1.AND.LEXP(3).EQ.1)THEN
      PWCORE=CPROD1*(2.D0*CVX*CVZ+CVX*BA(3)+CVZ*BA(1))+
     &       CPROD2*CVX*CVZ+
     &       CPROD3*(BA(1)+CVX)*(BA(3)+CVZ)
      END IF
C--------------------------------------------------------------------
C
C     YZ
C
C--------------------------------------------------------------------
      IF(LEXP(2).EQ.1.AND.LEXP(3).EQ.1)THEN
      PWCORE=CPROD1*(2.D0*CVY*CVZ+CVY*BA(3)+CVZ*BA(2))+
     &       CPROD2*CVY*CVZ+
     &       CPROD3*(BA(2)+CVY)*(BA(3)+CVZ)
      END IF
C--------------------------------------------------------------------
C
C     XX
C
C--------------------------------------------------------------------
      IF(LEXP(1).EQ.2)THEN
      PWCORE=CPROD+CPROD1*2.D0*(CVX*CVX+CVX*BA(1))+
     &       CPROD2*CVX*CVX+
     &       CPROD3*(BA(1)+CVX)*(BA(1)+CVX)
      END IF
C--------------------------------------------------------------------
C
C     YY
C
C--------------------------------------------------------------------
      IF(LEXP(2).EQ.2)THEN
      PWCORE=CPROD+CPROD1*2.D0*(CVY*CVY+CVY*BA(2))+
     &       CPROD2*CVY*CVY+
     &       CPROD3*(BA(2)+CVY)*(BA(2)+CVY)
      END IF
C--------------------------------------------------------------------
C
C     ZZ
C
C--------------------------------------------------------------------
      IF(LEXP(3).EQ.2)THEN
      PWCORE=CPROD+CPROD1*2.D0*(CVZ*CVZ+CVZ*BA(3))+
     &       CPROD2*CVZ*CVZ+
     &       CPROD3*(BA(3)+CVZ)*(BA(3)+CVZ)
      END IF
C
      PWCORE=CFAT0*PWCORE/4.D0/ALFA**2
ccc   RETURN
ccc   END
      end function pwcore
c----------------------------------------------------------------


c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION CF0(CARG,V1)
C
C     CALCULA A FUNCAO F
C
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16(C)
      CF0=(DEXP(V1)-CDERFM(CARG)*CDEXP(-CARG**2+V1))/CARG
ccc   RETURN
ccc   END
      end function cf0
c----------------------------------------------------------------


c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION CF1(CARG,V1)
C
C     CALCULA A DERIVADA PRIMEIRA DE F
C
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16(C)
      PI=4.D0*DATAN(1.D0)
      UME=2.D0/PI**0.5
      CF1=(UME*CDEXP(-CARG+V1)-CF0(CDSQRT(CARG),V1))/2.D0/CARG
ccc   RETURN
ccc   END
      end function cf1
c----------------------------------------------------------------


c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION CF2(CARG,V1)
C
C     CALCULA A DERIVADA SEGUNDA DE F
C
      IMPLICIT REAL*8 (A,B,D-H,O-Z),COMPLEX*16(C)
      PI=4.D0*DATAN(1.D0)
      CUME1=1.D0/2.D0/CARG
      UME2=2.D0/PI**0.5     
      CFEXP=CDEXP(-CARG+V1)
      CTERM1=CUME1*(CF0(CDSQRT(CARG),V1)/CARG-CF1(CARG,V1))
      CTERM2=CUME1*UME2*CFEXP*(1.D0/CARG+1.D0)
      CF2=CTERM1-CTERM2
ccc   RETURN
ccc   END
      end function cf2
c----------------------------------------------------------------


c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION PWION(AION,SIGMA,ALFA,CA,CK,LEXP,IZ)
C!====================================================================!
C!                                                                    !
C!    CALCULA O ELEMENTO DE MATRIZ DA PARTE IONICA DO PSEUDOPOTENCIAL !
C!    CENTRADO NA ORIGEM ENTRE UMA ONDA PLANA DADA POR                !
C!                                                                    !
C!    DEXP(-i*(KX*X+KY*Y+KZ*Z))                                       !
C!                                                                    !
C!    E UMA GAUSSIANA CARTESIANA CENTRADA EM CA DADA POR              !
C!                                                                    !
C!    N(L1,L2,L3)*(X-CA(1))**L1*(Y-CA(2))**L2*(Z-CA(3))**L3*          !
C!    DEXP(-ALFA*(R-CA)**2)                                           !
C!                                                                    !
C!    NO PROGRAMA FAZEMOS                                             !
C!                                                                    !
C!    KX=CK(1), KY=CK(2), KZ=CK(3)                                    !
C!    UIL = (-i)**LPQ                                                 !
C!    LPQ = MOMENTO ANGULAR (0, 1 E 2)                                !
C!    LGR = L1+L2+L3                                                  !
C!                                                                    !
C!====================================================================!
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 UIL
      COMMON/CHATO/PI,PI2,FAT,COMECO,UIL(3)
      DIMENSION AION(2,3,3),SIGMA(3,3),CA(3),CK(3),LEXP(3),CHAVE(3,5,3)
C
      IF(COMECO.EQ.0.123456789D0)GO TO 100
      PI=4.D0*DATAN(1.D0)
      PI2=PI*PI
      FAT=16.D0*PI2
      COMECO=0.123456789D0
      UIL(1)=(1.D0,0.D0)
      UIL(2)=(0.D0,-1.D0)
      UIL(3)=(-1.D0,0.D0)
C
C     CALCULA A PRIMEIRA CHAVE DO ELEMENTO DE MATRIZ
C
  100 LGR=LEXP(1)+LEXP(2)+LEXP(3)
      CALL KGINT(CHAVE,IZ,AION,SIGMA,CA,CK,ALFA)
      PWION=(0.D0,0.D0)
      DO 3 L=1,3
      LPQ=L-1
      QWION=0.D0
      DO 4 LAMB=1,L+LGR
      LAM=LAMB-1
      MIN=IABS(LAM-LPQ)
      IF(MIN.GT.LGR)GO TO 4
      DO 5 NU=MIN+1,LGR+1
      NNU=NU-1
C
C     TESTA A PARIDADE DE L, LAMBDA E NU
C
      IF(NNU+LPQ+LAM.NE.2*((NNU+LPQ+LAM)/2))GO TO 5
C
C     CALCULA A SEGUNDA CHAVE DO ELEMENTO DE MATRIZ
C
      CALL KGESTR(GIONK,CK,CA,LEXP,LPQ,LAM,NNU)
C
      QWION=QWION+CHAVE(NU,LAMB,L)*GIONK
    5 CONTINUE
    4 CONTINUE
      PWION=PWION+QWION*UIL(L)
    3 CONTINUE
      PWION=FAT*PWION
ccc   RETURN
ccc   END      
      end function pwion
c----------------------------------------------------------------


c----------------------------------------------------------------
      SUBROUTINE KGINT(CHAVE,IZ,AION,SIGMA,CA,CK,ALFA)
C!==============================================================!
C!                                                              !
C!    CALCULA A PRIMEIRA CHAVE DO ELEMENTO DE MATRIZ.           !
C!    CHAVE = RESULTADO DA DUPLA SOMA EM J E N                  !
C!    L = MOMENTO ANGULAR (0, 1 E 2)                            !
C!    LAM = LAMBDA                                              !
C!    NU = NU                                                   !
C!    IZ = NUMERO ATOMICO (APENAS PARA A IDENTIFICACAO DO ATOMO)!
C!    AION, SIGMA = PARAMETROS DO PSEUDOPOTENCIAL               !
C!    CA = CENTRO DA GAUSSIANA CARTESIANA                       !
C!    CK = VETOR DE ONDA                                        !
C!    ALFA = EXPOENTE DA GAUSSIANA CARTESIANA                   !
C!                                                              !
C!==============================================================!
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SALVAG/COMECO,CAANT,NORD(9720),PARAM(9720),
     * SALVA(15,9720),PI,NSALVA,INDERE(15)
      DIMENSION CA(3),CK(3),SIGMA(3,3),AION(2,3,3),CHAVE(45),
     * DENOSQ(3,3),DEXPO(3,3),GTIL(75,3,3)
C
      IF(COMECO.EQ.0.123456789D0)GO TO 100
      NSALVA=0
      PI=4.D0*DATAN(1.D0)
      CAANT=1.D20
      KONT=0
      DO 313 LL=1,3
      L=LL-1
      DO 313 LLAM=1,LL+2
      LAM=LLAM-1
      MIN=IABS(LAM-L)+1
      IF(MIN.GT.3)GO TO 313
      DO 333 NNU=MIN,3
      NU=NNU-1
      IAUX=NU+L+LAM
      IF(IAUX.NE.2*(IAUX/2))GO TO 333
      KONT=KONT+1
      INDERE(KONT)=NNU+3*LAM+15*L
  333 CONTINUE
  313 CONTINUE
      COMECO=0.123456789D0
  100 CKSQ=CK(1)**2+CK(2)**2+CK(3)**2
      IF(DABS(CKSQ-CAANT).LT.1.D-6)GO TO 200
      NSALVA=0
      CAANT=CKSQ
C
C   CALCULO DO PARAMETRO DE PROCURA
C  PARAMETRO DE PROCURA = PI*Z + ALFA*A
  200 CASQ=CA(1)**2+CA(2)**2+CA(3)**2
      PROCUR=IZ*PI**2+ALFA*PI+CASQ
      CALL REACH(NSALVA,PARAM(1),PROCUR,IGUAL,MENOR)
      IF(IGUAL.EQ.0)GO TO 300
C   NO CASO DE TER ACHADO NA MEMORIA SIMPLESMENTE TRANSFERE
C
      DO 13 KONT=1,15
   13 CHAVE(INDERE(KONT))=SALVA(KONT,NORD(IGUAL))
      RETURN
C
C   CASO EM QUE A CHAVE DEVE SER CALCULADA
  300 CKMOD=DSQRT(CKSQ)
      CAMOD=DSQRT(CASQ)
      NSALVA=NSALVA+1
      PARAM(NSALVA)=PROCUR
      NORD(NSALVA)=NSALVA
      DO 1 L=1,3
      DO 1 J=1,3
      DENO=SIGMA(J,L)+ALFA
      DENOSQ(J,L)=DSQRT(DENO)
      AA=2.D0*CAMOD*ALFA/DENOSQ(J,L)
      BB=CKMOD/DENOSQ(J,L)
      AUX=ALFA*SIGMA(J,L)*CASQ/DENO+CKSQ/4.D0/DENO
      DEXPO(J,L)=0.D0
      IF(AUX.LT.25.D0)DEXPO(J,L)=DEXP(-AUX)
      CALL INTGRL(AA,BB,GTIL(1,J,L))
    1 CONTINUE
      KONT=0
      DO 2 LL=1,3
      L=LL-1
      DO 2 LLAM=1,LL+2
      LAM=LLAM-1
      MIN=IABS(LAM-L)+1
      IF(MIN.GT.3)GO TO 2
      DO 22 NNU=MIN,3
      NU=NNU-1
      IAUX=NU+L+LAM
      IF(IAUX.NE.2*(IAUX/2))GO TO 22
      KONT=KONT+1
      AUX=0.D0
      DO 3 J=1,3
      BUX=0.D0
      DO 4 NN=1,2
      N=NN-1
      M=2*N+NU
    4 BUX=BUX+AION(NN,J,LL)*GTIL(M+1+LAM*5+L*25,J,LL)/DENOSQ(J,LL)**
     * (M+3)
    3 AUX=AUX+DEXPO(J,LL)*BUX
      CHAVE(INDERE(KONT))=AUX
      SALVA(KONT,NSALVA)=AUX
   22 CONTINUE
    2 CONTINUE
C  REORDENA OS PARAMETROS
      IF(MENOR+1.GT.NSALVA-1)RETURN
      DO 5 II=MENOR+1,NSALVA-1
      I=NSALVA+MENOR+1-II
      NORD(I)=NORD(I-1)
    5 PARAM(I)=PARAM(I-1)
      NORD(MENOR+1)=NSALVA
      PARAM(MENOR+1)=PROCUR
ccc   RETURN
ccc   END
      end subroutine kgint
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE KGESTR(GIONK,CA,A,LEXP,L,LAM,NU)
C   CALCULA O FATOR ESTRUTURAL PARA O ELEMENTO DE MATRIZ
C   DO TERMO IONICO  ENTRE ONDA PLANA E GAUSSIANA 
C   CARTESIANA
C  GIONK = ELEMENTO DE MATRIZ
C  CA = COMPONENTES DO VETOR DE ONDA K DA ONDA PLANA
C  A = VETOR DE POSICAO DO CENTRO DA GAUSSIANA RELATIVO
C      AO CENTRO DO POTENCIAL
C  LEXP = EXPOENTES DA PARTE CARTESIANA DA GAUSSIANA
C  L = MOMENTO ANGULAR DE UM DOS TERMOS IONICOS DO POTENCIAL
C  LAM = LAMBDA
C  NU = NU
C    SE LEXP, L, LAM, E NU NAO OBEDECEREM AS REGRAS DE SOMA,
C    GIONK E' FEITO IGUAL A ZERO.
C
C   ESTA' DIMENSIONADO PARA UM MAXIMO DE 24 ATOMOS     
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ESTR/AJA(277),YQ(381,277),NORD(277),NUMBER,
     * INDERE(9,10,3,5),PI,COMECO,CANT,YK(9)
C
      DIMENSION CA(3),A(3),LEXP(3),M3(27)
      DATA M3/1,4,7,3,8,0,6,0,0,2,9,0,10,0,0,0,0,0,5,0,0,0,0,0,0,0,0/
C
      MG=M3(9*LEXP(1)+3*LEXP(2)+LEXP(3)+1)
      LG=LEXP(1)+LEXP(2)+LEXP(3)
      GIONK=0.D0
      IF(L.GT.2)RETURN
      IF(LG.GT.2)RETURN
      IF(LAM.GT.L+LG)RETURN
      IF(NU.LT.IABS(LAM-L))RETURN
      IF(NU.GT.LG)RETURN
      IF(LAM+L+NU.NE.2*((LAM+L+NU)/2))RETURN
      IF(COMECO.EQ.0.123456789D0)GO TO 202
      KONT=0
      DO 215 LL=1,5
      LAMM=LL-1
      DO 215 KK=1,3
      DO 215 II=1,10
      LGG=0
      IF(II.GT.1)LGG=1
      IF(II.GT.4)LGG=2
      DO 215 JJ=1,9
      LP=0
      IF(JJ.GT.1)LP=1
      IF(JJ.GT.4)LP=2
      INDERE(JJ,II,KK,LL)=0
      IF(LAMM.GT.LP+LGG)GO TO 215
      IF(KK-1.LT.IABS(LAMM-LP))GO TO 215
      IF(KK.GT.LGG+1)GO TO 215
      IF(LAMM+LP+KK-1.NE.2*((LAMM+LP+KK-1)/2))GO TO 215
      KONT=KONT+1
      INDERE(JJ,II,KK,LL)=KONT
  215 CONTINUE
      NUMBER=0
      CANT=1.D20
      PI=4.D0*DATAN(1.D0)
      COMECO=0.123456789D0
  202 CONTINUE
C
C   CALCULA OU USA AS HARMONICAS ESFERICAS DO VETOR K CALCULADAS
C   ANTERIORMENTE
      CAMOD=DSQRT(CA(1)**2+CA(2)**2+CA(3)**2)
      IF(CAMOD.EQ.0.D0)THEN
        CAMOD=1.D0
        CA(1)=1.D0
        ENDIF
      CACT=(CA(1)+PI*CA(2)+PI**2*CA(3))/CAMOD
      IF(DABS(CACT-CANT).GT.1.D-8)CALL PEQHAR(CA(1),CA(2),CA(3),YK)
      CANT=CACT
C
C   PROCURA OU CALCULA A SOMA EM MU DO PRODUTO Y(A)*Q
      ACT=A(1)+PI*A(2)+PI**2*A(3)
      INAL=-1
      CALL REACH(NUMBER,AJA(1),-ACT,IGUAL,MENOR)
      IF(IGUAL.NE.0)GO TO 101
      INAL=1
      CALL REACH(NUMBER,AJA(1),ACT,IGUAL,MENOR)
C  SE A SOMA JA TINHA SIDO CALCULADA, DESVIA
      IF(IGUAL.NE.0)GO TO 101
      NUMBER=NUMBER+1
      IACT=NUMBER
      CALL ESTRUT(A,YQ(1,NUMBER),INDERE)
      AJA(NUMBER)=ACT
      NORD(NUMBER)=NUMBER
C   REORDENA 'AJA' E 'NORD'
      IF(MENOR+1.GT.NUMBER-1)GO TO 102
      DO 1 II=MENOR+1,NUMBER-1
      I=NUMBER+MENOR+1-II
      NORD(I)=NORD(I-1)
    1 AJA(I)=AJA(I-1)
      AJA(MENOR+1)=ACT
      NORD(MENOR+1)=NUMBER
      GO TO 102
C   CASO EM QUE A SOMA JA TINHA SIDO CALCULADA
  101 IACT=NORD(IGUAL)
  102 INICIO=INDERE(L**2+1,MG,NU+1,LAM+1)-L**2-1
      INAL=INAL**(LG+L)
      GIONK=0.D0
      DO 2 M=L**2+1,(L+1)**2
    2 GIONK=GIONK+YK(M)*YQ(INICIO+M,IACT)
      GIONK=GIONK*INAL
ccc   RETURN
ccc   END
      end subroutine kgestr
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE ESTRUT(A,YQ,INDERE)
C     CALCULA OS FATORES DE ESTRUTURA YQ(M,MG,NU,LAM) PARA O VETOR
C   (A(1) A(2) A(3))
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SALVA/AINIC,AINV(25,25),YP(9,25),RANV(3,25),PI
      DIMENSION A(3),Q(25,3),MREP(10),VET(3),QLM(10,3,25),
     1 AMAT(25,25),NL(25),NC(25),INDERE(9,10,3,5),YQ(381),
     2 YA(25)
C
      DATA MREP/0,100,10,1,200,20,2,11,101,110/
      IF(AINIC.EQ.0.123456789D0)GO TO 5
      PI=4.D0*DATAN(1.D0)
      NALEAT=0
      DO 1 IR=1,25
      RANV(3,IR)=-1.D0+2.D0*RAN(NALEAT)
      FI=-PI+2.D0*PI*RAN(NALEAT)
      AUX=DSQRT(1.D0-RANV(3,IR)**2)
      RANV(1,IR)=AUX*DCOS(FI)
      RANV(2,IR)=AUX*DSIN(FI)
      CALL HARMON(PI,4,RANV(1,IR),RANV(2,IR),RANV(3,IR),AMAT(1,IR))
      DO 1 I=1,9
    1 YP(I,IR)=AMAT(I,IR)
      CALL INVMAT(1,25,25,NL,NC,ISIG,AMAT,AINV)
      AINIC=0.123456789D0
    5 AMOD=DSQRT(A(1)**2+A(2)**2+A(3)**2)
      IF(AMOD.LT.1.D-8)AMOD=1.D0
      DO 2 MG=1,10
      LA=MREP(MG)/100
      NA=MREP(MG)-100*LA
      MA=NA/10
      NA=NA-10*MA
      DO 2 IR=1,25
      DO 2 N=1,3
      ZLM=1.D0
      IF(LA.GT.0)ZLM=ZLM*((N-1)*AMOD*RANV(1,IR)-A(1))
      IF(LA.GT.1)ZLM=ZLM*((N-1)*AMOD*RANV(1,IR)-A(1))
      IF(MA.GT.0)ZLM=ZLM*((N-1)*AMOD*RANV(2,IR)-A(2))
      IF(MA.GT.1)ZLM=ZLM*((N-1)*AMOD*RANV(2,IR)-A(2))
      IF(NA.GT.0)ZLM=ZLM*((N-1)*AMOD*RANV(3,IR)-A(3))
      IF(NA.GT.1)ZLM=ZLM*((N-1)*AMOD*RANV(3,IR)-A(3))
    2 QLM(MG,N,IR)=ZLM
      CALL HARMON(PI,4,A(1),A(2),A(3),YA(1))
      DO 3 MG=1,10
      DO 3 M=1,9
      DO 8 LAM=1,5
      DO 6 N=1,3
      IF(INDERE(M,MG,N,LAM).NE.0)GO TO 9
    6 CONTINUE
      GO TO 8
    9 DO 7 MU=(LAM-1)**2+1,LAM**2
      DO 4 N=1,3
      VET(N)=0.D0
      DO 4 IR=1,25
    4 VET(N)=VET(N)+AINV(IR,MU)*YP(M,IR)*QLM(MG,N,IR)
      VET(3)=0.5D0*(VET(3)+VET(1)-2.D0*VET(2))/AMOD**2
      VET(2)=(VET(2)-VET(1)-VET(3)*AMOD**2)/AMOD
      DO 7 N=1,3
    7 Q(MU,N)=VET(N)
      DO 10 N=1,3
      IZUX=INDERE(M,MG,N,LAM)
      IF(IZUX.EQ.0)GO TO 10
      ZUX=0.D0
      DO 11 MU=(LAM-1)**2+1,LAM**2
   11 ZUX=ZUX+YA(MU)*Q(MU,N)
      YQ(IZUX)=ZUX
   10 CONTINUE
    8 CONTINUE
    3 CONTINUE
ccc   RETURN
ccc   END
      end subroutine estrut
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE REACH(NUMBER,ARRAY,X,IGUAL,MENOR)
C   ACHA O NUMERO DE ORDEM DO ELEMENTO DO ARRAY ORDENADO 'ARRAY'
C   QUE COINCIDE COM 'X'.
C    IGUAL = NUMERO DE ORDEM ACHADO (ZERO SE NAO ACHAR)
C    MENOR = NUMERO DE ORDEM CUJO ARRAY E' IMEDIATEMENTE
C            MENOR QUE X. SE IGUAL NAO FOR NULO, MENOR
C            E' FEITO IGUAL A IGUAL.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARRAY(NUMBER+1)
      SMALL=1.D-8
      IGUAL=0
      MENOR=0
      IF(NUMBER.EQ.0)RETURN
      IF(X.LT.ARRAY(1)-SMALL)RETURN
      MENOR=NUMBER
      IF(X.GT.ARRAY(NUMBER)+SMALL)RETURN
      MIN=1
      MAX=NUMBER
      IF(DABS(X-ARRAY(MIN)).GT.SMALL)GO TO 1
      IGUAL=MIN
      MENOR=MIN
      RETURN
    1 IF(DABS(X-ARRAY(MAX)).GT.SMALL)GO TO 2
      IGUAL=MAX
      MENOR=MAX
      RETURN
    2 IAUX=(MIN+MAX)/2
      MENOR=MIN
      IF(IAUX.EQ.MIN)RETURN
      IF(DABS(X-ARRAY(IAUX)).GT.SMALL)GO TO 3
      IGUAL=IAUX
      MENOR=IAUX
      RETURN
    3 IF(X.GT.ARRAY(IAUX))GO TO 4
      MAX=IAUX
      GO TO 2
    4 MIN=IAUX
      GO TO 2
ccc   END
      end subroutine reach
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE HARMON(PI,LMAX,AX,AY,AZ,YLM)
CEQ*******************************************************************
CEQ*   CALCULA AS HARMONICAS ESFERICAS REAIS ('YLM').
CEQ*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YLM((LMAX+1)**2),FAT(16,16),P(17,16),C(31)
      LBOBO=LMAX+1
      IF(LBOBO.LT.3) LBOBO=3
      P(1,1)=1.D0
      C(16)=1.D0/DSQRT(2.D0)
      AUX=AX**2+AY**2+AZ**2
      IF(AUX.GT.1.D-8)GO TO 345
C   SE O VETOR FOR NULO, ARBITRA A DIRECAO (1 0 0)
      AUX=1.D0
      X=0.D0
      SQ=1.D0
      C(17)=1.D0
      C(15)=0.D0
      GO TO 346
  345 AUX=DSQRT(AUX)
      X=AZ/AUX
      SQ=DSQRT(1.D0-X**2)
      IF(SQ.GT.1.D-8)GO TO 347
C   NESTE CASO ARBITRA O ANGULO FI COMO ZERO
      C(17)=1.D0
      C(15)=0.D0
      SQ=1.D-16
      GO TO 346
  347 C(17)=AX/AUX/SQ
      C(15)=AY/AUX/SQ
  346 P(2,1)=X
      DO 1 L=3,LBOBO+1
 1    P(L,1)=((2*L-3)*X*P(L-1,1)-(L-2)*P(L-2,1))/(L-1)
      DO 2 M=2,LBOBO
      DO 2 L=M,LBOBO+1
 2    P(L,M)=((L-M+1)*X*P(L,M-1)-(L+M-3)*P(L-1,M-1))/SQ
      DO 3 M=18,15+LBOBO
      C(M)=C(17)*C(M-1)-C(15)*C(33-M)
 3    C(32-M)=C(15)*C(M-1)+C(17)*C(33-M)
      DO 4 L=1,LBOBO
 4    FAT(L,1)=1.D0
      DO 5 M=2,LBOBO
      DO 5 L=M,LBOBO
 5    FAT(L,M)=FAT(L,M-1)/(L-M+1)/(L+M-2)
      MM=0
      DO 6 L=1,LMAX+1
      DO 6 M=1,2*L-1
      MM=MM+1
      MMM=IABS(M-L)+1
      AUX=DSQRT((2*L-1)/2.D0/PI*FAT(L,MMM))
    6 YLM(MM)=AUX*P(L,MMM)*C(M-L+16)
ccc   RETURN
ccc   END
      end subroutine harmon
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE INVMAT(KEY,M,NDIM,NL,NC,ISIG,A,AINV)
CEQ*****************************************************************
CEQ*   INVERTE MATRIZES. VER EXPLICACAO NO PROGRAMA CELSIM
CEQ*****************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NL(M),NC(M),A(M,M),AINV(M,M)
      IF(KEY.EQ.0) GO TO 1
      ISIG=1
      AMAX=0.D0
      DO 2 I=1,NDIM
      NL(I)=I
      NC(I)=I
      DO 2 J=1,NDIM
      BMAX=DABS(A(I,J))
      IF(BMAX.LE.AMAX)GO TO 2
      AMAX=BMAX
      MI=I
      MJ=J
    2 CONTINUE
      IF(AMAX.NE.0.D0)GO TO 4
      STOP
CSTOP - PORQUE A MATRIZ EH IDENTICAMENTE NULA.
    4 IF(MI.EQ.1)GO TO 3
      NL(1)=MI
      NL(MI)=1
      ISIG=-ISIG
    3 IF(MJ.EQ.1)GO TO 1
      NC(1)=MJ
      NC(MJ)=1
      ISIG=-ISIG
    1 DO 6 I=1,NDIM
      IL=NL(I)
      IC=NC(I)
      AII=A(IL,IC)
      A(IL,IC)=1.D0
      AMAX=0.D0
      DO 7 J=1,NDIM
    7 A(IL,J)=A(IL,J)/AII
      DO 8 J=1,NDIM
      IF(J-I)9,8,9
    9 JL=NL(J)
      AUX=A(JL,IC)
      A(JL,IC)=0.D0
      DO 10 K=1,NDIM
      KC=NC(K)
      A(JL,KC)=A(JL,KC)-AUX*A(IL,KC)
      IF(KEY)11,10,11
   11 IF(J-I)10,10,12
   12 IF(K-I)10,10,13
   13 BMAX=DABS(A(JL,KC))
      IF(BMAX-AMAX)10,10,18
   18 MJ=J
      MK=K
      AMAX=BMAX
   10 CONTINUE
    8 CONTINUE
      IF(KEY.EQ.0.OR.I.EQ.NDIM)GO TO 6
      IF(AMAX.EQ.0.D0)STOP
CSTOP - PORQUE A MATRIZ EH SINGULAR E NAO PODE SER INVERTIDA
      IF(MJ.EQ.I+1)GO TO 14
      IAUX=NL(I+1)
      NL(I+1)=NL(MJ)
      NL(MJ)=IAUX
      ISIG=-ISIG
   14 IF(MK.EQ.I+1)GO TO 6
      IAUX=NC(I+1)
      NC(I+1)=NC(MK)
      NC(MK)=IAUX
      ISIG=-ISIG
    6 CONTINUE
      DO 19 I=1,NDIM
      IL=NL(I)
      IC=NC(I)
      DO 19 J=1,NDIM
      JL=NL(J)
      JC=NC(J)
   19 AINV(IC,JL)=A(IL,JC)
ccc   RETURN
ccc   END
      end subroutine invmat
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      REAL*8 FUNCTION RAN(I)
CEQ********************************************************************
CEQ*   SUBROTINA PARA GERAR NUMEROS ALEATORIOS.
CEQ********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      A=I
      B=(24298*A+99991)/199017
      IB=B
      RAN=B-IB
      I=RAN*199017
ccc   RETURN
ccc   END
      end function ran
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE PEQHAR(XX,YY,ZZ,YLM)
      IMPLICIT REAL*8 (A-H,O-Z)
C   DETERMINA AS HARMONICAS ESFERICAS ATE' L=2
      DIMENSION YLM(9)
      DATA R1O4P,R3O4P,R15O4P,R5O4P/ 2.820947917738781D-1,
     * 4.886025119029199D-1,1.092548430592079D0,6.307831305050400D-1/
      R=DSQRT(XX**2+YY**2+ZZ**2)
      IF(R.LT.1.D-10)THEN
        R=1.D0
        X=1.D0
        Y=0.D0
        Z=0.D0
        GO TO 1
        ENDIF
      X=XX/R
      Y=YY/R
      Z=ZZ/R
    1 YLM(1)=R1O4P
      YLM(2)=-R3O4P*Y
      YLM(3)=R3O4P*Z
      YLM(4)=-R3O4P*X
      YLM(5)=R15O4P*X*Y
      YLM(9)=R15O4P*0.5D0*(X**2-Y**2)
      YLM(6)=-R15O4P*Y*Z
      YLM(8)=-R15O4P*X*Z
      YLM(7)=R5O4P*0.5D0*(3.D0*Z**2-1.D0)
ccc   RETURN
ccc   END
      end subroutine peqhar
c----------------------------------------------------------------
 
 
c----------------------------------------------------------------
      SUBROUTINE INTGRL(AA,BB,GTIL)
C!===================================================================!
C!                                                                   !
C!    CALCULA AS INTEGRAIS GTIL(N,LAMBDA,L,A,B) DEFINIDAS POR        !
C!                                                                   !
C!    DEXP(-(A**2-B**2)/4.D0)/2.D0*INTEGRAL(-INFINITO,+INFINITO)     !
C!    r**N+2*DEXP(-r**2)*g(LAMBDA,A*r)*j(L,B*r) dr                   !
C!                                                                   !
C!    AS INTEGRAIS SAO OBTIDAS VIA RECORRENCIA.AS CONDICOES INICIAIS !
C!    SAO                                                            !
C!                                                                   !
C!    GTIL(2,-1,-1),GTIL(1,0,-1),GTIL(1,-1,0),GTIL(0,0,0)            !
C!                                                                   !
C!    COM AS QUAIS SE CALCULAM                                       !
C!                                                                   !
C!    GTIL(2,1,-1),GTIL(2,-1,1)                                      !
C!                                                                   !
C!    SENDO POSSIVEL COMPLETAR A TABELA (PAGINA 19).                 !
C!                                                                   !
C!    AQUI FAZEMOS                                                   !
C!                                                                   !
C!    GTIL(N,LAMBDA,L)=GTIL(I-1,J-1,K-1)                             !
C!                                                                   !
C!    NAO DEVEMOS ESQUECER DA PARIDADE NU+L+LAMBDA=PAR               !
C!                                                                   !
C!===================================================================!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GTIL(5,5,3)
      A=AA
      IF(A.LE.1.D-30)A=0.D0
      B=BB
      IF(B.LE.1.D-30)B=0.D0
      DO 1 I=1,5
      DO 1 J=1,5
      DO 1 K=1,3
    1 GTIL(I,J,K)=0.D0
      PI=4.D0*DATAN(1.D0)
      PIRT=DSQRT(PI)
cmb2guima
cbug  Mudanca de 0.1D0 para 0.2D0 para evitar overflow na GTIL3.
cbug  Isso parece resolver o problema.
      IF(A.LT.0.1D0.OR.B.LT.0.1D0)THEN
        IF(A.LT.0.1D0.AND.B.LT.0.1D0)THEN
          CALL GTIL2(A,B,GTIL)
        ELSE
          CALL GTIL3(A,B,GTIL)
        ENDIF
        RETURN
      ENDIF
      ARG=A*B/2.D0
      BRG=(A**2-B**2)/2.D0
C
C     CONDICOES INICIAIS 
C     GTIL(2,-1,-1)=G11
C     GTIL( 1,0,-1)=G12
C     GTIL( 1,-1,0)=G13
C     GTIL(  0,0,0)=GTIL(1,1,1)
C
      G11=PIRT*((1.D0+BRG)*DCOS(ARG)-A*B*DSIN(ARG))/4.D0/A/B
      G12=PIRT*(DCOS(ARG)/B-DSIN(ARG)/A)/4.D0
      G13=PIRT*(DCOS(ARG)/A+DSIN(ARG)/B)/4.D0
      GTIL(1,1,1)=PIRT*DSIN(ARG)/2.D0/A/B
C
C     AS DUAS PROXIMAS 
C     GTIL(2,1,-1)=G21
C     GTIL(2,-1,1)=G22
C
      G21=G11-G12/A
      G22=G13/B-G11
C
C     PARTIMOS PARA COMPLETAR A TABELA. PRIMEIRO FAZEMOS OS EIXOS
C     LAMBDA E L...
C
      GTIL(3,1,1)=(GTIL(1,1,1)+A*G13+B*G12)/2.D0
      GTIL(2,2,1)=G13-GTIL(1,1,1)/A
      GTIL(2,1,2)=GTIL(1,1,1)/B-G12
      GTIL(4,2,1)=(GTIL(2,2,1)+A*GTIL(3,1,1)+B*G21)/2.D0
      GTIL(4,1,2)=(GTIL(2,1,2)+A*G22+B*GTIL(3,1,1))/2.D0
      GTIL(5,1,1)=(5.D0*GTIL(3,1,1)+A*GTIL(4,2,1)-B*GTIL(4,1,2))/2.D0
      GTIL(3,3,1)=GTIL(3,1,1)-3.D0*GTIL(2,2,1)/A
      GTIL(3,1,3)=3.D0*GTIL(2,1,2)/B-GTIL(3,1,1)
      GTIL(5,3,1)=GTIL(5,1,1)-3.D0*GTIL(4,2,1)/A
      GTIL(5,1,3)=3.D0*GTIL(4,1,2)/B-GTIL(5,1,1)
C
C     ...E AGORA O RESTANTE.
C
      GTIL(3,2,2)=G22-GTIL(2,1,2)/A
      GTIL(1,2,2)=-2.D0*GTIL(3,2,2)+A*GTIL(2,1,2)+B*GTIL(2,2,1)
      GTIL(5,2,2)=(GTIL(3,2,2)+A*GTIL(4,1,2)+B*GTIL(4,2,1))/2.D0
      GTIL(2,3,2)=GTIL(2,1,2)-3.D0*GTIL(1,2,2)/A
      GTIL(4,3,2)=GTIL(4,1,2)-3.D0*GTIL(3,2,2)/A
      GTIL(2,2,3)=3.D0*GTIL(1,2,2)/B-GTIL(2,2,1)
      GTIL(4,2,3)=3.D0*GTIL(3,2,2)/B-GTIL(4,2,1)
      GTIL(3,4,2)=GTIL(3,2,2)-5.D0*GTIL(2,3,2)/A
      GTIL(5,4,2)=GTIL(5,2,2)-5.D0*GTIL(4,3,2)/A
      GTIL(3,3,3)=3.D0*GTIL(2,3,2)/B-GTIL(3,3,1)
      GTIL(1,3,3)=(-2.D0*GTIL(3,3,3)+A*GTIL(2,2,3)+B*GTIL(2,3,2))/3.D0
      GTIL(5,3,3)=(-GTIL(3,3,3)+A*GTIL(4,2,3)+B*GTIL(4,3,2))/2.D0
      GTIL(2,4,3)=GTIL(2,2,3)-5.D0*GTIL(1,3,3)/A
      GTIL(4,4,3)=GTIL(4,2,3)-5.D0*GTIL(3,3,3)/A
      GTIL(3,5,3)=GTIL(3,3,3)-7.D0*GTIL(2,4,3)/A
      GTIL(5,5,3)=GTIL(5,3,3)-7.D0*GTIL(4,4,3)/A
ccc   RETURN
ccc   END
      end subroutine INTGRL
c----------------------------------------------------------------


c----------------------------------------------------------------
      SUBROUTINE GTIL3(AA,BB,GTIL)
c  Calcula GTIL no caso de se usar regra de recorrencia para o maior
c  parametro AA ou BB, e serie de potencia para o menor BB ou AA.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GTIL(5,5,3),HTIL(5,5,5)
      DIMENSION ANT(2),POS(6),D(6),F(7,5),FIN(5),
     * KEY(10,6,6)
      DO 13 I=1,5
      DO 13 J=1,5
      DO 11 K=1,5
      KEY(I,J,K)=0
   11 HTIL(I,J,K)=0.D0
      DO 13 K=1,3
   13 GTIL(I,J,K)=0.D0
      PI=4.D0*DATAN(1.D0)
      IF(AA.EQ.0.D0.AND.BB.EQ.0.D0)THEN
        FATOR=DSQRT(PI)/4.D0
        DO 303 N=1,5,2
        GTIL(N,1,1)=FATOR
  303   FATOR=FATOR/2.D0*(N+2)
        RETURN
      ENDIF
      A=AA
      B=BB
      ISIG=1
      IF(B.GT.A)THEN
        B=AA
        A=BB
        ISIG=-1
      ENDIF
c  Inicializa os fatores F(M,N,m,B)
      FIN(1)=DEXP(ISIG*B**2/4.D0)
      DO 4 M=2,5
    4 FIN(M)=FIN(M-1)*B/(2*M-1)
      MIL=10000
      IF(B.EQ.0.D0)MIL=10
      DO 301 MG=1,MIL
c  Calcula os fatores D(M,l,A)
c  D(L)=D(MG-1,L-2,A)
c  POS(L)=D(MG,L-2,A)
c  ANT(L)=D(MG-2,L-2,A) somente para L=1 e L=2
      IF(MG.EQ.1)THEN
        ANT(2)=0.D0
        D(2)=DSQRT(PI)/4.D0
        ANT(1)=2.D0/A*D(2)
        D(1)=ISIG*A/2.D0*ANT(2)
        POS(1)=0.5D0*ANT(1)+ISIG*A/2.D0*D(2)
        POS(2)=0.5D0*ANT(2)+A/2.D0*D(1)
        DO 302 L=3,6
          POS(L)=ISIG*(POS(L-2)-(2*L-5)/A*D(L-1))
          IF(L.EQ.4)THEN
            D(4)=(ISIG-6.D0/A**2)*D(2)
            GO TO 302
          ENDIF
          D(L)=ISIG*((L-1)*D(L-2)/(4-L)-2*(2*L-5)/A/(4-L)*POS(L-1))
  302   CONTINUE
      ELSE
cmb2guima
cbug   Aparentemente o problema de overflow acontece no calculo
cbug   de POS(1) e POS(2), causado por ANT(1) e ANT(2).
        POS(1)=(MG-1)/2.D0*ANT(1)+ISIG*A/2.D0*D(2)
        POS(2)=(MG-1)/2.D0*ANT(2)+A/2.D0*D(1)
        DO 309 L=3,6
  309   POS(L)=ISIG*(POS(L-2)-(2*L-5)/A*D(L-1))
        ANT(1)=D(1)
        ANT(2)=D(2)
        DO 304 L=1,6
  304   D(L)=POS(L)
      ENDIF
c   Calcula os fatores F(M,N,m,B)
      DO 305 N=1,5
      DO 305 M=1,5
      MNM=MG-N-M+1
      IF(MNM.LT.0)GO TO 305
      IF(MNM.NE.2*(MNM/2))GO TO 305
      IF(MNM.EQ.0)THEN
        F(N,M)=FIN(M)
      ELSE
        F(N,M)=-ISIG*F(N,M)*B**2/MNM/(MNM+2*M-1)
      ENDIF      
      DO 306 L=1,5
      TERMO=F(N,M)*D(L+1)
      IF(N+L+M.EQ.2*((N+L+M)/2))KEY(N,L,M)=1
      IF(BB.LE.AA)THEN
        HTIL(N,L,M)=HTIL(N,L,M)+TERMO
        IF(DABS(TERMO).LT.1.D-12*DABS(HTIL(N,L,M)))KEY(N,L,M)=1
      ELSE
        HTIL(N,M,L)=HTIL(N,M,L)+TERMO
        IF(DABS(TERMO).LT.1.D-12*DABS(HTIL(N,M,L)))KEY(N,M,L)=1
      ENDIF
  306 CONTINUE
  305 CONTINUE
      DO 307 N=1,5
      DO 307 M=1,5
      DO 307 L=1,5
      IF(KEY(N,L,M).EQ.0)GO TO 301
  307 CONTINUE
      GO TO 310
  301 CONTINUE
  310 CONTINUE
      DO 308 N=1,5
      DO 308 L=1,5
      DO 308 M=1,3
      IF(N+L+M.EQ.2*((N+L+M)/2))THEN
        GTIL(N,L,M)=0.D0
      ELSE
        GTIL(N,L,M)=HTIL(N,L,M)
      ENDIF
  308 CONTINUE
      GTIL(5,5,1)=0.D0
      GTIL(4,5,2)=0.D0
      GTIL(4,4,1)=0.D0
ccc   RETURN
ccc   END
      end subroutine gtil3
c----------------------------------------------------------------


c----------------------------------------------------------------
      SUBROUTINE GTIL2(AA,BB,GTIL)
c  Calcula GTIL no caso de se expandir em serie de potencias de AA e BB.
c  Esta subrotina 'e 'util quando AA, BB < 1
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GTIL(5,5,3),HTIL(5,5,5)
      DO 13 I=1,5
      DO 13 J=1,5
      DO 11 K=1,5
   11 HTIL(I,J,K)=0.D0
      DO 13 K=1,3
   13 GTIL(I,J,K)=0.D0
      PI=4.D0*DATAN(1.D0)
      IF(AA.EQ.0.D0.AND.BB.EQ.0.D0)THEN
        FATOR=DSQRT(PI)/4.D0
        DO 303 N=1,5,2
        GTIL(N,1,1)=FATOR
  303   FATOR=FATOR/2.D0*(N+2)
        RETURN
      ENDIF
      A=AA
      B=BB
      IF(B.GT.A)THEN
        B=AA
        A=BB
      ENDIF
      DO 301 N=1,5
      DO 301 L=1,5
      DO 301 M=1,5
      IF(N+L+M.EQ.2*((N+L+M)/2))GO TO 301
      NLM=N+L+M-3
      BM=1.D0
      SOMA=0.D0
      IF(B.LT.1.D-10.AND.M.GT.1)GO TO 6
      IF(M.GT.1)BM=B**(M-1)
      FATOR=DSQRT(PI)/4.D0/2.D0**(NLM/2)*A**(L-1)*BM*
     * DEXP(-0.25D0*(AA**2-BB**2))
C  (N+l+m+1)!!
      DO 1 I=1,NLM/2+1
    1 FATOR=FATOR*(2*(I-1)+1)
      FATOR=FATOR/(NLM+1)
c  (2l+1)!!
      DO 2 I=1,L
    2 FATOR=FATOR/(2*I-1)
      FATOR=FATOR*(2*L-1)
c  (2m+1)!!
      DO 3 I=1,M
    3 FATOR=FATOR/(2*I-1)
      FATOR=FATOR*4.D0/A**2
      ISIG=1
      DO 4 I=1,1000
      FATOR=FATOR*(NLM+2*I-1)/(2*L+2*I-3)*A**2/4.D0
      SUM=1.D0
      IF(I.GT.1)THEN
        FATOR=FATOR/(I-1)
        FAT=1.D0
        IF(B.GT.1.D-10)THEN
          DO 5 J=2,I
          FAT=-FAT*(B/A)**2*(I-J+1)*(2*L+2*I-2*J+1)/(J-1)/(2*M+2*J-3)
    5     SUM=SUM+FAT
        ENDIF
      ENDIF
      TERMO=FATOR*ISIG*SUM
      SOMA=SOMA+TERMO
      IF(AA.LT.BB)ISIG=-ISIG
      IF(DABS(TERMO).LT.1.D-12)GO TO 6
    4 CONTINUE
    6 IF(BB.LE.AA)HTIL(N,L,M)=SOMA
      IF(BB.GT.AA)HTIL(N,M,L)=SOMA
  301 CONTINUE
      DO 308 N=1,5
      DO 308 L=1,5
      DO 308 M=1,3
      IF(N+L+M.EQ.2*((N+L+M)/2))THEN
        GTIL(N,L,M)=0.D0
      ELSE
        GTIL(N,L,M)=HTIL(N,L,M)
      ENDIF
  308 CONTINUE
      GTIL(5,5,1)=0.D0
      GTIL(4,5,2)=0.D0
      GTIL(4,4,1)=0.D0
ccc   RETURN
ccc   END
      end subroutine gtil2
c----------------------------------------------------------------

      end module cintg_pp
