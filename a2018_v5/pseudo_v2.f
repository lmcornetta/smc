      module pseudo_pp

      contains

c----------------------------------------------------------------
      REAL*8 FUNCTION PSEUDO(LMAX,ZVAL,CORE,ALFCOR,AION,ALFION,CENTRO,
     1 LEXP,ALFA)
C      CALCULA O ELEMENTO DE MATRIZ DO PSEUDO POTENCIAL ENTRE DUAS
C   GAUSSIANAS CARTESIANAS I (= 1 OU 2) DO TIPO
C       (X-CENTRO(1,I))**LEXP(1,I)
C      *(Y-CENTRO(2,I))**LEXP(2,I)
C      *(Z-CENTRO(3,I))**LEXP(3,I)
C      *EXP( - ALFA(I) * (R-CENTRO(I))**2 )
C
C     O POTENCIAL E' DADO NA PAG. 4208 DO ARTIGO 
C   G.B. BACHELET, D.R. HAMANN, AND M. SCHLUTER, PHYS. REV. B 26, 
C   4199 (1982).
C     A EXPRESSAO DO PSEUDO-POTENCIAL E'
C  PSEUDO(R) = VCORE(R) + SOMA(L=0,2) ( VION(R,L) * PROJETOR(L) )
C     ONDE
C VCORE(R) = - (ZVAL/R) * SOMA(I=1,2) (CORE(I) * ERF(R*ALFCOR(I)**0.5))
C
C VION(R,L) = SOMA(J=1,3) ( (AION(1,J,L)+AION(2,J,L)*R**2) *
C                            EXP(-ALFION*R**2) )
C
C     'LMAX' DEVE SER FORNECIDO. SIGNIFICA A MAXIMA SOMA DOS EXPOENTES
C  DO POLINOMIO DA GAUSSIANA CARTESIANA QUE JAMAIS SERA' USADA.
C
C
C                  ATENCAO
C                  ATENCAO
C 
C             YOUR ATTENTION, PLEASE
C             YOUR ATTENTION, PLEASE
C
C                  ATTENTION
C                  ATTENTION
C
C     DESTE PROGRAMA CONSTA A REAL*8 FUNCTION GERF (DE Guima-ERF ).
C  NADA MAIS E' QUE A VELHA FUNCAO ERF CALCULADA COM PRECISAO DE 1
C  PARTE EM 5 * 10**8. SE SEU COMPUTADOR POSSUI UMA FUNCAO ERF 
C  INTRINSICA DO FORTRAN COM PRECISAO MAIOR, USE-A. PARA ISTO, ELIMINA
C  A FUNCAO 'GERF' E , NAS CHAMADAS, SUBSTITUA 'GERF' POR 'ERF' OU
C  'DERF', OU PELO NOME COMO A FUNCAO INTRINSICA E' CHAMADA.
C
C  PRECIS = VALOR MINIMO DAS INTEGRAIS QUE SAO CALCULADAS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CORE(2),ALFCOR(2),AION(2,3,3),ALFION(3,3),
     1 CENTRO(3,2),LEXP(3,2),ALFA(2),LA(2),MA(2),NA(2)
C      EXPMIN=-DLOG(PRECIS)
      PRECIS=1.5D+7
      EXPMIN=PRECIS
      PSEUDO=VCORE(EXPMIN,LMAX,ZVAL,CORE,ALFCOR,CENTRO,ALFA,LEXP)      
      DO 1 I=1,2
      LL=0
      LA(I)=LEXP(1,I)
      MA(I)=LEXP(2,I)
      NA(I)=LEXP(3,I)
      LL=LL+LA(I)
      LL=LL+MA(I)
      LL=LL+NA(I)
C   SE A SOMA DOS EXPOENTES FOR SUPERIOR A 2, A CONTRIBUICAO DO
C   ION E' FEITA IGUAL A ZERO.
      IF(LL.GT.2)RETURN
    1 CONTINUE
      PSEUDO=PSEUDO+ELEMEN(EXPMIN,AION,ALFION,CENTRO,ALFA,LA,MA,NA)
      RETURN
ccc   END
      end function pseudo
c----------------------------------------------------------------

c----------------------------------------------------------------
      REAL*8 FUNCTION VCORE(EXPMIN,LMAX,ZVAL,C,RO,A,ALFA,LEXP)
C   ACHA O ELEMENTO DE MATRIZ DO V(CORE) ENTRE DUAS GUASSIANAS
C     CARTESIANAS.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(2),RO(2),A(3,2),ALFA(2),LEXP(3,2),LL(2),VET(3),
     1 YLM(81),Q(81,5),SOMAK(9,9,2)
C   LMAX = SOMA DOS EXPOENTES DAS COMPONENTES X, Y, Z MAXIMA PARA
C          TODOS OS CENTROS E TODAS AS GAUSSIANAS CARTESIANAS.
      VCORE=0.D0
      AMOD=0.D0
      VMOD=0.D0
      DO 301 I=1,3
      AMOD=AMOD+A(I,1)**2*ALFA(1)+A(I,2)**2*ALFA(2)
      VET(I)=ALFA(1)*A(I,1)+ALFA(2)*A(I,2)
  301 VMOD=VMOD+VET(I)**2
      AUX=AMOD-VMOD/(ALFA(1)+ALFA(2))
      IF(AUX.GT.EXPMIN)RETURN
      ABMOD=VMOD
      PI=4.D0*DATAN(1.D0)
      DO 1 I=1,2
      LL(I)=0
      DO 1 J=1,3
    1 LL(I)=LL(I)+LEXP(J,I)
      VMOD=2.D0*DSQRT(VMOD)
      CALL HARMON(PI,LL(1)+LL(2),VET(1),VET(2),VET(3),YLM(1))
      CALL FATORQ(LMAX,A,LEXP,Q)
      SAB=1.D0/DSQRT(ALFA(1)+ALFA(2))
      Y=SAB*VMOD
      DO 302 I=1,2
  302 CALL KTIL(LL(1)+LL(2),EXPMIN,RO(I)/(ALFA(1)+ALFA(2)),Y,
     * SOMAK(1,1,I))
      VCORE=0.D0
      DO 3 L=1,LL(1)+LL(2)+1
      LAM=L-1
      MAXT=(LL(1)+LL(2)+LAM)/2-LAM+1
      DO 6 IT=1,MAXT
      N=IT-1+LAM
      SOMA=0.D0
      DO 4 MU=LAM**2+1,(LAM+1)**2
    4 SOMA=SOMA+Q(MU,IT)*YLM(MU)
      TOMA=0.D0
      DO 5 I=1,2
    5 TOMA=TOMA+C(I)*SOMAK(N+1,LAM+1,I)
    6 VCORE=VCORE+SOMA*TOMA*SAB**(2*N-LAM+2)
    3 CONTINUE
      VCORE=-VCORE*ZVAL*DEXP(-AUX)
      RETURN
ccc   END
      end function vcore
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE FATORQ(LMAX,A,LEXP,Q)
C   ACHA O FATOR Q DO ELEMENTO DE MATRIZ DA PARTE LOCAL
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/JNVER/AINIC,EME(15333),ENE(55),RANV(3,81),RK(5),L2(5),
     1 L4(9),LDEMU(81)
      DIMENSION NL(81),NC(81),AMAT(6561),ANE(5,5),YMAT(81,81),A(3,2),
     1 LEXP(3,2),Q(81,5),LL(2)
      IF(AINIC.EQ.0.123456789D0)GO TO 100
      IF(LMAX.GT.4)STOP ' LMAX DEVE SER 4 NO MAXIMO'
C  ESTABELECE O INICIO DAS MATRIZES 'M' E 'N' NOS VETORES 'EME' E 'ENE'
      LS=0
      DO 2 L=1,5
      LS=LS+L**2
    2 L2(L)=LS-L**2+1
      LS=0
      DO 3 L=1,9
      LS=LS+L**4
    3 L4(L)=LS-L**4+1
C  ESTABELECE A MATRIZ DAS POTENCIAS (R(K)**2)**T E INVERTE
      DO 4 K=1,5
      RS=0.5D0*K
      RK(K)=DSQRT(RS)
      X=1.D0
      DO 4 IT=1,5
      ANE(IT,K)=X
    4 X=X*RS
      DO 5 L=1,LMAX+1
      DO 6 I=1,L
      DO 6 J=1,L
    6 AMAT(I+(J-1)*L)=ANE(I,J)
    5 CALL INVMAT(1,L,L,NL,NC,ISIG,AMAT,ENE(L2(L)))
C  ESTABELECE A MATRIZ DAS HARMONICAS ESFERICAS DE DIRECOES ALEATORIAS
C    E INVERTE
      PI=4.D0*DATAN(1.D0)
      NALEAT=0
      DO 9 IR=1,(2*LMAX+1)**2
      RANV(3,IR)=-1.D0+2.D0*RAN(NALEAT)
      FI=-PI+2.D0*PI*RAN(NALEAT)
      AUX=DSQRT(1.D0-RANV(3,IR)**2)
      RANV(1,IR)=AUX*DCOS(FI)
      RANV(2,IR)=AUX*DSIN(FI)
    9 CALL HARMON(PI,2*LMAX,RANV(1,IR),RANV(2,IR),RANV(3,IR),
     1 YMAT(1,IR))
      MU=0
      DO 7 L=1,2*LMAX+1
      NDIM=L**2
      DO 8 I=1,NDIM
      DO 8 J=1,NDIM
    8 AMAT(I+(J-1)*NDIM)=YMAT(I,J)
      CALL INVMAT(1,NDIM,NDIM,NL,NC,ISIG,AMAT,EME(L4(L)))
      DO 7 M=1,2*L-1
      MU=MU+1
    7 LDEMU(MU)=L-1
      AINIC=0.123456789D0
C  AQUI TERMINA A INICIALIZACAO E COMECA A EXECUCAO PROPRIAMENTE
C  ACHA A SOMA 'LL' DOS EXPOENTES
  100 DO 11 I=1,2
      LL(I)=0
      DO 11 J=1,3
   11 LL(I)=LL(I)+LEXP(J,I)
C  MAXJ = DIMENSAO DA MATRIZ 'M'
C  MAXT = DIMENSAO DA MATRIZ 'N'
      MAXJ=(LL(1)+LL(2)+1)**2
      LOINT=L4(LL(1)+LL(2)+1)-1
C   PARA CADA 'J' E CADA 'T' CALCULA O VALOR DE 'Q'
      DO 12 J=1,MAXJ
      MAXT=(LL(1)+LL(2)+LDEMU(J))/2-LDEMU(J)+1
      MOINT=LOINT+(J-1)*MAXJ
      IOINT=L2(MAXT)-1
      DO 12 IT=1,MAXT
      JOINT=IOINT+(IT-1)*MAXT
      Q(J,IT)=0.D0
      DO 15 K=1,MAXT
      RP=1.D0/RK(K)**LDEMU(J)
      AUX=0.D0
      DO 13 I=1,MAXJ
      BUX=1.D0
      DO 17 IB=1,2
      DO 14 IE=1,3
      IF(LEXP(IE,IB).EQ.0)GO TO 14
      DO 16 NADA=1,LEXP(IE,IB)
   16 BUX=BUX*(RK(K)*RANV(IE,I)-A(IE,IB))
   14 CONTINUE
   17 CONTINUE
   13 AUX=AUX+BUX*EME(MOINT+I)
   15 Q(J,IT)=Q(J,IT)+AUX*RP*ENE(JOINT+K)
   12 CONTINUE
      RETURN
ccc   END
      end subroutine fatorq
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE HARMON(PI,LMAX,AX,AY,AZ,YLM)
C   CALCULA AS HARMONICAS ESFERICAS REAIS ('YLM').
c  Versao de agosto de 1995
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YLM((LMAX+1)**2),FAT(16,16),P(17,16),C(31),Q(16,16,16)
      LBOBO=LMAX+1
      IF(LBOBO.LT.3) LBOBO=3
      DO 11 L=1,LBOBO
      DO 11 M=1,L
      DO 11 N=1,L+1-M
   11 Q(L,M,N)=0.D0
      Q(1,1,1)=1.D0
      Q(2,1,2)=1.D0
      DO 12 L=3,LBOBO
      DO 12 N=1,L
      IF(N.GT.1)Q(L,1,N)=(2*L-3)*Q(L-1,1,N-1)/(L-1)
   12 Q(L,1,N)=Q(L,1,N)-(L-2)*Q(L-2,1,N)/(L-1)
      DO 13 L=2,LBOBO
      DO 13 M=2,L
      DO 13 N=1,L+1-M
   13 Q(L,M,N)=Q(L,M-1,N+1)*N
      C(16)=1.D0/DSQRT(2.D0)
      AUX=AX**2+AY**2+AZ**2
      IF(AUX.NE.0.D0)GO TO 345
C   SE O VETOR FOR NULO, ARBITRA A DIRECAO (1 0 0)
      AUX=1.D0
      X=0.D0
      SQ=1.D0
      C(17)=1.D0
      C(15)=0.D0
      GO TO 346
  345 AUX=DSQRT(AUX)
      X=AZ/AUX
      SQ=1.D0-X**2
      IF(SQ.GT.0.D0)GO TO 347
C   NESTE CASO ARBITRA O ANGULO FI COMO ZERO
      C(17)=1.D0
      C(15)=0.D0
      SQ=0.D0
      GO TO 346
  347 SQ=DSQRT(SQ)
      C(17)=AX/AUX/SQ
      C(15)=AY/AUX/SQ
  346 DO 1 L=1,LBOBO
      SQN=1.D0
      DO 7 M=1,L
      P(L,M)=0.D0
      XN=1.D0
      DO 2 N=1,L+1-M
      P(L,M)=P(L,M)+Q(L,M,N)*XN
    2 XN=XN*X
      P(L,M)=P(L,M)*SQN
    7 SQN=-SQN*SQ
    1 CONTINUE
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
      RETURN
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
      RETURN
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
      RETURN
ccc   END
      end function ran
c----------------------------------------------------------------

c----------------------------------------------------------------
      REAL*8 FUNCTION ELEMEN(EXPMIN,C,SIGMA,A,ALFA,LA,MA,NA)
C      CALCULA O EMENTO DE MATRIZ DO PSEUDO-POTENCIAL DA ORIGEM ENTRE
C   AS FUNCOES
C     X**LA(1)*Y**MA(1)*Z**NA(1)*EXP(-ALFA(1)*R**2)
C     CENTRADA EM (A(1,1) A(2,1) A(3,1))
C   E
C     X**LA(2)*Y**MA(2)*Z**NA(2)*EXP(-ALFA(2)*R**2)
C     CENTRADA EM (A(1,2) A(2,2) A(3,2))
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ESTR/AJA(3,11),QJA(1655,11),YJA(25,11),PI,
     1 COMECO,MARQUI,NARQUI,INVIND(1655),NUQRED,INDERE(6750)
      DIMENSION A(3,2),ALFA(2),LA(2),MA(2),NA(2),Q(6750,2),
     1 B(3),C(2,3,3),SIGMA(3,3),AJOTA(360,3,3),Y(25,2),AEXP(3,3),
     2 AA(3,3),BB(3,3),MREP(10),LBRA(2),MBRA(2),XFACT(3)
      DATA MREP/0,100,10,1,200,20,2,11,101,110/
      ELEMEN=0.D0
      AMOD=DSQRT(A(1,1)**2+A(2,1)**2+A(3,1)**2)
      BMOD=DSQRT(A(1,2)**2+A(2,2)**2+A(3,2)**2)
      AUX=ALFA(1)*AMOD**2+ALFA(2)*BMOD**2
      BUX=(ALFA(1)*AMOD+ALFA(2)*BMOD)**2
      CUX=ALFA(1)+ALFA(2)
      DO 301 J=1,3
      DO 301 L=1,3
      IF(AUX-BUX/(CUX+SIGMA(J,L)).LT.EXPMIN)GO TO 302
  301 CONTINUE
      RETURN
  302 CONTINUE
      IF(COMECO.EQ.0.123456789D0)GO TO 202
      KONT=0
      DO 215 LL=1,25
      LAM=0
      IF(LL.GT.1)LAM=1
      IF(LL.GT.4)LAM=2
      IF(LL.GT.9)LAM=3
      IF(LL.GT.16)LAM=4
      DO 215 KK=1,3
      DO 215 JJ=1,9
      LP=0
      IF(JJ.GT.1)LP=1
      IF(JJ.GT.4)LP=2
      DO 215 II=1,10
      LG=0
      IF(II.GT.1)LG=1
      IF(II.GT.4)LG=2
      IAUX=II+(JJ-1)*10+(KK-1)*90+(LL-1)*270
      INDERE(IAUX)=0
      IF(KK-1.LT.IABS(LAM-LP))GO TO 215
      IF(KK.GT.LG+1)GO TO 215
      IF(LAM+LP+KK-1.NE.2*((LAM+LP+KK-1)/2))GO TO 215
      KONT=KONT+1
      INVIND(KONT)=IAUX
      INDERE(IAUX)=KONT
      IF(LAM+LP+LG.NE.2*((LAM+LP+LG)/2))INDERE(IAUX)=-INDERE(IAUX)
  215 CONTINUE
      NUQRED=KONT
      MARQUI=0
      NARQUI=0
      PI=4.D0*DATAN(1.D0)
      COMECO=0.123456789D0
  202 DO 1 IBRA=1,2
C   PROCURA OS FATORES DE ESTRUTURA 'Q' E AS HARMONICAS ESFERICAS 'Y' 
C     NO COMMON 
      CALL LESCQY(0,JACHOU,A(1,IBRA),Q(1,IBRA),Y(1,IBRA))
      IF(JACHOU.NE.0)GO TO 1
C   CALCULA O FATOR DE ESTRUTURA E AS HARMONICAS ESFERICAS
      CALL ESTRUT(A(1,IBRA),Q(1,IBRA),INDERE(1))
      CALL HARMON(PI,4,A(1,IBRA),A(2,IBRA),A(3,IBRA),Y(1,IBRA))
C   ESCREVE NO COMMON
      CALL LESCQY(1,JACHOU,A(1,IBRA),Q(1,IBRA),Y(1,IBRA))
    1 CONTINUE
C   ACHA A FUNCAO JOTA
      DO 17 L=1,3
      DO 17 IND=1,3
      IF(C(1,IND,L).EQ.0.D0.AND.C(2,IND,L).EQ.0.D0)GO TO 17
      AEXP(IND,L)=DSQRT(ALFA(1)+ALFA(2)+SIGMA(IND,L))
      AA(IND,L)=2.D0*ALFA(1)*DSQRT(A(1,1)**2+A(2,1)**2+A(3,1)**2)/
     1 AEXP(IND,L)
      IF(AA(IND,L).LT.1.D-8)AA(IND,L)=0.D0
      BB(IND,L)=2.D0*ALFA(2)*DSQRT(A(1,2)**2+A(2,2)**2+A(3,2)**2)/
     1 AEXP(IND,L)
      IF(BB(IND,L).LT.1.D-8)BB(IND,L)=0.D0
      CALL JOTA(AA(IND,L),BB(IND,L),AJOTA(1,IND,L))
   17 CONTINUE
C   FAZ O CALCULO DO ELEMENTO DE MATRIZ
      DO 21 IBRA=1,2
      LBRA(IBRA)=LA(IBRA)+MA(IBRA)+NA(IBRA)
      I=100*LA(IBRA)+10*MA(IBRA)+NA(IBRA)
      DO 221 J=1,10
      IF(MREP(J).NE.I)GO TO 221
      MBRA(IBRA)=J
      GO TO 21
  221 CONTINUE
   21 CONTINUE
      ELEMEN=0.D0
      DO 22 L=1,3
      NADA=0
      DO 32 IND=1,3
      XFACT(IND)=0.D0
      IF(C(1,IND,L).EQ.0.D0.AND.C(2,IND,L).EQ.0.D0)GO TO 32
      CUX=(AA(IND,L)+BB(IND,L))**2/4.D0-ALFA(1)*(A(1,1)**2+A(2,1)**2+
     1 A(3,1)**2)-ALFA(2)*(A(1,2)**2+A(2,2)**2+A(3,2)**2)
      IF(DABS(CUX).GT.EXPMIN)GO TO 32
      XFACT(IND)=4.D0*PI*DEXP(CUX)
      NADA=1
   32 CONTINUE
      IF(NADA.EQ.0)GO TO 22
      LAM1MX=LBRA(1)+L
      IF(AA(1,1).EQ.0.D0)LAM1MX=1
      DO 30 LAM1=1,LAM1MX
      IF(IABS(LAM1-L).GT.LBRA(1))GO TO 30
      LAM2MX=LBRA(2)+L
      IF(BB(1,1).EQ.0.D0)LAM2MX=1
      DO 31 LAM2=1,LAM2MX
      IF(IABS(LAM2-L).GT.LBRA(2))GO TO 31
      DO 24 NUSUM=IABS(LAM1-L)+IABS(LAM2-L)+1,LBRA(1)+LBRA(2)+1
      MAUX=NUSUM+LAM1+LAM2-3
      IF(MAUX.NE.2*(MAUX/2))GO TO 24
      MIN=IABS(LAM1-L)+1
      MAUX=NUSUM-LBRA(2)
      IF(MAUX.GT.MIN)MIN=MAUX
      MAX=LBRA(1)+1
      MAUX=NUSUM-IABS(LAM2-L)
      IF(MAUX.LT.MAX)MAX=MAUX
      IF(MIN.GT.MAX)GO TO 24
      PLEMEN=0.D0
      DO 33 IND=1,3
      IF(XFACT(IND).EQ.0.D0)GO TO 33
      DO 26 N=1,2
   26 PLEMEN=PLEMEN+C(N,IND,L)*AJOTA(2*N+NUSUM-2+(LAM1-1)*10+
     1 (LAM2-1)*60,IND,L)/AEXP(IND,L)**(2*N+NUSUM)*XFACT(IND)
   33 CONTINUE
      QLEMEN=0.D0
      DO 27 MU1=(LAM1-1)**2+1,LAM1**2
      IMU1=MBRA(1)+(MU1-1)*270
      RLEMEN=0.D0
      DO 28 MU2=(LAM2-1)**2+1,LAM2**2
      IMU2=MBRA(2)+(MU2-1)*270
      SLEMEN=0.D0
      DO 29 M=(L-1)**2+1,L**2
      IM1=IMU1+(M-1)*10
      IM2=IMU2+(M-1)*10
      DO 29 NU1=MIN,MAX
      NU2=NUSUM-NU1+1
      INU1=IM1+(NU1-1)*90
C   VERIFICA A REGRA DE PARIDADE PARA NU
      IF(INDERE(INU1).EQ.0)GO TO 29
      INU2=IM2+(NU2-1)*90
      SLEMEN=SLEMEN+Q(INU1,1)*Q(INU2,2)
   29 CONTINUE
   28 RLEMEN=RLEMEN+SLEMEN*Y(MU2,2)
   27 QLEMEN=QLEMEN+RLEMEN*Y(MU1,1)
      ELEMEN=ELEMEN+PLEMEN*QLEMEN
   24 CONTINUE
   31 CONTINUE
   30 CONTINUE
   22 CONTINUE
      RETURN
ccc   END
      end function elemen
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE LESCQY(KEY,JACHOU,A,Q,Y)
C   LE NO COMMON ( SE KEY=0 ) OU ESCREVE NO COMMON (SE KEY=1) OS 
C     FATORES 'Q' E 'Y' PARA DADO VETOR DE POSICAO 'A'.
C     MARQUI CONTA O NUMERO JA' ACHADO. SO' OS 11 ULTIMOS SAO ACUMULADOS
C     NO COMMON.
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ESTR/AJA(3,11),QJA(1655,11),YJA(25,11),PI,
     1 COMECO,MARQUI,NARQUI,INVIND(1655),NUQRED,INDERE(6750)
      DIMENSION A(3),Q(6750),Y(25)
      IF(KEY.NE.0)GO TO 100
      JACHOU=0
      IF(MARQUI.EQ.0)RETURN
C   VERIFICA SE O BRA OU O KET JA' ESTA' NA MEMORIA
      MAXIMO=MARQUI
      IF(MARQUI.GT.11)MAXIMO=11
      DO 3 I=1,MAXIMO
      DO 4 J=1,3
      IF(DABS(A(J)-AJA(J,I)).LT.0.1D-8)GO TO 4
      GO TO 5
    4 CONTINUE
C   JA' ESTA' NA MEMORIA
      DO 6 IAUX=1,NUQRED
      LL=INVIND(IAUX)
    6 Q(LL)=QJA(IAUX,I)
      DO 106 LL=1,25
  106 Y(LL)=YJA(LL,I)
      JACHOU=1
      RETURN
C   NAO ESTA' NA MEMORIA. ENTAO VERIFICA SE O VETOR NEGATIVO ESTA'
C   NA MEMORIA.
    5 DO 7 J=1,3
      IF(DABS(A(J)+AJA(J,I)).LT.0.1D-8)GO TO 7
      GO TO 3
    7 CONTINUE
C   O VETOR NEGATIVO ESTA' NA MEMORIA.
      DO 9 IAUX=1,NUQRED
      LL=INVIND(IAUX)
      Q(LL)=QJA(IAUX,I)
      IF(INDERE(LL).LT.0)Q(LL)=-Q(LL)
    9 CONTINUE
      LSIG=-1
      DO 109 LL=1,5
      LSIG=-LSIG
      DO 109 MM=(LL-1)**2+1,LL**2
  109 Y(MM)=LSIG*YJA(MM,I)
      JACHOU=1
      RETURN
    3 CONTINUE
      RETURN
C   ESCREVE NO COMMON
  100 CONTINUE
      MARQUI=MARQUI+1
      MAR=MARQUI-((MARQUI-1)/11)*11
      DO 116 LL=1,25
  116 YJA(LL,MAR)=Y(LL)
      DO 16 KK=1,3
  16  AJA(KK,MAR)=A(KK)
      DO 216 IAUX=1,NUQRED
      LL=INVIND(IAUX)
  216 QJA(IAUX,MAR)=Q(LL)
      RETURN
ccc   END
      end subroutine lescqy
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE ESTRUT(A,Q,INDERE)
C     CALCULA OS FATORES DE ESTRUTURA Q(MG,M,NU,MU) PARA O VETOR
C   (A(1) A(2) A(3))
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SALVA/AINIC,AINV(25,25),YP(9,25),RANV(3,25),PI
      DIMENSION A(3),Q(6750),MREP(10),VET(3),QLM(10,3,25),
     1 AMAT(25,25),NL(25),NC(25),INDERE(6750)
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
      DO 3 MG=1,10
      DO 3 M=1,9
      DO 3 MU=1,25
      IAUX=MG+(M-1)*10+(MU-1)*270
      JAUX=0
      DO 6 N=1,3
      Q(IAUX+(N-1)*90)=0.D0
    6 JAUX=JAUX+IABS(INDERE(IAUX+(N-1)*90))
      IF(JAUX.EQ.0)GO TO 3
      DO 4 N=1,3
      VET(N)=0.D0
      DO 4 IR=1,25
    4 VET(N)=VET(N)+AINV(IR,MU)*YP(M,IR)*QLM(MG,N,IR)
      VET(3)=0.5D0*(VET(3)+VET(1)-2.D0*VET(2))/AMOD**2
      VET(2)=(VET(2)-VET(1)-VET(3)*AMOD**2)/AMOD
      DO 7 N=1,3
    7 Q(IAUX+(N-1)*90)=VET(N)
    3 CONTINUE
      RETURN
ccc   END
      end subroutine estrut
c----------------------------------------------------------------

C     Marcio:
c     Voce estava certo. Havia o erro da exponencial sendo
c     dividida e nao multiplicada.
c     Para evitar que a exponencial exploda, cancelei juntamente
c     com a parte em F, conforme voce sugeriu.
c     O 1.D-30 que sugeri ser 1.D-10, voltei para 1.D-20. Acho
c     que qualquer VAX aguenta isto.
c
c
c     Marco e Alexandra:
c     O Marcio ainda esta' em plena forma. Descobriu um erro
c     na subrotina JOTA3. Favor corrigir.
c
c
c     Guima, 24 de dezembro de 1995
c
c     P.S. Feliz Natal.



c----------------------------------------------------------------
      SUBROUTINE JOTA3(AA,BB,YJ)
c Programa para argumentos AA>1 e BB<1
C   ACHA A TABELA DE FUNCOES J(N,L,L',A,B).
C     ONDE
C  J(N,L,L',A,B) = EXP(-(A+B)**2/4)*INTEGRAL(0,INFINITO)(r**(N+2)*
C       EXP(-r**2)*I(L,A*r)*I(L',B*r)*dr)
C     ONDE I(L,x) E' A FUNCAO DE BESSEL ESFERICA PARA ARGUMENTO
C     IMAGINARIO, ISTO E'
C  I(L,x) = i**L * j(L,-ix)
C  O programa calcula
c  YJ(N,L,M)=J(N-1,L-1,M-1,A,B)
C   DA REGRA DE PARIDADE, VE-SE QUE
C       N + L + M = IMPAR
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YJ(10,6,6),ANT(2),POS(6),D(6),F(7,5),FIN(5),
     * KEY(10,6,6)
      IF(AA.LE.1.D-20)AA=0.D0
      IF(BB.LE.1.D-20)BB=0.D0
      PI=4.D0*DATAN(1.D0)
      DO 2 N=1,7
      DO 2 L=1,5
      DO 2 M=1,5
      KEY(N,L,M)=0
    2 YJ(N,L,M)=0.D0
      IF(AA.EQ.0.D0.AND.BB.EQ.0.D0)THEN
        FATOR=DSQRT(PI)/4.D0
        DO 3 N=1,7,2
        YJ(N,1,1)=FATOR
    3   FATOR=FATOR/2.D0*(N+2)
        RETURN
      ENDIF
      A=AA
      B=BB
      IF(B.GT.A)THEN
        B=AA
        A=BB
      ENDIF
c  Inicializa os fatores F(M,N,m,B)
      FIN(1)=DEXP(-(2.D0*A*B+B**2)/4.D0)
      DO 4 M=2,5
    4 FIN(M)=FIN(M-1)*B/(2*M-1)
      MIL=1000
      IF(B.EQ.0.D0)MIL=10
      DO 301 MG=1,MIL
c  Calcula os fatores D(M,l,A)
c  D(L)=D(MG-1,L-2,A)
c  POS(L)=D(MG,L-2,A)
c  ANT(L)=D(MG-2,L-2,A) somente para L=1 e L=2
      IF(MG.EQ.1)THEN
!fk     ANT(2)=DSQRT(PI)/2.D0/A*GERF(A/2.D0)
        ANT(2)=DSQRT(PI)/2.D0/A*erf(A/2.D0)
        D(2)=DSQRT(PI)/4.D0
        ANT(1)=2.D0/A*D(2)
        D(1)=A/2.D0*ANT(2)
        POS(1)=0.5D0*ANT(1)+A/2.D0*D(2)
        POS(2)=0.5D0*ANT(2)+A/2.D0*D(1)
        DO 302 L=3,6
          POS(L)=POS(L-2)-(2*L-5)/A*D(L-1)
          IF(L.EQ.4)THEN
            D(4)=(1.D0-6.D0/A**2)*D(2)
            GO TO 302
          ENDIF
          D(L)=(L-1)*D(L-2)/(4-L)-2*(2*L-5)/A/(4-L)*POS(L-1)
  302   CONTINUE
      ELSE
        POS(1)=(MG-1)/2.D0*ANT(1)+A/2.D0*D(2)
        POS(2)=(MG-1)/2.D0*ANT(2)+A/2.D0*D(1)
        DO 303 L=3,6
  303   POS(L)=POS(L-2)-(2*L-5)/A*D(L-1)
        ANT(1)=D(1)
        ANT(2)=D(2)
        DO 304 L=1,6
  304   D(L)=POS(L)
      ENDIF
c   Calcula os fatores F(M,N,m,B)
      DO 305 N=1,7
      DO 305 M=1,5
      MNM=MG-N-M+1
      IF(MNM.LT.0)GO TO 305
      IF(MNM.NE.2*(MNM/2))GO TO 305
      IF(MNM.EQ.0)THEN
        F(N,M)=FIN(M)
      ELSE
        F(N,M)=F(N,M)*B**2/MNM/(MNM+2*M-1)
      ENDIF      
      DO 306 L=1,5
      TERMO=F(N,M)*D(L+1)
      IF(BB.LE.AA)THEN
        YJ(N,L,M)=YJ(N,L,M)+TERMO
        IF(TERMO.LT.1.D-10*YJ(N,L,M))KEY(N,L,M)=1
      ELSE
        YJ(N,M,L)=YJ(N,M,L)+TERMO
        IF(TERMO.LT.1.D-10*YJ(N,M,L))KEY(N,M,L)=1
      ENDIF
  306 CONTINUE
  305 CONTINUE
      DO 307 N=1,7
      DO 307 M=1,5
      DO 307 L=1,5
      IF(KEY(N,L,M).EQ.0)GO TO 301
  307 CONTINUE
      RETURN
  301 CONTINUE
ccc   END
      end subroutine jota3
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE JOTA2(AA,BB,YJ)
C   ACHA A TABELA DE FUNCOES J(N,L,L',A,B).
C     ONDE
C  J(N,L,L',A,B) = EXP(-(A+B)**2/4)*INTEGRAL(0,INFINITO)(r**(N+2)*
C       EXP(-r**2)*I(L,A*r)*I(L',B*r)*dr)
C     ONDE I(L,x) E' A FUNCAO DE BESSEL ESFERICA PARA ARGUMENTO
C     IMAGINARIO, ISTO E'
C  I(L,x) = i**L * j(L,-ix)
C  O programa calcula
c  YJ(N,L,M)=J(N-1,L-1,M-1,A,B)
C   DA REGRA DE PARIDADE, VE-SE QUE
C       N + L + M = IMPAR
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YJ(10,6,6)
      PI=4.D0*DATAN(1.D0)
      DO 302 N=1,7
      DO 302 L=1,5
      DO 302 M=1,5
  302 YJ(N,L,M)=0.D0
      IF(AA.EQ.0.D0.AND.BB.EQ.0.D0)THEN
        FATOR=DSQRT(PI)/4.D0
        DO 303 N=1,7,2
        YJ(N,1,1)=FATOR
  303   FATOR=FATOR/2.D0*(N+2)
        RETURN
      ENDIF
      A=AA
      B=BB
      IF(B.GT.A)THEN
        B=AA
        A=BB
      ENDIF
      DO 301 N=1,7
      DO 301 L=1,5
      DO 301 M=1,5
      IF(N+L+M.EQ.2*((N+L+M)/2))GO TO 301
      NLM=N+L+M-3
      BM=1.D0
      SOMA=0.D0
      IF(B.LT.1.D-10.AND.M.GT.1)GO TO 6
      IF(M.GT.1)BM=B**(M-1)
      FATOR=DSQRT(PI)/4.D0/2.D0**(NLM/2)*A**(L-1)*BM*
     * DEXP(-0.25D0*(A+B)**2)
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
      DO 4 I=1,1000
      FATOR=FATOR*(NLM+2*I-1)/(2*L+2*I-3)*A**2/4.D0
      SUM=1.D0
      IF(I.GT.1)THEN
        FATOR=FATOR/(I-1)
        FAT=1.D0
        IF(B.GT.1.D-10)THEN
          DO 5 J=2,I
          FAT=FAT*(B/A)**2*(I-J+1)*(2*L+2*I-2*J+1)/(J-1)/(2*M+2*J-3)
    5     SUM=SUM+FAT
        ENDIF
      ENDIF
      TERMO=FATOR*SUM
      SOMA=SOMA+TERMO
      IF(TERMO.LT.1.D-10*SOMA)GO TO 6
    4 CONTINUE
    6 IF(BB.LE.AA)YJ(N,L,M)=SOMA
      IF(BB.GT.AA)YJ(N,M,L)=SOMA
  301 CONTINUE
      RETURN
ccc   END
      end subroutine jota2
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE JOTA(A,B,YJ)
C   ACHA A TABELA DE FUNCOES J(N,L,L',A,B).
C     ONDE
C  J(N,L,L',A,B) = EXP(-(A+B)**2/4)*INTEGRAL(0,INFINITO)(r**(N+2)*
C       EXP(-r**2)*I(L,A*r)*I(L',B*r)*dr)
C     ONDE I(L,x) E' A FUNCAO DE BESSEL ESFERICA PARA ARGUMENTO
C     IMAGINARIO, ISTO E'
C  I(L,x) = i**L * j(L,-ix)
C   NO PROGRAMA TEMOS 
C     XJ(N,L,M) = J(N-3,L-2,M-2)
C   DA REGRA DE PARIDADE, VE-SE QUE
C       N + L + M = IMPAR
C   NO FINAL FAZEMOS
C      YJ(N,L,M) = J(N-1,L-1,M-1)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/IJOTA/PI,AP(9,9),AR(9,9),AKEY
      DIMENSION XJ(10,6,6),PSUM(9),RSUM(9),PDIF(9),RDIF(9),
     1 YJ(10,6,6)
C   INICIALIZA OS POLINOMIOS
c Atencao, a zeragem de XJ nao estava incluida no programa original
      DO 827 I=1,10
      DO 827 J=1,6
      DO 827 K=1,6
  827 XJ(I,J,K)=0.D0
      IF(AKEY.EQ.0.123456789D0)GO TO 200
      PI=4.D0*DATAN(1.D0)
      DO 11 I=1,9
      AP(I,1)=0.D0
   11 AR(I,1)=0.D0
      AP(1,1)=1.D0
      DO 12 L=2,9
      DO 13 I=1,9
      AP(I,L)=0.D0
   13 AR(I,L)=0.5D0*AP(I,L-1)
      DO 14 I=2,9
      AP(I-1,L)=(I-1)*AP(I,L-1)
   14 AR(I-1,L)=AR(I-1,L)+(I-1)*AR(I,L-1)
      DO 15 I=1,8
   15 AP(I+1,L)=AP(I+1,L)+0.5D0*AP(I,L-1)
   12 CONTINUE
      AKEY=0.123456789D0
  200 CONTINUE
C   SE ALGUM ARGUMENTO A OU B E' NULO, ADOTA OUTRO PROCEDIMENTO DE
C     CALCULO
      IF(A.LT.1.0D0.OR.B.LT.1.0D0)THEN
        IF(A.LT.1.0D0.AND.B.LT.1.0D0)THEN
          CALL JOTA2(A,B,YJ)
        ELSE
          CALL JOTA3(A,B,YJ)
        ENDIF
        RETURN
      ENDIF
C   ACHA J(N,0,0), J(N,-1,0), E J(N,0,-1) PARA N.GE.0
      CALL POLINO(A+B,PSUM,RSUM)
      CALL POLINO(A-B,PDIF,RDIF)
      SQPI=DSQRT(PI)/8.D0/A/B
      EP=0.D0
      IF(A*B.LT.40.D0)EP=DEXP(-A*B)
!fk   ERFSUM=GERF(0.5D0*(A+B))
      ERFSUM=erf(0.5D0*(A+B))
!fk   ERFDIF=GERF(0.5D0*(A-B))
      ERFDIF=erf(0.5D0*(A-B))
      EPSUM=0.D0
      IF((A+B)**2.LT.160.D0)EPSUM=DEXP(-0.25D0*(A+B)**2)/4.D0/A/B
      ISIG=1
      DO 1 I=1,7,2
      AUX=PSUM(I)+ISIG*PSUM(I)-EP*(PDIF(I)+ISIG*PDIF(I))
     1 +ERFSUM*(PSUM(I)-ISIG*PSUM(I))-EP*ERFDIF*(PDIF(I)-ISIG*PDIF(I))
      XJ(I+2,2,2)=SQPI*AUX+EPSUM*(RSUM(I)-ISIG*RSUM(I)-RDIF(I)+
     1 ISIG*RDIF(I))
    1 CONTINUE
      ISIG=-1
      DO 501 I=2,8,2
      AUX=SQPI*(PSUM(I)-ISIG*PSUM(I))
      BUX=SQPI*ERFSUM*(PSUM(I)+ISIG*PSUM(I))
      CUX=EPSUM*(RSUM(I)+ISIG*RSUM(I))
      DUX=SQPI*EP*(-PDIF(I)+ISIG*PDIF(I))
      EUX=SQPI*EP*ERFDIF*(-PDIF(I)-ISIG*PDIF(I))
      FUX=EPSUM*(-RDIF(I)-ISIG*RDIF(I))
      XJ(I+2,1,2)=AUX+BUX+CUX+DUX+EUX+FUX
  501 XJ(I+2,2,1)=AUX+BUX+CUX+ISIG*DUX-ISIG*EUX-ISIG*FUX
c***
ccc   AUX=EME(A+B)
ccc   BUX=EME(A-B)
      AUX=emme(A+B)
      BUX=emme(A-B)
c***
      XJ(2,1,2)=2.D0*SQPI*((A+B)*AUX-(A-B)*EP*BUX)
      XJ(2,2,1)=2.D0*SQPI*((A+B)*AUX+(A-B)*EP*BUX)
C   ACHA J(-2,0,0)
      XJ(1,2,2)=A*XJ(2,1,2)+B*XJ(2,2,1)-2.D0*XJ(3,2,2)
C   ACHA J(N,1,-1) E J(N,-1,1) PARA N.GE.0
      DO 2 I=3,9,2
      XJ(I,3,1)=2.D0/B*XJ(I+1,1,2)-(A/B+2.D0/A/B)*XJ(I,2,2)-
     1 (I-4)/B*XJ(I-1,1,2)+(I-4)/A/B*XJ(I-2,2,2)
    2 XJ(I,1,3)=2.D0/A*XJ(I+1,2,1)-(B/A+2.D0/A/B)*XJ(I,2,2)-
     1 (I-4)/A*XJ(I-1,2,1)+(I-4)/A/B*XJ(I-2,2,2)
C   USA AS REGRAS DE RECORRENCIA PARA COMPLETAR A TABELA
      DO 3 L=1,4
      DO 4 M=L,4
      DO 4 N=M+1,9
      IF(N+M+L.NE.2*((N+M+L)/2))GO TO 4
      XJ(N,M+2,L+1)=XJ(N,M,L+1)-(2*M-1)/A*XJ(N-1,M+1,L+1)
      XJ(N,L+1,M+2)=XJ(N,L+1,M)-(2*M-1)/B*XJ(N-1,L+1,M+1)
    4 CONTINUE
      DO 5 N=L+2,9
      IF(N.EQ.2*(N/2))GO TO 5
      XJ(N,L+2,L+2)=XJ(N,L,L+2)-(2*L-1)/A*XJ(N-1,L+1,L+2)
    5 CONTINUE
      IF(L.NE.2*(L/2))GO TO 3
      AUX=2.D0*XJ(L+3,L+2,L+2)-A*XJ(L+2,L+1,L+2)-B*XJ(L+2,L+2,L+1)
      XJ(L+1,L+2,L+2)=-AUX/(L+1)
    3 CONTINUE
      DO 301 N=1,7
      DO 301 L=1,5
      DO 301 M=1,5
      IF(N+L+M.EQ.2*((N+L+M)/2))GO TO 301
      YJ(N,L,M)=XJ(N+2,L+1,M+1)
  301 CONTINUE
      RETURN
ccc   END
      end subroutine jota
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE POLINO(X,P,R)
      IMPLICIT REAL*8 (A-H,O-Z)
C   ACHA OS POLINOMIOS P(N,X) E R(N,X)
      COMMON/IJOTA/PI,AP(9,9),AR(9,9),AKEY
      DIMENSION P(9),R(9)
      P(1)=1.D0
      R(1)=0.D0
      DO 1 I=2,9
      P(I)=AP(1,I)
      R(I)=AR(1,I)
      AUX=X
      DO 1 J=2,I
      P(I)=P(I)+AP(J,I)*AUX
      R(I)=R(I)+AR(J,I)*AUX
    1 AUX=AUX*X
      RETURN
ccc   END
      end subroutine polino
c----------------------------------------------------------------

c----------------------------------------------------------------
ccc   REAL*8 FUNCTION EME(X)
      REAL*8 FUNCTION emme(X)
      IMPLICIT REAL*8 (A-H,O-Z)
C  EME(X) = 2/X/SQRT(PI)*EXP(-X**2/4)*INTEGRAL(0,INFINITO)(EXP(-R**2)/R*
C           SINH(X*R))
ccc   EME=1.D0
      emme=1.D0
      IF(X.EQ.0.D0)RETURN
ccc   EME=0.D0
      emme=0.D0
      TERMO=-X**2/4.D0
ccc   IF(TERMO.GT.-40.D0)EME=EME+DEXP(TERMO)
      IF(TERMO.GT.-40.D0)emme=emme+DEXP(TERMO)
      DO 1 N=1,100000
      RAZAO=DLOG((2*N-1)*X**2/N/4/(2*N+1))
      TERMO=TERMO+RAZAO
      IF(TERMO.LT.-40.D0)GO TO 2
ccc   EME=EME+DEXP(TERMO)
      emme=emme+DEXP(TERMO)
      GO TO 1
    2 IF(RAZAO.LT.0.D0)GO TO 3
    1 CONTINUE
    3 RETURN
ccc   END
      end function emme
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE KTIL2(LSUM,EXPMIN,X,Y,SOMAK)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SOMAK(9,9),Q(1000),C(1001)
      RAZAO=1.D0+X
      NMAX=1
      Q(1)=1.D0/DSQRT(RAZAO)
      C(1)=1.D0
      DO 10 L=1,LSUM+1
      DO 10 N=1,LSUM+1
   10 SOMAK(N,L)=0.D0
      IF(X.EQ.0.D0)RETURN
      DO 1 L=1,LSUM+1
      YL=1.D0
      IF(L.GT.1.AND.Y.LT.1.D-10)RETURN
      IF(Y.GT.1.D-10)YL=Y**(L-1)
C  (2L+1)!!
        LDFAT=1
        DO 2 I=1,L
    2   LDFAT=LDFAT*(2*I-1)
      DO 1 N=L,LSUM+1
C  n!
        NFAT=1
        IF(N.GT.1)THEN
          DO 3 I=1,N-1
    3     NFAT=NFAT*I
        ENDIF
      FATOR=0.5D0*DSQRT(X)*DEXP(-0.25D0*Y**2)*YL*NFAT/LDFAT
      SOMA=0.D0
      MIL=1000
      IF(Y.EQ.0.D0)MIL=10
      DO 4 M=1,MIL
      IF(M.GT.1)THEN
        FATOR=FATOR*Y**2/2.D0/(M-1)/(2*L+2*M-3)*(N+M-2)
      ENDIF
c  Caso em que tem que recalcular o polinomio Q
      IF(N+M-1.GT.NMAX)THEN
        DO 5 JJ=1,NMAX+1
        J=NMAX+2-JJ
        AUX=0.D0
        IF(J.LE.NMAX)AUX=AUX+(NMAX+1-J)*C(J)/NMAX
        IF(J.GT.1)AUX=AUX+0.5D0*(2*J-3)*C(J-1)/NMAX
    5   C(J)=AUX
        NMAX=NMAX+1
        FAT=DSQRT(RAZAO)
        SUM=0.D0
        DO 6 K=1,NMAX
        FAT=FAT/RAZAO
    6   SUM=SUM+FAT*C(K)
        Q(NMAX)=SUM
      ENDIF
      TERMO=FATOR*Q(N+M-1)
      SOMA=SOMA+TERMO
      IF(TERMO.LT.1.D-10*SOMA)GO TO 7
    4 CONTINUE
    7 SOMAK(N,L)=SOMA
    1 CONTINUE
      RETURN
ccc   END
      end subroutine ktil2
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE KTIL(LSUM,EXPMIN,X,Y,SOMAK)
C   ACHA A TABELA DE FUNCOES KTIL PARA VALORES N E LAMBDA TAIS QUE
C              LSUM >= N >= LAMBDA
C              LSUM >= LAMBDA >= 0
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/KFUNC/H(16,16),DH(16,16),COMECO,SQPI
      DIMENSION SOMAK(9,9),HFAT(16),Z(17),GFAT(16),YJ(8,8)
      IF(LSUM.GT.8)STOP 'L(1)+L(2) MAIOR QUE 8'
      DO 1 L=1,LSUM+1
      DO 1 N=L,LSUM+1
   1  SOMAK(N,L)=0.D0
      IF(DABS(Y).LT.1.D0)THEN
        CALL KTIL2(LSUM,EXPMIN,X,Y,SOMAK)
        RETURN
      ENDIF
      IF(COMECO.EQ.0.123456789D0)GO TO 200
      SQPI=DSQRT(4.D0*DATAN(1.D0))
      DO 3 I=1,16
      DO 3 N=1,16
      H(I,N)=0.D0
    3 DH(I,N)=0.D0
      H(1,1)=1.D0
      DO 4 N=1,15
      DO 5 I=1,N
    5 DH(I,N)=I*H(I+1,N)
      DO 6 I=1,N
    6 H(I,N+1)=DH(I,N)
      DO 4 I=2,N+1
    4 H(I,N+1)=H(I,N+1)+2.D0*H(I-1,N)
      COMECO=0.123456789D0
  200 DO 7 N=1,2*LSUM
    7 HFAT(N)=0.D0
      SQX=DSQRT(1.D0+X)
      AUX=X*Y**2/4.D0/(1.D0+X)
      IF(AUX.GT.EXPMIN)GO TO 300
      FAT=Y/2.D0/SQX
      ARG=DEXP(-AUX)
      DO 8 I=1,2*LSUM
      DO 9 N=I,2*LSUM
    9 HFAT(N)=HFAT(N)+H(I,N)*ARG
    8 ARG=ARG*FAT
      FAT=DSQRT(X)/2.D0/SQX
      DO 10 N=1,2*LSUM
      GFAT(N)=HFAT(N)*FAT
   10 FAT=FAT/2.D0/SQX
!fk 300 Z(1)=SQPI/2.D0*GERF(Y*DSQRT(X)/2.D0/SQX)
  300 Z(1)=SQPI/2.D0*erf(Y*DSQRT(X)/2.D0/SQX)
      IF(LSUM.EQ.0)GO TO 400
      Z(2)=Y*Z(1)/2.D0+GFAT(1)
      DO 11 N=3,2*LSUM+1
   11 Z(N)=(N-2)/2.D0*Z(N-2)+Y/2.D0*Z(N-1)+GFAT(N-1)
  400 DO 12 N=1,LSUM+1
   12 SOMAK(N,1)=Z(2*N-1)/Y
      IF(LSUM.EQ.0)RETURN
      FAT=2.D0*SQPI/Y
      DO 14 N=1,LSUM
      FAT=FAT/4.D0/(1.D0+X)
   14 YJ(N,1)=HFAT(2*N)*FAT
      IF(LSUM.EQ.1)GO TO 17
      DO 15 L=2,LSUM
      DO 16 N=L,LSUM
   16 YJ(N,L)=2.D0*(1.D0+X)/Y*YJ(N,L-1)-(2*N-1)/Y*YJ(N-1,L-1)
   15 CONTINUE
   17 CONTINUE
      DO 13 L=2,LSUM+1
      DO 13 N=L,LSUM+1
   13 SOMAK(N,L)=2.D0*(SOMAK(N,L-1)-(N-1)*SOMAK(N-1,L-1))/Y-
     * DSQRT(X)/SQPI/Y*YJ(N-1,L-1)
      RETURN
ccc   END
      end subroutine ktil
c----------------------------------------------------------------

c----------------------------------------------------------------
      real*8 function gerf (x)
c
c  error function erf(x)
c july 1977 edition.  w. fullerton, c3, los alamos sci. lab.
c modified by j.l.martins dec 84
c 16 significant figures
c
      implicit real*8 (a-h,o-z)
      dimension erfcs(13), erfccs(24), erc2cs(23)
      data spi/1.772453850905516d0/
      data erfcs( 1) /   -.0490461212 34691808d0 /
      data erfcs( 2) /   -.1422612051 03713640d0 /
      data erfcs( 3) /    .0100355821 87599796d0 /
      data erfcs( 4) /   -.0005768764 69976748d0 /
      data erfcs( 5) /    .0000274199 31252196d0 /
      data erfcs( 6) /   -.0000011043 17550734d0 /
      data erfcs( 7) /    .0000000384 88755420d0 /
      data erfcs( 8) /   -.0000000011 80858253d0 /
      data erfcs( 9) /    .0000000000 32334215d0 /
      data erfcs(10) /   -.0000000000 00799101d0 /
      data erfcs(11) /    .0000000000 00017990d0 /
      data erfcs(12) /   -.0000000000 00000371d0 /
      data erfcs(13) /    .0000000000 00000007d0 /
      data erc2cs( 1) /   -.0696013466 02309501d0 /
      data erc2cs( 2) /   -.0411013393 62620893d0 /
      data erc2cs( 3) /    .0039144958 66689626d0 /
      data erc2cs( 4) /   -.0004906395 65054897d0 /
      data erc2cs( 5) /    .0000715747 90013770d0 /
      data erc2cs( 6) /   -.0000115307 16341312d0 /
      data erc2cs( 7) /    .0000019946 70590201d0 /
      data erc2cs( 8) /   -.0000003642 66647159d0 /
      data erc2cs( 9) /    .0000000694 43726100d0 /
      data erc2cs(10) /   -.0000000137 12209021d0 /
      data erc2cs(11) /    .0000000027 88389661d0 /
      data erc2cs(12) /   -.0000000005 81416472d0 /
      data erc2cs(13) /    .0000000001 23892049d0 /
      data erc2cs(14) /   -.0000000000 26906391d0 /
      data erc2cs(15) /    .0000000000 05942614d0 /
      data erc2cs(16) /   -.0000000000 01332386d0 /
      data erc2cs(17) /    .0000000000 00302804d0 /
      data erc2cs(18) /   -.0000000000 00069666d0 /
      data erc2cs(19) /    .0000000000 00016208d0 /
      data erc2cs(20) /   -.0000000000 00003809d0 /
      data erc2cs(21) /    .0000000000 00000904d0 /
      data erc2cs(22) /   -.0000000000 00000216d0 /
      data erc2cs(23) /    .0000000000 00000052d0 /
      data erfccs( 1) /   0.0715179310 20292500d0 /
      data erfccs( 2) /   -.0265324343 37606719d0 /
      data erfccs( 3) /    .0017111539 77920853d0 /
      data erfccs( 4) /   -.0001637516 63458512d0 /
      data erfccs( 5) /    .0000198712 93500549d0 /
      data erfccs( 6) /   -.0000028437 12412769d0 /
      data erfccs( 7) /    .0000004606 16130901d0 /
      data erfccs( 8) /   -.0000000822 77530261d0 /
      data erfccs( 9) /    .0000000159 21418724d0 /
      data erfccs(10) /   -.0000000032 95071356d0 /
      data erfccs(11) /    .0000000007 22343973d0 /
      data erfccs(12) /   -.0000000001 66485584d0 /
      data erfccs(13) /    .0000000000 40103931d0 /
      data erfccs(14) /   -.0000000000 10048164d0 /
      data erfccs(15) /    .0000000000 02608272d0 /
      data erfccs(16) /   -.0000000000 00699105d0 /
      data erfccs(17) /    .0000000000 00192946d0 /
      data erfccs(18) /   -.0000000000 00054704d0 /
      data erfccs(19) /    .0000000000 00015901d0 /
      data erfccs(20) /   -.0000000000 00004729d0 /
      data erfccs(21) /    .0000000000 00001432d0 /
      data erfccs(22) /   -.0000000000 00000439d0 /
      data erfccs(23) /    .0000000000 00000138d0 /
      data erfccs(24) /   -.0000000000 00000048d0 /
c     smallest positive number on computer
c     small=6.8d-8618
c     relative precision
c     rel=3.6d-15
c     xsml=-sqrt(-log(spi*rel))
c     xmax=sqrt(-log(spi*small))
c     xmax=xmax-0.5d0*log(xmax)/xmax-0.01d0
c     sqeps=sqrt(2.0d0*rel)
      xsml=-5.7d0
      xmax=140.0d0
      sqeps=8.5d-8
      if(x.lt.xsml) then
        gerf=-1.0d0
        return
      endif
      if(x.gt.xmax) then
        gerf=1.0d0
        return
      endif
      y=abs(x)
      if(y.lt.sqeps) then
        gerf=2.0d0*x/spi
        return
      endif
      if(y.lt.1.0d0) then
        z=2.0d0*x*x-1.0d0
        twoz=2.0d0*z
        b1=0.0d0
        b0=0.0d0
        do 10 i=1,13
          b2=b1
          b1=b0
          b0=twoz*b1-b2+erfcs(14-i)
10      continue
        gerf=x*(1.0d0+0.5d0*(b0-b2))
        return
      endif
      yt=y*y
      if(yt.le.4.0d0) then
        z=(8.0d0/yt-5.0d0)/3.0d0
        b1=0.0d0
        b0=0.0d0
        twoz=2.0d0*z
        do 20 i=1,23
          b2=b1
          b1=b0
          b0=twoz*b1-b2+erc2cs(24-i)
20      continue
        erfc=(exp(-yt)/y)*0.5d0*(1.0d0+b0-b2)
      else
        z=8.0d0/yt-1.0d0
        b1=0.0d0
        b0=0.0d0
        twoz=2.0d0*z
        do 30 i=1,24
          b2=b1
          b1=b0
          b0=twoz*b1-b2+erfccs(25-i)
30      continue
        erfc=(exp(-yt)/y)*0.5d0*(1.0d0+b0-b2)
      endif
      if(x.lt.0.0d0) erfc=2.0d0-erfc
      gerf=1.d0-erfc
      return
ccc   end
      end function gerf
c----------------------------------------------------------------

      end module pseudo_pp
