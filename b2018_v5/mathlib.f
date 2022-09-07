      module mathlib

      contains

c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION CDERF(ZP)
CRVX$ END
CRVX$ CRAY
C     COMPLEX FUNCTION CDERF(ZP)
CRVX$ END
C=======================================================================
C
C      ERROR FUNCTION OF COMPLEX ARGUMENT.
C
C=======================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAMBDA
      COMPLEX*16 DCERFC,ZP,AI,Z,ZA,CERF
      LOGICAL B
      DATA AI/(0.D0,1.D0)/
      Z=AI*ZP
CRVX$ VAX
      X=DREAL(Z)
      Y=DIMAG(Z)
CRVX$ END
CRVX$ CRAY
C     X=REAL(Z)
C     Y=AIMAG(Z)
CRVX$ END
      IQ=1
      IF (X .LT. 0.D0) IQ=2
      IF (Y .LT. 0.D0) IQ=3
      IF (Y .LT. 0.D0.AND. X  .GE. 0.D0) IQ=4
      X = DABS(X)
      Y = DABS(Y)
CRVX$ VAX
      ZA=DCMPLX(X,Y)
CRVX$ END
CRVX$ CRAY
C     ZA=CMPLX(X,Y)
CRVX$ END
      IF (Y .LT.4.29D0.AND.X .LT. 5.33D0) GO TO 10
      GO TO 15
   10 S=(1.D0-Y/4.29D0)*DSQRT(1.D0-X*X/28.41D0)
      H=1.6D0*S
      H2=2.D0*H
      K= 6+ 23*S
      LAMBDA=H2**K
      NU=9+ 21*S
      GO TO 20
   15 H=0.D0
      K=0
      NU=8
C=======================================================================
C
C  20 IF(H .GT. 0.D0) LAMBDA=H2**K            THIS CARD IS DELETED
C
C=======================================================================
   20 CONTINUE
      B=.FALSE.
      IF(H.EQ.0.D0 .OR. LAMBDA .EQ.0.D0) B=.TRUE.
      R1=0.D0
      R2=0.D0
      S1=0.D0
      S2=0.D0
      N=NU+1
   30 N=N-1
      IF(N .LT. 0) GO TO 50
      NP1=N+1
      T1=Y+H+NP1*R1
      T2= X-NP1*R2
      C=.5D0/(T1*T1+T2*T2)
      R1=C*T1
      R2=C*T2
      IF(H.LE.0.D0.OR.  N .GT. K) GO TO 30
      T1=LAMBDA+S1
      S1=R1* T1-R2*S2
      S2=R2*T1+R1*S2
      LAMBDA=LAMBDA/H2
      GO TO 30
   50 CONTINUE
      IF (Y .NE.0.D0) GO TO 60
      RE=DEXP(-X*X)/1.12837916709551D0
      GO TO 70
   60 RE=S1
      IF (B) RE=R1
   70 AIM=S2
      IF (B) AIM=R2
CRVX$ VAX
      CERF=1.12837916709551D0*DCMPLX(RE,AIM)
CRVX$ END
CRVX$ CRAY
C     CERF=1.12837916709551D0*CMPLX(RE,AIM)
CRVX$ END
      IF (IQ .EQ. 1) GO TO 90
      IF (IQ .NE. 2) GO TO 80
CRVX$ VAX
      CERF=DCONJG(CERF)
      GO TO 90
   80 CERF=2.D0 *CDEXP(-ZA*ZA)-CERF
CRVX$ END
CRVX$ CRAY
C     CERF=CONJG(CERF)
C     GO TO 90
C  80 CERF=2.D0 *CEXP(-ZA*ZA)-CERF
CRVX$ END
      IF (IQ .EQ. 3) GO TO 90
CRVX$ VAX
      CERF=DCONJG(CERF)
 90   DCERFC = CERF*CDEXP(Z*Z)
CRVX$ END
CRVX$ CRAY
C     CERF=CONJG(CERF)
C90   DCERFC = CERF*CEXP(Z*Z)
CRVX$ END
      CDERF=1.0D0-DCERFC
ccc   RETURN
ccc   END
      end function cderf
c----------------------------------------------------------------

c----------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION CDERFM(ZP)
CRVX$ END
CRVX$ CRAY
C     COMPLEX FUNCTION CDERFM(ZP)
CRVX$ END
C=======================================================================
C
C      ERROR FUNCTION OF COMPLEX ARGUMENT.  MODIFIED TO RETURN ONLY THE
C      VALUE OF CERF(ZP) FOR USE IN SUBROUTINE INTG.  TLG:9/5/84.
C
C=======================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAMBDA
      COMPLEX*16 DCERFC,ZP,AI,Z,ZA,CERF
      LOGICAL B
      DATA AI/(0.D0,1.D0)/
      Z=AI*ZP
CRVX$ VAX
      X=DREAL(Z)
      Y=DIMAG(Z)
CRVX$ END
CRVX$ CRAY
C     X=REAL(Z)
C     Y=AIMAG(Z)
CRVX$ END
      IQ=1
      IF (X .LT. 0.D0) IQ=2
      IF (Y .LT. 0.D0) IQ=3
      IF (Y .LT. 0.D0.AND. X  .GE. 0.D0) IQ=4
      X = DABS(X)
      Y = DABS(Y)
CRVX$ VAX
      ZA=DCMPLX(X,Y)
CRVX$ END
CRVX$ CRAY
C     ZA=CMPLX(X,Y)
CRVX$ END
      IF (Y .LT.4.29D0.AND.X .LT. 5.33D0) GO TO 10
      GO TO 15
   10 S=(1.D0-Y/4.29D0)*DSQRT(1.D0-X*X/28.41D0)
      H=1.6D0*S
      H2=2.D0*H
      K= 6+ 23*S
      LAMBDA=H2**K
      NU=9+ 21*S
      GO TO 20
   15 H=0.D0
      K=0
      NU=8
C=======================================================================
C
C  20 IF(H .GT. 0.D0) LAMBDA=H2**K            THIS CARD IS DELETED
C
C=======================================================================
   20 CONTINUE
      B=.FALSE.
      IF(H.EQ.0.D0 .OR. LAMBDA .EQ.0.D0) B=.TRUE.
      R1=0.D0
      R2=0.D0
      S1=0.D0
      S2=0.D0
      N=NU+1
   30 N=N-1
      IF(N .LT. 0) GO TO 50
      NP1=N+1
      T1=Y+H+NP1*R1
      T2= X-NP1*R2
      C=.5D0/(T1*T1+T2*T2)
      R1=C*T1
      R2=C*T2
      IF(H.LE.0.D0.OR.  N .GT. K) GO TO 30
      T1=LAMBDA+S1
      S1=R1* T1-R2*S2
      S2=R2*T1+R1*S2
      LAMBDA=LAMBDA/H2
      GO TO 30
   50 CONTINUE
      IF (Y .NE.0.D0) GO TO 60
      RE=DEXP(-X*X)/1.12837916709551D0
      GO TO 70
   60 RE=S1
      IF (B) RE=R1
   70 AIM=S2
      IF (B) AIM=R2
CRVX$ VAX
      CERF=1.12837916709551D0*DCMPLX(RE,AIM)
CRVX$ END
CRVX$ CRAY
C     CERF=1.12837916709551D0*CMPLX(RE,AIM)
CRVX$ END
      IF (IQ .EQ. 1) GO TO 90
      IF (IQ .NE. 2) GO TO 80
CRVX$ VAX
      CERF=DCONJG(CERF)
      GO TO 90
   80 CERF=2.D0 *CDEXP(-ZA*ZA)-CERF
CRVX$ END
CRVX$ CRAY
C     CERF=CONJG(CERF)
C     GO TO 90
C  80 CERF=2.D0 *CEXP(-ZA*ZA)-CERF
CRVX$ END
      IF (IQ .EQ. 3) GO TO 90
CRVX$ VAX
      CERF=DCONJG(CERF)
 90   DCERFC = CERF
CRVX$ END
CRVX$ CRAY
C     CERF=CONJG(CERF)
C90   DCERFC = CERF
CRVX$ END
      CDERFM=DCERFC
ccc   RETURN
ccc   END
      end function cderfm
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE GAUSS(NPQ,XX,WW,NQP)

      implicit none

      integer :: NPQ, NQP
      real*8  :: XX(NQP), WW(NQP)

      integer :: i, ii, NPQH
      real*8  :: X2(1),W2(1),X4(2),W4(2),
     1  X6(3),W6(3),X8(4),W8(4),X10(5),W10(5),
     1  X12(6),W12(6),X14(7),W14(7),X16(8),W16(8),X18(9),
     1   W18(9),X20(10),W20(10),X22(11),W22(11),X24(12),W24(12),X26(13),
     2      W26(13),X32(16),W32(16),
     3      X48(24),W48(24)
      DATA X2/0.57735026918962576450D0/
      DATA W2/1.0D0/
      DATA X4/0.86113631159405257522D0,0.33998104358485626480D0/
      DATA W4/0.34785484513745385737D0,0.65214515486254614262D0/
      DATA X6/
     1      0.93246951420315202781D00,0.66120938646626451366D00,
     2      0.23861918608319690863D00/
      DATA W6/
     1      0.17132449237917034504D00,0.36076157304813860756D00,
     2      0.46791393457269104738D00/
      DATA X8/
     1      0.96028985649753623168D00,0.79666647741362673959D00,
     2      0.52553240991632898581D00,0.18343464249564980493D00/
        DATA W8/
     1      0.10122853629037625915D00,0.22238103445337447054D00,
     2      0.31370664587788728733D00,0.36268378337836198296D00/
      DATA X10/
     1      0.97390652851717172007D00,0.86506336668898451073D00,
     2      0.67940956829902440623D00,0.43339539412924719079D00,
     3      0.14887433898163121088D00/
         DATA W10/
     1      0.66671344308688137593D-1,0.14945134915058059314D00,
     2      0.21908636251598204399D00,0.26926671930999635509D00,
     3      0.29552422471475287017D00/
      DATA X12/
     1  0.98156063424671925069D00,0.90411725637047485667D00,
     2      0.76990267419430468703D00,0.58731795428661744729D00,
     3      0.36783149899818019375D00,0.12523340851146891547D00/
      DATA W12/
     1      0.47175336386511827194D-1,0.10693932599531843096D00,
     2      0.16007832854334622633D00,0.20316742672306592174D00,
     3      0.23349253653835480876D00,0.24914704581340278500D00/
      DATA X14/
     1      0.98628380869681233884D00,0.92843488366357351733D00,
     2      0.82720131506976499318D00,0.68729290481168547014D00,
     3      0.51524863635815409196D00,0.31911236892788976043D00,
     4      0.10805494870734366206D00/
      DATA W14/
     1      0.35119460331751863031D-1,0.80158087159760209805D-1,
     2      0.12151857068790318468D00,0.15720316715819353456D00,
     3      0.18553839747793781374D00,0.20519846372129560396D00,
     4      0.21526385346315779019D00/
           DATA X16/
     1      0.98940093499164993259D00,0.94457502307323257607D00,
     2      0.86563120238783174388D00,0.75540440835500303389D00,
     3      0.61787624440264374844D00,0.45801677765722738634D00,
     4      0.28160355077925891323D00,0.95012509837637440185D-1/
      DATA W16/
     1      0.27152459411754094851D-1,0.62253523938647892862D-1,
     2      0.95158511682492784809D-1,0.12462897125553387205D00,
     3      0.14959598881657673208D00,0.16915651939500253818D00,
     4      0.18260341504492358886D00,0.18945061045506849628D00/
      DATA X18/
     1      0.99156516842093094673D00,0.95582394957139775518D00,
     2      0.89260246649755573920D00,0.80370495897252311568D00,
     3      0.69168704306035320787D00,0.55977083107394753460D00,
     4      0.41175116146284264603D00,0.25188622569150550958D00,
     5      0.84775013041735301242D-1/
      DATA W18/
     1      0.21616013526483310313D-1,0.49714548894969796453D-1,
     2      0.76425730254889056529D-1,0.10094204410628716556D00,
     3      0.12255520671147846018D00,0.14064291467065065120D00,
     4      0.15468467512626524492D00,0.16427648374583272298D00,
     5      0.16914238296314359184D00/
           DATA X20/
     1      0.99312859918509492478D00,0.96397192727791379126D00,
     2      0.91223442825132590586D00,0.83911697182221882339D00,
     3      0.74633190646015079261D00,0.63605368072651502545D00,
     4      0.51086700195082709800D00,0.37370608871541956067D00,
     5      0.22778585114164507808D00,0.76526521133497333754D-1/
      DATA W20/
     1      0.17614007139152118311D-1,0.40601429800386941331D-1,
     2      0.62672048334109063569D-1,0.83276741576704748724D-1,
     3      0.10193011981724043503D00,0.11819453196151841731D00,
     4      0.13168863844917662689D00,0.14209610931838205132D00,
     5      0.14917298647260374678D00,0.15275338713072585069D00/
      DATA X22/
     1      0.99429458548239929207D00,0.97006049783542872712D00,
     2      0.92695677218717400052D00,0.86581257772030013653D00,
     3      0.78781680597920816200D00,0.69448726318668278005D00,
     4      0.58764040350691159295D00,0.46935583798675702640D00,
     5      0.34193582089208422515D00,0.20786042668822128547D00,
     6      0.69739273319722221213D-1/
      DATA W22/
     1      0.14627995298272200684D-1,0.33774901584814154793D-1,
     2      0.52293335152683285940D-1,0.69796468424520488094D-1,
     3      0.85941606217067727414D-1,0.10041414444288096493D00,
     4      0.11293229608053921839D00,0.12325237681051242428D00,
     5      0.13117350478706237073D00,0.13654149834601517135D00,
     6      0.13925187285563199337D00/
      DATA X24/
     1      0.99518721999702136017D00,0.97472855597130949819D00,
     2      0.93827455200273275852D00,0.88641552700440103421D00,
     3      0.82000198597390292195D00,0.74012419157855436424D00,
     4      0.64809365193697556925D00,0.54542147138883953565D00,
     5      0.43379350762604513848D00,0.31504267969616337438D00,
     6      0.19111886747361630915D00,0.64056892862605626085D-1/
      DATA W24/
     1      0.12341229799987199546D-1,0.28531388628933663181D-1,
     2      0.44277438817419806168D-1,0.59298584915436780746D-1,
     3      0.73346481411080305734D-1,0.86190161531953275917D-1,
     4      0.97618652104113888269D-1,0.10744427011596563478D00,
     5      0.11550566805372560135D00,0.12167047292780339120D00,
     6      0.12583745634682829612D00,0.12793819534675215697D00/
      DATA X26/
     1      0.99588570114561692900D00,0.97838544595647099110D00,
     2      0.94715906666171425013D00,0.90263786198430707421D00,
     3      0.84544594278849801879D00,0.77638594882067885619D00,
     4      0.69642726041995726486D00,0.60669229301761806323D00,
     5      0.50844071482450571769D00,0.40305175512348630648D00,
     6      0.29200483948595689514D00,0.17685882035689018396D00,
     7      0.59230093429313207093D-1/
      DATA W26/
     1      0.10551372617343007155D-1,0.24417851092631908789D-1,
     2      0.37962383294362763950D-1,0.50975825297147811998D-1,
     3      0.63274046329574835539D-1,0.74684149765659745887D-1,
     4      0.85045894313485239210D-1,0.94213800355914148463D-1,
     5      0.10205916109442542323D00,0.10847184052857659065D00,
     6      0.11336181654631966654D00,0.11666044348529658204D00,
     7      0.11832141527926227651D00/
      DATA X32/
     1      0.99726386184948156354D0,0.98561151154526833540D0,
     2      0.96476225558750643077D0,0.93490607593773968917D0,
     3      0.89632115576605212397D0,0.84936761373256997013D0,
     4      0.79448379596794240696D0,0.73218211874028968039D0,
     5      0.66304426693021520098D0,0.58771575724076232904D0,
     6      0.50689990893222939002D0,0.42135127613063534536D0,
     7      0.33186860228212764978D0,0.23928736225213707454D0,
     8      0.14447196158279649349D0,0.04830766568773831623D0/
      DATA W32/
     1      0.00701861000947009660D0,0.01627439473090567061D0,
     2      0.02539206530926205945D0,0.03427386291302143310D0,
     3      0.04283589802222668066D0,0.05099805926237617620D0,
     4      0.05868409347853554714D0,0.06582222277636184684D0,
     5      0.07234579410884850623D0,0.07819389578707030647D0,
     6      0.08331192422694675522D0,0.08765209300440381114D0,
     7      0.09117387869576388471D0,0.09384439908080456564D0,
     8      0.09563872007927485942D0,0.09654008851472780057D0/
      DATA X48/
     1      0.99877100725242611860D0,0.99353017226635075755D0,
     2      0.98412458372282685774D0,0.97059159254624725046D0,
     3      0.95298770316043086072D0,0.93138669070655433311D0,
     4      0.90587913671556967282D0,0.87657202027424788591D0,
     5      0.84358826162439353071D0,0.80706620402944262708D0,
     6      0.76715903251574033925D0,0.72403413092381465467D0,
     7      0.67787237963266390521D0,0.62886739677651362400D0,
     8      0.57722472608397270382D0,0.52316097472223303368D0,
     9      0.46690290475095840454D0,0.40868648199071672992D0,
     1      0.34875588629216073816D0,0.28736248735545557674D0,
     2      0.22476379039468906122D0,0.16122235606889171806D0,
     3      0.09700469920946269893D0,0.03238017096286936203D0/
      DATA W48/
     1      0.00315334605230583862D0,0.00732755390127626210D0,
     2      0.01147723457923453948D0,0.01557931572294384873D0,
     3      0.01961616045735552781D0,0.02357076083932437914D0,
     4      0.02742650970835694820D0,0.03116722783279808890D0,
     5      0.03477722256477043889D0,0.03824135106583070632D0,
     6      0.04154508294346474921D0,0.04467456085669428042D0,
     7      0.04761665849249047482D0,0.05035903555385447496D0,
     8      0.05289018948519366710D0,0.05519950369998416287D0,
     9      0.05727729210040321570D0,0.05911483969839563575D0,
     1  0.06070443916589388005D0,0.06203942315989266390D0,
     2      0.06311419228625402566D0,0.06392423858464818662D0,
     3      0.06446616443595008221D0,0.06473769681268392250D0/
      NPQH=NPQ*2+1
      GO TO (2,4,6,8,10,12,14,16,18,20,22,24,26,99,99,32),NPQ
      IF(NPQ.EQ.24) GO TO 48
C 
    2  CONTINUE
      DO 200 I=1,NPQ
      II=NPQH-I
      XX(II)=X2(I)
      XX(I)=-XX(II)
      WW(I)=W2(I)
  200      WW(II)=W2(I)
      RETURN
    4  CONTINUE
      DO 400 I=1,NPQ
      II=NPQH-I
      XX(II)=X4(I)
      XX(I)=-XX(II)
      WW(I)=W4(I)
  400      WW(II)=W4(I)
      RETURN
    6  CONTINUE
      DO 600 I=1,NPQ
      II=NPQH-I
      XX(II)=X6(I)
      XX(I)=-XX(II)
      WW(I)=W6(I)
  600      WW(II)=W6(I)
      RETURN
    8  CONTINUE
      DO 800 I=1,NPQ
      II=NPQH-I
      XX(II)=X8(I)
      XX(I)=-XX(II)
      WW(I)=W8(I)
  800      WW(II)=W8(I)
      RETURN
   10   CONTINUE
      DO 1000 I=1,NPQ
      II=NPQH-I
      XX(II)=X10(I)
      XX(I)=-XX(II)
      WW(I)=W10(I)
 1000      WW(II)=W10(I)
      RETURN
   12   CONTINUE
      DO 1200 I=1,NPQ
      II=NPQH-I
      XX(II)=X12(I)
      XX(I)=-XX(II)
            WW(I)=W12(I)
 1200      WW(II)=W12(I)
      RETURN
   14  CONTINUE
      DO 1400 I=1,NPQ
      II=NPQH-I
      XX(II)=X14(I)
      XX(I)=-XX(II)
      WW(I)=W14(I)
 1400      WW(II)=W14(I)
      RETURN
   16   CONTINUE
      DO 1600 I=1,NPQ
      II=NPQH-I
      XX(II)=X16(I)
      XX(I)=-XX(II)
      WW(I)=W16(I)
 1600      WW(II)=W16(I)
      RETURN
   18   CONTINUE
      DO 1800 I=1,NPQ
      II=NPQH-I
      XX(II)=X18(I)
      XX(I)=-XX(II)
      WW(I)=W18(I)
 1800      WW(II)=W18(I)
      RETURN
   20   CONTINUE
      DO 2000 I=1,NPQ
      II=NPQH-I
      XX(II)=X20(I)
      XX(I)=-XX(II)
      WW(I)=W20(I)
 2000      WW(II)=W20(I)
      RETURN
   22   CONTINUE
      DO 2200 I=1,NPQ
      II=NPQH-I
      XX(II)=X22(I)
      XX(I)=-XX(II)
      WW(I)=W22(I)
 2200      WW(II)=W22(I)
      RETURN
   24   CONTINUE
      DO 2400 I=1,NPQ
      II=NPQH-I
      XX(II)=X24(I)
      XX(I)=-XX(II)
      WW(I)=W24(I)
 2400      WW(II)=W24(I)
      RETURN
   26   CONTINUE
      DO 2600 I=1,NPQ
      II=NPQH-I
      XX(II)=X26(I)
      XX(I)=-XX(II)
      WW(I)=W26(I)
 2600      WW(II)=W26(I)
      RETURN
   32   CONTINUE
      DO 3200 I=1,NPQ
      II=NPQH-I
      XX(II)=X32(I)
      XX(I)=-XX(II)
      WW(I)=W32(I)
 3200      WW(II)=W32(I)
      RETURN
   48 CONTINUE
      DO 4800 I=1,NPQ
      II=NPQH-I
      XX(II)=X48(I)
      XX(I)=-XX(II)
      WW(I)=W48(I)
 4800      WW(II)=W48(I)
      RETURN
   99 WRITE(6,300)
  300 FORMAT(//'??? 28 OR 30 POINT-QUADRATURE IS NOT DEFINED ??')

      end subroutine gauss
c----------------------------------------------------------------

c----------------------------------------------------------------
      subroutine lebedev_driver(x, y, z, w, n)

      implicit none

      integer, intent(in) :: n
      real*8, intent(out) :: x(*), y(*), z(*), w(*)

c***
c  This is a driver to the Lebedev-Laikov subroutine downloaded from the CCL homepage:
c
c  http://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/index.shtml
c
c  1) As requested, the following paper should be cited:
c
c      [1] V.I. Lebedev, and D.N. Laikov
c          "A quadrature formula for the sphere of the 131st
c           algebraic order of accuracy"
c          Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
c
c  2) The HUGE quadratures below were removed:
c     N = 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
c***

      if( n .eq. 6 ) then
         call LD0006(X,Y,Z,W,N)
      else if( n .eq. 14 ) then
         call LD0014(X,Y,Z,W,N)
      else if( n .eq. 26 ) then
         call LD0026(X,Y,Z,W,N)
      else if( n .eq. 38 ) then
         call LD0038(X,Y,Z,W,N)
      else if( n .eq. 50 ) then
         call LD0050(X,Y,Z,W,N)
      else if( n .eq. 74 ) then
         call LD0074(X,Y,Z,W,N)
      else if( n .eq. 86 ) then
         call LD0086(X,Y,Z,W,N)
      else if( n .eq. 110 ) then
         call LD0110(X,Y,Z,W,N)
      else if( n .eq. 146 ) then
         call LD0146(X,Y,Z,W,N)
      else if( n .eq. 170 ) then
         call LD0170(X,Y,Z,W,N)
      else if( n .eq. 194 ) then
         call LD0194(X,Y,Z,W,N)
      else if( n .eq. 230 ) then
         call LD0230(X,Y,Z,W,N)
      else if( n .eq. 266 ) then
         call LD0266(X,Y,Z,W,N)
      else if( n .eq. 302 ) then
         call LD0302(X,Y,Z,W,N)
      else if( n .eq. 350 ) then
         call LD0350(X,Y,Z,W,N)
      else if( n .eq. 434 ) then
         call LD0434(X,Y,Z,W,N)
      else if( n .eq. 590 ) then
         call LD0590(X,Y,Z,W,N)
      else if( n .eq. 770 ) then
         call LD0770(X,Y,Z,W,N)
      else if( n .eq. 974 ) then
         call LD0974(X,Y,Z,W,N)
      else if( n .eq. 1202) then
         call LD1202(X,Y,Z,W,N)
      else if( n .eq. 1454) then
         call LD1454(X,Y,Z,W,N)
      else if( n .eq. 1730) then
         call LD1730(X,Y,Z,W,N)
      else if( n .eq. 2030) then
         call LD2030(X,Y,Z,W,N)
      else if( n .eq. 2354) then
         call LD2354(X,Y,Z,W,N)
ccc   else if( n .eq. 2702) then
ccc      call LD2702(X,Y,Z,W,N)
ccc   else if( n .eq. 3074) then
ccc      call LD3074(X,Y,Z,W,N)
ccc   else if( n .eq. 3470) then
ccc      call LD3470(X,Y,Z,W,N)
ccc   else if( n .eq. 3890) then
ccc      call LD3890(X,Y,Z,W,N)
ccc   else if( n .eq. 4334) then
ccc      call LD4334(X,Y,Z,W,N)
ccc   else if( n .eq. 4802) then
ccc      call LD4802(X,Y,Z,W,N)
ccc   else if( n .eq. 5294) then
ccc      call LD5294(X,Y,Z,W,N)
ccc   else if( n .eq. 5810) then
ccc      call LD5810(X,Y,Z,W,N)
      else
         write(6,"('INVALID NUMBER OF LEBEDEV-LAIKOV',
     .             ' QUADRATURE POINTS')")
         write(6,"('THE ALLOWED VALUES ARE:')")
         write(6,"('0006 0014 0026 0038 0050 0074 0086 0110')")
         write(6,"('0146 0170 0194 0230 0266 0302 0350 0434')")
         write(6,"('0590 0770 0974 1202 1454 1730 2030 2354')")

         write(6,"(/,'ABORT...')")
         stop
      end if

      end subroutine lebedev_driver
c----------------------------------------------------------------


c----------------------------------------------------------------
       subroutine gen_oh(code, num, x, y, z, w, a, b, v)
       implicit logical(a-z)
       double precision x(*),y(*),z(*),w(*)
       double precision a,b,v
       integer code
       integer num
       double precision c
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated from C to fortran77 by hand.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
cvw
cvw    Given a point on a sphere (specified by a and b), generate all
cvw    the equivalent points under Oh symmetry, making grid points with
cvw    weight v.
cvw    The variable num is increased by the number of different points
cvw    generated.
cvw
cvw    Depending on code, there are 6...48 different but equivalent
cvw    points.
cvw
cvw    code=1:   (0,0,1) etc                                (  6 points)
cvw    code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
cvw    code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
cvw    code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
cvw    code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
cvw    code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
cvw
       goto (1,2,3,4,5,6) code
       write (6,*) 'Gen_Oh: Invalid Code'
       stop 
    1  continue
       a=1.0d0
       x(1) =  a
       y(1) =  0.0d0
       z(1) =  0.0d0
       w(1) =  v
       x(2) = -a
       y(2) =  0.0d0
       z(2) =  0.0d0
       w(2) =  v
       x(3) =  0.0d0
       y(3) =  a
       z(3) =  0.0d0
       w(3) =  v
       x(4) =  0.0d0
       y(4) = -a
       z(4) =  0.0d0
       w(4) =  v
       x(5) =  0.0d0
       y(5) =  0.0d0
       z(5) =  a
       w(5) =  v
       x(6) =  0.0d0
       y(6) =  0.0d0
       z(6) = -a
       w(6) =  v
       num=num+6
       return
cvw
    2  continue
       a=sqrt(0.5d0)
       x( 1) =  0d0
       y( 1) =  a
       z( 1) =  a
       w( 1) =  v
       x( 2) =  0d0
       y( 2) = -a
       z( 2) =  a
       w( 2) =  v
       x( 3) =  0d0
       y( 3) =  a
       z( 3) = -a
       w( 3) =  v
       x( 4) =  0d0
       y( 4) = -a
       z( 4) = -a
       w( 4) =  v
       x( 5) =  a
       y( 5) =  0d0
       z( 5) =  a
       w( 5) =  v
       x( 6) = -a
       y( 6) =  0d0
       z( 6) =  a
       w( 6) =  v
       x( 7) =  a
       y( 7) =  0d0
       z( 7) = -a
       w( 7) =  v
       x( 8) = -a
       y( 8) =  0d0
       z( 8) = -a
       w( 8) =  v
       x( 9) =  a
       y( 9) =  a
       z( 9) =  0d0
       w( 9) =  v
       x(10) = -a
       y(10) =  a
       z(10) =  0d0
       w(10) =  v
       x(11) =  a
       y(11) = -a
       z(11) =  0d0
       w(11) =  v
       x(12) = -a
       y(12) = -a
       z(12) =  0d0
       w(12) =  v
       num=num+12
       return
cvw
    3  continue
       a = sqrt(1d0/3d0)
       x(1) =  a
       y(1) =  a
       z(1) =  a
       w(1) =  v
       x(2) = -a
       y(2) =  a
       z(2) =  a
       w(2) =  v
       x(3) =  a
       y(3) = -a
       z(3) =  a
       w(3) =  v
       x(4) = -a
       y(4) = -a
       z(4) =  a
       w(4) =  v
       x(5) =  a
       y(5) =  a
       z(5) = -a
       w(5) =  v
       x(6) = -a
       y(6) =  a
       z(6) = -a
       w(6) =  v
       x(7) =  a
       y(7) = -a
       z(7) = -a
       w(7) =  v
       x(8) = -a
       y(8) = -a
       z(8) = -a
       w(8) =  v
       num=num+8
       return
cvw
    4  continue
       b = sqrt(1d0 - 2d0*a*a)
       x( 1) =  a
       y( 1) =  a
       z( 1) =  b
       w( 1) =  v
       x( 2) = -a
       y( 2) =  a
       z( 2) =  b
       w( 2) =  v
       x( 3) =  a
       y( 3) = -a
       z( 3) =  b
       w( 3) =  v
       x( 4) = -a
       y( 4) = -a
       z( 4) =  b
       w( 4) =  v
       x( 5) =  a
       y( 5) =  a
       z( 5) = -b
       w( 5) =  v
       x( 6) = -a
       y( 6) =  a
       z( 6) = -b
       w( 6) =  v
       x( 7) =  a
       y( 7) = -a
       z( 7) = -b
       w( 7) =  v
       x( 8) = -a
       y( 8) = -a
       z( 8) = -b
       w( 8) =  v
       x( 9) =  a
       y( 9) =  b
       z( 9) =  a
       w( 9) =  v
       x(10) = -a
       y(10) =  b
       z(10) =  a
       w(10) =  v
       x(11) =  a
       y(11) = -b
       z(11) =  a
       w(11) =  v
       x(12) = -a
       y(12) = -b
       z(12) =  a
       w(12) =  v
       x(13) =  a
       y(13) =  b
       z(13) = -a
       w(13) =  v
       x(14) = -a
       y(14) =  b
       z(14) = -a
       w(14) =  v
       x(15) =  a
       y(15) = -b
       z(15) = -a
       w(15) =  v
       x(16) = -a
       y(16) = -b
       z(16) = -a
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  a
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  a
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  a
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  a
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -a
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -a
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -a
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -a
       w(24) =  v
       num=num+24
       return
cvw
    5  continue
       b=sqrt(1d0-a*a)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  0d0
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  0d0
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  0d0
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  0d0
       w( 4) =  v
       x( 5) =  b
       y( 5) =  a
       z( 5) =  0d0
       w( 5) =  v
       x( 6) = -b
       y( 6) =  a
       z( 6) =  0d0
       w( 6) =  v
       x( 7) =  b
       y( 7) = -a
       z( 7) =  0d0
       w( 7) =  v
       x( 8) = -b
       y( 8) = -a
       z( 8) =  0d0
       w( 8) =  v
       x( 9) =  a
       y( 9) =  0d0
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  0d0
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) =  0d0
       z(11) = -b
       w(11) =  v
       x(12) = -a
       y(12) =  0d0
       z(12) = -b
       w(12) =  v
       x(13) =  b
       y(13) =  0d0
       z(13) =  a
       w(13) =  v
       x(14) = -b
       y(14) =  0d0
       z(14) =  a
       w(14) =  v
       x(15) =  b
       y(15) =  0d0
       z(15) = -a
       w(15) =  v
       x(16) = -b
       y(16) =  0d0
       z(16) = -a
       w(16) =  v
       x(17) =  0d0
       y(17) =  a
       z(17) =  b
       w(17) =  v
       x(18) =  0d0
       y(18) = -a
       z(18) =  b
       w(18) =  v
       x(19) =  0d0
       y(19) =  a
       z(19) = -b
       w(19) =  v
       x(20) =  0d0
       y(20) = -a
       z(20) = -b
       w(20) =  v
       x(21) =  0d0
       y(21) =  b
       z(21) =  a
       w(21) =  v
       x(22) =  0d0
       y(22) = -b
       z(22) =  a
       w(22) =  v
       x(23) =  0d0
       y(23) =  b
       z(23) = -a
       w(23) =  v
       x(24) =  0d0
       y(24) = -b
       z(24) = -a
       w(24) =  v
       num=num+24
       return
cvw
    6  continue
       c=sqrt(1d0 - a*a - b*b)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  c
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  c
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  c
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  c
       w( 4) =  v
       x( 5) =  a
       y( 5) =  b
       z( 5) = -c
       w( 5) =  v
       x( 6) = -a
       y( 6) =  b
       z( 6) = -c
       w( 6) =  v
       x( 7) =  a
       y( 7) = -b
       z( 7) = -c
       w( 7) =  v
       x( 8) = -a
       y( 8) = -b
       z( 8) = -c
       w( 8) =  v
       x( 9) =  a
       y( 9) =  c
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  c
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) = -c
       z(11) =  b
       w(11) =  v
       x(12) = -a
       y(12) = -c
       z(12) =  b
       w(12) =  v
       x(13) =  a
       y(13) =  c
       z(13) = -b
       w(13) =  v
       x(14) = -a
       y(14) =  c
       z(14) = -b
       w(14) =  v
       x(15) =  a
       y(15) = -c
       z(15) = -b
       w(15) =  v
       x(16) = -a
       y(16) = -c
       z(16) = -b
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  c
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  c
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  c
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  c
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -c
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -c
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -c
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -c
       w(24) =  v
       x(25) =  b
       y(25) =  c
       z(25) =  a
       w(25) =  v
       x(26) = -b
       y(26) =  c
       z(26) =  a
       w(26) =  v
       x(27) =  b
       y(27) = -c
       z(27) =  a
       w(27) =  v
       x(28) = -b
       y(28) = -c
       z(28) =  a
       w(28) =  v
       x(29) =  b
       y(29) =  c
       z(29) = -a
       w(29) =  v
       x(30) = -b
       y(30) =  c
       z(30) = -a
       w(30) =  v
       x(31) =  b
       y(31) = -c
       z(31) = -a
       w(31) =  v
       x(32) = -b
       y(32) = -c
       z(32) = -a
       w(32) =  v
       x(33) =  c
       y(33) =  a
       z(33) =  b
       w(33) =  v
       x(34) = -c
       y(34) =  a
       z(34) =  b
       w(34) =  v
       x(35) =  c
       y(35) = -a
       z(35) =  b
       w(35) =  v
       x(36) = -c
       y(36) = -a
       z(36) =  b
       w(36) =  v
       x(37) =  c
       y(37) =  a
       z(37) = -b
       w(37) =  v
       x(38) = -c
       y(38) =  a
       z(38) = -b
       w(38) =  v
       x(39) =  c
       y(39) = -a
       z(39) = -b
       w(39) =  v
       x(40) = -c
       y(40) = -a
       z(40) = -b
       w(40) =  v
       x(41) =  c
       y(41) =  b
       z(41) =  a
       w(41) =  v
       x(42) = -c
       y(42) =  b
       z(42) =  a
       w(42) =  v
       x(43) =  c
       y(43) = -b
       z(43) =  a
       w(43) =  v
       x(44) = -c
       y(44) = -b
       z(44) =  a
       w(44) =  v
       x(45) =  c
       y(45) =  b
       z(45) = -a
       w(45) =  v
       x(46) = -c
       y(46) =  b
       z(46) = -a
       w(46) =  v
       x(47) =  c
       y(47) = -b
       z(47) = -a
       w(47) =  v
       x(48) = -c
       y(48) = -b
       z(48) = -a
       w(48) =  v
       num=num+48
       return
       END subroutine gen_oh
       SUBROUTINE LD0006(X,Y,Z,W,N)
       DOUBLE PRECISION X(   6)
       DOUBLE PRECISION Y(   6)
       DOUBLE PRECISION Z(   6)
       DOUBLE PRECISION W(   6)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV    6-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.1666666666666667D+0
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0006
       SUBROUTINE LD0014(X,Y,Z,W,N)
       DOUBLE PRECISION X(  14)
       DOUBLE PRECISION Y(  14)
       DOUBLE PRECISION Z(  14)
       DOUBLE PRECISION W(  14)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   14-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.6666666666666667D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7500000000000000D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0014
       SUBROUTINE LD0026(X,Y,Z,W,N)
       DOUBLE PRECISION X(  26)
       DOUBLE PRECISION Y(  26)
       DOUBLE PRECISION Z(  26)
       DOUBLE PRECISION W(  26)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   26-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.4761904761904762D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3809523809523810D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0026
       SUBROUTINE LD0038(X,Y,Z,W,N)
       DOUBLE PRECISION X(  38)
       DOUBLE PRECISION Y(  38)
       DOUBLE PRECISION Z(  38)
       DOUBLE PRECISION W(  38)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   38-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.9523809523809524D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4597008433809831D+0
       V=0.2857142857142857D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0038
       SUBROUTINE LD0050(X,Y,Z,W,N)
       DOUBLE PRECISION X(  50)
       DOUBLE PRECISION Y(  50)
       DOUBLE PRECISION Z(  50)
       DOUBLE PRECISION W(  50)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   50-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.1269841269841270D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2257495590828924D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2109375000000000D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3015113445777636D+0
       V=0.2017333553791887D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0050
       SUBROUTINE LD0074(X,Y,Z,W,N)
       DOUBLE PRECISION X(  74)
       DOUBLE PRECISION Y(  74)
       DOUBLE PRECISION Z(  74)
       DOUBLE PRECISION W(  74)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   74-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.5130671797338464D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1660406956574204D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2958603896103896D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4803844614152614D+0
       V=0.2657620708215946D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3207726489807764D+0
       V=0.1652217099371571D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0074
       SUBROUTINE LD0086(X,Y,Z,W,N)
       DOUBLE PRECISION X(  86)
       DOUBLE PRECISION Y(  86)
       DOUBLE PRECISION Z(  86)
       DOUBLE PRECISION W(  86)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   86-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.1154401154401154D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1194390908585628D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3696028464541502D+0
       V=0.1111055571060340D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6943540066026664D+0
       V=0.1187650129453714D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3742430390903412D+0
       V=0.1181230374690448D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0086
       SUBROUTINE LD0110(X,Y,Z,W,N)
       DOUBLE PRECISION X( 110)
       DOUBLE PRECISION Y( 110)
       DOUBLE PRECISION Z( 110)
       DOUBLE PRECISION W( 110)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  110-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.3828270494937162D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9793737512487512D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1851156353447362D+0
       V=0.8211737283191111D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6904210483822922D+0
       V=0.9942814891178103D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3956894730559419D+0
       V=0.9595471336070963D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4783690288121502D+0
       V=0.9694996361663028D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0110
       SUBROUTINE LD0146(X,Y,Z,W,N)
       DOUBLE PRECISION X( 146)
       DOUBLE PRECISION Y( 146)
       DOUBLE PRECISION Z( 146)
       DOUBLE PRECISION W( 146)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  146-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.5996313688621381D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7372999718620756D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7210515360144488D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6764410400114264D+0
       V=0.7116355493117555D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174961227965453D+0
       V=0.6753829486314477D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1574676672039082D+0
       V=0.7574394159054034D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1403553811713183D+0
       B=0.4493328323269557D+0
       V=0.6991087353303262D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0146
       SUBROUTINE LD0170(X,Y,Z,W,N)
       DOUBLE PRECISION X( 170)
       DOUBLE PRECISION Y( 170)
       DOUBLE PRECISION Z( 170)
       DOUBLE PRECISION W( 170)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  170-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.5544842902037365D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6071332770670752D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6383674773515093D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2551252621114134D+0
       V=0.5183387587747790D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6743601460362766D+0
       V=0.6317929009813725D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318910696719410D+0
       V=0.6201670006589077D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2613931360335988D+0
       V=0.5477143385137348D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4990453161796037D+0
       B=0.1446630744325115D+0
       V=0.5968383987681156D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0170
       SUBROUTINE LD0194(X,Y,Z,W,N)
       DOUBLE PRECISION X( 194)
       DOUBLE PRECISION Y( 194)
       DOUBLE PRECISION Z( 194)
       DOUBLE PRECISION W( 194)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  194-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.1782340447244611D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5716905949977102D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5573383178848738D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6712973442695226D+0
       V=0.5608704082587997D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2892465627575439D+0
       V=0.5158237711805383D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4446933178717437D+0
       V=0.5518771467273614D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1299335447650067D+0
       V=0.4106777028169394D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3457702197611283D+0
       V=0.5051846064614808D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1590417105383530D+0
       B=0.8360360154824589D+0
       V=0.5530248916233094D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0194
       SUBROUTINE LD0230(X,Y,Z,W,N)
       DOUBLE PRECISION X( 230)
       DOUBLE PRECISION Y( 230)
       DOUBLE PRECISION Z( 230)
       DOUBLE PRECISION W( 230)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  230-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=-0.5522639919727325D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4450274607445226D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4492044687397611D+0
       V=0.4496841067921404D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2520419490210201D+0
       V=0.5049153450478750D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6981906658447242D+0
       V=0.3976408018051883D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6587405243460960D+0
       V=0.4401400650381014D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4038544050097660D-1
       V=0.1724544350544401D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5823842309715585D+0
       V=0.4231083095357343D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3545877390518688D+0
       V=0.5198069864064399D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2272181808998187D+0
       B=0.4864661535886647D+0
       V=0.4695720972568883D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0230
       SUBROUTINE LD0266(X,Y,Z,W,N)
       DOUBLE PRECISION X( 266)
       DOUBLE PRECISION Y( 266)
       DOUBLE PRECISION Z( 266)
       DOUBLE PRECISION W( 266)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  266-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=-0.1313769127326952D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2522728704859336D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4186853881700583D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7039373391585475D+0
       V=0.5315167977810885D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1012526248572414D+0
       V=0.4047142377086219D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4647448726420539D+0
       V=0.4112482394406990D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3277420654971629D+0
       V=0.3595584899758782D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6620338663699974D+0
       V=0.4256131351428158D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8506508083520399D+0
       V=0.4229582700647240D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3233484542692899D+0
       B=0.1153112011009701D+0
       V=0.4080914225780505D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2314790158712601D+0
       B=0.5244939240922365D+0
       V=0.4071467593830964D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0266
       SUBROUTINE LD0302(X,Y,Z,W,N)
       DOUBLE PRECISION X( 302)
       DOUBLE PRECISION Y( 302)
       DOUBLE PRECISION Z( 302)
       DOUBLE PRECISION W( 302)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  302-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.8545911725128148D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3599119285025571D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3515640345570105D+0
       V=0.3449788424305883D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6566329410219612D+0
       V=0.3604822601419882D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4729054132581005D+0
       V=0.3576729661743367D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9618308522614784D-1
       V=0.2352101413689164D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2219645236294178D+0
       V=0.3108953122413675D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7011766416089545D+0
       V=0.3650045807677255D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644152887060663D+0
       V=0.2982344963171804D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718955891878961D+0
       V=0.3600820932216460D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2510034751770465D+0
       B=0.8000727494073952D+0
       V=0.3571540554273387D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1233548532583327D+0
       B=0.4127724083168531D+0
       V=0.3392312205006170D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0302
       SUBROUTINE LD0350(X,Y,Z,W,N)
       DOUBLE PRECISION X( 350)
       DOUBLE PRECISION Y( 350)
       DOUBLE PRECISION Z( 350)
       DOUBLE PRECISION W( 350)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  350-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.3006796749453936D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3050627745650771D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7068965463912316D+0
       V=0.1621104600288991D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4794682625712025D+0
       V=0.3005701484901752D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1927533154878019D+0
       V=0.2990992529653774D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6930357961327123D+0
       V=0.2982170644107595D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3608302115520091D+0
       V=0.2721564237310992D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6498486161496169D+0
       V=0.3033513795811141D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1932945013230339D+0
       V=0.3007949555218533D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3800494919899303D+0
       V=0.2881964603055307D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2899558825499574D+0
       B=0.7934537856582316D+0
       V=0.2958357626535696D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9684121455103957D-1
       B=0.8280801506686862D+0
       V=0.3036020026407088D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1833434647041659D+0
       B=0.9074658265305127D+0
       V=0.2832187403926303D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0350
       SUBROUTINE LD0434(X,Y,Z,W,N)
       DOUBLE PRECISION X( 434)
       DOUBLE PRECISION Y( 434)
       DOUBLE PRECISION Z( 434)
       DOUBLE PRECISION W( 434)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  434-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.5265897968224436D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2548219972002607D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2512317418927307D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6909346307509111D+0
       V=0.2530403801186355D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1774836054609158D+0
       V=0.2014279020918528D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4914342637784746D+0
       V=0.2501725168402936D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6456664707424256D+0
       V=0.2513267174597564D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2861289010307638D+0
       V=0.2302694782227416D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7568084367178018D-1
       V=0.1462495621594614D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927259763368002D+0
       V=0.2445373437312980D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8818132877794288D+0
       V=0.2417442375638981D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9776428111182649D+0
       V=0.1910951282179532D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2054823696403044D+0
       B=0.8689460322872412D+0
       V=0.2416930044324775D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5905157048925271D+0
       B=0.7999278543857286D+0
       V=0.2512236854563495D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5550152361076807D+0
       B=0.7717462626915901D+0
       V=0.2496644054553086D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9371809858553722D+0
       B=0.3344363145343455D+0
       V=0.2236607760437849D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0434
       SUBROUTINE LD0590(X,Y,Z,W,N)
       DOUBLE PRECISION X( 590)
       DOUBLE PRECISION Y( 590)
       DOUBLE PRECISION Z( 590)
       DOUBLE PRECISION W( 590)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  590-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.3095121295306187D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1852379698597489D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7040954938227469D+0
       V=0.1871790639277744D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6807744066455243D+0
       V=0.1858812585438317D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6372546939258752D+0
       V=0.1852028828296213D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044419707800358D+0
       V=0.1846715956151242D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4215761784010967D+0
       V=0.1818471778162769D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3317920736472123D+0
       V=0.1749564657281154D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2384736701421887D+0
       V=0.1617210647254411D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1459036449157763D+0
       V=0.1384737234851692D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095034115507196D-1
       V=0.9764331165051050D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6116843442009876D+0
       V=0.1857161196774078D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3964755348199858D+0
       V=0.1705153996395864D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1724782009907724D+0
       V=0.1300321685886048D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5610263808622060D+0
       B=0.3518280927733519D+0
       V=0.1842866472905286D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4742392842551980D+0
       B=0.2634716655937950D+0
       V=0.1802658934377451D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5984126497885380D+0
       B=0.1816640840360209D+0
       V=0.1849830560443660D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3791035407695563D+0
       B=0.1720795225656878D+0
       V=0.1713904507106709D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778673190586244D+0
       B=0.8213021581932511D-1
       V=0.1555213603396808D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5033564271075117D+0
       B=0.8999205842074875D-1
       V=0.1802239128008525D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0590
       SUBROUTINE LD0770(X,Y,Z,W,N)
       DOUBLE PRECISION X( 770)
       DOUBLE PRECISION Y( 770)
       DOUBLE PRECISION Z( 770)
       DOUBLE PRECISION W( 770)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  770-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.2192942088181184D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1436433617319080D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1421940344335877D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5087204410502360D-1
       V=0.6798123511050502D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1228198790178831D+0
       V=0.9913184235294912D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2026890814408786D+0
       V=0.1180207833238949D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2847745156464294D+0
       V=0.1296599602080921D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3656719078978026D+0
       V=0.1365871427428316D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4428264886713469D+0
       V=0.1402988604775325D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140619627249735D+0
       V=0.1418645563595609D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6306401219166803D+0
       V=0.1421376741851662D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6716883332022612D+0
       V=0.1423996475490962D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6979792685336881D+0
       V=0.1431554042178567D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1446865674195309D+0
       V=0.9254401499865368D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3390263475411216D+0
       V=0.1250239995053509D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5335804651263506D+0
       V=0.1394365843329230D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6944024393349413D-1
       B=0.2355187894242326D+0
       V=0.1127089094671749D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2269004109529460D+0
       B=0.4102182474045730D+0
       V=0.1345753760910670D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8025574607775339D-1
       B=0.6214302417481605D+0
       V=0.1424957283316783D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1467999527896572D+0
       B=0.3245284345717394D+0
       V=0.1261523341237750D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1571507769824727D+0
       B=0.5224482189696630D+0
       V=0.1392547106052696D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2365702993157246D+0
       B=0.6017546634089558D+0
       V=0.1418761677877656D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7714815866765732D-1
       B=0.4346575516141163D+0
       V=0.1338366684479554D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062936666210730D+0
       B=0.4908826589037616D+0
       V=0.1393700862676131D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3822477379524787D+0
       B=0.5648768149099500D+0
       V=0.1415914757466932D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0770
       SUBROUTINE LD0974(X,Y,Z,W,N)
       DOUBLE PRECISION X( 974)
       DOUBLE PRECISION Y( 974)
       DOUBLE PRECISION Z( 974)
       DOUBLE PRECISION W( 974)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  974-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.1438294190527431D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1125772288287004D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4292963545341347D-1
       V=0.4948029341949241D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1051426854086404D+0
       V=0.7357990109125470D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1750024867623087D+0
       V=0.8889132771304384D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2477653379650257D+0
       V=0.9888347838921435D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3206567123955957D+0
       V=0.1053299681709471D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3916520749849983D+0
       V=0.1092778807014578D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4590825874187624D+0
       V=0.1114389394063227D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5214563888415861D+0
       V=0.1123724788051555D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6253170244654199D+0
       V=0.1125239325243814D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6637926744523170D+0
       V=0.1126153271815905D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910410398498301D+0
       V=0.1130286931123841D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7052907007457760D+0
       V=0.1134986534363955D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1236686762657990D+0
       V=0.6823367927109931D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2940777114468387D+0
       V=0.9454158160447096D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4697753849207649D+0
       V=0.1074429975385679D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6334563241139567D+0
       V=0.1129300086569132D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974048614181342D-1
       B=0.2029128752777523D+0
       V=0.8436884500901954D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1375760408473636D+0
       B=0.4602621942484054D+0
       V=0.1075255720448885D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3391016526336286D+0
       B=0.5030673999662036D+0
       V=0.1108577236864462D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1271675191439820D+0
       B=0.2817606422442134D+0
       V=0.9566475323783357D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2693120740413512D+0
       B=0.4331561291720157D+0
       V=0.1080663250717391D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1419786452601918D+0
       B=0.6256167358580814D+0
       V=0.1126797131196295D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6709284600738255D-1
       B=0.3798395216859157D+0
       V=0.1022568715358061D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7057738183256172D-1
       B=0.5517505421423520D+0
       V=0.1108960267713108D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2783888477882155D+0
       B=0.6029619156159187D+0
       V=0.1122790653435766D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1979578938917407D+0
       B=0.3589606329589096D+0
       V=0.1032401847117460D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2087307061103274D+0
       B=0.5348666438135476D+0
       V=0.1107249382283854D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4055122137872836D+0
       B=0.5674997546074373D+0
       V=0.1121780048519972D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD0974
       SUBROUTINE LD1202(X,Y,Z,W,N)
       DOUBLE PRECISION X(1202)
       DOUBLE PRECISION Y(1202)
       DOUBLE PRECISION Z(1202)
       DOUBLE PRECISION W(1202)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 1202-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.1105189233267572D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9205232738090741D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9133159786443561D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712636449657089D-1
       V=0.3690421898017899D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9140060412262223D-1
       V=0.5603990928680660D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1531077852469906D+0
       V=0.6865297629282609D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2180928891660612D+0
       V=0.7720338551145630D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2839874532200175D+0
       V=0.8301545958894795D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3491177600963764D+0
       V=0.8686692550179628D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4121431461444309D+0
       V=0.8927076285846890D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4718993627149127D+0
       V=0.9060820238568219D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5273145452842337D+0
       V=0.9119777254940867D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6209475332444019D+0
       V=0.9128720138604181D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6569722711857291D+0
       V=0.9130714935691735D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6841788309070143D+0
       V=0.9152873784554116D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012604330123631D+0
       V=0.9187436274321654D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1072382215478166D+0
       V=0.5176977312965694D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2582068959496968D+0
       V=0.7331143682101417D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4172752955306717D+0
       V=0.8463232836379928D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5700366911792503D+0
       V=0.9031122694253992D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9827986018263947D+0
       B=0.1771774022615325D+0
       V=0.6485778453163257D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9624249230326228D+0
       B=0.2475716463426288D+0
       V=0.7435030910982369D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9402007994128811D+0
       B=0.3354616289066489D+0
       V=0.7998527891839054D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9320822040143202D+0
       B=0.3173615246611977D+0
       V=0.8101731497468018D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9043674199393299D+0
       B=0.4090268427085357D+0
       V=0.8483389574594331D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8912407560074747D+0
       B=0.3854291150669224D+0
       V=0.8556299257311812D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8676435628462708D+0
       B=0.4932221184851285D+0
       V=0.8803208679738260D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8581979986041619D+0
       B=0.4785320675922435D+0
       V=0.8811048182425720D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8396753624049856D+0
       B=0.4507422593157064D+0
       V=0.8850282341265444D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8165288564022188D+0
       B=0.5632123020762100D+0
       V=0.9021342299040653D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8015469370783529D+0
       B=0.5434303569693900D+0
       V=0.9010091677105086D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7773563069070351D+0
       B=0.5123518486419871D+0
       V=0.9022692938426915D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7661621213900394D+0
       B=0.6394279634749102D+0
       V=0.9158016174693465D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7553584143533510D+0
       B=0.6269805509024392D+0
       V=0.9131578003189435D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7344305757559503D+0
       B=0.6031161693096310D+0
       V=0.9107813579482705D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7043837184021765D+0
       B=0.5693702498468441D+0
       V=0.9105760258970126D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD1202
       SUBROUTINE LD1454(X,Y,Z,W,N)
       DOUBLE PRECISION X(1454)
       DOUBLE PRECISION Y(1454)
       DOUBLE PRECISION Z(1454)
       DOUBLE PRECISION W(1454)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 1454-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.7777160743261247D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7557646413004701D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3229290663413854D-1
       V=0.2841633806090617D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8036733271462222D-1
       V=0.4374419127053555D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1354289960531653D+0
       V=0.5417174740872172D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1938963861114426D+0
       V=0.6148000891358593D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2537343715011275D+0
       V=0.6664394485800705D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3135251434752570D+0
       V=0.7025039356923220D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3721558339375338D+0
       V=0.7268511789249627D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4286809575195696D+0
       V=0.7422637534208629D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4822510128282994D+0
       V=0.7509545035841214D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5320679333566263D+0
       V=0.7548535057718401D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6172998195394274D+0
       V=0.7554088969774001D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6510679849127481D+0
       V=0.7553147174442808D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6777315251687360D+0
       V=0.7564767653292297D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963109410648741D+0
       V=0.7587991808518730D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7058935009831749D+0
       V=0.7608261832033027D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9955546194091857D+0
       V=0.4021680447874916D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9734115901794209D+0
       V=0.5804871793945964D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9275693732388626D+0
       V=0.6792151955945159D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8568022422795103D+0
       V=0.7336741211286294D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7623495553719372D+0
       V=0.7581866300989608D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707522908892223D+0
       B=0.4387028039889501D+0
       V=0.7538257859800743D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5196463388403083D+0
       B=0.3858908414762617D+0
       V=0.7483517247053123D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4646337531215351D+0
       B=0.3301937372343854D+0
       V=0.7371763661112059D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4063901697557691D+0
       B=0.2725423573563777D+0
       V=0.7183448895756934D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3456329466643087D+0
       B=0.2139510237495250D+0
       V=0.6895815529822191D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831395121050332D+0
       B=0.1555922309786647D+0
       V=0.6480105801792886D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197682022925330D+0
       B=0.9892878979686097D-1
       V=0.5897558896594636D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1564696098650355D+0
       B=0.4598642910675510D-1
       V=0.5095708849247346D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6027356673721295D+0
       B=0.3376625140173426D+0
       V=0.7536906428909755D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496032320255096D+0
       B=0.2822301309727988D+0
       V=0.7472505965575118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4921707755234567D+0
       B=0.2248632342592540D+0
       V=0.7343017132279698D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4309422998598483D+0
       B=0.1666224723456479D+0
       V=0.7130871582177445D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3664108182313672D+0
       B=0.1086964901822169D+0
       V=0.6817022032112776D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2990189057758436D+0
       B=0.5251989784120085D-1
       V=0.6380941145604121D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6268724013144998D+0
       B=0.2297523657550023D+0
       V=0.7550381377920310D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707324144834607D+0
       B=0.1723080607093800D+0
       V=0.7478646640144802D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5096360901960365D+0
       B=0.1140238465390513D+0
       V=0.7335918720601220D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4438729938312456D+0
       B=0.5611522095882537D-1
       V=0.7110120527658118D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6419978471082389D+0
       B=0.1164174423140873D+0
       V=0.7571363978689501D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5817218061802611D+0
       B=0.5797589531445219D-1
       V=0.7489908329079234D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD1454
       SUBROUTINE LD1730(X,Y,Z,W,N)
       DOUBLE PRECISION X(1730)
       DOUBLE PRECISION Y(1730)
       DOUBLE PRECISION Z(1730)
       DOUBLE PRECISION W(1730)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 1730-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.6309049437420976D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6398287705571748D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6357185073530720D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2860923126194662D-1
       V=0.2221207162188168D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7142556767711522D-1
       V=0.3475784022286848D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1209199540995559D+0
       V=0.4350742443589804D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1738673106594379D+0
       V=0.4978569136522127D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2284645438467734D+0
       V=0.5435036221998053D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2834807671701512D+0
       V=0.5765913388219542D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3379680145467339D+0
       V=0.6001200359226003D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3911355454819537D+0
       V=0.6162178172717512D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4422860353001403D+0
       V=0.6265218152438485D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4907781568726057D+0
       V=0.6323987160974212D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5360006153211468D+0
       V=0.6350767851540569D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6142105973596603D+0
       V=0.6354362775297107D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6459300387977504D+0
       V=0.6352302462706235D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718056125089225D+0
       V=0.6358117881417972D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910888533186254D+0
       V=0.6373101590310117D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7030467416823252D+0
       V=0.6390428961368665D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8354951166354646D-1
       V=0.3186913449946576D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2050143009099486D+0
       V=0.4678028558591711D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3370208290706637D+0
       V=0.5538829697598626D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4689051484233963D+0
       V=0.6044475907190476D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5939400424557334D+0
       V=0.6313575103509012D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1394983311832261D+0
       B=0.4097581162050343D-1
       V=0.4078626431855630D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1967999180485014D+0
       B=0.8851987391293348D-1
       V=0.4759933057812725D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2546183732548967D+0
       B=0.1397680182969819D+0
       V=0.5268151186413440D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3121281074713875D+0
       B=0.1929452542226526D+0
       V=0.5643048560507316D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3685981078502492D+0
       B=0.2467898337061562D+0
       V=0.5914501076613073D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4233760321547856D+0
       B=0.3003104124785409D+0
       V=0.6104561257874195D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4758671236059246D+0
       B=0.3526684328175033D+0
       V=0.6230252860707806D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5255178579796463D+0
       B=0.4031134861145713D+0
       V=0.6305618761760796D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718025633734589D+0
       B=0.4509426448342351D+0
       V=0.6343092767597889D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2686927772723415D+0
       B=0.4711322502423248D-1
       V=0.5176268945737826D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3306006819904809D+0
       B=0.9784487303942695D-1
       V=0.5564840313313692D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3904906850594983D+0
       B=0.1505395810025273D+0
       V=0.5856426671038980D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4479957951904390D+0
       B=0.2039728156296050D+0
       V=0.6066386925777091D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027076848919780D+0
       B=0.2571529941121107D+0
       V=0.6208824962234458D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5542087392260217D+0
       B=0.3092191375815670D+0
       V=0.6296314297822907D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020850887375187D+0
       B=0.3593807506130276D+0
       V=0.6340423756791859D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4019851409179594D+0
       B=0.5063389934378671D-1
       V=0.5829627677107342D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4635614567449800D+0
       B=0.1032422269160612D+0
       V=0.6048693376081110D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215860931591575D+0
       B=0.1566322094006254D+0
       V=0.6202362317732461D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5758202499099271D+0
       B=0.2098082827491099D+0
       V=0.6299005328403779D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6259893683876795D+0
       B=0.2618824114553391D+0
       V=0.6347722390609353D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5313795124811891D+0
       B=0.5263245019338556D-1
       V=0.6203778981238834D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5893317955931995D+0
       B=0.1061059730982005D+0
       V=0.6308414671239979D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6426246321215801D+0
       B=0.1594171564034221D+0
       V=0.6362706466959498D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6511904367376113D+0
       B=0.5354789536565540D-1
       V=0.6375414170333233D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD1730
       SUBROUTINE LD2030(X,Y,Z,W,N)
       DOUBLE PRECISION X(2030)
       DOUBLE PRECISION Y(2030)
       DOUBLE PRECISION Z(2030)
       DOUBLE PRECISION W(2030)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 2030-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.4656031899197431D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5421549195295507D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2540835336814348D-1
       V=0.1778522133346553D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6399322800504915D-1
       V=0.2811325405682796D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1088269469804125D+0
       V=0.3548896312631459D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1570670798818287D+0
       V=0.4090310897173364D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2071163932282514D+0
       V=0.4493286134169965D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2578914044450844D+0
       V=0.4793728447962723D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3085687558169623D+0
       V=0.5015415319164265D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3584719706267024D+0
       V=0.5175127372677937D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4070135594428709D+0
       V=0.5285522262081019D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4536618626222638D+0
       V=0.5356832703713962D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4979195686463577D+0
       V=0.5397914736175170D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5393075111126999D+0
       V=0.5416899441599930D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6115617676843916D+0
       V=0.5419308476889938D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6414308435160159D+0
       V=0.5416936902030596D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6664099412721607D+0
       V=0.5419544338703164D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6859161771214913D+0
       V=0.5428983656630975D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6993625593503890D+0
       V=0.5442286500098193D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7062393387719380D+0
       V=0.5452250345057301D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7479028168349763D-1
       V=0.2568002497728530D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1848951153969366D+0
       V=0.3827211700292145D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3059529066581305D+0
       V=0.4579491561917824D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4285556101021362D+0
       V=0.5042003969083574D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5468758653496526D+0
       V=0.5312708889976025D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6565821978343439D+0
       V=0.5438401790747117D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1253901572367117D+0
       B=0.3681917226439641D-1
       V=0.3316041873197344D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1775721510383941D+0
       B=0.7982487607213301D-1
       V=0.3899113567153771D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2305693358216114D+0
       B=0.1264640966592335D+0
       V=0.4343343327201309D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2836502845992063D+0
       B=0.1751585683418957D+0
       V=0.4679415262318919D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3361794746232590D+0
       B=0.2247995907632670D+0
       V=0.4930847981631031D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3875979172264824D+0
       B=0.2745299257422246D+0
       V=0.5115031867540091D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4374019316999074D+0
       B=0.3236373482441118D+0
       V=0.5245217148457367D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4851275843340022D+0
       B=0.3714967859436741D+0
       V=0.5332041499895321D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5303391803806868D+0
       B=0.4175353646321745D+0
       V=0.5384583126021542D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726197380596287D+0
       B=0.4612084406355461D+0
       V=0.5411067210798852D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2431520732564863D+0
       B=0.4258040133043952D-1
       V=0.4259797391468714D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3002096800895869D+0
       B=0.8869424306722721D-1
       V=0.4604931368460021D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3558554457457432D+0
       B=0.1368811706510655D+0
       V=0.4871814878255202D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097782537048887D+0
       B=0.1860739985015033D+0
       V=0.5072242910074885D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4616337666067458D+0
       B=0.2354235077395853D+0
       V=0.5217069845235350D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5110707008417874D+0
       B=0.2842074921347011D+0
       V=0.5315785966280310D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577415286163795D+0
       B=0.3317784414984102D+0
       V=0.5376833708758905D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6013060431366950D+0
       B=0.3775299002040700D+0
       V=0.5408032092069521D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3661596767261781D+0
       B=0.4599367887164592D-1
       V=0.4842744917904866D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237633153506581D+0
       B=0.9404893773654421D-1
       V=0.5048926076188130D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4786328454658452D+0
       B=0.1431377109091971D+0
       V=0.5202607980478373D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5305702076789774D+0
       B=0.1924186388843570D+0
       V=0.5309932388325743D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5793436224231788D+0
       B=0.2411590944775190D+0
       V=0.5377419770895208D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6247069017094747D+0
       B=0.2886871491583605D+0
       V=0.5411696331677717D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4874315552535204D+0
       B=0.4804978774953206D-1
       V=0.5197996293282420D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5427337322059053D+0
       B=0.9716857199366665D-1
       V=0.5311120836622945D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5943493747246700D+0
       B=0.1465205839795055D+0
       V=0.5384309319956951D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6421314033564943D+0
       B=0.1953579449803574D+0
       V=0.5421859504051886D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020628374713980D+0
       B=0.4916375015738108D-1
       V=0.5390948355046314D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529222529856881D+0
       B=0.9861621540127005D-1
       V=0.5433312705027845D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD2030
       SUBROUTINE LD2354(X,Y,Z,W,N)
       DOUBLE PRECISION X(2354)
       DOUBLE PRECISION Y(2354)
       DOUBLE PRECISION Z(2354)
       DOUBLE PRECISION W(2354)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 2354-POINT ANGULAR GRID
CVW
chvd
chvd   This subroutine is part of a set of subroutines that generate
chvd   Lebedev grids [1-6] for integration on a sphere. The original 
chvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
chvd   translated into fortran by Dr. Christoph van Wuellen.
chvd   This subroutine was translated using a C to fortran77 conversion
chvd   tool written by Dr. Christoph van Wuellen.
chvd
chvd   Users of this code are asked to include reference [1] in their
chvd   publications, and in the user- and programmers-manuals 
chvd   describing their codes.
chvd
chvd   This code was distributed through CCL (http://www.ccl.net/).
chvd
chvd   [1] V.I. Lebedev, and D.N. Laikov
chvd       "A quadrature formula for the sphere of the 131st
chvd        algebraic order of accuracy"
chvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
chvd
chvd   [2] V.I. Lebedev
chvd       "A quadrature formula for the sphere of 59th algebraic
chvd        order of accuracy"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
chvd
chvd   [3] V.I. Lebedev, and A.L. Skorokhodov
chvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
chvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
chvd
chvd   [4] V.I. Lebedev
chvd       "Spherical quadrature formulas exact to orders 25-29"
chvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
chvd
chvd   [5] V.I. Lebedev
chvd       "Quadratures on a sphere"
chvd       Computational Mathematics and Mathematical Physics, Vol. 16,
chvd       1976, pp. 10-24. 
chvd
chvd   [6] V.I. Lebedev
chvd       "Values of the nodes and weights of ninth to seventeenth 
chvd        order Gauss-Markov quadrature formulae invariant under the
chvd        octahedron group with inversion"
chvd       Computational Mathematics and Mathematical Physics, Vol. 15,
chvd       1975, pp. 44-51.
chvd
       N=1
       V=0.3922616270665292D-4
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4703831750854424D-3
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4678202801282136D-3
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2290024646530589D-1
       V=0.1437832228979900D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5779086652271284D-1
       V=0.2303572493577644D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9863103576375984D-1
       V=0.2933110752447454D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1428155792982185D+0
       V=0.3402905998359838D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1888978116601463D+0
       V=0.3759138466870372D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2359091682970210D+0
       V=0.4030638447899798D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831228833706171D+0
       V=0.4236591432242211D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3299495857966693D+0
       V=0.4390522656946746D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758840802660796D+0
       V=0.4502523466626247D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4204751831009480D+0
       V=0.4580577727783541D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4633068518751051D+0
       V=0.4631391616615899D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5039849474507313D+0
       V=0.4660928953698676D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5421265793440747D+0
       V=0.4674751807936953D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6092660230557310D+0
       V=0.4676414903932920D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6374654204984869D+0
       V=0.4674086492347870D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6615136472609892D+0
       V=0.4674928539483207D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6809487285958127D+0
       V=0.4680748979686447D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6952980021665196D+0
       V=0.4690449806389040D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041245497695400D+0
       V=0.4699877075860818D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6744033088306065D-1
       V=0.2099942281069176D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1678684485334166D+0
       V=0.3172269150712804D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2793559049539613D+0
       V=0.3832051358546523D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935264218057639D+0
       V=0.4252193818146985D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5052629268232558D+0
       V=0.4513807963755000D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6107905315437531D+0
       V=0.4657797469114178D-3
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1135081039843524D+0
       B=0.3331954884662588D-1
       V=0.2733362800522836D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1612866626099378D+0
       B=0.7247167465436538D-1
       V=0.3235485368463559D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2100786550168205D+0
       B=0.1151539110849745D+0
       V=0.3624908726013453D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2592282009459942D+0
       B=0.1599491097143677D+0
       V=0.3925540070712828D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3081740561320203D+0
       B=0.2058699956028027D+0
       V=0.4156129781116235D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3564289781578164D+0
       B=0.2521624953502911D+0
       V=0.4330644984623263D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035587288240703D+0
       B=0.2982090785797674D+0
       V=0.4459677725921312D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491671196373903D+0
       B=0.3434762087235733D+0
       V=0.4551593004456795D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4928854782917489D+0
       B=0.3874831357203437D+0
       V=0.4613341462749918D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5343646791958988D+0
       B=0.4297814821746926D+0
       V=0.4651019618269806D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5732683216530990D+0
       B=0.4699402260943537D+0
       V=0.4670249536100625D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2214131583218986D+0
       B=0.3873602040643895D-1
       V=0.3549555576441708D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2741796504750071D+0
       B=0.8089496256902013D-1
       V=0.3856108245249010D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3259797439149485D+0
       B=0.1251732177620872D+0
       V=0.4098622845756882D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3765441148826891D+0
       B=0.1706260286403185D+0
       V=0.4286328604268950D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4255773574530558D+0
       B=0.2165115147300408D+0
       V=0.4427802198993945D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4727795117058430D+0
       B=0.2622089812225259D+0
       V=0.4530473511488561D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5178546895819012D+0
       B=0.3071721431296201D+0
       V=0.4600805475703138D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5605141192097460D+0
       B=0.3508998998801138D+0
       V=0.4644599059958017D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6004763319352512D+0
       B=0.3929160876166931D+0
       V=0.4667274455712508D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3352842634946949D+0
       B=0.4202563457288019D-1
       V=0.4069360518020356D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3891971629814670D+0
       B=0.8614309758870850D-1
       V=0.4260442819919195D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4409875565542281D+0
       B=0.1314500879380001D+0
       V=0.4408678508029063D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4904893058592484D+0
       B=0.1772189657383859D+0
       V=0.4518748115548597D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5375056138769549D+0
       B=0.2228277110050294D+0
       V=0.4595564875375116D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5818255708669969D+0
       B=0.2677179935014386D+0
       V=0.4643988774315846D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232334858144959D+0
       B=0.3113675035544165D+0
       V=0.4668827491646946D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4489485354492058D+0
       B=0.4409162378368174D-1
       V=0.4400541823741973D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5015136875933150D+0
       B=0.8939009917748489D-1
       V=0.4514512890193797D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5511300550512623D+0
       B=0.1351806029383365D+0
       V=0.4596198627347549D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5976720409858000D+0
       B=0.1808370355053196D+0
       V=0.4648659016801781D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6409956378989354D+0
       B=0.2257852192301602D+0
       V=0.4675502017157673D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5581222330827514D+0
       B=0.4532173421637160D-1
       V=0.4598494476455523D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6074705984161695D+0
       B=0.9117488031840314D-1
       V=0.4654916955152048D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6532272537379033D+0
       B=0.1369294213140155D+0
       V=0.4684709779505137D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6594761494500487D+0
       B=0.4589901487275583D-1
       V=0.4691445539106986D-3
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END subroutine LD2354
c----------------------------------------------------------------

c----------------------------------------------------------------
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8  x1,x2,x(n),w(n)
      REAL*8  EPS
      PARAMETER (EPS=3.d-15) 

      INTEGER i,j,m
      REAL*8  p1,p2,p3,pp,xl,xm,z,z1
      real*8  pi
ccc
c EPS is the relative precision. Given the lower and upper limits 
c of integration x1 and x2, and given n, this routine returns
c arrays x(1:n) and w(1:n) of length n, containing the abscissas 
c and weights of the Gauss-Legendre n-point quadrature formula.
c
c The roots are symmetric in the interval, so we only 
c have to find half of them.
ccc

      m=(n+1)/2 
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)

      pi = acos(-1.0d0)

      do i=1,m
        z=cos(pi*(i-.25d0)/(n+.5d0))
 1      continue
        p1=1.d0
        p2=0.d0

        do j=1,n  
          p3=p2 
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do

        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
      end do

      END SUBROUTINE gauleg
c----------------------------------------------------------------

      end module mathlib
