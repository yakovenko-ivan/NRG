!*DECK DDSCL
      SUBROUTINE DDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
!C***BEGIN PROLOGUE  DDSCL
!C***SUBSIDIARY
!C***PURPOSE  Subroutine DDSCL rescales the YH array whenever the step
!C            size is changed.
!C***LIBRARY   SLATEC (SDRIVE)
!C***TYPE      DOUBLE PRECISION (SDSCL-S, DDSCL-D, CDSCL-C)
!C***AUTHOR  Kahaner, D. K., (NIST)
!C             National Institute of Standards and Technology
!C             Gaithersburg, MD  20899
!C           Sutherland, C. D., (LANL)
!C             Mail Stop D466
!C             Los Alamos National Laboratory
!C             Los Alamos, NM  87545
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   790601  DATE WRITTEN
!C   900329  Initial submission to SLATEC.
!C***END PROLOGUE  DDSCL
      INTEGER I, J, N, NQ
      DOUBLE PRECISION H, HMAX, RC, RH, RMAX, R1, YH(N,*)
!C***FIRST EXECUTABLE STATEMENT  DDSCL
      IF (H .LT. 1.D0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.D0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
      RETURN
      END
