      SUBROUTINE GLine(Line,At,Carts)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 Carts(3)
      CHARACTER*(*) Line
      CHARACTER*30  Tmp
      CHARACTER*2   At
      L=LEN(Line)
      DO I=1,L
         IF(Line(I:I).NE.' ')THEN
            J=I
            GOTO 101
         ENDIF
      ENDDO
      WRITE(*,*)' Line = ',Line
      STOP ' All blanks in GLine! '
 101  CONTINUE
      At=Line(J:J+1)
      J=J+2
      DO 200 N=1,3
         DO I=J,L
            IF(Line(I:I).NE.' ')THEN
               K1=I
               GOTO 201
            ENDIF
         ENDDO
 201     CONTINUE
         DO I=K1,L
            IF(Line(I:I).EQ.' ')THEN
               K2=I
               GOTO 202
            ENDIF
         ENDDO
 202     CONTINUE
         Tmp=Line(K1:K2)
         CALL Squish(Tmp)
         READ(Tmp,'(D26.16)')Carts(N)
         J=K2+1
 200  CONTINUE
      RETURN
      END

      SUBROUTINE Squish(String)
      CHARACTER*(*) String
      L=LEN(String)
      J=0
      DO I=1,L
         IF (String(I:I).NE.' ') THEN
             J=J+1
             String(J:J)=String(I:I)
         ENDIF
      ENDDO
      J=J+1
      IF(J.LE.L) String(J:L)=' '
      RETURN
      END

