      SUBROUTINE LowCase(String)
C
      IMPLICIT INTEGER (a-h,o-z)
C
      CHARACTER*(*) String
C
      CHARACTER*26 Lower
      CHARACTER*26 Upper
C
      Lower='abcdefghijklmnopqrstuvwxyz'
      Upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
c
      DO I=1,LEN(String)
         J=Index(Upper,String(I:I))
         IF(J.NE.0)String(I:I)=Lower(J:J)
      ENDDO
C
      RETURN
      END

