! This file is part of SusHi.
! 
! It includes a conversion functions for strings.
! 
C
      CHARACTER*10 FUNCTION TOSTRING(X)
C
C  Converts an arbitrary integer with at most a length of 10
C  to a string.
C
C  Example for usage:
C
C    FILE = 'file.'//TOSTRING(29389)
C
C  produces the string file.29389
C
C  written by rh, January 1997
C
      CHARACTER*10 STR,DUMMY
      INTEGER X,XLOC,X1,LENG

      IF (X*10.E-10.GT.1.) THEN
         write(6,*) 'Error in TOSTRING. Type <return> to continue.'
         read(*,*)
         WRITE(6,*) 'Number to be converted?'
         READ(5,900) STR
         GOTO 102
      ENDIF

      LENG = 10
      STR = ''
      
      XLOC = X
      DO 101 II=1,LENG
         X1 = MOD(XLOC,10)
         DUMMY = CHAR(48+X1)
         STR = DUMMY(1:1)//STR(1:II)
         IF (XLOC.LT.10) GOTO 102
         XLOC = (XLOC-X1)/10
 101  CONTINUE
 102  CONTINUE
      TOSTRING = STR
 900  FORMAT(A10)
      END
