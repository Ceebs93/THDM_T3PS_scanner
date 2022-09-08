ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_67 = (MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_69 = (MDL_COSBMA*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
      GC_87 = (MDL_EE*MDL_COMPLEXI*MDL_MW*MDL_SINBMA)/MDL_SW
      GC_127 = -((MDL_COMPLEXI*MDL_SINBMA*MDL_YB)/MDL_SQRT__2)
      GC_141 = (MDL_EE*MDL_COMPLEXI*MDL_YMC)/(MDL_MW*MDL_SW*MDL_TB
     $ *MDL_SQRT__2)
      GC_142 = -(MDL_COSBMA*MDL_EE*MDL_COMPLEXI*MDL_YMC)/(2.000000D+00
     $ *MDL_MW*MDL_SW*MDL_TB)-(MDL_COMPLEXI*MDL_SINBMA*MDL_YC)
     $ /MDL_SQRT__2
      GC_145 = -((MDL_EE*MDL_COMPLEXI*MDL_YMDO)/(MDL_MW*MDL_SW*MDL_TB
     $ *MDL_SQRT__2))
      GC_146 = -(MDL_COSBMA*MDL_EE*MDL_COMPLEXI*MDL_YMDO)/(2.000000D
     $ +00*MDL_MW*MDL_SW*MDL_TB)-(MDL_COMPLEXI*MDL_SINBMA*MDL_YDO)
     $ /MDL_SQRT__2
      GC_157 = -((MDL_EE*MDL_COMPLEXI*MDL_YMS)/(MDL_MW*MDL_SW*MDL_TB
     $ *MDL_SQRT__2))
      GC_165 = (MDL_EE*MDL_COMPLEXI*MDL_YMUP)/(MDL_MW*MDL_SW*MDL_TB
     $ *MDL_SQRT__2)
      GC_169 = -(MDL_COSBMA*MDL_EE*MDL_COMPLEXI*MDL_YMS)/(2.000000D+00
     $ *MDL_MW*MDL_SW*MDL_TB)-(MDL_COMPLEXI*MDL_SINBMA*MDL_YS)
     $ /MDL_SQRT__2
      GC_173 = -((MDL_COMPLEXI*MDL_SINBMA*MDL_YT)/MDL_SQRT__2)
      GC_181 = -(MDL_COSBMA*MDL_EE*MDL_COMPLEXI*MDL_YMUP)/(2.000000D
     $ +00*MDL_MW*MDL_SW*MDL_TB)-(MDL_COMPLEXI*MDL_SINBMA*MDL_YUP)
     $ /MDL_SQRT__2
      END
