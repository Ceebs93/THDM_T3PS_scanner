C-{{{ function sgg3(yy)

      real*8 function sgg3(xt)
C
      implicit none
      include '../commons/common-expand.f'
      real*8 sgg3exact,sgg3exp,xt
C
      if (nexpand(3).eq.0) then
         sgg3 = sgg3exact(xt)
      else
         sgg3 = sgg3exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function sqg3(yy)

      real*8 function sqg3(xt)
C
      implicit none
      include '../commons/common-expand.f'
      real*8 sqg3exact,sqg3exp,xt
C
      if (nexpand(3).eq.0) then
         sqg3 = sqg3exact(xt)
      else
         sqg3 = sqg3exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function sqq3(yy)

      real*8 function sqq3(xt)
C
      implicit none
      include '../commons/common-expand.f'
      real*8 sqq3exact,sqq3exp,xt
C
      if (nexpand(3).eq.0) then
         sqq3 = sqq3exact(xt)
      else
         sqq3 = sqq3exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function squ3(yy)

      real*8 function squ3(xt)
C
      implicit none
      include '../commons/common-expand.f'
      real*8 squ3exact,squ3exp,xt
C
      if (nexpand(3).eq.0) then
         squ3 = squ3exact(xt)
      else
         squ3 = squ3exp(xt)
      endif

      return
      end

C-}}}
C-{{{ function sqqb3(yy)

      real*8 function sqqb3(xt)
C
      implicit none
      include '../commons/common-expand.f'
      real*8 sqqb3exact,sqqb3exp,xt
C
      if (nexpand(3).eq.0) then
         sqqb3 = sqqb3exact(xt)
      else
         sqqb3 = sqqb3exp(xt)
      endif

      return
      end

C-}}}

C-{{{ function sgg3exact(yy)

      real*8 function sgg3exact(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sgg3exact not defined yet.')

      return
      end

C-}}}
C-{{{ function sqg3exact(yy)

      real*8 function sqg3exact(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sqg3exact not defined yet.')

      return
      end

C-}}}
C-{{{ function sqq3exact(yy)

      real*8 function sqq3exact(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sqq3exact not defined yet.')

      return
      end

C-}}}
C-{{{ function squ3exact(yy)

      real*8 function squ3exact(xt)
C
      implicit none
      real*8 xt

      call printdieggh('squ3exact not defined yet.')

      return
      end

C-}}}
C-{{{ function sqqb3exact(yy)

      real*8 function sqqb3exact(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sqqb3exact not defined yet.')

      return
      end

C-}}}


C-{{{ function sgg3mt(yy)

      real*8 function sgg3mt(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sgg3mt not defined yet.')

      return
      end

C-}}}
C-{{{ function sqg3mt(yy)

      real*8 function sqg3mt(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sqg3mt not defined yet.')

      return
      end

C-}}}
C-{{{ function sqq3mt(yy)

      real*8 function sqq3mt(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sqq3mt not defined yet.')

      return
      end

C-}}}
C-{{{ function squ3mt(yy)

      real*8 function squ3mt(xt)
C
      implicit none
      real*8 xt

      call printdieggh('squ3mt not defined yet.')

      return
      end

C-}}}
C-{{{ function sqqb3mt(yy)

      real*8 function sqqb3mt(xt)
C
      implicit none
      real*8 xt

      call printdieggh('sqqb3mt not defined yet.')

      return
      end

C-}}}


