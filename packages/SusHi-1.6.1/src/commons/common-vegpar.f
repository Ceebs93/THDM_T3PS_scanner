
      real*8 acc
      integer itmx1sushi,itmx2sushi,ncall1sushi,ncall2sushi,nprnvsushi
      integer itmx1ggh,itmx2ggh,ncall1ggh,ncall2ggh,nprnvggh
      integer itmx1bbh,itmx2bbh,ncall1bbh,ncall2bbh,nprnvbbh
      logical lveg1sushi,lveg1ggh,lveg1bbh
      common/vegpar/ acc
      common/vegparsushi/ itmx1sushi,itmx2sushi,ncall1sushi,ncall2sushi
     &     ,nprnvsushi,lveg1sushi
      common/vegparbbh/ itmx1bbh,itmx2bbh,ncall1bbh,ncall2bbh
     &     ,nprnvbbh,lveg1bbh
      common/vegparggh/ itmx1ggh,itmx2ggh,ncall1ggh,ncall2ggh
     &     ,nprnvggh,lveg1ggh
