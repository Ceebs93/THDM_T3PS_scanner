      double precision z,muF,muFfacggh,mufggh,mh2,minpt,maxpt,
     &minrap,maxrap,pt,y
      integer iset,dist,SUiset(4),pdforder
      character SUpdfs(4)*50
      logical ptcut,rapcut,pseudorap,pty,subtr,juser,scalespt
      
      common /integration/ z,muF,muFfacggh,mufggh,
     &     minpt,maxpt,minrap,maxrap,pt,y,iset,dist,
     &     SUiset,pdforder,ptcut,rapcut,pseudorap,pty,subtr,
     &     juser,scalespt,SUpdfs

      common /mh2_sushi/ mh2
