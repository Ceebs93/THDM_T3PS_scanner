           subroutine lh_readin(param_name)
c----------------------------------------------------------------------
c Read the parameters from the lh file
c
c 1. Input values for the EW sector
c 2. Higgs mass and width
c 3. Fermion masses (pole and MSbar) and widths
c----------------------------------------------------------------------
      implicit none
c
c     parameters
c
      integer maxpara
      parameter (maxpara=400)
c
c     local
c
	  character*(*) param_name	
          integer npara,l1,l2,id
          character*132 buff
          real*8 real_value
          real*8 value(maxpara)
          integer ivalue(maxpara),n
          character*20 name(maxpara),bn
          logical block_found,done,fopened
      integer iunit,i,name_length,idum
c
c       block info
c
      character*20 block_name
c
c   Common
c
          include 'coupl.inc'
	  include 'input.inc'
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb !CKM matrix elements
      common/values/    alpha,gfermi,alfas,
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb
c
c----------
c     start
c----------
c
c     open file
c
      iunit=14
      call open_file_mdl(iunit,param_name,fopened)
           done=.false.

       n=0
           do while(.not.done)
           block_found=.false.
c
c looks for the blocks or for decay
c
      do while(.not.block_found)
       read(iunit,'(a132)',end=99,err=99) buff
c--     change to lower case
           call case_trap(buff,20)
       if(buff(1:5).eq.'block') then
         if(index(buff,"#").ne.0) l1=index(buff,"#")-1
         block_name=buff(6:min(l1,26))
         call no_spaces(block_name,name_length)
c        write(*,*) block_name(1:name_length)
         block_found=.true.
           elseif(buff(1:5).eq.'decay') then
               n=n+1
               l1=30
               if(index(buff,"#").ne.0) l1=index(buff,"#")-1 ! ignore comments
               read(buff(6:l1),*) ivalue(n),value(n)
               name(n)="decay"
       endif
      end do
c
c

      if(block_found) then
 	  do while(.true.)
          read(iunit,'(a132)',end=99,err=99) buff
     	  call case_trap(buff,20)
          if(buff(1:1).eq.'b'.or.buff(1:1).eq.'d') then
          	backspace iunit
          	exit
          endif	         	
            if(buff(1:1).ne.'#'.and.buff.ne.'') then  !if it not a comment
              n=n+1	       
              l1=30
              if(index(buff,"#").ne.0) l1=index(buff,"#")-1 ! ignore comments       
c
c  WARNING:... not all blocks have the same sintax!! You need to change it
c  depending on the block you are reading
c
             if(block_name(1:5).eq."mgckm") then
                read(buff(1:l1),*) ivalue(n),idum,value(n)
                ivalue(n)=ivalue(n)*10+idum
             elseif(block_name(1:5).eq."alpha") then
                read(buff(1:l1),*) value(n)
   	     else
                read(buff(1:l1),*) ivalue(n),value(n)
              endif  
              name(n)=block_name(1:name_length)
              write(*,"(1x,i2,2x,e16.8,1x,a)") 
     &        ivalue(n),value(n),name(n)
           	  endif
      end do ! do while in the block
      else
        done=.true.
      endif
      end do ! do while the entire file
	  

 99    continue      
	
       bn="sminputs"
       call set_it(n,ivalue,value,name,1,bn,alpha,128.9d0)
       alpha=1d0/alpha
       call set_it(n,ivalue,value,name,2,bn,gfermi,0.1166d-4)
       call set_it(n,ivalue,value,name,3,bn,alfas,0.119d0)
       call set_it(n,ivalue,value,name,4,bn,zmass,91.188d0)
       call set_it(n,ivalue,value,name,6,bn,tmass,174.3d0)
       call set_it(n,ivalue,value,name,7,bn,lmass,1.777d0)
       bn="mgyukawa"
       call set_it(n,ivalue,value,name,4,bn,mcMS,1.25d0)
       call set_it(n,ivalue,value,name,5,bn,mbMS,4.2d0)
       call set_it(n,ivalue,value,name,6,bn,mtMS,174d0)
       call set_it(n,ivalue,value,name,15,bn,mtaMS,1.777d0)
       bn="mgckm"
       call set_it(n,ivalue,value,name,11,bn,vud,1d0)
       call set_it(n,ivalue,value,name,12,bn,vus,0d0)
       call set_it(n,ivalue,value,name,13,bn,vub,0d0)
       call set_it(n,ivalue,value,name,21,bn,vcd,0d0)
       call set_it(n,ivalue,value,name,22,bn,vcs,1d0)
       call set_it(n,ivalue,value,name,23,bn,vcb,0d0)
       call set_it(n,ivalue,value,name,31,bn,vtd,0d0)
       call set_it(n,ivalue,value,name,32,bn,vts,0d0)
       call set_it(n,ivalue,value,name,33,bn,vtb,1d0)
       bn="mass"
       call set_it(n,ivalue,value,name,4,bn,cmass,1.4d0)
       call set_it(n,ivalue,value,name,5,bn,bmass,4.7d0)
       call set_it(n,ivalue,value,name,6,bn,tmass,tmass*1d0)
       call set_it(n,ivalue,value,name,15,bn,lmass,lmass*1d0)
       call set_it(n,ivalue,value,name,25,bn,hmass,120d0)
       call set_it(n,ivalue,value,name,23,bn,zmass,zmass*1d0)
       call set_it(n,ivalue,value,name,24,bn,wmass,80.419d0)
       call set_it(n,ivalue,value,name,25,bn,h1mass,100d0)
       call set_it(n,ivalue,value,name,35,bn,h2mass,100d0)
       call set_it(n,ivalue,value,name,36,bn,h3mass,100d0)
       call set_it(n,ivalue,value,name,37,bn,hcmass,100d0)
       bn="decay"
       call set_it(n,ivalue,value,name,6,bn,twidth,1.5083d0)
       call set_it(n,ivalue,value,name,25,bn,hwidth,0.0037d0)
       call set_it(n,ivalue,value,name,23,bn,zwidth,2.441d0)
       call set_it(n,ivalue,value,name,24,bn,wwidth,2.0476d0)

       call set_it(n,ivalue,value,name,25,bn,h1width,1d0)
       call set_it(n,ivalue,value,name,35,bn,h2width,1d0)
       call set_it(n,ivalue,value,name,36,bn,h3width,1d0)
       call set_it(n,ivalue,value,name,37,bn,hcwidth,1d0)
       bn='mguser'
       call set_it(n,ivalue,value,name,1,bn,resgh1ee   ,0d0)
       call set_it(n,ivalue,value,name,2,bn,imsgh1ee   ,0d0)
       call set_it(n,ivalue,value,name,3,bn,repgh1ee   ,0d0)
       call set_it(n,ivalue,value,name,4,bn,impgh1ee   ,0d0)
       call set_it(n,ivalue,value,name,5,bn,resgh2ee   ,0d0)
       call set_it(n,ivalue,value,name,6,bn,imsgh2ee   ,0d0)
       call set_it(n,ivalue,value,name,7,bn,repgh2ee   ,0d0)
       call set_it(n,ivalue,value,name,8,bn,impgh2ee   ,0d0)
       call set_it(n,ivalue,value,name,9,bn,resgh3ee   ,0d0)
       call set_it(n,ivalue,value,name,10,bn,imsgh3ee   ,0d0)
       call set_it(n,ivalue,value,name,11,bn,repgh3ee   ,0d0)
       call set_it(n,ivalue,value,name,12,bn,impgh3ee   ,0d0)
       call set_it(n,ivalue,value,name,13,bn,resgh1mumu ,0d0)
       call set_it(n,ivalue,value,name,14,bn,imsgh1mumu ,0d0)
       call set_it(n,ivalue,value,name,15,bn,repgh1mumu ,0d0)
       call set_it(n,ivalue,value,name,16,bn,impgh1mumu ,0d0)
       call set_it(n,ivalue,value,name,17,bn,resgh2mumu ,0d0)
       call set_it(n,ivalue,value,name,18,bn,imsgh2mumu ,0d0)
       call set_it(n,ivalue,value,name,19,bn,repgh2mumu ,0d0)
       call set_it(n,ivalue,value,name,20,bn,impgh2mumu ,0d0)
       call set_it(n,ivalue,value,name,21,bn,resgh3mumu ,0d0)
       call set_it(n,ivalue,value,name,22,bn,imsgh3mumu ,0d0)
       call set_it(n,ivalue,value,name,23,bn,repgh3mumu ,0d0)
       call set_it(n,ivalue,value,name,24,bn,impgh3mumu ,0d0)
       call set_it(n,ivalue,value,name,25,bn,resgh1tata ,0d0)
       call set_it(n,ivalue,value,name,26,bn,imsgh1tata ,0d0)
       call set_it(n,ivalue,value,name,27,bn,repgh1tata ,0d0)
       call set_it(n,ivalue,value,name,28,bn,impgh1tata ,0d0)
       call set_it(n,ivalue,value,name,29,bn,resgh2tata ,0d0)
       call set_it(n,ivalue,value,name,30,bn,imsgh2tata ,0d0)
       call set_it(n,ivalue,value,name,31,bn,repgh2tata ,0d0)
       call set_it(n,ivalue,value,name,32,bn,impgh2tata ,0d0)
       call set_it(n,ivalue,value,name,33,bn,resgh3tata ,0d0)
       call set_it(n,ivalue,value,name,34,bn,imsgh3tata ,0d0)
       call set_it(n,ivalue,value,name,35,bn,repgh3tata ,0d0)
       call set_it(n,ivalue,value,name,36,bn,impgh3tata ,0d0)
       call set_it(n,ivalue,value,name,37,bn,resgh1uu   ,0d0)
       call set_it(n,ivalue,value,name,38,bn,imsgh1uu   ,0d0)
       call set_it(n,ivalue,value,name,39,bn,repgh1uu   ,0d0)
       call set_it(n,ivalue,value,name,40,bn,impgh1uu   ,0d0)
       call set_it(n,ivalue,value,name,41,bn,resgh2uu   ,0d0)
       call set_it(n,ivalue,value,name,42,bn,imsgh2uu   ,0d0)
       call set_it(n,ivalue,value,name,43,bn,repgh2uu   ,0d0)
       call set_it(n,ivalue,value,name,44,bn,impgh2uu   ,0d0)
       call set_it(n,ivalue,value,name,45,bn,resgh3uu   ,0d0)
       call set_it(n,ivalue,value,name,46,bn,imsgh3uu   ,0d0)
       call set_it(n,ivalue,value,name,47,bn,repgh3uu   ,0d0)
       call set_it(n,ivalue,value,name,48,bn,impgh3uu   ,0d0)
       call set_it(n,ivalue,value,name,49,bn,resgh1cc   ,0d0)
       call set_it(n,ivalue,value,name,50,bn,imsgh1cc   ,0d0)
       call set_it(n,ivalue,value,name,51,bn,repgh1cc   ,0d0)
       call set_it(n,ivalue,value,name,52,bn,impgh1cc   ,0d0)
       call set_it(n,ivalue,value,name,53,bn,resgh2cc   ,0d0)
       call set_it(n,ivalue,value,name,54,bn,imsgh2cc   ,0d0)
       call set_it(n,ivalue,value,name,55,bn,repgh2cc   ,0d0)
       call set_it(n,ivalue,value,name,56,bn,impgh2cc   ,0d0)
       call set_it(n,ivalue,value,name,57,bn,resgh3cc   ,0d0)
       call set_it(n,ivalue,value,name,58,bn,imsgh3cc   ,0d0)
       call set_it(n,ivalue,value,name,59,bn,repgh3cc   ,0d0)
       call set_it(n,ivalue,value,name,60,bn,impgh3cc   ,0d0)
       call set_it(n,ivalue,value,name,61,bn,resgh1tt   ,0d0)
       call set_it(n,ivalue,value,name,62,bn,imsgh1tt   ,0d0)
       call set_it(n,ivalue,value,name,63,bn,repgh1tt   ,0d0)
       call set_it(n,ivalue,value,name,64,bn,impgh1tt   ,0d0)
       call set_it(n,ivalue,value,name,65,bn,resgh2tt   ,0d0)
       call set_it(n,ivalue,value,name,66,bn,imsgh2tt   ,0d0)
       call set_it(n,ivalue,value,name,67,bn,repgh2tt   ,0d0)
       call set_it(n,ivalue,value,name,68,bn,impgh2tt   ,0d0)
       call set_it(n,ivalue,value,name,69,bn,resgh3tt   ,0d0)
       call set_it(n,ivalue,value,name,70,bn,imsgh3tt   ,0d0)
       call set_it(n,ivalue,value,name,71,bn,repgh3tt   ,0d0)
       call set_it(n,ivalue,value,name,72,bn,impgh3tt   ,0d0)
       call set_it(n,ivalue,value,name,73,bn,resgh1dd   ,0d0)
       call set_it(n,ivalue,value,name,74,bn,imsgh1dd   ,0d0)
       call set_it(n,ivalue,value,name,75,bn,repgh1dd   ,0d0)
       call set_it(n,ivalue,value,name,76,bn,impgh1dd   ,0d0)
       call set_it(n,ivalue,value,name,77,bn,resgh2dd   ,0d0)
       call set_it(n,ivalue,value,name,78,bn,imsgh2dd   ,0d0)
       call set_it(n,ivalue,value,name,79,bn,repgh2dd   ,0d0)
       call set_it(n,ivalue,value,name,80,bn,impgh2dd   ,0d0)
       call set_it(n,ivalue,value,name,81,bn,resgh3dd   ,0d0)
       call set_it(n,ivalue,value,name,82,bn,imsgh3dd   ,0d0)
       call set_it(n,ivalue,value,name,83,bn,repgh3dd   ,0d0)
       call set_it(n,ivalue,value,name,84,bn,impgh3dd   ,0d0)
       call set_it(n,ivalue,value,name,85,bn,resgh1ss   ,0d0)
       call set_it(n,ivalue,value,name,86,bn,imsgh1ss   ,0d0)
       call set_it(n,ivalue,value,name,87,bn,repgh1ss   ,0d0)
       call set_it(n,ivalue,value,name,88,bn,impgh1ss   ,0d0)
       call set_it(n,ivalue,value,name,89,bn,resgh2ss   ,0d0)
       call set_it(n,ivalue,value,name,90,bn,imsgh2ss   ,0d0)
       call set_it(n,ivalue,value,name,91,bn,repgh2ss   ,0d0)
       call set_it(n,ivalue,value,name,92,bn,impgh2ss   ,0d0)
       call set_it(n,ivalue,value,name,93,bn,resgh3ss   ,0d0)
       call set_it(n,ivalue,value,name,94,bn,imsgh3ss   ,0d0)
       call set_it(n,ivalue,value,name,95,bn,repgh3ss   ,0d0)
       call set_it(n,ivalue,value,name,96,bn,impgh3ss   ,0d0)
       call set_it(n,ivalue,value,name,97,bn,resgh1bb   ,0d0)
       call set_it(n,ivalue,value,name,98,bn,imsgh1bb   ,0d0)
       call set_it(n,ivalue,value,name,99,bn,repgh1bb   ,0d0)
       call set_it(n,ivalue,value,name,100,bn,impgh1bb   ,0d0)
       call set_it(n,ivalue,value,name,101,bn,resgh2bb   ,0d0)
       call set_it(n,ivalue,value,name,102,bn,imsgh2bb   ,0d0)
       call set_it(n,ivalue,value,name,103,bn,repgh2bb   ,0d0)
       call set_it(n,ivalue,value,name,104,bn,impgh2bb   ,0d0)
       call set_it(n,ivalue,value,name,105,bn,resgh3bb   ,0d0)
       call set_it(n,ivalue,value,name,106,bn,imsgh3bb   ,0d0)
       call set_it(n,ivalue,value,name,107,bn,repgh3bb   ,0d0)
       call set_it(n,ivalue,value,name,108,bn,impgh3bb   ,0d0)
       call set_it(n,ivalue,value,name,109,bn,resghcud   ,0d0)
       call set_it(n,ivalue,value,name,110,bn,imsghcud   ,0d0)
       call set_it(n,ivalue,value,name,111,bn,repghcud   ,0d0)
       call set_it(n,ivalue,value,name,112,bn,impghcud   ,0d0)
       call set_it(n,ivalue,value,name,113,bn,resghcus   ,0d0)
       call set_it(n,ivalue,value,name,114,bn,imsghcus   ,0d0)
       call set_it(n,ivalue,value,name,115,bn,repghcus   ,0d0)
       call set_it(n,ivalue,value,name,116,bn,impghcus   ,0d0)
       call set_it(n,ivalue,value,name,117,bn,resghcub   ,0d0)
       call set_it(n,ivalue,value,name,118,bn,imsghcub   ,0d0)
       call set_it(n,ivalue,value,name,119,bn,repghcub   ,0d0)
       call set_it(n,ivalue,value,name,120,bn,impghcub   ,0d0)
       call set_it(n,ivalue,value,name,121,bn,resghccd   ,0d0)
       call set_it(n,ivalue,value,name,122,bn,imsghccd   ,0d0)
       call set_it(n,ivalue,value,name,123,bn,repghccd   ,0d0)
       call set_it(n,ivalue,value,name,124,bn,impghccd   ,0d0)
       call set_it(n,ivalue,value,name,125,bn,resghccs   ,0d0)
       call set_it(n,ivalue,value,name,126,bn,imsghccs   ,0d0)
       call set_it(n,ivalue,value,name,127,bn,repghccs   ,0d0)
       call set_it(n,ivalue,value,name,128,bn,impghccs   ,0d0)
       call set_it(n,ivalue,value,name,129,bn,resghccb   ,0d0)
       call set_it(n,ivalue,value,name,130,bn,imsghccb   ,0d0)
       call set_it(n,ivalue,value,name,131,bn,repghccb   ,0d0)
       call set_it(n,ivalue,value,name,132,bn,impghccb   ,0d0)
       call set_it(n,ivalue,value,name,133,bn,resghctd   ,0d0)
       call set_it(n,ivalue,value,name,134,bn,imsghctd   ,0d0)
       call set_it(n,ivalue,value,name,135,bn,repghctd   ,0d0)
       call set_it(n,ivalue,value,name,136,bn,impghctd   ,0d0)
       call set_it(n,ivalue,value,name,137,bn,resghcts   ,0d0)
       call set_it(n,ivalue,value,name,138,bn,imsghcts   ,0d0)
       call set_it(n,ivalue,value,name,139,bn,repghcts   ,0d0)
       call set_it(n,ivalue,value,name,140,bn,impghcts   ,0d0)
       call set_it(n,ivalue,value,name,141,bn,resghctb   ,0d0)
       call set_it(n,ivalue,value,name,142,bn,imsghctb   ,0d0)
       call set_it(n,ivalue,value,name,143,bn,repghctb   ,0d0)
       call set_it(n,ivalue,value,name,144,bn,impghctb   ,0d0)
       call set_it(n,ivalue,value,name,145,bn,resghcvee  ,0d0)
       call set_it(n,ivalue,value,name,146,bn,imsghcvee  ,0d0)
       call set_it(n,ivalue,value,name,147,bn,repghcvee  ,0d0)
       call set_it(n,ivalue,value,name,148,bn,impghcvee  ,0d0)
       call set_it(n,ivalue,value,name,149,bn,resghcvmmu ,0d0)
       call set_it(n,ivalue,value,name,150,bn,imsghcvmmu ,0d0)
       call set_it(n,ivalue,value,name,151,bn,repghcvmmu ,0d0)
       call set_it(n,ivalue,value,name,152,bn,impghcvmmu ,0d0)
       call set_it(n,ivalue,value,name,153,bn,resghcvtta ,0d0)
       call set_it(n,ivalue,value,name,154,bn,imsghcvtta ,0d0)
       call set_it(n,ivalue,value,name,155,bn,repghcvtta ,0d0)
       call set_it(n,ivalue,value,name,156,bn,impghcvtta ,0d0)
       call set_it(n,ivalue,value,name,157,bn,regwwh1     ,0d0)
       call set_it(n,ivalue,value,name,158,bn,imgwwh1     ,0d0)
       call set_it(n,ivalue,value,name,159,bn,regwwh2     ,0d0)
       call set_it(n,ivalue,value,name,160,bn,imgwwh2     ,0d0)
       call set_it(n,ivalue,value,name,161,bn,regzzh1     ,0d0)
       call set_it(n,ivalue,value,name,162,bn,imgzzh1     ,0d0)
       call set_it(n,ivalue,value,name,163,bn,regzzh2     ,0d0)
       call set_it(n,ivalue,value,name,164,bn,imgzzh2     ,0d0)
       call set_it(n,ivalue,value,name,165,bn,regahchc    ,0d0)
       call set_it(n,ivalue,value,name,166,bn,imgahchc    ,0d0)
       call set_it(n,ivalue,value,name,167,bn,regzh1h2    ,0d0)
       call set_it(n,ivalue,value,name,168,bn,imgzh1h2    ,0d0)
       call set_it(n,ivalue,value,name,169,bn,regzh1h3    ,0d0)
       call set_it(n,ivalue,value,name,170,bn,imgzh1h3    ,0d0)
       call set_it(n,ivalue,value,name,171,bn,regzh2h3    ,0d0)
       call set_it(n,ivalue,value,name,172,bn,imgzh2h3    ,0d0)
       call set_it(n,ivalue,value,name,173,bn,regzhchc    ,0d0)
       call set_it(n,ivalue,value,name,174,bn,imgzhchc    ,0d0)
       call set_it(n,ivalue,value,name,175,bn,regwphch1   ,0d0)
       call set_it(n,ivalue,value,name,176,bn,imgwphch1   ,0d0)
       call set_it(n,ivalue,value,name,177,bn,regwphch2   ,0d0)
       call set_it(n,ivalue,value,name,178,bn,imgwphch2   ,0d0)
       call set_it(n,ivalue,value,name,179,bn,regwphch3   ,0d0)
       call set_it(n,ivalue,value,name,180,bn,imgwphch3   ,0d0)
       call set_it(n,ivalue,value,name,181,bn,regh1h1h1   ,0d0)
       call set_it(n,ivalue,value,name,182,bn,imgh1h1h1   ,0d0)
       call set_it(n,ivalue,value,name,183,bn,regh1h1h2   ,0d0)
       call set_it(n,ivalue,value,name,184,bn,imgh1h1h2   ,0d0)
       call set_it(n,ivalue,value,name,185,bn,regh1h2h2   ,0d0)
       call set_it(n,ivalue,value,name,186,bn,imgh1h2h2   ,0d0)
       call set_it(n,ivalue,value,name,187,bn,regh1h3h3   ,0d0)
       call set_it(n,ivalue,value,name,188,bn,imgh1h3h3   ,0d0)
       call set_it(n,ivalue,value,name,189,bn,regh2h2h2   ,0d0)
       call set_it(n,ivalue,value,name,190,bn,imgh2h2h2   ,0d0)
       call set_it(n,ivalue,value,name,191,bn,regh2h3h3   ,0d0)
       call set_it(n,ivalue,value,name,192,bn,imgh2h3h3   ,0d0)
       call set_it(n,ivalue,value,name,193,bn,regh1hchc   ,0d0)
       call set_it(n,ivalue,value,name,194,bn,imgh1hchc   ,0d0)
       call set_it(n,ivalue,value,name,195,bn,regh2hchc   ,0d0)
       call set_it(n,ivalue,value,name,196,bn,imgh2hchc   ,0d0)
       call set_it(n,ivalue,value,name,197,bn,regaahchc   ,0d0)
       call set_it(n,ivalue,value,name,198,bn,imgaahchc   ,0d0)
       call set_it(n,ivalue,value,name,199,bn,regazhchc   ,0d0)
       call set_it(n,ivalue,value,name,200,bn,imgazhchc   ,0d0)
       call set_it(n,ivalue,value,name,201,bn,regawphch1  ,0d0)
       call set_it(n,ivalue,value,name,202,bn,imgawphch1  ,0d0)
       call set_it(n,ivalue,value,name,203,bn,regawphch2  ,0d0)
       call set_it(n,ivalue,value,name,204,bn,imgawphch2  ,0d0)
       call set_it(n,ivalue,value,name,205,bn,regawphch3  ,0d0)
       call set_it(n,ivalue,value,name,206,bn,imgawphch3  ,0d0)
       call set_it(n,ivalue,value,name,207,bn,regzzh1h1   ,0d0)
       call set_it(n,ivalue,value,name,208,bn,imgzzh1h1   ,0d0)
       call set_it(n,ivalue,value,name,209,bn,regzzh2h2   ,0d0)
       call set_it(n,ivalue,value,name,210,bn,imgzzh2h2   ,0d0)
       call set_it(n,ivalue,value,name,211,bn,regzzh3h3   ,0d0)
       call set_it(n,ivalue,value,name,212,bn,imgzzh3h3   ,0d0)
       call set_it(n,ivalue,value,name,213,bn,regzzhchc   ,0d0)
       call set_it(n,ivalue,value,name,214,bn,imgzzhchc   ,0d0)
       call set_it(n,ivalue,value,name,215,bn,regzwphch1  ,0d0)
       call set_it(n,ivalue,value,name,216,bn,imgzwphch1  ,0d0)
       call set_it(n,ivalue,value,name,217,bn,regzwphch2  ,0d0)
       call set_it(n,ivalue,value,name,218,bn,imgzwphch2  ,0d0)
       call set_it(n,ivalue,value,name,219,bn,regzwphch3  ,0d0)
       call set_it(n,ivalue,value,name,220,bn,imgzwphch3  ,0d0)
       call set_it(n,ivalue,value,name,221,bn,regwwh1h1   ,0d0)
       call set_it(n,ivalue,value,name,222,bn,imgwwh1h1   ,0d0)
       call set_it(n,ivalue,value,name,223,bn,regwwh2h2   ,0d0)
       call set_it(n,ivalue,value,name,224,bn,imgwwh2h2   ,0d0)
       call set_it(n,ivalue,value,name,225,bn,regwwh3h3   ,0d0)
       call set_it(n,ivalue,value,name,226,bn,imgwwh3h3   ,0d0)
       call set_it(n,ivalue,value,name,227,bn,regwwhchc   ,0d0)
       call set_it(n,ivalue,value,name,228,bn,imgwwhchc   ,0d0)
      return
      end

      
      subroutine set_it(npara,ivalue,value,name,id,
     &                  block_name,var,def_value)
c----------------------------------------------------------------------------------
c     finds the parameter value  in block_name and associate var to it.
c     If it is not found a default is given.
c----------------------------------------------------------------------------------
      implicit none

c
c     parameters
c
      integer maxpara
      parameter (maxpara=400)
c
c     arguments
c
      integer npara,ivalue(maxpara),id
      character*20  block_name,name(maxpara)
      real*8 var,def_value,value(maxpara)
c
c     local
c
      logical found
      integer i
c
c     start
c
	  found=.false.
      do i=1,npara
         found = (id.eq.ivalue(i)).and.(name(i).eq.block_name)
 	               if(found) then
         	var=value(i)
            exit
          endif	
      enddo
      
      if (.not.found) then
c         write (*,*) "Warning: parameter not found"
c         write (*,*) "         setting it to default value ",def_value
         var=def_value
      endif
      return

      end
      
      
      subroutine case_trap(string,length)
c**********************************************************    
c change string to lowercase if the input is not
c**********************************************************
      implicit none
c
c     ARGUMENT
c      
      character*(*) string
      integer length
c
c     LOCAL
c
      integer i,k

      do i=1,length
         k=ichar(string(i:i))
         if(k.ge.65.and.k.le.90) then  !upper case A-Z
            k=ichar(string(i:i))+32   
            string(i:i)=char(k)        
         endif
      enddo

      return
      end

! c********************************************************************
      subroutine open_file_mdl(lun,filename,fopened)
c***********************************************************************
c     opens file input-card.dat in current directory or above
c***********************************************************************
      implicit none
c
c     Arguments
c
      integer lun
      logical fopened
      character*(*) filename
      character*30  tempname
      integer fine
      integer dirup

c-----
c  Begin Code
c-----
c
c     first check that we will end in the main directory
c
      open(unit=lun,file=filename,status='old',err=10)
      fopened=.true.
      return      
 10   close(lun)

      open(unit=lun,file="Source/makefile",status='old',err=20)
      dirup=0
      goto 100
 20   close(lun)

      open(unit=lun,file="../Source/makefile",status='old',err=30)
      dirup=1
      goto 100
 30   close(lun)

      open(unit=lun,file="../../Source/makefile",status='old',err=40)
      dirup=2
      goto 100
 40   close(lun)

      open(unit=lun,file="../../../Source/makefile",status='old',err=50)
      dirup=3
      goto 100
 50   close(lun)

      open(unit=lun,file="../../../../Source/makefile",status='old',err=60)
      dirup=4
      goto 100
 60   close(lun)

 100  continue
      close(lun)

      fopened=.true.
      tempname=filename
      fine=index(tempname,' ')
      if(fine.eq.0) fine=len(tempname)
c
c         if I have to read a card
c
          if(index(filename,"_card").gt.0) then
             tempname='/Cards/'//tempname(1:fine)
             fine=fine+7
      endif

      if(dirup.eq.0) open(unit=lun,file=tempname(1:fine),status='old',err=110)
      if(dirup.eq.1) open(unit=lun,file='../'//tempname(1:fine),status='old',err=110)
      if(dirup.eq.2) open(unit=lun,file='../../'//tempname(1:fine),status='old',err=110)
      if(dirup.eq.3) open(unit=lun,file='../../../'//tempname(1:fine),status='old',err=110)
      if(dirup.eq.4) open(unit=lun,file='../../../../'//tempname(1:fine),status='old',err=110)
      return

 110  fopened=.false.
      close(lun)
      write (*,*) 'Warning: file ',tempname(1:fine),' is not in the main directory'

      return
      end

      subroutine no_spaces(buff,nchars)
c**********************************************************************
c     Given buff a buffer of words separated by spaces
c     returns it where all space are moved to the right
c     returns also the length of the single word.
c     maxlength is the length of the buffer
c**********************************************************************
      implicit none
c
c     Constants
c
      integer    maxline
      parameter (maxline=20)
      character*1 null
      parameter  (null=' ')
c
c     Arguments
c
      character*(maxline) buff
      integer nchars,maxlength
c
c     Local
c
      integer i,j
      character*(maxline) temp
c-----
c  Begin Code
c-----
      nchars=0
c      write (*,*) "buff=",buff(1:maxlength)
      do i=1,maxline
         if(buff(i:i).ne.null) then
            nchars=nchars+1
            temp(nchars:nchars)=buff(i:i)
         endif
c         write(*,*) i,":",buff(1:maxlength),":",temp(1:nchars),":"
      enddo
      buff=temp      
      end

      
