! This file is part of SusHi.
! 
! The routines of this file are used by the input routine
! to be found in inputoutput.f .
! 
C-{{{ RCS:

c..
c..   $Id: readslha.f,v 1.9 2004/03/17 17:31:15 rharland Exp $
c..   $Log: readslha.f,v $
c..   Revision 1.9  2004/03/17 17:31:15  rharland
c..   Input strings are now treated case sensitive (not keywords like BLOCK
c..   etc)
c..
c..   Revision 1.8  2004/02/19 18:18:18  rharland
c..   bug in readblocks fixed.
c..
c..   Revision 1.7  2004/02/19 17:58:07  rharland
c..   *** empty log message ***
c..
c..   Revision 1.6  2004/02/19 00:00:21  rharland
c..   BLOCK CREIN implemented.
c..
c..   Revision 1.5  2004/02/10 18:50:32  rharland
c..   small change that makes g77 -pedantic run through.
c..
c..   Revision 1.4  2004/02/07 15:08:47  rharland
c..   *** empty log message ***
c..
c..   Revision 1.3  2004/02/07 15:07:44  rharland
c..   READSLHA got additional argument IFOUND.
c..
c..   Revision 1.2  2004/02/06 15:55:18  rharland
c..   Id and Log introduced.
c..
c..

C-}}}
C----------------------------------------------------------------------
c..
c..   readslha.f
c..
c..   needs 
c..   common-slha.f
c..
c..   NOTE: maximum length of input lines is currently 
c..         set to 200 characters.
c..
c..   To add a block:
c..   (1) add an array for its value in common-slha.f: NEWARRAY(100)
c..   (2) possibly add an array for the names of its entries in 
c..       common-slha.f:  CNEWARRAY(100)
c..   (3) add a section in subroutine READDATA that fills NEWARRAY.
c..   (4) add a section in subroutine SLHAMAP that fills CNEWARRAY.
c..
C----------------------------------------------------------------------

C-{{{ subroutine readblocks:

      subroutine readblocks(iunit,blocks)
c..
c..   This is just an interface to readslha to read multiple blocks.
c..   
c..   Example: 
c..   BLOCKS = ('mass','sminputs','minpar')
c..
      implicit real*8(a-h,o-z)
      character blocks(20)*15
      character blocktype*15
      
c..   if blocks(i) is shorter than *15, this will fill up the string
c..   with whitespace:
      blocktype = '               '

      i=1
      do while (blocks(i).ne.' ')
         blocktype = blocks(i)
         call readslha(iunit,blocktype,ifound)
         if (ifound.eq.0) then
            write(6,*) 'No Block ',blocktype,' found.'
         endif
         i=i+1
      enddo

      end

C-}}}
C-{{{ subroutine readslha:

      subroutine readslha(iunit,blocktypin,ifound)
c..
c..   Read block BLOCKTYPIN from SLHA input file
c..   and fill the corresponding COMMON block.
c..   SLHA input file must be open in unit IUNIT.
c..
c..   Example:
c..   BLOCKTYPIN = 'mass'
c..   
c..   If BLOCK is found in input file, IFOUND is set to 1.
c..   IFOUND = 0 otherwise.
c..   
      integer klen
      character cline*200,clineout*200,emptyline*200
      character blocktype*15,blocktypin*15
      character uppercase*1

      rewind(iunit)

      write(emptyline,1002)
      cline = emptyline
      blocktype=blocktypin
      lenblock = lnblnk(blocktype)

c..   change input to upper case:
      do j=1,lenblock
         blocktype(j:j)=uppercase(blocktype(j:j))
      enddo

      lcount=0
      iread=0
      ifound=0
c..   I use a GOTO loop, because I don't know any other way
c..   how to read to the end of a file without knowing its length.
 100  lcount = lcount+1
      read(iunit,1001,END=101) cline
      call normalline(cline,clineout)
c..   
c..   cline:    original input with comments and multiple blanks removed
c..   clineout: uppercase of cline
c..   
      cline = clineout
      do i=1,200
         clineout(i:i) = uppercase(clineout(i:i))
      enddo
c..   check if line is valid input line:
      if ((clineout(1:1).ne.' ').and.(clineout(1:1).ne.'#').and.
     &     (clineout(1:1).ne.'B').and.(clineout(1:1).ne.'D')) then
         write(6,*) '*** Error reading SLHA input file (UNIT=',iunit
     &        ,'):'
         write(6,*) 'Line ',lcount,': ',
     &        clineout(1:1),' = char(',
     &        ichar(clineout(1:1)),') not allowed in first column'
         stop
      endif
      if (clineout.eq.emptyline) goto 100
c..   If we are looking at a keyword (BLOCK, DECAY), then either
c..   - start reading data if it's the correct BLOCK
c..   - stop reading if we are already reading data
      if ((clineout(1:5).eq.'BLOCK').or.(clineout(1:5).eq.'DECAY')) then
         klen = 7
         do while (clineout(klen:klen).ne.' ') 
            klen=klen+1
         enddo
         klen=klen-1
         iread = 0
         if ((klen.eq.lenblock+6).and.
     &        (clineout(1:klen).eq.'BLOCK '//blocktype(1:lenblock)))
     &        then
            ifound = 1
            iread = 1
            goto 100
         endif
      endif
      if (iread.eq.1) then
         call slhablocks(blocktype,cline,ierr)
         if (ierr.eq.1) then
            write(6,*) '*** Error reading SLHA input file (UNIT
     &           =',iunit,'):'
            write(6,*) 'Line ',lcount,':>',cline(1:60)
            stop
         endif
      endif
      goto 100

 1001 format(200a)
 1002 format(200(' '))
 101  continue
      end

C-}}}
C-{{{ subroutine normalline:

      subroutine normalline(clinein,clineout)
c..
c..   Remove multiple blanks, remove comments,
c..   and change everything to upper case
c..
c..   Example:
c..   clinein = '  some    string with   tabs and  #  comments'
c..   
c..   Returns:
c..   clineout = ' some string with tabs and     '
c..   
c..   The output string is always 200 chars long (trailing blanks)
c..   
      implicit real*8 (a-h,o-z)
      character cline*200,clinein*200,clineinv*200,clineout*200
      character uppercase*1
      logical nonempty

      cline = clinein

c..   remove comments, change tabs to spaces, and change to upper case:
      do i=1,200
         if (cline(i:i).eq.char(9)) cline(i:i) = ' '
         if (cline(i:i).eq.char(13)) cline(i:i) = ' '
         if (cline(i:i).eq.'#') then
            do j=i,200
               cline(j:j) = ' '
            enddo
         endif
      enddo

c..   remove multiple whitespace:
      do i=1,200
         clineinv(i:i) = cline(200+1-i:200+1-i)
      enddo
      clineout=' '
      i=1
      nonempty = .false.
      do while (i.le.200)
         do while ((clineinv(i:i).ne.' ').and.(i.le.200))
            nonempty = .true.
            clineout = clineinv(i:i)//clineout
            i=i+1
         enddo
         if (i.le.200) then
            clineout = ' '//clineout
            i=i+1
         endif
         do while ((clineinv(i:i).eq.' ').and.(i.le.200))
            i=i+1
         enddo
      enddo

      end

C-}}}
C-{{{ function uppercase:

      function uppercase(str)

      character*1 str,uppercase

      do ich=97,122
         if (char(ich).eq.str) then
            uppercase=char(ich-32)
            return
         else
            uppercase=str
         endif
      enddo

      end

C-}}}
