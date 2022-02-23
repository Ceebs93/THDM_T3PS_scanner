
!******************************************************
program HBSLHAinputblocksfromFH
!
! Run with
!   ./HBSLHAinputblocksfromFH <path to SLHA file>
! where <path to SLHA file> is the path to the SLHA file providing the
! input for FeynHiggs
!
! The program will create a new file, called
!    <path to SLHA file>.fh
! which will contain the blocks from the original SLHA file and
! the results from FeynHiggs to the SLHA file
! in the standard SLHA format, but also the blocks
!     HiggsBoundsInputHiggsCouplingsFermions
!       and
!     HiggsBoundsInputHiggsCouplingsBosons
! which are defined in the HiggsBounds manual.
!
!******************************************************
    use extra_bits_for_SLHA, only: addcouplingsblocktoSLHAfile
    use usefulbits, only: couplratio, np, Hneut, Hplus, allocate_couplratio_parts
#ifdef NAGf90Fortran
    use F90_UNIX_ENV, only: iargc, getarg
#endif

    implicit none

    character(len=300) :: infile, outfile, outfile2
    type(couplratio) :: effC(1)

    character(LEN=300) :: temp
    integer :: number_args
#ifndef NAGf90Fortran
    integer :: iargc
#endif

    integer :: i

    number_args = IARGC()

    if (number_args .ne. 1) then
        write (*, *) '!******************************************************'
        write (*, *) ' '
        write (*, *) ' Run with'
        write (*, *) '     ./HBSLHAinputblocksfromFH <path to SLHA file>'
        write (*, *) ' where <path to SLHA file> is the path to the SLHA file providing the'
        write (*, *) ' input for FeynHiggs'
        write (*, *) ''
        write (*, *) ' The program will create a new file, called'
        write (*, *) '     <path to SLHA file>.fh'
        write (*, *) ' which will contain the blocks from the original SLHA file and'
        write (*, *) ' the results from FeynHiggs to the SLHA file '
        write (*, *) ' in the standard SLHA format, but also the blocks'
        write (*, *) '     HiggsBoundsInputHiggsCouplingsFermions'
        write (*, *) '       and'
        write (*, *) '     HiggsBoundsInputHiggsCouplingsBosons'
        write (*, *) ' which are defined in the HiggsBounds manual.'
        write (*, *) ''
        write (*, *) '!******************************************************'

        stop "Incorrect number of arguments given to HBSLHAinputblocksfromFH"
    endif

    ! Read arguments into text strings.
    i = 1
    temp = ""
    call GETARG(i, temp)
    infile = ""
    infile = trim(temp)

    outfile = trim(adjustl(infile))//'.fh'
    outfile2 = trim(adjustl(infile))//'.in'

    call system("cp "//trim(infile)//" "//trim(outfile2))

    np = 0
    np(Hneut) = 3
    np(Hplus) = 1
    call allocate_couplratio_parts(effC)

    call createSLHAfilewithFHwithoutHBinputblocks(infile, outfile, &
                                                  effC(1)%hjbb_s, effC(1)%hjbb_p, effC(1)%hjtt_s, effC(1)%hjtt_p, &
                                                  effC(1)%hjtautau_s, effC(1)%hjtautau_p, &
                                                  effC(1)%hjWW, effC(1)%hjZZ, &
                                                  effC(1)%hjgg, effC(1)%hjhiZ)

    call addcouplingsblocktoSLHAfile(outfile, effC(1))

    call addcouplingsblocktoSLHAfile(outfile2, effC(1))

end program HBSLHAinputblocksfromFH

