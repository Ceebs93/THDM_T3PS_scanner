! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables_type1
!******************************************************************

    implicit none

    !table type 1-----------------------------
    type table1
        integer :: id, nx, particle_x !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
        character(LEN=45) :: label
        character(LEN=35) :: citekey
        character(LEN=100) :: desc
        character(LEN=3) :: expt
        double precision :: lumi, energy
        double precision :: xmax, xmin, sep, deltax
        integer :: SMlike
        integer :: llh
        double precision, allocatable :: dat(:, :) !in dat(a,b), a=row, b=1,2 for obs,pred
    end type
    !------------------------------------------

    integer, parameter :: file_id_1 = 10 !same as file_id_common in usefulbits.f90

contains

    !************************************************************
    subroutine initializetables_type1_blank(tablet1)
        !***********************************************************
        ! still leaves dat unallocated
        integer:: i
        type(table1) :: tablet1(:)

        do i = lbound(tablet1, dim=1), ubound(tablet1, dim=1)
            tablet1(i)%id = -1
            tablet1(i)%nx = -1
            tablet1(i)%particle_x = -1
            tablet1(i)%label = ''
            tablet1(i)%citekey = ''
            tablet1(i)%desc = ''
            tablet1(i)%expt = ''
            tablet1(i)%lumi = -1.0D0
            tablet1(i)%energy = -1.0D0
            tablet1(i)%xmax = -1.0D0
            tablet1(i)%xmin = -1.0D0
            tablet1(i)%sep = -1.0D0
            tablet1(i)%deltax = -1.0D0
            tablet1(i)%SMlike = 0
            tablet1(i)%llh = 0
        end do

    end subroutine initializetables_type1_blank

    !************************************************************
    subroutine copy_type1(tablet1_orig, tablet1_copy)
        !***********************************************************
        ! note tablet1_1,tablet1_2 are not arrays
        ! still leaves dat uncopied
        type(table1) :: tablet1_orig
        type(table1) :: tablet1_copy

        tablet1_copy%id = tablet1_orig%id
        tablet1_copy%nx = tablet1_orig%nx
        tablet1_copy%particle_x = tablet1_orig%particle_x
        tablet1_copy%label = tablet1_orig%label
        tablet1_copy%citekey = tablet1_orig%citekey
        tablet1_copy%expt = tablet1_orig%expt
        tablet1_copy%xmax = tablet1_orig%xmax
        tablet1_copy%xmin = tablet1_orig%xmin
        tablet1_copy%sep = tablet1_orig%sep
        tablet1_copy%deltax = tablet1_orig%deltax

    end subroutine copy_type1
    !***********************************************************
    function t1elementnumberfromid(t1, id)
        !--------------------------------------input
        type(table1), intent(in) :: t1(:)
        integer, intent(in) :: id
        !-----------------------------------function
        integer :: t1elementnumberfromid
        !-----------------------------------internal
        integer :: n, x
        !-------------------------------------------

        n = 0
        do x = lbound(t1, dim=1), ubound(t1, dim=1)
            if (t1(x)%id .eq. id) then
                n = n + 1
                t1elementnumberfromid = x
            end if
        end do

        if (n .ne. 1) stop 'problem in function t3elementnumberfromid 1'

    end function t1elementnumberfromid

    !************************************************************
    subroutine initializetables1(S95_t1)
        !***********************************************************
        ! fills S95_t1
        !***********************************************************
        use install_data
        use usefulbits, only: Hneut, Hplus, file_id_common2
        implicit none

        !--------------------------------------input
        type(table1) :: S95_t1(:)
        !-----------------------------------internal
        integer :: x, xbeg, xend
        character(len=100), allocatable :: filename(:)
        character(LEN=pathname_length + 150) :: fullfilename
        integer :: col
        integer :: ios
        !-------------------------------------------

        xbeg = lbound(S95_t1, dim=1)
        xend = ubound(S95_t1, dim=1)

        allocate (filename(xbeg:xend))
        x = xbeg - 1

        !instead, could read in the values of xmin,xmax,sep from the
        !files, but it's kinda nice having them all here to refer to

        x = x + 1
        S95_t1(x)%id = 142
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0602042, table 14b (LEP)' ! table 14b
        S95_t1(x)%citekey = 'Schael:2006cr'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 12.0D0
        S95_t1(x)%xmax = 120.0D0
        S95_t1(x)%sep = 0.5D0
        filename(x) = 'lep210_hbb'

        x = x + 1
        S95_t1(x)%id = 143
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0602042, table 14c (LEP)' ! table 14c
        S95_t1(x)%citekey = 'Schael:2006cr'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 120.0D0
        S95_t1(x)%sep = 0.5D0
        filename(x) = 'lep210_htt_interpol'

        x = x + 1
        S95_t1(x)%id = 300
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0206022 (OPAL)'
        S95_t1(x)%citekey = 'Abbiendi:2002qp'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 1.0D0
        S95_t1(x)%xmax = 100.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'lep_decaymodeindep'

        x = x + 1
        S95_t1(x)%id = 400
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0107032 (LEP)'
        S95_t1(x)%citekey = 'Searches:2001ab'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 118.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'LEP_h-invisible'

        x = x + 1
        S95_t1(x)%id = 500
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'LHWG Note 2002-02'
        S95_t1(x)%citekey = 'ALEPH:2002gcw'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 20.0D0
        S95_t1(x)%xmax = 116.0D0
        S95_t1(x)%sep = 2.0D0
        filename(x) = 'LEP_h-gammagamma'

        x = x + 1
        S95_t1(x)%id = 600
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'LHWG (unpublished)'!uses hep-ex/0510022,hep-ex/0205055,hep-ex/0312042,hep-ex/0408097
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 20.0D0
        S95_t1(x)%xmax = 128.6D0
        S95_t1(x)%sep = 0.1D0
        filename(x) = 'LEP_h-2jets'

        x = x + 1
        S95_t1(x)%id = 601
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0107034 (LHWG)'!uses hep-ex/0510022,hep-ex/0205055,hep-ex/0312042,hep-ex/0408097
        S95_t1(x)%citekey = 'Searches:2001aa'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 60.0D0
        S95_t1(x)%xmax = 114.5D0
        S95_t1(x)%sep = 0.1D0
        filename(x) = 'LEP_h-2jets_0107034'

        x = x + 1
        S95_t1(x)%id = 711
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0410017 (DELPHI)'
        S95_t1(x)%citekey = 'Abdallah:2004wy'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 12.0D0
        S95_t1(x)%xmax = 50.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'Delphi_yuk_h_bbbb'

        x = x + 1
        S95_t1(x)%id = 713
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0410017 (DELPHI)'
        S95_t1(x)%citekey = 'Abdallah:2004wy'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 12.0D0
        S95_t1(x)%xmax = 50.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'Delphi_yuk_a_bbbb'

        x = x + 1
        S95_t1(x)%id = 721
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0410017 (DELPHI)'
        S95_t1(x)%citekey = 'Abdallah:2004wy'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 50.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'Delphi_yuk_h_bbtautau'

        x = x + 1
        S95_t1(x)%id = 741
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0111010 (OPAL)'
        S95_t1(x)%citekey = 'Abbiendi:2001kp'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 12.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'OPAL_yuk_h_bbtautau'

        x = x + 1
        S95_t1(x)%id = 723
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0410017 (DELPHI)'
        S95_t1(x)%citekey = 'Abdallah:2004wy'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 50.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'Delphi_yuk_a_bbtautau'

        x = x + 1
        S95_t1(x)%id = 743
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0111010 (OPAL)'
        S95_t1(x)%citekey = 'Abbiendi:2001kp'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 12.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'OPAL_yuk_a_bbtautau'

        x = x + 1
        S95_t1(x)%id = 731
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0410017 (DELPHI)'
        S95_t1(x)%citekey = 'Abdallah:2004wy'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 27.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'Delphi_yuk_h_tautautautau'

        x = x + 1
        S95_t1(x)%id = 733
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0410017 (DELPHI)'
        S95_t1(x)%citekey = 'Abdallah:2004wy'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 26.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'Delphi_yuk_a_tautautautau'

        x = x + 1
        S95_t1(x)%id = 402
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0401022 (Delphi)'
        S95_t1(x)%citekey = 'Abdallah:2003ry'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 40.0D0
        S95_t1(x)%xmax = 114.0D0
        S95_t1(x)%sep = 2.0D0
        filename(x) = 'Delphi_h-invisible'

        x = x + 1
        S95_t1(x)%id = 403
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0501033 (L3)'
        S95_t1(x)%citekey = 'Achard:2004cf'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 50.0D0
        S95_t1(x)%xmax = 110.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'L3_h-invisible'

        x = x + 1
        S95_t1(x)%id = 401
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = '[hep-ex] arXiv:0707.0373 (OPAL)'
        S95_t1(x)%citekey = 'Abbiendi:2007ac'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 5.0D0
        S95_t1(x)%xmax = 115.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'OPAL_h-invisible'

        x = x + 1
        S95_t1(x)%id = 821
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0107031 (LHWG)'
        S95_t1(x)%citekey = 'Searches:2001ac'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 60.0D0
        S95_t1(x)%xmax = 90.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = 'LEP_HpHm_qqqq'

        x = x + 1
        S95_t1(x)%id = 811
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0404012 (Delphi)'
        S95_t1(x)%citekey = 'Abdallah:2003wd'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 52.0D0
        S95_t1(x)%xmax = 94.0D0
        S95_t1(x)%sep = 2.0D0
        filename(x) = 'Delphi_HpHm_qqqq'

        x = x + 1
        S95_t1(x)%id = 813
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'LEP'
        S95_t1(x)%label = 'hep-ex/0404012 (Delphi)'
        S95_t1(x)%citekey = 'Abdallah:2003wd'
        S95_t1(x)%energy = 0.208D0
        S95_t1(x)%xmin = 52.0D0
        S95_t1(x)%xmax = 94.0D0
        S95_t1(x)%sep = 2.0D0
        filename(x) = 'Delphi_HpHm_taunutaunu'

!----------------------- Z H -> l l b b -------------------------

        x = x + 1
        S95_t1(x)%id = 10799
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10799'
        S95_t1(x)%citekey = 'TEVNPH:2012ab' ! cite combination
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.45D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'CDF_ZH_llbb_9.45fb_10799'

        x = x + 1
        S95_t1(x)%id = 6296
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6296'
        S95_t1(x)%citekey = 'D0:2012osa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'D0_ZH_llbb_9.7fb_6296'

        x = x + 1
        S95_t1(x)%id = 3564
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:1008.3564 (D0)'
        S95_t1(x)%citekey = 'Abazov:2010zk'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4.2D0
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'D0_ZH_llbb_4.2fb_3564'

!----------------------- V H -> b b Etmiss -------------------------

        x = x + 1
        S95_t1(x)%id = 10798
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10798'
        S95_t1(x)%citekey = 'TEVNPH:2012ab'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.45D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_VH_Metbb_9.45fb_10798'

        x = x + 1
        S95_t1(x)%id = 6299
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6340' ! 'D0 Note 6299'
        S95_t1(x)%citekey = 'J.-F.Grivaz:2012gba'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.5D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_VH_bb_9.5fb_6299'

        x = x + 1
        S95_t1(x)%id = 2012161
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-161'
        S95_t1(x)%citekey = 'ATLAS:2012aha'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 17.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 130.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012161_Atlas_VH-Vbb_17.7fb-1'

        x = x + 1
        S95_t1(x)%id = 13012
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1310.3687 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2013zna'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 135.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13012_CMS_VH_bb_24fb-1'

!----------------------- VBF(H), H -> b b -------------------------

        x = x + 1
        S95_t1(x)%id = 13011
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-13-011'
        S95_t1(x)%citekey = 'CMS:2013jda'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 135.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13011_CMS_VBF_bb_19fb-1'

!----------------------- W H -> b b -------------------------

        x = x + 1
        S95_t1(x)%id = 6309
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6309'
        S95_t1(x)%citekey = 'D0:2012psa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_WH_lnubb_9.7_6309'

        x = x + 1
        S95_t1(x)%id = 10796
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10796'
        S95_t1(x)%citekey = 'TEVNPH:2012ab' ! cite combination
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.45D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'CDF_WH_lnubb_9.45fb_10796'

!----------------------- V H, H -> invisible -------------------------

        x = x + 1
        S95_t1(x)%id = 3244
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1402.3244 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2014iia'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.8D0
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 400.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '3244_Atlas_H-inv_24.8fb-1'

        x = x + 1
        S95_t1(x)%id = 13442
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1404.1344 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2014tja'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.60
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 105.0D0
        S95_t1(x)%xmax = 145.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1344_CMS_ZH-invisible_24.6fb-1'

!----------------------- VBF, H -> invisible -------------------------

        x = x + 1
        S95_t1(x)%id = 13441
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1404.1344 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2014tja'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.5D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 400.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1344_CMS_VBF-invisible_19.5fb-1'

        x = x + 1
        S95_t1(x)%id = 13443
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1404.1344 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2014tja'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%lumi = 24.60
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 145.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1344_CMS_ZH_VBF-invisible_24.6fb-1'

        x = x + 1
        S95_t1(x)%id = 6682
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1809.06682 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018sfi'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 75.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 5.0D0
        filename(x) = '6682_Atlas_VBF_h-invisible_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 5105
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ph] 1904.05105 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2019rtt'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.610000E+01
        S95_t1(x)%xmin = 1.200000E+02
        S95_t1(x)%xmax = 1.300000E+02
        S95_t1(x)%sep = 1.000000E+01
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '5105_ATLAS_H-inv_combined_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 59371
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ph] 1809.05937 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018owy'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.820000E+01
        S95_t1(x)%xmin = 1.200000E+02
        S95_t1(x)%xmax = 1.300000E+02
        S95_t1(x)%sep = 1.000000E+01
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '59371_CMS_H-inv_combined_38.2fb-1'

        x = x + 1
        S95_t1(x)%id = 59372
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ph] 1809.05937 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018owy'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.820000E+01
        S95_t1(x)%xmin = 1.100000E+02
        S95_t1(x)%xmax = 1.000000E+03
        S95_t1(x)%sep = 1.000000E+01
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '59372_CMS_H-inv_38.2fb-1'

!----------------------- H -> W W -------------------------

        x = x + 1
        S95_t1(x)%id = 5757
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 5757'
        S95_t1(x)%citekey = 'HaraldFox:2008yta'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 3.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_ppH_WW_ll_3.0fb_5757'

        x = x + 1
        S95_t1(x)%id = 3930
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = '[hep-ex] arXiv:0809.3930 (CDF)'
        S95_t1(x)%citekey = 'Aaltonen:2008ec'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 3.0D0
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'CDF_ggH_WW_3.0fb_3930'

        x = x + 1
        S95_t1(x)%id = 6276
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6276'
        S95_t1(x)%citekey = 'D0:2012ssa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_H-VV_9.7fb_6276'

        x = x + 1
        S95_t1(x)%id = 6301
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6301'
        S95_t1(x)%citekey = 'D0:2012jsa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'D0_VH_VWW_9.7fb_6301'

        x = x + 1
        S95_t1(x)%id = 10600
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10599'
        S95_t1(x)%citekey = 'Benjamin:2011sv' ! cite combination
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 8.2D0
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'CDF_ggH-WW_8.2fb_10600_interpol'

        x = x + 1
        S95_t1(x)%id = 10599
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10599'
        S95_t1(x)%citekey = 'Benjamin:2011sv' ! cite combination
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 8.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_H-WW_8.2fb_10599'

        x = x + 1
        S95_t1(x)%id = 4468
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = '[hep-ex] arXiv:1001.4468 (CDF)'
        S95_t1(x)%citekey = 'Aaltonen:2010cm'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_H-WW_4.8fb_4468_interpol'

        x = x + 1
        S95_t1(x)%id = 6302
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6302'
        S95_t1(x)%citekey = 'D0:2012rsa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_H-WW_9.7fb_6302'

        x = x + 1
        S95_t1(x)%id = 6183
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6183'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 8.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 130.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_SM_combined_6183'

        x = x + 1
        S95_t1(x)%id = 4481
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:1001.4481 (D0)'
        S95_t1(x)%citekey = 'Abazov:2010ct'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 5.4D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_H-WW_5.4fb_4481'

        x = x + 1
        S95_t1(x)%id = 3331
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'TCB'
        S95_t1(x)%label = '[hep-ex] arXiv:1108.3331 (TEVNPHWG)'!CDF note 10608, D0 Note 6230
        S95_t1(x)%citekey = 'Benjamin:2011sv'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 8.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_D0_combined_gg-H-WW_8.2fb_3331_bayesian_interpol'

        x = x + 1
        S95_t1(x)%id = 2012012
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-012'
        S95_t1(x)%citekey = 'ATLAS:2015dka'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012012_Atlas_H-WW-lnulnu_4.7fb-1'

        x = x + 1
        S95_t1(x)%id = 2013030
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2013-030'
        S95_t1(x)%citekey = 'ATLAS:2013wla'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2013030_Atlas_H-WW-lnulnu_25fb-1'

        x = x + 1
        S95_t1(x)%id = 2577
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arxiv:1112.2577'
        S95_t1(x)%citekey = 'ATLAS:2011aa'
        S95_t1(x)%desc = 'pp->h + X->W W* + X ->l l nu nu'
        S95_t1(x)%energy = 7.D0
        S95_t1(x)%lumi = 2.05D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2577_Atlas_H-WW-lnulnu_2.05fb-1'

        x = x + 1
        S95_t1(x)%id = 13003
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-13-003'
        S95_t1(x)%citekey = 'CMS:bxa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13003_CMS_H-WW-lnulnu_25fb-1'

        x = x + 1
        S95_t1(x)%id = 13027
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-13-027'
        S95_t1(x)%citekey = 'CMS:2012bea'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 170.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13027_CMS_H-WW-lnujj_24.3fb-1'

        x = x + 1
        S95_t1(x)%id = 13022
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-13-022'
        S95_t1(x)%citekey = 'CMS:2013yea'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.4D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13022_CMS_VBF-WW_25.4fb-1'

        x = x + 1
        S95_t1(x)%id = 2016062
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1710.07235 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2017fgj'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.2D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 500.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 50.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016062_Atlas_gg-H-WW-lnuqq_13.2fb-1'

!--------------- H -> VV -> l nu l nu ----------------

        x = x + 1
        S95_t1(x)%id = 3357
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1109.3357 (ATLAS)'
        S95_t1(x)%citekey = 'ATLAS:2011af'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 1.04D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 20.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '3357_Atlas_H-ZZ-llnunu_1.04fb-1'

        x = x + 1
        S95_t1(x)%id = 2012016
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-016'
        S95_t1(x)%citekey = 'ATLAS:2012nja'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012016_Atlas_H-ZZ-llnunu_4.7fb-1'

        x = x + 1
        S95_t1(x)%id = 3478
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arxiv:1202.3478 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2012ft'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 250.0D0
        S95_t1(x)%xmax = 590.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '11026_CMS_H-ZZ-llnunu_4.6fb-1'

!----------------------- H -> W W -----------------------------

        x = x + 1
        S95_t1(x)%id = 12046
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-12-046'
        S95_t1(x)%citekey = 'CMS:wxa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 17D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 170.0D0
        S95_t1(x)%xmax = 580.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '12046_CMS_H-WW-lnuqq_17fb-1'

        x = x + 1
        S95_t1(x)%id = 20160741
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-074'
        S95_t1(x)%citekey = 'ATLAS:2016kjy'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 50.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016074_Atlas_gg-H-WW-lnulnu_13.2fb-1'

        x = x + 1
        S95_t1(x)%id = 20160742
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-074'
        S95_t1(x)%citekey = 'ATLAS:2016kjy'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 50.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016074_Atlas_vbf-H-WW-lnulnu_13.2fb-1'

        x = x + 1
        S95_t1(x)%id = 170331
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1912.01594 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019pqw'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '170331_CMS_gg_VBF-H-WW_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 170332
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1912.01594 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019pqw'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '170332_CMS_gg-H-WW_13TeV_35.9fb-1'

! H -> V V combination

        x = x + 1
        S95_t1(x)%id = 23801
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '(hep-ex) arXiv:1808.02380 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018bun'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 100.0D0
        filename(x) = '23801_Atlas_gg-phi-VV_36.1fb-1_13TeV'

        x = x + 1
        S95_t1(x)%id = 23802
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '(hep-ex) arXiv:1808.02380 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018bun'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 500.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 100.0D0
        filename(x) = '23802_Atlas_vbf-phi-VV_36.1fb-1_13TeV'

!------------------------- H -> Z Z --------------------------

        x = x + 1
        S95_t1(x)%id = 5064
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1108.5064 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2011ec'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 1.04D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '5064_Atlas_H-ZZ-llqq_1.04fb-1'

        x = x + 1
        S95_t1(x)%id = 2012017
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-017'
        S95_t1(x)%citekey = 'ATLAS:2012mja'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012017_Atlas_H-ZZ-llqq_4.7fb-1'

        x = x + 1
        S95_t1(x)%id = 14161
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1202.1416 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2012sn'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 130.0D0
        S95_t1(x)%xmax = 164.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '110271_CMS_H-ZZ-llqq_4.6fb-1'

        x = x + 1
        S95_t1(x)%id = 14162
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1202.1416 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2012sn'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '110272_CMS_H-ZZ-llqq_4.6fb-1'

        x = x + 1
        S95_t1(x)%id = 1415
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1202.1415 (ATLAS)'
        S95_t1(x)%citekey = 'ATLAS:2012ac'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1415_Atlas_H-ZZ-4l_4.8fb-1'

        x = x + 1
        S95_t1(x)%id = 2012092
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-092'
        S95_t1(x)%citekey = 'ATLAS:2012foa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 10.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012092_Atlas_H-ZZ-4l_10.6fb-1'

        x = x + 1
        S95_t1(x)%id = 20130131
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2013-013'
        S95_t1(x)%citekey = 'ATLAS:2013nma'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 180.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2013013-1_Atlas_H-ZZ-4l_incl_25fb-1'

        x = x + 1
        S95_t1(x)%id = 20130132
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2013-013'
        S95_t1(x)%citekey = 'ATLAS:2013nma'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 21.0D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = '2013013-2_Atlas_H-ZZ-4l_ggF_21fb-1'

        x = x + 1
        S95_t1(x)%id = 20130133
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2013-013'
        S95_t1(x)%citekey = 'ATLAS:2013nma'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 21.0D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = '2013013-3_Atlas_H-ZZ-4l_VBFVH_21fb-1'

        x = x + 1
        S95_t1(x)%id = 20160792
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-079' !superseded by 1712.06386
        S95_t1(x)%citekey = 'ATLAS:2016oum'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 14.8D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 400.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = '20160792_gg-H-ZZ-4l_14.8fb-1_lowmass'

        x = x + 1
        S95_t1(x)%id = 20160793
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-079' !superseded by 1712.06386
        S95_t1(x)%citekey = 'ATLAS:2016oum'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 14.8D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = '20160793_VBF-H-ZZ-4l_14.8fb-1'

        x = x + 1
        S95_t1(x)%id = 1997
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arxiv:1202.1997 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2012dg'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '11025_CMS_H-ZZ-4l_4.7fb-1'

        x = x + 1
        S95_t1(x)%id = 130021
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1312.5353 (CMS)' !'CMS-PAS-HIG-13-002'
        S95_t1(x)%citekey = 'Chatrchyan:2013mxa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13002-1_CMS_H-ZZ-4l_lowm_25fb-1'

        x = x + 1
        S95_t1(x)%id = 130022
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1312.5353 (CMS)' !'CMS-PAS-HIG-13-002'
        S95_t1(x)%citekey = 'Chatrchyan:2013mxa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 150.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13002-2_CMS_H-ZZ-4l_highm_25fb-1'

        x = x + 1
        S95_t1(x)%id = 009361
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1504.00936 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2015cwa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 145.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '009361_CMS_H-VV_5.1fb-1_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 063861
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'arXiv:1712.06386 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2017rel'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 400.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = '063861_ATLAS_gg-H-ZZ-4l+2l2nu_NWA_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 20160821
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-082' ! superseded by 1708.09638
        S95_t1(x)%citekey = 'ATLAS:2016npe'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.2D0
        S95_t1(x)%xmin = 500.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = '20160821_Atlas_gg-H-ZZ-llnunu_13.2fb-1'

        x = x + 1
        S95_t1(x)%id = 20160822
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-082' ! superseded by 1708.09638
        S95_t1(x)%citekey = 'ATLAS:2016npe'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.2D0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = '20160822_Atlas_gg-H-ZZ-llqq_13.2fb-1'

        x = x + 1
        S95_t1(x)%id = 20160823
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-082' ! superseded by 1708.09638
        S95_t1(x)%citekey = 'ATLAS:2016npe'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.2D0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 1.0D0
        filename(x) = '20160823_Atlas_qq-Hqq-ZZ-llqq_13.2fb-1'

        x = x + 1
        S95_t1(x)%id = 16034
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-16-034'
        S95_t1(x)%citekey = 'CMS:2017sbi'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 12.9D0
        S95_t1(x)%xmin = 550.0D0
        S95_t1(x)%xmax = 2000.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = '16034_CMS_H-ZZ-llqq_12.9fb-1'

!------------------------- SM combined -----------------------

        x = x + 1
        S95_t1(x)%id = 10439
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10439'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 6.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_H-tautau_6.0fb_10439'

        x = x + 1
        S95_t1(x)%id = 5845
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 5845'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 105.0D0
        S95_t1(x)%xmax = 145.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_tautauqq_4.9fb_5845'

        x = x + 1
        S95_t1(x)%id = 9999
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 9999'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_SM_combined_2.0-4.8fb_9999_interpol'

        x = x + 1
        S95_t1(x)%id = 6305
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6305'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 7.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 105.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_H-tautau_7.3fb_6305'

        x = x + 1
        S95_t1(x)%id = 10010
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10010'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_hadr_4fb_10010'

        x = x + 1
        S95_t1(x)%id = 6171
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6171'
        S95_t1(x)%citekey = 'SubhenduChakrabarti:2011gfa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 105.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_H-tautaujj-4.3fb_6171'

        x = x + 1
        S95_t1(x)%id = 6304
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6304'
        S95_t1(x)%citekey = 'D0:2012asa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_SM_combined_6304'

        x = x + 1
        S95_t1(x)%id = 10500
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10500'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 6.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_VH_tautau_6.2fb_10500'

        x = x + 1
        S95_t1(x)%id = 10573
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10573'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 8.2D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 120.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_H-WW_8.2fb_10573'

        x = x + 1
        S95_t1(x)%id = 6286
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6286'
        S95_t1(x)%citekey = 'D0:2012msa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 7.0D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_tautaumu_7.0fb_6286'

        x = x + 1
        S95_t1(x)%id = 10884
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'TCB'
        S95_t1(x)%label = '[hep-ex] arXiv:1207.0449 (TEVNPHWG)'
        S95_t1(x)%citekey = 'Group:2012zca'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 10.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '10884_CDF_D0_SM_combined_10fb-1'

        x = x + 1
        S95_t1(x)%id = 6436
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'TCB'
        S95_t1(x)%label = '[hep-ex] arXiv:1207.6436 (TEVNPHWG)'
        S95_t1(x)%citekey = 'Aaltonen:2012qt'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1207.6436_CDF_D0_combined_h-bb_9.7fb-1'

        x = x + 1
        S95_t1(x)%id = 1408
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '(hep-ex) arxiv:1202.1408 (ATLAS)'
        S95_t1(x)%citekey = 'ATLAS:2012ae'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1408_Atlas_SMcombined_4.9fb-1'

        x = x + 1
        S95_t1(x)%id = 2012019
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-019'
        S95_t1(x)%citekey = 'ATLAS:2012kja'
        S95_t1(x)%desc = '(p p)->h+...'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012019_Atlas_SMcombined_4.9fb-1'

        x = x + 1
        S95_t1(x)%id = 7214
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '(hep-ex) arXiv:1207.7214 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2012tfa'
        S95_t1(x)%desc = '(p p)->h+...'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 10.5D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '7214_Atlas_SMcombined_10.5fb-1'

        x = x + 1
        S95_t1(x)%id = 2011157
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2011-157, CMS-PAS-HIG-11-023'
        S95_t1(x)%citekey = 'ATLAS:2011jka'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 2.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2011157_Atlas_CMS_SMcombined_2.7fb-1'

        x = x + 1
        S95_t1(x)%id = 1488
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arxiv:1202.1488 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2012tx'
        S95_t1(x)%energy = 7.D0
        S95_t1(x)%lumi = 4.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1488_CMS_SMcombined_4.8fb-1'

        x = x + 1
        S95_t1(x)%id = 12045
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-12-045'
        S95_t1(x)%citekey = 'CMS:aya'
        S95_t1(x)%desc = '(p p)->h+...'
        S95_t1(x)%energy = 8.D0
        S95_t1(x)%lumi = 17.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '12045_CMS_SMcombined_17.3fb-1'

!----------------- b H -> 3 b jets -------------------

        x = x + 1
        S95_t1(x)%id = 1931
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:1011.1931 (D0)'
        S95_t1(x)%citekey = 'Abazov:2010ci'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 5.2D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'D0_bbH_bb_5.2fb_1931'

        x = x + 1
        S95_t1(x)%id = 4782
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'arXiv:1106.4782 (CDF)'
        S95_t1(x)%citekey = 'Abazov:2010ci'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 2.6D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 350.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'CDF_bbH_bb_2.6fb_4782'

        x = x + 1
        S95_t1(x)%id = 1508329
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'arXiv:1506.08329 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2015tra'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 900.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = '1508329_CMS_Hb-bbb_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 2015080
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'CERN-THESIS-2015-080 (ATLAS)'
        S95_t1(x)%citekey = 'Malone:2015mia'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.5D0
        S95_t1(x)%xmin = 450.0D0
        S95_t1(x)%xmax = 800.0D0
        S95_t1(x)%sep = 50.0D0
        filename(x) = '2015080_Atlas_hb-bbb_19.5fb-1'

        x = x + 1
        S95_t1(x)%id = 2749
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1907.02749 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2019zwb'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 2.780000E+01
        S95_t1(x)%xmin = 4.500000E+02
        S95_t1(x)%xmax = 1.400000E+03
        S95_t1(x)%sep = 5.000000E+01
        S95_t1(x)%deltax = 2.000000E+01
        filename(x) = '2749_ATL_Hb-bbbb_27.8fb-1'

        x = x + 1
        S95_t1(x)%id = 12191
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1805.12191 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018taj'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.570000E+01
        S95_t1(x)%xmin = 3.000000E+02
        S95_t1(x)%xmax = 1.300000E+03
        S95_t1(x)%sep = 2.000000E+01
        filename(x) = '12191_CMS_bH-bbb_35.7fb-1'
        S95_t1(x)%deltax = 1.000000E+01

!----------------------- H -> b b (inclusive) -------------------------
        x = x + 1
        S95_t1(x)%id = 180509299
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1805.09299 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018tqo'
        S95_t1(x)%energy = 13D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 1.25D3
        S95_t1(x)%xmax = 3D3
        S95_t1(x)%sep = 5D1
        filename(x) = '180509299_Atlas_bb_36.1fb-1'
        S95_t1(x)%deltax = 2.5D1

! ------------------- H -> bb (boosted) --------------------

        x = x + 1
        S95_t1(x)%id = 1810118221
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1810.11822 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018ikr'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 5.000000E+01
        S95_t1(x)%xmax = 3.500000E+02
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '181011822_CMS_H_bb_35.9fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 1810118222
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1810.11822 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018ikr'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 5.000000E+01
        S95_t1(x)%xmax = 3.500000E+02
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '181011822_CMS_A_bb_35.9fb-1'
        S95_t1(x)%deltax = 0.000000E+00

!-------------------- H -> tau tau -------------------------

        x = x + 1
        S95_t1(x)%id = 4555
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:1106.4555 (D0)'
        S95_t1(x)%citekey = 'Abazov:2011jh'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 5.4D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'D0_H-tautau_5.4fb_4555_interpol'

        x = x + 1
        S95_t1(x)%id = 1014
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = '[hep-ex] arXiv:0906.1014 (CDF)'
        S95_t1(x)%citekey = 'Aaltonen:2009vf'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 1.8D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 250.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'CDF_H-tautau_1.8fb_1014_interpol'

        x = x + 1
        S95_t1(x)%id = 3363
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'TCB'
        S95_t1(x)%label = '[hep-ex] arXiv:1003.3363 (TEVNPHWG)'
        S95_t1(x)%citekey = 'Benjamin:2010xb'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 2.2D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'CDF_D0_combinedH-tautau_2.2fb_3363_CLs'

        x = x + 1
        S95_t1(x)%id = 2012160
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-160'
        S95_t1(x)%citekey = 'ATLAS:2012dsy'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 17.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012160_Atlas_H-tautau_17.6fb-1'

! Analysis is disabled when using CMS tau tau likelihood
        x = x + 1
        S95_t1(x)%id = 2014049
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'arXiv:1409.6064' ! ATLAS-CONF-2014-049
        S95_t1(x)%citekey = 'Aad:2014vgg'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%llh = 1
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 20.0D0
        filename(x) = '2014049_Atlas_bbh_h-tautau_20.3fb-1'

! Analysis is disabled when using CMS tau tau likelihood
        x = x + 1
        S95_t1(x)%id = 20140492
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'arXiv:1409.6064' ! ATLAS-CONF-2014-049
        S95_t1(x)%citekey = 'Aad:2014vgg'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%llh = 1
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 20.0D0
        filename(x) = '20140492_Atlas_ggh_h-tautau_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 12043
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-12-043'
        S95_t1(x)%citekey = 'CMS:lng'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 17.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 145.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '12043_CMS_H_tautau_17fb-1'

!-------------------- Higgs to Higgs decays -------------------------
!
! Requires hi -> hj Z branching ratio in input!
!
        x = x + 1
        S95_t1(x)%id = 14011
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'arXiv:1504.04710 (CMS)' !CMS-PAS-HIG-14-011
        S95_t1(x)%citekey = 'Khachatryan:2015lba'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 225.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '14011_CMS_A-Zh-llbb_19.7fb-1'
!
        x = x + 1
        S95_t1(x)%id = 04670
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.04670 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015xja'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '04670_Atlas_H-hh-combined_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 046701
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.04670 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015xja'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 500.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '04670_Atlas_H-hh-gagaWW_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 8567
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1807.08567 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018ewm'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 500.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '8567_Atlas_X-HH-WWgaga_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 4873
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1807.04873 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018ewm'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '4873_Atlas_X-HH-gagabb_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 046702
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.04670 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015xja'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '04670_Atlas_H-hh-bbtautau_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 336
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1808.00336 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018sfw'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '336_Atlas_X-HH-bbtautau_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 2016049
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-049' ! superseded by 1804.06174
        S95_t1(x)%citekey = 'ATLAS:2016ixk'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 50.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016049_Atlas_H-hh-bbbb_13.3fb-1'

        x = x + 1
        S95_t1(x)%id = 011811
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1510.01181 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2015tha'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 350.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '140341_CMS_H-hh-bbtautau_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 17020321
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1701.02032 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2017mnf'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 5.0D0
        S95_t1(x)%xmax = 15.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '170102032_CMS_HSM-aa-tautau_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 18006
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1907.07235 (CMS)' ! 'CMS-PAS-HIG-18-006'
        S95_t1(x)%citekey = 'Sirunyan:2019gou'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 15.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.5D0
        filename(x) = '18006_CMS_pp-H-hh-tautautautau_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 7355
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1806.07355 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018iil'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 20.0D0
        S95_t1(x)%xmax = 60.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '7355_Atlas_VH-aa-bbbb_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 17020322
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1701.02032 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2017mnf'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 25.0D0
        S95_t1(x)%xmax = 62.5D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '170102032_CMS_HSM-aa-mumubb_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 539
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1807.00539 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018esj'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 20.0D0
        S95_t1(x)%xmax = 60.0D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '539_Atlas_H-aa-bbmumu_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 6359
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1812.06359 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018mot'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 20.50D0
        S95_t1(x)%xmax = 61.5D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '6359_CMS_pp-H-hh-mumubb_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 17020323
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1701.02032 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2017mnf'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 15.0D0
        S95_t1(x)%xmax = 62.5D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '170102032_CMS_HSM-aa-mumutautau_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 4865
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1805.04865 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018mbx'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 15.0D0
        S95_t1(x)%xmax = 61.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '4865_CMS_H125-aa-tautaumumu_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 10191
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1805.10191 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2018pzn'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 15.0D0
        S95_t1(x)%xmax = 60.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 2.5D0
        filename(x) = '10191_CMS_H-aa-bbtautau_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 5051
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.05051 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015bua'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 10.0D0
        S95_t1(x)%xmax = 61.75D0
        S95_t1(x)%sep = 0.05D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '5051_Atlas_gg-H-aa-gagagaga_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 1506534
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1510.06534 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2015nba'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 4.0D0
        S95_t1(x)%xmax = 8.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1506534_CMS_gg-H-hh-tautautautau_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 17006
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1708.04188 (CMS)' !'CMS-PAS-HIG 17-006'
        S95_t1(x)%citekey = 'Sirunyan:2017guj'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.0D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 900.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '17006_CMS_H-hh-bbVV-bblnulnu_36fb-1'

        x = x + 1
        S95_t1(x)%id = 011812
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1510.01181 (CMS)'
        S95_t1(x)%citekey = 'Khachatryan:2015tha'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 350.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '140342_CMS_A-hZ-tautaull_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 044781
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1502.04478 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015wra'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '044781_ATLAS_gg-A-hZ-tautaull_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 044782
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1502.04478 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015wra'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '044782_ATLAS_gg-A-hZ-bbll_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 180051
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1903.00941 (CMS)' ! 'CMS-PAS-HIG-18-005'
        S95_t1(x)%citekey = 'Sirunyan:2019xls'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '180051_CMS_gg-A-Zh-Zbb_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 180052
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1903.00941 (CMS)' ! 'CMS-PAS-HIG-18-005'
        S95_t1(x)%citekey = 'Sirunyan:2019xls'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '180052_CMS_bb-A-Zh-Zbb_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 1712065181
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1712.06518 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2017cxo'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 2000.0D0
        S95_t1(x)%sep = 2.0D1
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '171206518_ATLAS_gg-A-Zh-Zbb_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 1712065182
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1712.06518 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2017cxo'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 220.0D0
        S95_t1(x)%xmax = 2000.0D0
        S95_t1(x)%sep = 2.0D1
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '171206518_ATLAS_bb-A-Zh-Zbb_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 2020043
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2020-043'
        S95_t1(x)%citekey = 'ATLAS:2020pgp'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.390000E+02
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 3.000000E+02
        S95_t1(x)%xmax = 2.000000E+03
        S95_t1(x)%sep = 2.000000E+01
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2020043_ATLAS_gg-A-Zh-Zbb_13TeV_139fb-1'

        x = x + 1
        S95_t1(x)%id = 16002
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-16-002' ! superseded by 1806.03548
        S95_t1(x)%citekey = 'CMS:2016tlj'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 2.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 1200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '16002_CMS_H-hh-bbbb_2.3fb-1'

        x = x + 1
        S95_t1(x)%id = 17002
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1707.02909 (CMS)' ! 'CMS-PAS-HIG-17-002'
        S95_t1(x)%citekey = 'Sirunyan:2017djm'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 250.0D0
        S95_t1(x)%xmax = 900.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '17002_CMS_H-hh-tautaubb_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 14013
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'arXiv:1503.04114 (CMS)' ! CMS-PAS-HIG-14-013
        S95_t1(x)%citekey = 'Khachatryan:2015yea'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 17.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 270.0D0
        S95_t1(x)%xmax = 1097.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '14013_H-hh-bbbb_17.9fb-1'

        x = x + 1
        S95_t1(x)%id = 17030
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1811.09689 (CMS)' ! 'CMS-PAS-HIG-17-030'
        S95_t1(x)%citekey = 'Sirunyan:2018two'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 270.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '17030_CMS_pp-X-HH_comb_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 2025
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'arXiv:1906.02025 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2019uzh'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2025_ATLAS_pp-X-HH_comb_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 110281
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'arXiv:1811.11028 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018ksn'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 260.0D0
        S95_t1(x)%xmax = 500.0D0
        S95_t1(x)%sep = 20.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '11028_ATLAS_X-HH-WWWW_13TeV_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 18013
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-18-013'
        S95_t1(x)%citekey = 'CMS:2019vgr'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 2.600000E+02
        S95_t1(x)%xmax = 1.000000E+03
        S95_t1(x)%sep = 2.000000E+01
        filename(x) = '18013_CMS_HH-bbZZ_35.9fb-1'
        S95_t1(x)%deltax = 1.000000E+01

!-------------------- H -> mu mu -------------------------

        x = x + 1
        S95_t1(x)%id = 2013010
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2013-010'
        S95_t1(x)%citekey = 'ATLAS:2013qma'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 21D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2013010_Atlas_H-mumu_21fb-1'

        x = x + 1
        S95_t1(x)%id = 7663
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1406.7663 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2014xva'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 120.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '7663_Atlas_H-mumu_24.8fb-1'

        x = x + 1
        S95_t1(x)%id = 31521
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1907.03152 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019tkw'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 1.300000E+02
        S95_t1(x)%xmax = 1.000000E+03
        S95_t1(x)%sep = 5.000000E+00
        S95_t1(x)%deltax = 0.000000E+00
        filename(x) = '3152_gg_CMS_H-mumu_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 31522
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1907.03152 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019tkw'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 1.300000E+02
        S95_t1(x)%xmax = 1.000000E+03
        S95_t1(x)%sep = 5.000000E+00
        S95_t1(x)%deltax = 0.000000E+00
        filename(x) = '3152_bb_CMS_H-mumu_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 81441
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1901.08144 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2019sgt'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.610000E+01
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 1.000000E+03
        S95_t1(x)%sep = 1.000000E+01
        S95_t1(x)%deltax = 0.000000E+00
        filename(x) = '8144_gg_ATL_H-mumu_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 81442
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1901.08144 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2019sgt'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.610000E+01
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 1.000000E+03
        S95_t1(x)%sep = 1.000000E+01
        S95_t1(x)%deltax = 0.000000E+00
        filename(x) = '8144_bb_ATL_H-mumu_36.1fb-1'

!----------------------- H -> mu tau ----------------------------
        x = x + 1
        S95_t1(x)%id = 180171
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1911.10267 (CMS)' ! 'CMS-PAS-HIG-18-017'
        S95_t1(x)%citekey = 'Sirunyan:2019shc'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 9.000000E+02
        S95_t1(x)%sep = 2.000000E+01
        S95_t1(x)%deltax = 0.000000E+00
        filename(x) = '18017_CMS_H-mutau_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 61312
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1907.06131 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2019ugc'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.610000E+01
        S95_t1(x)%xmin = 1.200000E+02
        S95_t1(x)%xmax = 1.300000E+02
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '6131_Atlas_H125-mutau_36.1fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 1807065734
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1807.06573 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018jff'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 5.000000E+02
        S95_t1(x)%xmax = 3.00000E+03
        S95_t1(x)%sep = 1.000000E+02
        filename(x) = '180706573_Atlas_H-mutau_36.1fb-1'
        S95_t1(x)%deltax = 0.000000E+00

!----------------------- H -> e tau ----------------------------
        x = x + 1
        S95_t1(x)%id = 180172
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1911.10267 (CMS)' ! 'CMS-PAS-HIG-18-017'
        S95_t1(x)%citekey = 'Sirunyan:2019shc'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 9.000000E+02
        S95_t1(x)%sep = 2.000000E+01
        S95_t1(x)%deltax = 0.000000E+00
        filename(x) = '18017_CMS_H-etau_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 61311
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1907.06131 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2019ugc'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.610000E+01
        S95_t1(x)%xmin = 1.200000E+02
        S95_t1(x)%xmax = 1.300000E+02
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '6131_Atlas_H125-etau_36.1fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 1807065733
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1807.06573 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018jff'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 5.000000E+02
        S95_t1(x)%xmax = 3.00000E+03
        S95_t1(x)%sep = 1.000000E+02
        filename(x) = '180706573_Atlas_H-etau_36.1fb-1'
        S95_t1(x)%deltax = 0.000000E+00

!----------------------- H -> e mu ----------------------------

        x = x + 1
        S95_t1(x)%id = 102351
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1909.10235 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2019ojw'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.390000E+02
        S95_t1(x)%xmin = 1.200000E+02
        S95_t1(x)%xmax = 1.300000E+02
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '10235_Atlas_H125-emu_139fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 180201122
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1802.01122 (ATLAS)'
        S95_t1(x)%citekey = 'Sirunyan:2018zhy'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 3.00000E+03
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '180201122_CMS_H_emu_35.9fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 1807065731
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1807.06573 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018jff'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 5.000000E+02
        S95_t1(x)%xmax = 3.00000E+03
        S95_t1(x)%sep = 1.000000E+02
        filename(x) = '180706573_Atlas_H-emu_36.1fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 1807065732
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1807.06573 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018jff'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 5.000000E+02
        S95_t1(x)%xmax = 3.00000E+03
        S95_t1(x)%sep = 1.000000E+02
        filename(x) = '180706573_Atlas_ggH-emu_36.1fb-1'
        S95_t1(x)%deltax = 0.000000E+00

!------------------------ H W -> W W W --------------------------

        x = x + 1
        S95_t1(x)%id = 5873
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 5873'
        S95_t1(x)%citekey = 'MaikoTakahashi:2009mxa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 3.6D0
        S95_t1(x)%xmin = 120.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 20.0D0
        filename(x) = 'D0_WH_WWW_llnunu_3.6fb_5873'

        x = x + 1
        S95_t1(x)%id = 1268
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:1107.1268 (D0)'
        S95_t1(x)%citekey = 'Abazov:2011ed'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 5.3D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 115.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_VH_ll_5.3fb_1268'

        x = x + 1
        S95_t1(x)%id = 13009
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-13-009'
        S95_t1(x)%citekey = 'CMS:zwa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13009_CMS_WH-WWW_25fb-1'

        x = x + 1
        S95_t1(x)%id = 12006
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-12-006'
        S95_t1(x)%citekey = 'CMS:2012yta'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 140.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '12006_CMS_WH-Wtautau_4.7fb-1'

        x = x + 1
        S95_t1(x)%id = 12051
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-12-051'
        S95_t1(x)%citekey = 'CMS:xxa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 17.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 145.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '12051_CMS_VH_Vtautau_17fb-1'

        x = x + 1
        S95_t1(x)%id = 2012078
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-078'
        S95_t1(x)%citekey = 'ATLAS:2012toa'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 300.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012078_Atlas_WH_WWW_4.7fb-1'

!-------------------- H -> gamma gamma --------------------

        x = x + 1
        S95_t1(x)%id = 6295
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6295'
        S95_t1(x)%citekey = 'GuoChen:2012yga'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 9.7D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 2.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_gaga_9.7fb_6295'

        x = x + 1
        S95_t1(x)%id = 4960
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'TCB'
        S95_t1(x)%label = '[hep-ex] arXiv:1107.4960 (TEVNPHWG)'!CDF Note 10510, D0 Note 6203
        S95_t1(x)%citekey = 'TEVNPHWorking:2011aa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 8.2D0
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_D0_combined_SM_Hgaga_8.2fb_4960'

        x = x + 1
        S95_t1(x)%id = 1414
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1202.1414 (ATLAS)'
        S95_t1(x)%citekey = 'ATLAS:2012ad'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '1414_Atlas_H-gaga_4.9fb-1'

        x = x + 1
        S95_t1(x)%id = 2012168
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-168'
        S95_t1(x)%citekey = 'ATLAS:2012znl'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 17.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012168_Atlas_H-gaga_17.8fb-1'

        x = x + 1
        S95_t1(x)%id = 6583
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1407.6583'
        S95_t1(x)%citekey = 'Aad:2014ioa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 65.0D0
        S95_t1(x)%xmax = 600.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '6583_Atlas_H-gaga_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 059301
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1507.05930 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015kna'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 140.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '05930_Atlas_gg-H-ZZ-20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 059302
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1507.05930 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015kna'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 140.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '05930_Atlas_VBF-H-ZZ-20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 0038911
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.00389 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015agg'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 100.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '0038911_Atlas_ggH-WW-NWA_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 0038912
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.00389 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015agg'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 100.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '0038912_Atlas_VBF-WW-NWA_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 0038913
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.00389 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015agg'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 20.0D0
        S95_t1(x)%deltax = 20.0D0
        filename(x) = '0038913_Atlas_ggH-WW-CPS_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 0038914
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1509.00389 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015agg'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 20.0D0
        S95_t1(x)%deltax = 20.0D0
        filename(x) = '0038914_Atlas_VBF-WW-CPS_20.3fb-1'

        x = x + 1
        S95_t1(x)%id = 17013
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1811.08459 (CMS)' ! 'CMS-PAS-HIG-17-013'
        S95_t1(x)%citekey = 'Sirunyan:2018aui'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 110.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '17013_CMS_H-gaga_lowmass_8TeV_20fb-1_and_13TeV_36fb-1'

        x = x + 1
        S95_t1(x)%id = 13001
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1407.0558 (CMS)' ! 'CMS-PAS-HIG-13-001'
        S95_t1(x)%citekey = 'Khachatryan:2014ira'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 25.D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13001_CMS_H-gaga_MVA_25fb-1'

        x = x + 1
        S95_t1(x)%id = 14037
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-14-037'
        S95_t1(x)%citekey = 'CMS:2015ocq'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 110.0D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '14037_CMS_H-gaga_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 14031
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-14-031'
        S95_t1(x)%citekey = 'CMS:2015lza'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 500.0D0
        S95_t1(x)%sep = 1.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '14031_CMS_H-Zgamma-19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 2018025
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2018-025'
        S95_t1(x)%citekey = 'ATLAS:2018xad'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 80.0D0
        S95_t1(x)%SMlike = 0
        S95_t1(x)%xmin = 65.0D0
        S95_t1(x)%xmax = 110.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2018025_Atlas_h-gaga_13TeV_80fb-1'

!-------------------- H -> gamma Z --------------------

        x = x + 1
        S95_t1(x)%id = 13075515
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] arXiv:1307.5515 (CMS)'
        S95_t1(x)%citekey = 'Chatrchyan:2013vaa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.6D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 120.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13075515_CMS_H-Zgamma_24.6fb-1'

        x = x + 1
        S95_t1(x)%id = 3051
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1402.3051 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2014fia'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 24.8D0
        S95_t1(x)%SMlike = 1
        S95_t1(x)%xmin = 120.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '3051_Atlas_H-Zga_24.8fb-1'

!--------------------- b H -> b tau tau ---------------------------

        x = x + 1
        S95_t1(x)%id = 4885
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:1106.4885 (D0)'
        S95_t1(x)%citekey = 'Abazov:2011qz'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 7.3D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 320.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'D0_Hb_tautaub_7.3fb_4885'

        x = x + 1
        S95_t1(x)%id = 6083
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 6083'
        S95_t1(x)%citekey = 'FabriceCouderc:2010yca'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 4.3D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 200.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'D0_Hb_tautaub_4.3fb_6083'

!-------------------- t t H -> t t b b --------------------

        x = x + 1
        S95_t1(x)%id = 5739
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = 'D0 Note 5739'
        S95_t1(x)%citekey = 'ChristianSchwanenberger:2008ita'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 2.1D0
        S95_t1(x)%xmin = 105.0D0
        S95_t1(x)%xmax = 155.0D0
        S95_t1(x)%sep = 10.0D0
        filename(x) = 'D0_ttH_ttbb_2.1fb_5739'

        x = x + 1
        S95_t1(x)%id = 10574
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = 'CDF Note 10574'
        S95_t1(x)%citekey = ''
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 7.5D0
        S95_t1(x)%xmin = 100.0D0
        S95_t1(x)%xmax = 170.0D0
        S95_t1(x)%sep = 5.0D0
        filename(x) = 'CDF_ttH_ttbb_7.5fb_10574_interpol'

        x = x + 1
        S95_t1(x)%id = 2012135
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2012-135'
        S95_t1(x)%citekey = 'ATLAS:2012cpa'
        S95_t1(x)%energy = 7D0
        S95_t1(x)%lumi = 4.7D0
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 140.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2012135_Atlas_ttH_Hbb_4.7fb-1'

        x = x + 1
        S95_t1(x)%id = 12025
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-12-025'
        S95_t1(x)%citekey = 'CMS:2012ywa'
        S95_t1(x)%energy = 7D0
        S95_t1(x)%lumi = 5.D0
        S95_t1(x)%xmin = 110.0D0
        S95_t1(x)%xmax = 140.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '12025_CMS_ttH-ttbb_5fb-1'

!-------------------- t t H -> t t mu mu --------------------
        x = x + 1
        S95_t1(x)%id = 4968
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1911.04968 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019bgz'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.370000E+02
        S95_t1(x)%xmin = 1.500000E+01
        S95_t1(x)%xmax = 3.400000E+02
        S95_t1(x)%sep = 5.000000E+00
        filename(x) = '4968_CMS_ttH-mumu_137fb-1'
        S95_t1(x)%deltax = 0.000000E+00

!-------------------- ttH -> tttt ---------------------------

        x = x + 1
        S95_t1(x)%id = 1908064631
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1908.06463 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019wxt'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.370000E+02
        S95_t1(x)%xmin = 3.500000E+02
        S95_t1(x)%xmax = 6.500000E+02
        S95_t1(x)%sep = 2.000000E+01
        filename(x) = '190806463_CMS_ttH-tttt_13TeV_137fb-1'
        S95_t1(x)%deltax = 0D0

        x = x + 1
        S95_t1(x)%id = 1908064632
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1908.06463 (CMS)'
        S95_t1(x)%citekey = 'Sirunyan:2019wxt'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.370000E+02
        S95_t1(x)%xmin = 3.500000E+02
        S95_t1(x)%xmax = 6.500000E+02
        S95_t1(x)%sep = 2.000000E+01
        filename(x) = '190806463_CMS_ttA-tttt_13TeV_137fb-1'
        S95_t1(x)%deltax = 0D0

!-------------------- H -> Z gamma --------------------------

        x = x + 1
        S95_t1(x)%id = 0611
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:0806.0611 (D0)'
        S95_t1(x)%citekey = 'Abazov:2008wg'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 1.1D0
        S95_t1(x)%xmin = 120.0D0
        S95_t1(x)%xmax = 320.0D0
        S95_t1(x)%sep = 20.0D0
        filename(x) = 'D0_H-Zgamma_1.0-1.1fb_0611'

!---------------- Daniel's attempts --------------
        x = x + 1
        S95_t1(x)%id = 2016025
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-HIG-16-025'
        S95_t1(x)%citekey = 'CMS:2016ncz'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 2.69D0
        S95_t1(x)%xmin = 550.0D0
        S95_t1(x)%xmax = 1200.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016025_CMS_H-bb_2.69fb-1'

        x = x + 1
        S95_t1(x)%id = 201608391
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1606.08391 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2016oyb'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 3.2D0
        S95_t1(x)%xmin = 20.0D0
        S95_t1(x)%xmax = 60.0D0
        S95_t1(x)%sep = 0.5D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '160608391_Atlas_HW_Waa_Wbbbb_3.2fb-1'

        x = x + 1
        S95_t1(x)%id = 2016044
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-044' ! superseded by 1708.00212
        S95_t1(x)%citekey = 'ATLAS:2016lri'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.3D0
        S95_t1(x)%xmin = 250.0D0
        S95_t1(x)%xmax = 2500.0D0
        S95_t1(x)%sep = 2D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016044_Atlas_H_ZA_13.3fb-1'

        x = x + 1
        S95_t1(x)%id = 1604833
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 1606.04833 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2016okv'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 3.2D0
        S95_t1(x)%xmin = 500.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 10D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '1604833_Atlas_H-VV_3.2fb-1'

        x = x + 1
        S95_t1(x)%id = 2016056
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-056' ! superseded by 1708.09624
        S95_t1(x)%citekey = 'ATLAS:2016bza'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 13.3D0
        S95_t1(x)%xmin = 300.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 5D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '2016056_Atlas_H-ZZ-llnunu_13.3fb-1'

        x = x + 1
        S95_t1(x)%id = 20160551
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-055' ! superseded by 1708.04445
        S95_t1(x)%citekey = 'ATLAS:2016yqq'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 15.5D0
        S95_t1(x)%xmin = 1200.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 100D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '20160551_Atlas_pp-H-ZZ-qqqq_15.5fb-1'

        x = x + 1
        S95_t1(x)%id = 20160552
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2016-055' ! superseded by 1708.04445
        S95_t1(x)%citekey = 'ATLAS:2016yqq'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 15.5D0
        S95_t1(x)%xmin = 1200.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 100D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '20160552_Atlas_pp-H-WW-qqqq_15.5fb-1'

        x = x + 1
        S95_t1(x)%id = 15009
        S95_t1(x)%particle_x = Hneut
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1707.07283 (CMS)' ! 'CMS-PAS-HIG-2015-009'
        S95_t1(x)%citekey = 'Sirunyan:2017uvf'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%xmin = 25.0D0
        S95_t1(x)%xmax = 60.0D0
        S95_t1(x)%sep = 0.1D0
        S95_t1(x)%deltax = 1.0D0
        filename(x) = '15009_CMS_bbA_mumu_19.7fb-1'
!---------------- charged Higgs ------------------

        x = x + 1
        S95_t1(x)%id = 1811
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:0908.1811 (D0)'
        S95_t1(x)%citekey = 'Abazov:2009aa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 1.0D0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 155.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_Hp_qq_1.0fb_1811_interpol'

        x = x + 1
        S95_t1(x)%id = 1269
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = '[hep-ex] arXiv:0907.1269 (CDF) lower mass'
        S95_t1(x)%citekey = 'Aaltonen:2009ke'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 2.2D0
        S95_t1(x)%xmin = 60.0D0
        S95_t1(x)%xmax = 70.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_t-Hplusb_1269_lowmass'

        x = x + 1
        S95_t1(x)%id = 1270
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CDF'
        S95_t1(x)%label = '[hep-ex] arXiv:0907.1269 (CDF) higher mass'
        S95_t1(x)%citekey = 'Aaltonen:2009ke'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 2.2D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'CDF_t-Hplusb_1269_highmass'

        x = x + 1
        S95_t1(x)%id = 1812
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = ' D0'
        S95_t1(x)%label = '[hep-ex] arXiv:0908.1811 (D0)'
        S95_t1(x)%citekey = 'Abazov:2009aa'
        S95_t1(x)%energy = 1.96D0
        S95_t1(x)%lumi = 1.0D0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 155.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = 'D0_Hp_taunu_1.0fb_1811_interpol'

        x = x + 1
        S95_t1(x)%id = 2011094
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2011-094' ! superseded by 1302.3694
        S95_t1(x)%citekey = 'Aad:2013hla'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 0.035D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 130.0D0
        S95_t1(x)%sep = 20.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2011094_Atlas_Hplus-cs_35pb-1'

        x = x + 1
        S95_t1(x)%id = 2760
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1204.2760 (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2012tj'
        S95_t1(x)%energy = 7.0D0
        S95_t1(x)%lumi = 4.6D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2760_Atlas_Hp_taunu_4.6fb-1'

        x = x + 1
        S95_t1(x)%id = 2014050
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2014-050'
        S95_t1(x)%citekey = 'ATLAS:2014zha'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.5D0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '2014050_Atlas_Hp_taunu_19.5fb-1'

        x = x + 1
        S95_t1(x)%id = 79151
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1807.07915 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018gjj'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 2000.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '079151_Atlas_tbHp-taunu_36fb-1'

        x = x + 1
        S95_t1(x)%id = 79152
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1807.07915 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018gjj'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '079152_Atlas_t-Hpb-taunu_36fb-1'

        x = x + 1
        S95_t1(x)%id = 14020
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1508.07774 (CMS)' ! 'CMS-PAS-HIG-14-020'
        S95_t1(x)%citekey = 'Khachatryan:2015qxa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '14020_CMS_Hp_taunu_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 13035
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1510.04252 (CMS)' ! 'CMS-PAS-HIG-13-035'
        S95_t1(x)%citekey = 'Khachatryan:2015uua'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '13035_CMS_chargedH-cs_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 18021
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-18-021'
        S95_t1(x)%citekey = 'CMS:2019sxh'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 8.000000E+01
        S95_t1(x)%xmax = 1.600000E+02
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '18021_CMS_t-H+b-csb_35.9fb-1'
        S95_t1(x)%deltax = 0.000000E+00

        x = x + 1
        S95_t1(x)%id = 16030
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1808.06575 (CMS)' ! 'CMS-PAS-HIG-16-030'
        S95_t1(x)%citekey = 'Sirunyan:2018dvm'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 19.7D0
        S95_t1(x)%xmin = 90.0D0
        S95_t1(x)%xmax = 150.0D0
        S95_t1(x)%sep = 10.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '16030_CMS_t_Hplusb_cb_19.7fb-1'

        x = x + 1
        S95_t1(x)%id = 1504233
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'arXiv:1503.04233 [hep-ex] (ATLAS)'
        S95_t1(x)%citekey = 'Aad:2015nfa'
        S95_t1(x)%energy = 8.0D0
        S95_t1(x)%lumi = 20.3D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 1000.0D0
        S95_t1(x)%sep = 20.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '1504233_Atlas_Hplus_VBF-WZ_20.3fb-1'

! Still leave this in, although strictly speaking superseded by 18-014, because different implementation.
        x = x + 1
        S95_t1(x)%id = 160311
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = 'CMS-PAS-HIG-16-031'
        S95_t1(x)%citekey = 'CMS:2016szv'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 12.9D0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 160.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '160311_CMS_t_Hplusb_taunu_12.9fb-1'

        x = x + 1
        S95_t1(x)%id = 18014
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 1903.04560 (CMS)' ! 'CMS-PAS-HIG-18-014'
        S95_t1(x)%citekey = 'Sirunyan:2019hkq'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 35.9D0
        S95_t1(x)%xmin = 80.0D0
        S95_t1(x)%xmax = 3000.0D0
        S95_t1(x)%sep = 5.0D0
        S95_t1(x)%deltax = 0.0D0
        filename(x) = '18014_CMS_pp-tbHpm-taunu_13TeV_35.9fb-1'

        x = x + 1
        S95_t1(x)%id = 3599
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] arXiv:1808.03599 (ATLAS)'
        S95_t1(x)%citekey = 'Aaboud:2018cwk'
        S95_t1(x)%energy = 13.0D0
        S95_t1(x)%lumi = 36.1D0
        S95_t1(x)%xmin = 200.0D0
        S95_t1(x)%xmax = 2000.0D0
        S95_t1(x)%sep = 25.0D0
        S95_t1(x)%deltax = 10.0D0
        filename(x) = '3599_Atlas_pp-tbHpm-tbtb_13Tev_36.1fb-1'

        x = x + 1
        S95_t1(x)%id = 18015
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'CMS'
        S95_t1(x)%label = '[hep-ex] 2001.07763 (CMS)' ! 'CMS-PAS-HIG-18-015'
        S95_t1(x)%citekey = 'Sirunyan:2020hwv'
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 3.590000E+01
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 3.000000E+03
        S95_t1(x)%sep = 1.000000E+01
        filename(x) = '18015_CMS_tbHpm-tbtb_35.9fb-1'
        S95_t1(x)%deltax = 5.000000E+00

        x = x + 1
        S95_t1(x)%id = 2020039
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = 'ATLAS-CONF-2020-039'
        S95_t1(x)%citekey = "ATLAS:2020jqj"
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.390000E+02
        S95_t1(x)%xmin = 2.000000E+02
        S95_t1(x)%xmax = 2.000000E+03
        S95_t1(x)%sep = 2.500000E+01
        filename(x) = '2020039_ATLAS_Hp_tb_139fb-1'

        x = x + 1
        S95_t1(x)%id = 2002113251
        S95_t1(x)%particle_x = Hplus
        S95_t1(x)%expt = 'ATL'
        S95_t1(x)%label = '[hep-ex] 2002.11325 (ATLAS)'
        S95_t1(x)%citekey = "Aad:2020kep"
        S95_t1(x)%energy = 1.300000E+01
        S95_t1(x)%lumi = 1.390000E+02
        S95_t1(x)%xmin = 6D+02
        S95_t1(x)%xmax = 2D+03
        S95_t1(x)%sep = 2D2
        filename(x) = '200211325_ATLAS_dijet+lep_139fb-1_tbHp'
!----------------------------------------------

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        ! checks we've filled the whole array
        if (x .ne. xend) then
            print *, "found", x, "but expected", xend, "type 1 limits."
            stop 'error in initializetables1 (a)'
        end if

        ! do loop to read in S95 tables
        col = 3
        do x = xbeg, xend
            S95_t1(x)%nx = nint((S95_t1(x)%xmax - S95_t1(x)%xmin)/S95_t1(x)%sep) + 1
            allocate (S95_t1(x)%dat(S95_t1(x)%nx, col - 1))
        end do

        open (file_id_common2, file=trim(adjustl(pathname))//'Expt_tables/'// &
              'S95_t1.binary', form='unformatted')

        read (file_id_common2, iostat=ios) S95_t1(xbeg)%dat

        if (ios .eq. 0) then

            do x = xbeg + 1, xend
                read (file_id_common2) S95_t1(x)%dat
            end do

        else
            rewind (file_id_common2)
            do x = xbeg, xend
                fullfilename = trim(adjustl(pathname))//'Expt_tables/' &
                               //trim(adjustl(S95_t1(x)%expt))//'tables/' &
                               //trim(filename(x))//'.txt'

                call read_tabletype1(S95_t1(x), 5, col, fullfilename)
#ifndef WEBVERSION
                write (file_id_common2) S95_t1(x)%dat
#endif
            end do
        end if

        close (file_id_common2)

        deallocate (filename)

    end subroutine initializetables1

    !************************************************************
    subroutine read_tabletype1(t1, skip, col, fullfilename)
        !************************************************************
        !fills t1%dat
        !--------------------------------------input
        type(table1) :: t1
        integer :: skip, col
        character(LEN=*) :: fullfilename
        !-----------------------------------internal
        integer :: i, n
        double precision :: xdummy, xdummy_store
        !-------------------------------------------

        t1%dat = 0.0D0

        open (file_id_1, file=(trim(fullfilename)))

        do i = 1, skip
            read (file_id_1, *) !skip lines
        end do

        xdummy_store = t1%xmin - t1%sep
        do i = 1, t1%nx
            read (file_id_1, *) xdummy, (t1%dat(i, n), n=1, col - 1)
!    if(minval(t1%dat(i,:)).lt.0) then
!     write(*,*) xdummy,minval(t1%dat(i,:))
!    endif
            ! checks that x are evenly spaced as expected
            if ((abs(xdummy - xdummy_store - t1%sep) .gt. 1.0D-7) &
                .or. (abs(xdummy - (t1%xmin + dble(i - 1)*t1%sep)) .gt. 1.0D-7)) then
                write (*, *) i, t1%id, xdummy, t1%xmin + dble(i - 1)*t1%sep
                stop 'error in read_tabletype1 (a1)'
            end if

            xdummy_store = xdummy

        end do

        if (abs(xdummy - t1%xmax) .gt. 1.0D-7) stop 'error in read_tabletype1 (a2)'

        close (file_id_1)

    end subroutine read_tabletype1
    !************************************************************
end module S95tables_type1
!************************************************************
