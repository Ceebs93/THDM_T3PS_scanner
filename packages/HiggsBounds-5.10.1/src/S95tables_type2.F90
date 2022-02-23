! This file is part of HiggsBounds
!  -KW
!******************************************************************
module S95tables_type2
!******************************************************************
    use S95tables_type1
    implicit none
!#define fast

    !table type 2------------------------------
    type table2
        logical :: needs_M2_gt_2M1
        integer :: id, nx1, nx2, particle_x1, particle_x2 !see usefulbits.f90 for key to particle codes n.b. they're NOT pdg
        character(LEN=45) :: label
        character(LEN=35) :: citekey
        character(LEN=3) :: expt
        double precision :: lumi, energy
        double precision :: xmax1, xmin1, xmax2, xmin2, sep1, sep2, deltax
        double precision, allocatable :: dat(:, :, :) !in dat(a,b,1:2)...obs=1,pred=2. 1st component of dat is y, 2nd is x
        double precision :: maxdatval
        double precision :: z !only used in slices_t2
    end type

    integer, parameter :: file_id_2_exp = 10  !same as file_id_common in usefulbits.f90
    integer, parameter :: file_id_2_obs = 11

    !------------------------------------------

contains

    !************************************************************
    subroutine initializetables_type2_blank(tablet2)
        !***********************************************************
        ! still leaves dat unallocated
        integer:: i
        type(table2) :: tablet2(:)

        do i = lbound(tablet2, dim=1), ubound(tablet2, dim=1)
            tablet2(i)%id = -1
            tablet2(i)%nx1 = -1
            tablet2(i)%nx2 = -1
            tablet2(i)%particle_x1 = -1
            tablet2(i)%particle_x2 = -1
            tablet2(i)%label = ''
            tablet2(i)%citekey = ''
            tablet2(i)%expt = ''
            tablet2(i)%lumi = -1.0D0
            tablet2(i)%energy = -1.0D0
            tablet2(i)%xmax1 = -1.0D0
            tablet2(i)%xmax2 = -1.0D0
            tablet2(i)%xmin1 = -1.0D0
            tablet2(i)%xmin2 = -1.0D0
            tablet2(i)%sep1 = -1.0D0
            tablet2(i)%sep2 = -1.0D0
            tablet2(i)%deltax = -1.0D0
            tablet2(i)%maxdatval = -1.0D0

            tablet2(i)%z = -1.0D9 !only used in slices_t2

            tablet2(i)%needs_M2_gt_2M1 = .False.
        end do

    end subroutine initializetables_type2_blank

    !***********************************************************
    subroutine initializetables2(S95_t2)
        !***********************************************************
        ! fills S95_t2
        !***********************************************************
        use install_data, only: pathname
        use usefulbits, only: Hneut, Hplus, Chineut, Chiplus, anyH, &
                              small, file_id_common2, not_a_particle
        implicit none

        !--------------------------------------input
        type(table2) :: S95_t2(:)
        !-----------------------------------internal
        integer :: i, tno, j, x, xbeg, xend, k, ios
        character(LEN=2) :: tableno
        character(len=100), allocatable :: filename(:)
        double precision :: dummy
        double precision, allocatable :: testrow(:)
        integer :: file_id_arr(2)
        double precision :: maxdatval
        !-------------------------------------------
        file_id_arr(1) = file_id_2_exp
        file_id_arr(2) = file_id_2_obs

        xbeg = lbound(S95_t2, dim=1)
        xend = ubound(S95_t2, dim=1)

        allocate (filename(xbeg:xend))
        x = xbeg - 1

        tno = 14
        do i = 1, 8
            x = x + 1
            tno = tno + 1
            if ((x .eq. 3) .or. (x .eq. 7)) tno = tno + 1
            write (tableno, '(I2)') tno

            S95_t2(x)%id = tno*10
            S95_t2(x)%expt = 'LEP'
            S95_t2(x)%energy = 0.208D0
            S95_t2(x)%deltax = 0.0D0
            S95_t2(x)%particle_x1 = Hneut
            S95_t2(x)%particle_x2 = Hneut

            select case (S95_t2(x)%id)
            case (220, 230, 240)
                S95_t2(x)%label = 'hep-ex/0602042 (LEP)'
            case default
                S95_t2(x)%label = 'hep-ex/0602042, table '//tableno//' (LEP)'
            end select
            S95_t2(x)%citekey = 'Schael:2006cr'

            S95_t2(x)%sep1 = 1.0D0
            S95_t2(x)%sep2 = 1.0D0
            S95_t2(x)%maxdatval = 1.0D2
            !S95_t2(x)%OBid=x+2

            select case (S95_t2(x)%id)
            case (150, 160, 220)
                S95_t2(x)%xmin1 = 1.0D0
                S95_t2(x)%xmax1 = 60.0D0
                S95_t2(x)%xmin2 = 2.0D0
                S95_t2(x)%xmax2 = 120.0D0
                S95_t2(x)%needs_M2_gt_2M1 = .True.
            case (180, 190, 230, 240)
                S95_t2(x)%xmin1 = 1.0D0
                S95_t2(x)%xmax1 = 180.0D0
                S95_t2(x)%xmin2 = 1.0D0
                S95_t2(x)%xmax2 = 180.0D0
                S95_t2(x)%needs_M2_gt_2M1 = .False.
            case (200, 210)
                S95_t2(x)%xmin1 = 1.0D0
                S95_t2(x)%xmax1 = 90.0D0
                S95_t2(x)%xmin2 = 2.0D0
                S95_t2(x)%xmax2 = 180.0D0
                S95_t2(x)%needs_M2_gt_2M1 = .True.
            case default
                write (*, *) 'error in initializetables2 (a)'
                stop
            end select

            filename(x) = 'table'//tableno//'full'
        end do

        do i = 5, 10
            x = x + 1
            tno = i
            write (tableno, '(I2)') tno

            S95_t2(x)%id = 900 + tno
            S95_t2(x)%expt = 'LEP'
            S95_t2(x)%energy = 0.208D0
            S95_t2(x)%deltax = 0.0D0
            S95_t2(x)%label = 'hep-ex/0401026, fig '//trim(adjustl(tableno))//' (OPAL)'
            S95_t2(x)%citekey = 'Abbiendi:2003sc'
            S95_t2(x)%sep1 = 1.0D0
            S95_t2(x)%sep2 = 1.0D0
            S95_t2(x)%maxdatval = 1.0D6 !these tables are in fb

            select case (tno)
            case (5, 6, 7, 8)
                S95_t2(x)%xmin1 = 0.0D0
                S95_t2(x)%xmax1 = 100.0D0
                S95_t2(x)%xmin2 = 75.0D0
                S95_t2(x)%xmax2 = 120.0D0
                S95_t2(x)%needs_M2_gt_2M1 = .False.
                S95_t2(x)%particle_x1 = Chineut
                S95_t2(x)%particle_x2 = Chiplus
            case (9, 10)
                S95_t2(x)%xmin1 = 0.0D0
                S95_t2(x)%xmax1 = 100.0D0
                S95_t2(x)%xmin2 = 50.0D0
                S95_t2(x)%xmax2 = 200.0D0
                S95_t2(x)%needs_M2_gt_2M1 = .False.
                S95_t2(x)%particle_x1 = Chineut
                S95_t2(x)%particle_x2 = Chineut
            case default
                write (*, *) 'error in initializetables2 (b)'
                stop
            end select

            filename(x) = '1026_fig'//trim(adjustl(tableno))
        end do

        x = x + 1
        S95_t2(x)%id = 6065
        S95_t2(x)%expt = 'LEP'
        S95_t2(x)%energy = 0.208D0
        S95_t2(x)%deltax = 0.0D0
        S95_t2(x)%label = '[hep-ex] arXiv:1301.6065 (LEP)'
        S95_t2(x)%citekey = 'Abbiendi:2013hk'
        S95_t2(x)%sep1 = 0.1D0
        S95_t2(x)%sep2 = 0.5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = 1.0D0
        S95_t2(x)%xmin2 = 43.0D0
        S95_t2(x)%xmax2 = 95.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hplus
        filename(x) = "6065_LEP_HpHm_fig4"

        x = x + 1
        S95_t2(x)%id = 02671
        S95_t2(x)%expt = 'LEP'
        S95_t2(x)%energy = 0.208D0
        S95_t2(x)%deltax = 0.0D0
        S95_t2(x)%label = '[hep-ex] arXiv:0812.0267, Fig.10a (OPAL,LEP)'
        S95_t2(x)%citekey = 'Abbiendi:2008aa'
        S95_t2(x)%sep1 = 0.5D0
        S95_t2(x)%sep2 = 0.5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 12.0D0
        S95_t2(x)%xmax1 = 90.0D0
        S95_t2(x)%xmin2 = 40.0D0
        S95_t2(x)%xmax2 = 93.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hplus
        filename(x) = "OPAL_H+H-_AWAW_0267"

        x = x + 1
        S95_t2(x)%id = 02672
        S95_t2(x)%expt = 'LEP'
        S95_t2(x)%energy = 0.208D0
        S95_t2(x)%deltax = 0.0D0
        S95_t2(x)%label = '[hep-ex] arXiv:0812.0267, Fig.10b (OPAL,LEP)'
        S95_t2(x)%citekey = 'Abbiendi:2008aa'
        S95_t2(x)%sep1 = 0.5D0
        S95_t2(x)%sep2 = 0.5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 12.0D0
        S95_t2(x)%xmax1 = 77.0D0
        S95_t2(x)%xmin2 = 40.0D0
        S95_t2(x)%xmax2 = 80.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hplus
        filename(x) = "OPAL_H+H-_AWtaunu_0267"

        x = x + 1
        S95_t2(x)%id = 3381
        S95_t2(x)%expt = ' D0'
        S95_t2(x)%energy = 1.96D0
        S95_t2(x)%lumi = 4.2
        S95_t2(x)%deltax = 0.0D0
        S95_t2(x)%label = '[hep-ex] arXiv:0905.3381, table I (D0)'
        S95_t2(x)%citekey = 'Abazov:2009yi'
        S95_t2(x)%sep1 = 0.1D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.2D0
        S95_t2(x)%xmax1 = 3.0D0
        S95_t2(x)%xmin2 = 80.0D0
        S95_t2(x)%xmax2 = 200.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = "D0_h-aa-mumumumu_3381"

        x = x + 1
        S95_t2(x)%id = 3382
        S95_t2(x)%expt = ' D0'
        S95_t2(x)%energy = 1.96D0
        S95_t2(x)%lumi = 4.2
        S95_t2(x)%deltax = 0.0D0
        S95_t2(x)%label = '[hep-ex] arXiv:0905.3381, table II (D0)'
        S95_t2(x)%citekey = 'Abazov:2009yi'
        S95_t2(x)%sep1 = 0.2D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 3.6D0
        S95_t2(x)%xmax1 = 19.0D0
        S95_t2(x)%xmin2 = 85.0D0
        S95_t2(x)%xmax2 = 200.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = "D0_h-aa-tautaumumu_3381"

        x = x + 1
        S95_t2(x)%id = 6227
        S95_t2(x)%expt = ' D0'
        S95_t2(x)%label = 'D0 Note 6227'
        S95_t2(x)%citekey = ''
        S95_t2(x)%energy = 1.96D0
        S95_t2(x)%lumi = 7.3
        S95_t2(x)%sep1 = 0.04D0
        S95_t2(x)%sep2 = 10.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.06D0
        S95_t2(x)%xmax1 = 0.18D0
        S95_t2(x)%xmin2 = 90.0D0
        S95_t2(x)%xmax2 = 300.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = 'D0_h-bb_h-tautau_comb_5.2-7.3fb_6227'

        x = x + 1
        S95_t2(x)%id = 5053
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:1406.5053 (ATLAS)'
        S95_t2(x)%citekey = 'Aad:2014yja'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 20.3
        S95_t2(x)%sep1 = 5.0D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 120.0D0
        S95_t2(x)%xmax1 = 130.0D0
        S95_t2(x)%xmin2 = 260.0D0
        S95_t2(x)%xmax2 = 500.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '5053_H-hh-gagabb_20fb-1'

        x = x + 1
        S95_t2(x)%id = 06896
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = '[hep-ex] arXiv:1603.06896 (CMS)'
        S95_t2(x)%citekey = 'Khachatryan:2016sey'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 19.7
        S95_t2(x)%sep1 = 5.0D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 120.0D0
        S95_t2(x)%xmax1 = 130.0D0
        S95_t2(x)%xmin2 = 260.0D0
        S95_t2(x)%xmax2 = 1100.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '06896_H-hh-gagabb_19.7fb-1'

        x = x + 1
        S95_t2(x)%id = 16032
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'CMS-PAS-HIG-16-032'
        S95_t2(x)%citekey = 'CMS:2016vpz'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 2.70
        S95_t2(x)%sep1 = 2.5D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 122.5D0
        S95_t2(x)%xmax1 = 127.5D0
        S95_t2(x)%xmin2 = 250.0D0
        S95_t2(x)%xmax2 = 900.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '16032_CMS_H-hh-gagabb_2.70fb-1'

        x = x + 1
        S95_t2(x)%id = 180406174
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:1804.06174 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2018knk'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 36.1D0
        S95_t2(x)%sep1 = 30.D0
        S95_t2(x)%sep2 = 20.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 110.0D0
        S95_t2(x)%xmax1 = 140.0D0
        S95_t2(x)%xmin2 = 260.0D0
        S95_t2(x)%xmax2 = 3000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '180406174_ATLAS_H-hh-bbbb_36.1fb-1'

        x = x + 1
        S95_t2(x)%id = 200105178
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:2001.05178 (ATLAS)'
        S95_t2(x)%citekey = 'Aad:2020kub'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 120D0
        S95_t2(x)%sep1 = 30.D0
        S95_t2(x)%sep2 = 20.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 110.0D0
        S95_t2(x)%xmax1 = 140.0D0
        S95_t2(x)%xmin2 = 260.0D0
        S95_t2(x)%xmax2 = 1000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '200105178_ATLAS_vbf_H-hh-bbbb_136fb-1'

        x = x + 1
        S95_t2(x)%id = 14022
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'CMS-PAS-HIG-14-022'
        S95_t2(x)%citekey = 'CMS:2015iga'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 19.7
        S95_t2(x)%sep1 = 1.0D0
        S95_t2(x)%sep2 = 2.5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 5.0D0
        S95_t2(x)%xmax1 = 15.0D0
        S95_t2(x)%xmin2 = 122.5D0
        S95_t2(x)%xmax2 = 127.5D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '14022_H-hh-tautau_19.7fb-1'

        x = x + 1
        S95_t2(x)%id = 150600424
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'arXiv:1506.00424 (CMS)'
        S95_t2(x)%citekey = 'Khachatryan:2015wka'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 20.7
        S95_t2(x)%sep1 = 0.05D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.25D0
        S95_t2(x)%xmax1 = 3.55D0
        S95_t2(x)%xmin2 = 86.0D0
        S95_t2(x)%xmax2 = 150.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '150600424_CMS_H-aa-mumu'

        x = x + 1
        S95_t2(x)%id = 16035
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'CMS-PAS-HIG-16-035' ! superseded by 1812.00380
        S95_t2(x)%citekey = 'CMS:2016tgd'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 2.8
        S95_t2(x)%sep1 = 0.05D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.25D0
        S95_t2(x)%xmax1 = 3.55D0
        S95_t2(x)%xmin2 = 86.0D0
        S95_t2(x)%xmax2 = 150.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .True.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '16035_CMS_H-aa-mumumu_2.8fb-1'

        x = x + 1
        S95_t2(x)%id = 02301
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = '[hep-ex] arXiv:1506.02301 (CMS)'
        S95_t2(x)%citekey = 'Khachatryan:2015qba'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 19.7
        S95_t2(x)%sep1 = 0.01D0
        S95_t2(x)%sep2 = 10.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.00D0
        S95_t2(x)%xmax1 = 0.10D0
        S95_t2(x)%xmin2 = 150.0D0
        S95_t2(x)%xmax2 = 840.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '14006_CMS_h-gaga_2D'

        x = x + 1
        S95_t2(x)%id = 003892
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:1509.00389 (ATLAS)'
        S95_t2(x)%citekey = 'Aad:2015agg'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 20.3
        S95_t2(x)%sep1 = 0.05D0
        S95_t2(x)%sep2 = 100.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.20D0
        S95_t2(x)%xmax1 = 0.80D0
        S95_t2(x)%xmin2 = 300.0D0
        S95_t2(x)%xmax2 = 1000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '00389_H-WW_20.3fb-1_2D'

        x = x + 1
        S95_t2(x)%id = 20160791
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = 'ATLAS-CONF-2016-079' ! superseded by 1712.06386
        S95_t2(x)%citekey = 'ATLAS:2016oum'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 14.8
        S95_t2(x)%sep1 = 0.01D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = 0.10D0
        S95_t2(x)%xmin2 = 400.0D0
        S95_t2(x)%xmax2 = 1000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '20160791_gg-H-ZZ-4l_14.8fb-1'

        x = x + 1
        S95_t2(x)%id = 01123
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = 'arXiv:1710.01123 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2017gsl'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 36.1
        S95_t2(x)%sep1 = 5.0D0
        S95_t2(x)%sep2 = 100.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 5.0D0
        S95_t2(x)%xmax1 = 15.0D0
        S95_t2(x)%xmin2 = 200.0D0
        S95_t2(x)%xmax2 = 4000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '01123_ATLAS_gg-H-WW-e_nu_mu_nu_36.1fb-1'

        x = x + 1
        S95_t2(x)%id = 06386
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = 'arXiv:1712.06386 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2017rel'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 36.1
        S95_t2(x)%sep1 = 0.01D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = .10D0
        S95_t2(x)%xmin2 = 400.0D0
        S95_t2(x)%xmax2 = 1000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '06386_ATLAS_gg-H-ZZ-4l+2l2nu_36.1fb-1'

        x = x + 1
        S95_t2(x)%id = 170121
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = '[hep-ex] 1804.01939 (CMS)' ! 'CMS-PAS-HIG 17-012'
        S95_t2(x)%citekey = 'Sirunyan:2018qlb'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 35.9
        S95_t2(x)%sep1 = 10.0D0
        S95_t2(x)%sep2 = 10.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = 100.0D0
        S95_t2(x)%xmin2 = 130.0D0
        S95_t2(x)%xmax2 = 3000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '17012_CMS_gg-H-ZZ_35.9fb-1'

        x = x + 1
        S95_t2(x)%id = 160332
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'CMS-PAS-HIG-16-033'
        S95_t2(x)%citekey = 'CMS:2016ilx'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 12.9
        S95_t2(x)%sep1 = 2.0D0
        S95_t2(x)%sep2 = 0.1D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = 40.0D0
        S95_t2(x)%xmin2 = 130.0D0
        S95_t2(x)%xmax2 = 2520.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '160332_CMS_H-ZZ-4l_VBF_12.9fb-1'

        x = x + 1
        S95_t2(x)%id = 110282
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = 'arXiv:1811.11028 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2018ksn'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 36.1
        S95_t2(x)%sep1 = 5.0D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 135.0D0
        S95_t2(x)%xmax1 = 165.0D0
        S95_t2(x)%xmin2 = 280.0D0
        S95_t2(x)%xmax2 = 340.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '11028_ATLAS_X-SS-WWWW_13TeV_36.1fb-1'

        x = x + 1
        S95_t2(x)%id = 150011
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'arXiv:1603.02991 (CMS)' ! CMS-PAS-HIG-15-001
        S95_t2(x)%citekey = 'Khachatryan:2016are'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 19.8
        S95_t2(x)%sep1 = 5.0D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 15.0D0
        S95_t2(x)%xmax1 = 1000.0D0
        S95_t2(x)%xmin2 = 15.0D0
        S95_t2(x)%xmax2 = 1000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '150011_H2-ZH1-lltautau_19.8fb-1'

        x = x + 1
        S95_t2(x)%id = 150012
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = 'arXiv:1603.02991 (CMS)' ! CMS-PAS-HIG-15-001
        S95_t2(x)%citekey = 'Khachatryan:2016are'
        S95_t2(x)%energy = 8.0D0
        S95_t2(x)%lumi = 19.8
        S95_t2(x)%sep1 = 2.4D0
        S95_t2(x)%sep2 = 2.4D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 24.0D0
        S95_t2(x)%xmax1 = 1200.0D0
        S95_t2(x)%xmin2 = 24.0D0
        S95_t2(x)%xmax2 = 1200.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '150012_H2-ZH1-llbb_19.8fb-1'

        x = x + 1
        S95_t2(x)%id = 18012
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = '[hep-ex] 1911.03781 (CMS)' ! 'CMS-PAS-HIG-18-012'
        S95_t2(x)%citekey = 'Sirunyan:2019wrn'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 35.9
        S95_t2(x)%sep1 = 5.0D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 30.0D0
        S95_t2(x)%xmax1 = 770.0D0
        S95_t2(x)%xmin2 = 135.0D0
        S95_t2(x)%xmax2 = 995.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '18012_CMS_H2-ZH1-llbb_35.9fb-1'

        x = x + 1
        S95_t2(x)%id = 170271
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = '[hep-ex] 1908.01115' ! 'CMS-PAS-HIG-17-027'
        S95_t2(x)%citekey = 'Sirunyan:2019wph'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 35.9
        S95_t2(x)%sep1 = 0.5D0
        S95_t2(x)%sep2 = 5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.50D0
        S95_t2(x)%xmax1 = 25.0D0
        S95_t2(x)%xmin2 = 400.0D0
        S95_t2(x)%xmax2 = 750.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '170271_CMS_A-tt_35.9fb-1'

        x = x + 1
        S95_t2(x)%id = 170272
        S95_t2(x)%expt = 'CMS'
        S95_t2(x)%label = '[hep-ex] 1908.01115' ! 'CMS-PAS-HIG-17-027'
        S95_t2(x)%citekey = 'Sirunyan:2019wph'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 35.9
        S95_t2(x)%sep1 = 0.5D0
        S95_t2(x)%sep2 = 5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.50D0
        S95_t2(x)%xmax1 = 25.0D0
        S95_t2(x)%xmin2 = 400.0D0
        S95_t2(x)%xmax2 = 750.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '170272_CMS_H-tt_35.9fb-1'

        x = x + 1
        S95_t2(x)%id = 4147
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:1707.04147 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2017yyg'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 36.7
        S95_t2(x)%sep1 = 0.02D0
        S95_t2(x)%sep2 = 0.5D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.00D0
        S95_t2(x)%xmax1 = 0.10D0
        S95_t2(x)%xmin2 = 200.0D0
        S95_t2(x)%xmax2 = 2700.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '4147_Atlas_pp-H-gaga_36.7fb-1'

! TS HB beyond Higgs implementation:
        x = x + 1
        S95_t2(x)%id = 109171
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:1901.10917 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2019zxd'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 79.8D0
        S95_t2(x)%sep1 = 0.01D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.00D0
        S95_t2(x)%xmax1 = 0.15D0
        S95_t2(x)%xmin2 = 250.0D0
        S95_t2(x)%xmax2 = 1050.0D0
        S95_t2(x)%deltax = 20.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = anyH
        filename(x) = '10917_ATL_dijet_photon_FI_76p6fb-1'

        x = x + 1
        S95_t2(x)%id = 109172
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] arXiv:1901.10917 (ATLAS)'
        S95_t2(x)%citekey = 'Aaboud:2019zxd'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 79.8D0
        S95_t2(x)%sep1 = 0.01D0
        S95_t2(x)%sep2 = 1.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.00D0
        S95_t2(x)%xmax1 = 0.15D0
        S95_t2(x)%xmin2 = 250.0D0
        S95_t2(x)%xmax2 = 1050.0D0
        S95_t2(x)%deltax = 20.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '10917_ATL_dijet_photon_2b_76p6fb-1'

! By DDD, Nov 8 2018
        x = x + 1
        S95_t2(x)%id = 20160341
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] 1804.01126 (ATLAS), bbh'
        S95_t2(x)%citekey = 'Aaboud:2018eoy'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 26.1
        S95_t2(x)%sep1 = 10.0D0
        S95_t2(x)%sep2 = 10.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 130.0D0
        S95_t2(x)%xmax1 = 700.0D0
        S95_t2(x)%xmin2 = 230.0D0
        S95_t2(x)%xmax2 = 800.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '20160341_ATLAS_A_HZ_llbb_36.1fb-1_bbh'

! By DDD, Nov 8 2018
        x = x + 1
        S95_t2(x)%id = 20160342
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] 1804.01126 (ATLAS), ggf'
        S95_t2(x)%citekey = 'Aaboud:2018eoy'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 36.1
        S95_t2(x)%sep1 = 10.0D0
        S95_t2(x)%sep2 = 10.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 130.0D0
        S95_t2(x)%xmax1 = 700.0D0
        S95_t2(x)%xmin2 = 230.0D0
        S95_t2(x)%xmax2 = 800.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = Hneut
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = '20160342_ATLAS_A_HZ_llbb_36.1fb-1_ggf'

        x = x + 1
        S95_t2(x)%id = 1903062481
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] 1903.06248 (ATLAS) mumu'
        S95_t2(x)%citekey = 'Aad:2019fac'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 139.D0
        S95_t2(x)%sep1 = 2.5D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = 10.D0
        S95_t2(x)%xmin2 = 250.0D0
        S95_t2(x)%xmax2 = 3000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = "190306248_ATLAS_H-mumu_139fb-1"

        x = x + 1
        S95_t2(x)%id = 1903062482
        S95_t2(x)%expt = 'ATL'
        S95_t2(x)%label = '[hep-ex] 1903.06248 (ATLAS) ee'
        S95_t2(x)%citekey = 'Aad:2019fac'
        S95_t2(x)%energy = 13.0D0
        S95_t2(x)%lumi = 139.D0
        S95_t2(x)%sep1 = 2.5D0
        S95_t2(x)%sep2 = 5.0D0
        S95_t2(x)%maxdatval = 1.0D6
        S95_t2(x)%xmin1 = 0.0D0
        S95_t2(x)%xmax1 = 10.D0
        S95_t2(x)%xmin2 = 250.0D0
        S95_t2(x)%xmax2 = 3000.0D0
        S95_t2(x)%needs_M2_gt_2M1 = .False.
        S95_t2(x)%particle_x1 = not_a_particle
        S95_t2(x)%particle_x2 = Hneut
        filename(x) = "190306248_ATLAS_H-ee_139fb-1"

        ! checks we've filled the whole array
        if (x .ne. xend) then
            write (*, *) 'error in initializetables2 (c)', x, xend
            stop
        end if

        ! read in the tables
        do x = xbeg, xend

            S95_t2(x)%nx2 = nint((S95_t2(x)%xmax2 - S95_t2(x)%xmin2)/S95_t2(x)%sep2) + 1
            S95_t2(x)%nx1 = nint((S95_t2(x)%xmax1 - S95_t2(x)%xmin1)/S95_t2(x)%sep1) + 1

            allocate (S95_t2(x)%dat(S95_t2(x)%nx2, S95_t2(x)%nx1, 2))
        end do

        ! read in the tables
        open (file_id_common2, file=trim(adjustl(pathname))//'Expt_tables/'// &
              'S95_t2.binary', form='unformatted')

        read (file_id_common2, iostat=ios) S95_t2(xbeg)%dat
        if (ios .eq. 0) then
            do x = xbeg + 1, xend
                read (file_id_common2) S95_t2(x)%dat
            end do

        else

            do x = xbeg, xend
                open (file_id_2_exp, file=trim(adjustl(pathname))//('Expt_tables/' &
                                                                    //trim(adjustl(S95_t2(x)%expt))//'tables/' &
                                                                    //trim(adjustl(S95_t2(x)%expt))//'tables2/' &
                                                                    //trim(adjustl(filename(x)))//'_pred.txt'))
                open (file_id_2_obs, file=trim(adjustl(pathname))//('Expt_tables/' &
                                                                    //trim(adjustl(S95_t2(x)%expt))//'tables/' &
                                                                    //trim(adjustl(S95_t2(x)%expt))//'tables2/' &
                                                                    //trim(adjustl(filename(x)))//'_obs.txt'))

                ! fill S95 from file
                ! row 0 and column 0 in LEP file contain higgs masses
                ! and (0,0) ie top left set to -100
                ! so avoid them
                allocate (testrow(0:S95_t2(x)%nx1))

                do k = lbound(file_id_arr, dim=1), ubound(file_id_arr, dim=1)
                    read (file_id_arr(k), *) (testrow(i), i=0, S95_t2(x)%nx1)
                    if ((testrow(0) + 100.0D0) .gt. small) stop 'error in initializetables2 (d)' !top left number should be -100
                    do i = 1, S95_t2(x)%nx1
                        if (abs(testrow(i) - (S95_t2(x)%xmin1 + dble(i - 1)*S95_t2(x)%sep1)) .gt. small*S95_t2(x)%sep1) then
                            write (*, *) S95_t2(x)%id, testrow(i), (S95_t2(x)%xmin1 + dble(i - 1)*S95_t2(x)%sep1)
                            stop 'error in initializetables2 (e)'
                        end if
                    end do
                end do
                deallocate (testrow)

                do j = 1, S95_t2(x)%nx2
                    read (file_id_2_exp, *) dummy, (S95_t2(x)%dat(j, i, 2), i=1, S95_t2(x)%nx1)
                    if (abs(dummy - (S95_t2(x)%xmin2 + dble(j - 1)*S95_t2(x)%sep2)) .gt. small*S95_t2(x)%sep2) then
                        write (*, *) S95_t2(x)%nx1, dummy, S95_t2(x)%dat(j, :, 2)
                        write (*, *) j, S95_t2(x)%nx2, S95_t2(x)%xmin2 + dble(j - 1)*S95_t2(x)%sep2
                        stop 'error in initializetables2 (f)'
                    end if
                    read (file_id_2_obs, *) dummy, (S95_t2(x)%dat(j, i, 1), i=1, S95_t2(x)%nx1)
                    if (abs(dummy - (S95_t2(x)%xmin2 + dble(j - 1)*S95_t2(x)%sep2)) .gt. small*S95_t2(x)%sep2) then
                        write (*, *) "Problematic analysis: ", S95_t2(x)%id
                        stop 'error in initializetables2 (g)'
                    end if
                end do

                maxdatval = S95_t2(x)%maxdatval
                if (maxdatval .gt. 0.0D0) then
                    ! set entries .ge. S95_t2(x)%maxdatval to (-4): they will not be relevent
                    where (S95_t2(x)%dat .ge. maxdatval) S95_t2(x)%dat = -4.0D0
                end if

                close (file_id_2_exp)
                close (file_id_2_obs)
            end do

            rewind (file_id_common2)
#ifndef WEBVERSION
            do x = xbeg, xend
                write (file_id_common2) S95_t2(x)%dat
            end do
#endif

        end if

        close (file_id_common2)
        deallocate (filename)

    end subroutine initializetables2
    !***********************************************************
    function t2elementnumberfromid(t2, id)
        !--------------------------------------input
        type(table2), intent(in) :: t2(:)
        integer, intent(in) :: id
        !-----------------------------------function
        integer :: t2elementnumberfromid
        !-----------------------------------internal
        integer :: n, x
        !-------------------------------------------

        n = 0
        do x = lbound(t2, dim=1), ubound(t2, dim=1)
            if (t2(x)%id .eq. id) then
                n = n + 1
                t2elementnumberfromid = x
            end if
        end do

        if (n .ne. 1) stop 'problem in function t2elementnumberfromid 1'

    end function t2elementnumberfromid
    !***********************************************************
    subroutine fill_slices_t1_from_slices_of_t2(t2, v1orv2, xy_selection, ftype_selection, slices_t1)
        ! if this subroutine is used,
        ! don't forget to deallocate slices_t1(x)%dat at some point
        !***********************************************************
        implicit none
        !--------------------------------------input
        type(table2), intent(in) :: t2
        integer, intent(in) :: v1orv2
        integer, intent(in) :: xy_selection(:)
        integer, intent(in) :: ftype_selection(:)
        !-------------------------------------output
        type(table1) :: slices_t1(:)  !i.e. 2 slices
        !-----------------------------------internal
        integer :: i, j, k, n
        integer :: n_ftype_selection
        !-------------------------------------------
        n_ftype_selection = ubound(ftype_selection, dim=1)

        do n = lbound(ftype_selection, dim=1), n_ftype_selection
            if (ftype_selection(n) .lt. lbound(t2%dat, dim=3)) stop 'problem in fill_slices_t1_from_slices_of_t2 3a'
            if (ftype_selection(n) .gt. ubound(t2%dat, dim=3)) stop 'problem in fill_slices_t1_from_slices_of_t2 3b'
        end do

        if (lbound(xy_selection, dim=1) .ne. lbound(slices_t1, dim=1)) then
            stop 'problem in fill_slices_t1_from_slices_of_t2 1a'
        end if
        if (ubound(xy_selection, dim=1) .ne. ubound(slices_t1, dim=1)) then
            stop 'problem in fill_slices_t1_from_slices_of_t2 1b'
        end if

        select case (v1orv2)
        case (1)

            do n = lbound(slices_t1, dim=1), ubound(slices_t1, dim=1)

                if (xy_selection(n) .lt. lbound(t2%dat, dim=1)) stop 'problem in fill_slices_t1_from_slices_of_t2 4a'
                if (xy_selection(n) .gt. ubound(t2%dat, dim=1)) stop 'problem in fill_slices_t1_from_slices_of_t2 4b'

                slices_t1(n)%id = t2%id
                slices_t1(n)%nx = t2%nx1
                slices_t1(n)%xmax = t2%xmax1
                slices_t1(n)%xmin = t2%xmin1
                slices_t1(n)%sep = t2%sep1
                slices_t1(n)%deltax = t2%deltax

                allocate (slices_t1(n)%dat(slices_t1(n)%nx, n_ftype_selection))
                slices_t1(n)%dat = -1.0D0

                do i = 1, slices_t1(n)%nx
                    do k = 1, n_ftype_selection
                        slices_t1(n)%dat(i, k) = t2%dat(xy_selection(n), i, ftype_selection(k))
                    end do
                end do

            end do
        case (2)

            do n = lbound(slices_t1, dim=1), ubound(slices_t1, dim=1)

                if (xy_selection(n) .lt. lbound(t2%dat, dim=2)) stop 'problem in fill_slices_t1_from_slices_of_t2 4aa'
                if (xy_selection(n) .gt. ubound(t2%dat, dim=2)) stop 'problem in fill_slices_t1_from_slices_of_t2 4bb'

                slices_t1(n)%id = t2%id
                slices_t1(n)%nx = t2%nx2
                slices_t1(n)%xmax = t2%xmax2
                slices_t1(n)%xmin = t2%xmin2
                slices_t1(n)%sep = t2%sep2
                slices_t1(n)%deltax = t2%deltax

                allocate (slices_t1(n)%dat(slices_t1(n)%nx, n_ftype_selection))
                slices_t1(n)%dat = -1.0D0

                do j = 1, slices_t1(n)%nx
                    do k = 1, n_ftype_selection
                        slices_t1(n)%dat(j, k) = t2%dat(j, xy_selection(n), ftype_selection(k))
                    end do
                end do

            end do
        case default
            stop 'problem in fill_slices_t1_from_slices_of_t2 5'
        end select

    end subroutine fill_slices_t1_from_slices_of_t2
    !***********************************************************

    !***********************************************************
    subroutine fill_t1_from_t2(t2, v1orv2, xy_selection, ftype_selection, t1)
        ! if this subroutine is used,
        ! don't forget to deallocate slices_t1(x)%dat at some point
        !***********************************************************
        implicit none
        !--------------------------------------input
        type(table2), intent(in) :: t2
        integer, intent(in) :: v1orv2
        integer, intent(in) :: xy_selection
        integer, intent(in) :: ftype_selection(:)
        !-------------------------------------output
        type(table1) :: t1
        !-----------------------------------internal
        integer :: i, j, k, n
        integer :: n_ftype_selection
        !-------------------------------------------
        n_ftype_selection = ubound(ftype_selection, dim=1)

        do n = lbound(ftype_selection, dim=1), n_ftype_selection
            if (ftype_selection(n) .lt. lbound(t2%dat, dim=3)) stop 'problem in fill_t1_from_t2 3a'
            if (ftype_selection(n) .gt. ubound(t2%dat, dim=3)) stop 'problem in fill_t1_from_t2 3b'
        end do

        t1%id = t2%id
        t1%deltax = t2%deltax

        select case (v1orv2)
        case (1)

            if (xy_selection .lt. lbound(t2%dat, dim=1)) stop 'problem in fill_t1_from_t2 4a'
            if (xy_selection .gt. ubound(t2%dat, dim=1)) stop 'problem in fill_t1_from_t2 4b'

            t1%nx = t2%nx1
            t1%xmax = t2%xmax1
            t1%xmin = t2%xmin1
            t1%sep = t2%sep1

            allocate (t1%dat(t1%nx, n_ftype_selection))
            t1%dat = -1.0D0

            do i = 1, t1%nx
                do k = 1, n_ftype_selection
                    t1%dat(i, k) = t2%dat(xy_selection, i, ftype_selection(k))
                end do
            end do

        case (2)

            if (xy_selection .lt. lbound(t2%dat, dim=2)) stop 'problem in fill_t1_from_t2 4aa'
            if (xy_selection .gt. ubound(t2%dat, dim=2)) stop 'problem in fill_t1_from_t2 4bb'

            t1%nx = t2%nx2
            t1%xmax = t2%xmax2
            t1%xmin = t2%xmin2
            t1%sep = t2%sep2

            allocate (t1%dat(t1%nx, n_ftype_selection))
            t1%dat = -1.0D0

            do j = 1, t1%nx
                do k = 1, n_ftype_selection
                    t1%dat(j, k) = t2%dat(j, xy_selection, ftype_selection(k))
                end do
            end do

        case default
            stop 'problem in fill_t1_from_t2 5'
        end select

    end subroutine fill_t1_from_t2
    !***********************************************************
    function efficiency_type2(table_id, x1, x2, particle1, particle2)
        use usefulbits, only: not_a_particle
        implicit none

        integer, intent(in) :: table_id, particle1, particle2
        double precision, intent(in) :: x1, x2
        integer :: i
        double precision :: efficiency_type2

        double precision, allocatable :: eff_table(:, :)

        ! silence unused variable warnings
        IF (x2 >= 0) CONTINUE
        IF (particle2 > 0) CONTINUE

        select case (table_id)

        case (109171, 109172) ! Implemented values for g=0.2 for Z' model.

            allocate (eff_table(2, 7))
            select case (table_id)

            case (109171)
                eff_table(1, :) = (/250.0D0, 350.0D0, 449.999999D0, 450.0D0, 550.0D0, 750.0D0, 950.0D0/)
                eff_table(2, :) = (/0.785D0, 0.806D0, 0.799D0, 0.762D0, 0.758D0, 0.765D0, 0.734D0/)

            case (109172) ! Implemented values for g=0.2 for Z' model.
                eff_table(1, :) = (/250.0D0, 350.0D0, 449.999999D0, 450.0D0, 550.0D0, 750.0D0, 950.0D0/)
                eff_table(2, :) = (/0.445D0, 0.440D0, 0.449D0, 0.417D0, 0.446D0, 0.420D0, 0.397D0/)
            end select

            if (particle1 .ne. not_a_particle) then
                do i = lbound(eff_table, dim=2), ubound(eff_table, dim=2)
                    if (x1 .le. eff_table(1, i)) then
                        efficiency_type2 = eff_table(2, i) ! Constant value  outside of range.
                        exit
                    else
                        if (i .ne. ubound(eff_table, dim=2)) then
                            if (x1 .le. eff_table(1, i + 1)) then
                                ! Linear interpolation
                                efficiency_type2 = (eff_table(2, i + 1) - eff_table(2, i))/(eff_table(1, i + 1) - eff_table(1, i)) &
                                                   *(x1 - eff_table(1, i)) + eff_table(2, i)
                                exit
                            end if
                        else
                            efficiency_type2 = eff_table(2, i)
                        end if
                    end if
                end do
            else
                stop 'Error: Efficiency requested for unknown particle type!'
            end if

            deallocate (eff_table)
        end select

    end function efficiency_type2

    function acceptance_type2(table_id, x1, x2, particle1, particle2, selection, flag, specificCharge)
        use usefulbits, only: Hneut, Hplus
        implicit none

        integer, intent(in) :: table_id, particle1, particle2, selection, flag
        integer, intent(in), optional::specificCharge
        double precision, intent(in) :: x1, x2
        double precision :: acceptance_type2
        integer :: charge

        double precision, allocatable :: fit_coeff(:, :)

        ! silence unused variable warnings
        IF (x2 >= 0) CONTINUE
        IF (particle2 > 0) CONTINUE

        if (present(specificCharge)) then
            charge = specificCharge
        else
            charge = 0
        end if

        acceptance_type2 = 1.0D0 ! initial value. If error occurs, it will be set to zero.

        if (flag .le. 0) then
            stop 'Error: integer flag has to be greater than 0.'
        end if

        select case (table_id)
        case (109171) ! ATLAS qq->phi + photon, flavor-inclusive
            ! x1 is assumed to be the mass of particle1,
            ! x2 and particle2 are not used here.
            ! selection destinguishes between the acceptance before (0) and after (1) the mjj cut.
            ! flag selects the relevant quark flavor combination.
            if ((flag .gt. 9) .and. (particle1 .eq. Hneut)) then
                stop "Error: Unknown flag for the qq-> neutral scalar + photon process!"
            end if
            if ((flag .gt. 6) .and. (particle1 .eq. Hplus)) then
                stop "Error: Unknown flag for the qq-> charged scalar + photon process!"
            end if

            if (x1 .lt. 450.0D0 .and. x1 .ge. 225.0D0) then
                ! single photon trigger
                if (selection .eq. 0) then ! preselecton
                    if (particle1 .eq. Hneut) then
                        allocate (fit_coeff(9, 3))
                        fit_coeff(1, :) = (/-0.07453818172065914D0, 3.824816674227349D-05, 0.016692363267610276D0/) !uu
                        fit_coeff(2, :) = (/-0.08365502599272057D0, 2.472432083149923D-05, 0.018160770096650507D0/) !dd
                        fit_coeff(3, :) = (/-0.03036038738244472D0, 4.221832611704063D-05, 0.006949082950767061D0/) !cc
                        fit_coeff(4, :) = (/-0.04243315863987981D0, 3.866314211122518D-05, 0.009348295231225559D0/) !ss
                        fit_coeff(5, :) = (/-0.02914949139483562D0, 3.419393921510185D-05, 0.006807103486670389D0/) !bb
                        fit_coeff(6, :) = (/-0.0689894294152424D0, 3.8902627867545764D-05, 0.015043024297551865D0/) !uc
                        fit_coeff(7, :) = (/-0.0453815009102044D0, 3.744108513321609D-05, 0.010614770882040572D0/) !ds
                        fit_coeff(8, :) = (/-0.048765522647684335D0, 3.456462028873385D-05, 0.01114951745281092D0/) !db
                        fit_coeff(9, :) = (/-0.04967090498121169D0, 3.350800124436956D-05, 0.010467748237181393D0/) !sb
                    else if (particle1 .eq. Hplus) then
                        allocate (fit_coeff(6, 3))
                        if (charge .eq. -1) then
                            fit_coeff(1, :) = (/-0.04484663445244751D0, 2.917809069238506D-05, 0.009845381304258653D0/) !ud
                            fit_coeff(2, :) = (/-0.02438929995556706D0, 3.289201683319401D-05, 0.005487672382653153D0/) !cs
                            fit_coeff(3, :) = (/-0.008308403859409175D0, 4.0866452739329755D-05, 0.0026511859056011886D0/) !us
                            fit_coeff(4, :) = (/-0.0620312387951807D0, 2.267517735452839D-05, 0.012612883806140607D0/) !cd
                            fit_coeff(5, :) = (/-0.023580025808217338D0, 3.0657222293821745D-05, 0.005352044311661076D0/) !ub
                            fit_coeff(6, :) = (/-0.012801225226009098D0, 3.236168398393879D-05, 0.0032221044092733972D0/) !cb
                        elseif (charge .eq. +1) then
                            fit_coeff(1, :) = (/-0.07335163824881155D0, 2.520750232361978D-05, 0.01592191407492654D0/) !ud
                            fit_coeff(2, :) = (/-0.029160671828213205D0, 2.8358926187903368D-05, 0.006174353616870782D0/) !cs
                            fit_coeff(3, :) = (/-0.06261671128216949D0, 2.7693711065568606D-05, 0.01351445163517528D0/) !us
                            fit_coeff(4, :) = (/-0.014512733956164226D0, 3.2090057309420846D-05, 0.0037939580782129917D0/) !cd
                            fit_coeff(5, :) = (/-0.07217431189967927D0, 2.5143326400233142D-05, 0.015194175446230258D0/) !ub
                            fit_coeff(6, :) = (/-0.032923747335909186D0, 2.375420069088394D-05, 0.006937870226189841D0/) !cb
                        else
                            stop "Charge of charged particle in acceptance_type2 has to be +1 or -1."
                        end if
                    else
                        stop "Error: Unknown particle in function acceptance_type2!"
                    end if
                elseif (selection .eq. 1) then ! after mjj cut
                    if (particle1 .eq. Hneut) then
                        allocate (fit_coeff(9, 3))
                        fit_coeff(1, :) = (/0.38525532656239947D0, -2.3628125387441066D-05, 0.036583038825438785D0/) !uu
                        fit_coeff(2, :) = (/0.39943561504923303D0, 1.784454676602849D-05, 0.03677608469507186D0/) !dd
                        fit_coeff(3, :) = (/-0.08698019570012273D0, -0.00017882404715921708D0, 0.11167279264415338D0/) !cc
                        fit_coeff(4, :) = (/0.5248353167619432D0, 4.6061903293067066D-05, 0.01282638879632233D0/) !ss
                        fit_coeff(5, :) = (/0.14312123365961413D0, -3.5483692688954893D-05, 0.060523486229233724D0/) !bb
                        fit_coeff(6, :) = (/0.4373782503387766D0, -2.1099726560508866D-05, 0.023776375433901743D0/) !uc
                        fit_coeff(7, :) = (/0.09950722594177391D0, -5.9714373125369575D-05, 0.08759342621226496D0/) !ds
                        fit_coeff(8, :) = (/0.2548839475478191D0, 1.2711069457280673D-05, 0.04845909024311625D0/) !db
                        fit_coeff(9, :) = (/0.6817853117541732D0, 9.064290919677945D-05, -0.023320207176590838D0/) !sb
                    elseif (particle1 .eq. Hplus) then
                        allocate (fit_coeff(6, 3))
                        if (charge .eq. -1) then
                            fit_coeff(1, :) = (/0.5044416517236326D0, 3.185865787709505D-05, 0.017651476316145345D0/) !ud
                            fit_coeff(2, :) = (/-0.058784400414788184D0, -0.00018154856031669657D0, 0.11332390252462604D0/) !cs
                            fit_coeff(3, :) = (/0.16119824265350222D0, -7.56869322478831D-05, 0.07597218513158863D0/) ! us
                            fit_coeff(4, :) = (/0.710682341234595D0, 0.00011764645546976406D0, -0.025532282301762386D0/) !cd
                            fit_coeff(5, :) = (/0.5598552848055856D0, 6.35468637968493D-05, -0.002042606314717672D0/) !ub
                            fit_coeff(6, :) = (/0.1535119999958045D0, -7.130706112897657D-05, 0.06527192110263923D0/) !cb
                        elseif (charge .eq. +1) then
                            fit_coeff(1, :) = (/0.28920582367011693D0, -5.392197698683376D-05, 0.059170447528441354D0/) !ud
                            fit_coeff(2, :) = (/0.1360543536099816D0, -8.8595728151652D-05, 0.07775522896899756D0/) !cs
                            fit_coeff(3, :) = (/0.36403229644970786D0, -2.0313192545953587D-05, 0.04273987985113689D0/) ! us
                            fit_coeff(4, :) = (/0.03753114628032014D0, -0.00010725799982331615D0, 0.09494959625977027D0/) !cd
                            fit_coeff(5, :) = (/-0.058835087810597805D0, -0.00014543013924738662D0, 0.10940912532011073D0/) !ub
                            fit_coeff(6, :) = (/-0.12533421966904126D0, -0.00015829179090177497D0, 0.11481519723740238D0/) !cb
                        else
                            stop "Charge of charged particle in acceptance_type2 has to be +1 or -1."
                        end if
                    else
                        stop "Error: Unknown particle in function acceptance_type2!"
                    end if
                else
                    stop "Error: Unknown selection in function acceptance_type2!"
                end if
            elseif (x1 .ge. 450.0D0 .and. x1 .le. 1100.0D0) then
                ! combined trigger
                if (selection .eq. 0) then ! preselecton
                    if (particle1 .eq. Hneut) then
                        allocate (fit_coeff(9, 3))
                        !fit function: c[0] + c[1]*(M - 225) + c[2]*np.log(M)
                        fit_coeff(1, :) = (/-0.3495397605469713D0, -1.2642459929452903D-05, 0.06945421537488733D0/) !uu
                        fit_coeff(2, :) = (/-0.3471775210703532D0, -2.5919492355969983D-05, 0.06893383402208757D0/) !dd
                        fit_coeff(3, :) = (/-0.25771965185274837D0, 5.2868594897241616D-06, 0.05060836725764018D0/) !cc
                        fit_coeff(4, :) = (/-0.2161349881915616D0, 2.0265589959806386D-05, 0.04342174072764107D0/) !ss
                        fit_coeff(5, :) = (/-0.2613043396316116D0, -3.894324500526707D-06, 0.051132326713863516D0/) !bb
                        fit_coeff(6, :) = (/-0.3485345843201749D0, -1.223497402338667D-05, 0.06839930939545714D0/) !uc
                        fit_coeff(7, :) = (/-0.25385483138026416D0, 6.894729090512683D-06, 0.05128069341398131D0/) !ds
                        fit_coeff(8, :) = (/-0.2726412110713955D0, 7.846779718611548D-07, 0.054373087775165D0/) !db
                        fit_coeff(9, :) = (/-0.2784463150279004D0, -2.3232981890858514D-06, 0.054315795840985275D0/) !sb
                    else if (particle1 .eq. Hplus) then
                        !fit function: c[0] + c[1]*(M - 225) + c[2]*np.log(M)
                        allocate (fit_coeff(6, 3))
                        if (charge .eq. -1) then
                            fit_coeff(1, :) = (/-0.25710823612127615D0, -5.034398780010226D-06, 0.05058736718447029D0/) !ud
                            fit_coeff(2, :) = (/-0.18967205251280128D0, 1.3486331956279841D-05, 0.037541628609521455D0/) !cs
                            fit_coeff(3, :) = (/-0.16896121633864306D0, 2.2335584697813463D-05, 0.033933534398936774D0/) !us
                            fit_coeff(4, :) = (/-0.25812436133618843D0, -4.261253080856257D-06, 0.050327388108494284D0/) !cd
                            fit_coeff(5, :) = (/-0.19271439731198028D0, 8.241909640369068D-06, 0.03809804952102698D0/) !ub
                            fit_coeff(6, :) = (/-0.1848252094694455D0, 6.595742525862169D-06, 0.036528211472614765D0/) !cb
                        elseif (charge .eq. +1) then
                            fit_coeff(1, :) = (/-0.32901275609576014D0, -2.329295062067327D-05, 0.06498087881664942D0/) !ud
                            fit_coeff(2, :) = (/-0.2138055657620705D0, -5.777785130482657D-07, 0.04162969122251766D0/) !cs
                            fit_coeff(3, :) = (/-0.29351696001599786D0, -1.115565824728008D-05, 0.057883413664741076D0/) !us
                            fit_coeff(4, :) = (/-0.19235684407120268D0, 6.230537342199852D-06, 0.038157270574518884D0/) !cd
                            fit_coeff(5, :) = (/-0.3087878807658839D0, -1.4223220997063607D-05, 0.060621421727487405D0/) !ub
                            fit_coeff(6, :) = (/-0.17933973462875433D0, 1.0897748362352643D-05, 0.035417916916325634D0/) !cb
                        else
                            stop "Charge of charged particle in acceptance_type2 has to be +1 or -1."
                        end if
                    else
                        stop "Error: Unknown particle in function acceptance_type2!"
                    end if
                elseif (selection .eq. 1) then ! after mjj cut
                    if (particle1 .eq. Hneut) then
                        allocate (fit_coeff(9, 3))
                        fit_coeff(1, :) = (/1.4179009959807805D0, 0.00019318402078541907D0, -0.13213062021217128D0/) !uu
                        fit_coeff(2, :) = (/1.5431432938623932D0, 0.0002603304761813639D0, -0.15281513143882092D0/) !dd
                        fit_coeff(3, :) = (/1.2978470088122862D0, 0.00015943440694614355D0, -0.11976935313211168D0/) !cc
                        fit_coeff(4, :) = (/1.4554725336733858D0, 0.00022367903663272575D0, -0.1380001570888785D0/) !ss
                        fit_coeff(5, :) = (/1.4594279432344206D0, 0.0002497796055038875D0, -0.15908339475374172D0/) !bb
                        fit_coeff(6, :) = (/1.352570526510758D0, 0.00017897228740841698D0, -0.12577105956572918D0/) !uc
                        fit_coeff(7, :) = (/1.2528303938956713D0, 0.0001800583268401668D0, -0.10282440526626438D0/) !ds
                        fit_coeff(8, :) = (/1.3204031551147875D0, 0.00022174772906246543D0, -0.12586100686977048D0/) !db
                        fit_coeff(9, :) = (/1.6547033810359737D0, 0.00029549146128929273D0, -0.1829824410338565D0/) !sb
                    elseif (particle1 .eq. Hplus) then
                        allocate (fit_coeff(6, 3))
                        if (charge .eq. -1) then
                            fit_coeff(1, :) = (/1.3818977638832035D0, 0.00020755445448888415D0, -0.12381726801188014D0/) !ud
                            fit_coeff(2, :) = (/0.8143455081368138D0, 2.4380225860519337D-05, -0.02858203998980347D0/) !cs
                            fit_coeff(3, :) = (/1.24000868172595D0, 0.00015736056718772585D0, -0.09943870123404218D0/) !us
                            fit_coeff(4, :) = (/1.039022396261481D0, 0.0001343541213216677D0, -0.06984861338158219D0/) !cd
                            fit_coeff(5, :) = (/1.6368519487358948D0, 0.00029132603866931765D0, -0.1791466833427859D0/) !ub
                            fit_coeff(6, :) = (/0.9315940135101826D0, 9.266106762835439D-05, -0.059585540188055924D0/) !cb
                        elseif (charge .eq. +1) then
                            fit_coeff(1, :) = (/1.288081714061273D0, 0.00016647838199530424D0, -0.10550214628762097D0/) !ud
                            fit_coeff(2, :) = (/1.1585741495362076D0, 0.0001499435859582795D0, -0.09041290604805757D0/) !cs
                            fit_coeff(3, :) = (/1.2135143760634988D0, 0.00016756223554864104D0, -0.09590287762626909D0/) !us
                            fit_coeff(4, :) = (/1.2909510263482682D0, 0.00018627160475987738D0, -0.11257857644665563D0/) !cd
                            fit_coeff(5, :) = (/1.2017501461720643D0, 0.00016579167960714954D0, -0.10232834425581878D0/) !ub
                            fit_coeff(6, :) = (/1.140050828111984D0, 0.00013845098165656946D0, -0.09551852523990538D0/) !cb
                        else
                            stop "Charge of charged particle in acceptance_type2 has to be +1 or -1."
                        end if
                    else
                        stop "Error: Unknown particle in function acceptance_type2!"
                    end if
                else
                    stop "Error: Unknown selection in function acceptance_type2!"
                end if
            else
                ! write(*,*) "WARNING: Acceptance function used outside range of validity!"
                acceptance_type2 = 0.0D0
            end if

            if (acceptance_type2 .gt. 0.0D0) then ! True if mass value is within range of validity
                acceptance_type2 = fit_coeff(flag, 1) + fit_coeff(flag, 2)*(x1 - 225.0D0) + fit_coeff(flag, 3)*dlog(x1)
            end if
        case (109172) ! ATLAS qq->phi + photon, 2b category
            if (particle1 .ne. Hneut) then
                stop "Error: Acceptance of 2b-tagged category only available for neutral scalars!"
            end if
            if (x1 .lt. 450.0D0 .and. x1 .ge. 225.0D0) then
                ! single photon trigger
                acceptance_type2 = 0.7349093272718935D0 - 3.918931699816906D-06*(x1 - 225.0D0) + 0.028483738485628442D0*dlog(x1)
            elseif (x1 .ge. 450.0D0 .and. x1 .le. 1100.0D0) then
                ! combined trigger
                acceptance_type2 = 0.301508632706373D0 - 0.00011289893741976084D0*(x1 - 225.0D0) + 0.10288525863298086D0*dlog(x1)
            else
                ! write(*,*) "WARNING: Acceptance function used outside range of validity!"
                acceptance_type2 = 0.0D0
            end if

        case default
            acceptance_type2 = 1.0
        end select

        ! write(*,*) "Acceptance fct.:", table_id, x1, particle2, selection, flag, acceptance_type2

        if (allocated(fit_coeff)) deallocate (fit_coeff)
    end function acceptance_type2

    function shifted_mass(table_id, x, particle, flag)
        use usefulbits, only: Hneut, Hplus
        implicit none

        integer, intent(in) :: table_id, particle, flag
        double precision, intent(in) :: x
        double precision :: shifted_mass

        double precision, allocatable :: fit_coeff(:, :)

        shifted_mass = x ! initial value. If error occurs, it will be set to zero.

        if (flag .le. 0) then
            stop 'Error: integer flag has to be greater than 0.'
        end if

        select case (table_id)
        case (109171) ! ATLAS qq->phi + photon, flavor-inclusive
            ! x is assumed to be the true mass of particle
            ! flag selects the relevant quark flavor combination.
            if ((flag .gt. 9) .and. (particle .eq. Hneut)) then
                stop "Error: Unknown flag for the qq-> neutral scalar + photon process!"
            end if
            if ((flag .gt. 6) .and. (particle .eq. Hplus)) then
                stop "Error: Unknown flag for the qq-> charged scalar + photon process!"
            end if

            if (x .lt. 450.0D0 .and. x .ge. 225.0D0) then
                ! single photon trigger
                if (particle .eq. Hneut) then
                    allocate (fit_coeff(9, 3))
                    fit_coeff(1, :) = (/281.4319562313968D0, 0.9673914266961734D0, -12.235462255102286D0/)
                    fit_coeff(2, :) = (/258.81946927233946D0, 0.9633838763717543D0, -8.394837719212434D0/)
                    fit_coeff(3, :) = (/224.11771622447955D0, 0.9529383908273734D0, -2.1334395373697994D0/)
                    fit_coeff(4, :) = (/284.55872258550187D0, 0.9525467718329332D0, -12.822919404917098D0/)
                    fit_coeff(5, :) = (/220.9499460459883D0, 0.9266918231874117D0, -2.069497910863455D0/)
                    fit_coeff(6, :) = (/289.3880456860291D0, 0.9629350202691872D0, -13.78476867531609D0/)
                    fit_coeff(7, :) = (/193.21640722258093D0, 0.9349332966650078D0, 3.7982633636725636D0/)
                    fit_coeff(8, :) = (/216.08451117396797D0, 0.9327140780790565D0, -0.8275101427359572D0/)
                    fit_coeff(9, :) = (/238.35298494195436D0, 0.9419221069111677D0, -4.847633989956753D0/)
                else if (particle .eq. Hplus) then
                    !fit function: c[0] + c[1]*(M - 225) + c[2]*np.log(M)
                    allocate (fit_coeff(6, 3))
                    fit_coeff(1, :) = (/272.117378980746D0, 0.9685260914441841D0, -10.599655101495404D0/) ! ud
                    fit_coeff(2, :) = (/258.9945452730261D0, 0.9564116331738203D0, -8.177737850503428D0/) ! cs
                    fit_coeff(3, :) = (/154.2868145807532D0, 0.923627818671771D0, 10.522652458614687D0/) ! us
                    fit_coeff(4, :) = (/243.18626334632324D0, 0.9531014690805794D0, -5.539394174464264D0/) ! cd
                    fit_coeff(5, :) = (/273.31450166305825D0, 0.954571643579162D0, -11.148508503812018D0/) ! ub
                    fit_coeff(6, :) = (/209.0255570171927D0, 0.9211762223473134D0, 0.6774235947175775D0/) ! cb
                else
                    stop "Error: Unknown particle in function shifted_mass!"
                end if
            elseif (x .ge. 450.0D0 .and. x .le. 1100.0D0) then
                ! combined trigger
                if (particle .eq. Hneut) then
                    allocate (fit_coeff(9, 3))
                    !fit function: c[0] + c[1]*(M - 225) + c[2]*np.log(M)
                    fit_coeff(1, :) = (/268.6419095201171D0, 0.9600671239644399D0, -9.859514713821117D0/)
                    fit_coeff(2, :) = (/272.21432619456056D0, 0.9643979603066067D0, -10.626478039631243D0/)
                    fit_coeff(3, :) = (/254.35644581068794D0, 0.9611958334854936D0, -7.505911098316428D0/)
                    fit_coeff(4, :) = (/254.71289010244473D0, 0.9425188137399285D0, -7.543214214890992D0/)
                    fit_coeff(5, :) = (/267.9237702016445D0, 0.9388069130521747D0, -10.263128783184035D0/)
                    fit_coeff(6, :) = (/280.28817594827103D0, 0.9587600118446112D0, -12.115839738882553D0/)
                    fit_coeff(7, :) = (/223.71851983237093D0, 0.9454616308413796D0, -1.7582329141468864D0/)
                    fit_coeff(8, :) = (/233.5373585585733D0, 0.9367802595332158D0, -3.953300640830744D0/)
                    fit_coeff(9, :) = (/229.9860946965195D0, 0.9387178657637238D0, -3.373496593320827D0/)
                else if (particle .eq. Hplus) then
                    !fit function: c[0] + c[1]*(M - 225) + c[2]*np.log(M)
                    allocate (fit_coeff(6, 3))
                    fit_coeff(1, :) = (/276.0391663421334D0, 0.9661644772641924D0, -11.14797232157632D0/)
                    fit_coeff(2, :) = (/242.3824874015853D0, 0.9481249846138476D0, -5.194289746546004D0/)
                    fit_coeff(3, :) = (/221.30860410941688D0, 0.9471664169027666D0, -1.4268813388883332D0/)
                    fit_coeff(4, :) = (/225.17451363750908D0, 0.9431376149033346D0, -2.241990843366313D0/)
                    fit_coeff(5, :) = (/264.80736830157275D0, 0.9503271038762466D0, -9.478877550413694D0/)
                    fit_coeff(6, :) = (/253.24791198178133D0, 0.9382083805622718D0, -7.390505995052736D0/)
                else
                    stop "Error: Unknown particle in function shifted_mass!"
                end if
            else
                ! write(*,*) "WARNING: shifted_mass function used outside range of validity!"
                shifted_mass = 0.0D0
            end if

            if (shifted_mass .gt. 0.0D0) then ! True if mass value is within range of validity
                shifted_mass = fit_coeff(flag, 1) + fit_coeff(flag, 2)*(x - 225.0D0) + fit_coeff(flag, 3)*dlog(x)
            end if

            if (allocated(fit_coeff)) deallocate (fit_coeff)

        case (109172) ! ATLAS qq->phi + photon, 2b category
            if (particle .ne. Hneut) then
                stop "Error: Acceptance of 2b-tagged category only available for neutral scalars!"
            end if
            if (x .lt. 450.0D0 .and. x .ge. 225.0D0) then
                ! single photon trigger
                shifted_mass = 243.16680381218623D0 + 0.9358995435370389D0*(x - 225.0D0) - 6.0740001337184495D0*dlog(x)
            elseif (x .ge. 450.0D0 .and. x .le. 1100.0D0) then
                ! combined trigger
                shifted_mass = 278.3215572624196D0 + 0.9453027258944698D0*(x - 225.0D0) - 12.232501554891819D0*dlog(x)
            else
                ! write(*,*) "WARNING: shifted_mass function used outside range of validity!"
                shifted_mass = 0.0D0
            end if

        case default
            shifted_mass = x
        end select

    end function shifted_mass

end module S95tables_type2
!************************************************************
