!> @file
!> module 'likelihoods'.
!! Calculates the exclusion likelihood for LHC searches. Currently, this
!! is done only for BSM Higgs searches with \f$ \tau^+ \tau^-\f$ final states.
!! The module is mostly used internally. The likelihood values can be accessed via
!! #higgsbounds_get_likelihood() (general result).
module likelihoods
!******************************************************************
! Some options:
    logical :: rescale = .False.
    logical :: fudgelimit = .False. ! For CMS-14-029

    !! use processed likelihood tables to prevent overexclusion in case of an
    !! observed excesses (see manual)
    logical, parameter :: preventOverexclusion = .True.
!--
    type likelihood2D
        integer :: analysisID
        logical :: pub
        integer :: particle_x
        character(LEN=45) :: label
        character(LEN=35) :: citekey
        character(LEN=3) :: expt
        double precision :: lumi, energy
        character(LEN=100) :: description
        double precision :: mass_sep
        double precision :: mass_range_tolerance
        integer :: n2Dslices
        double precision :: mass
        double precision :: xstep
        double precision :: ystep
        double precision :: xmin
        double precision :: xmax
        double precision :: ymin
        double precision :: ymax
        integer :: size
        integer :: xsize
        integer :: ysize
        double precision, dimension(200, 200) :: llhgrid_pred, llhgrid_obs
        double precision, dimension(1:3) :: bestfit_pred, bestfit_obs
    end type

    type likelihood3D
        integer :: n2Dslices
        type(likelihood2D), dimension(38) :: D2llhdata ! <- Set to the highest number of n2Dslices available!
    end type

    integer, parameter :: nllhs = 4
    integer, dimension(1:nllhs) :: n2Dslices = (/38, 29, 15, 14/)

    type(likelihood3D), dimension(1:nllhs) :: llhdata

    type threetuple
        double precision, dimension(1:3) :: val
    end type

! If more likelihoods become available (nllhs > 1), have to generalize this
! to a list of likelihoods

! CMS related stuff
!  integer, parameter :: n2Dslices = 38

contains
!--------------------------------------------------------------
    subroutine setup_likelihoods
!--------------------------------------------------------------
        use usefulbits, only: file_id_common4, Hneut
        use install_data, only: pathname
!  implicit none

        integer :: ios, i, j, status

! for CMS-3316 likelihoods
!   integer :: analysisID = 3316
!   character(LEN=100) :: label = '[hep-ex] arXiv:1408.3316 (CMS)'
!   character(LEN=100) :: stem = 'Expt_tables/CMStables/CMS_tautau_llh_1408.3316/'
!  character(LEN=100) :: stem = 'Expt_tables/CMStables/CMS_tautau_llh_1408.3316_withSMHinBG/'
        integer :: analysisID
        logical :: pub
        character(LEN=45) :: label
        character(LEN=35) :: citekey
        character(LEN=100) :: stem
        character(LEN=100) :: prefix
        character(LEN=100) :: description
        character(LEN=3) :: expt
        integer :: particle_type
        double precision :: lumi
        double precision :: energy
        double precision :: mass_sep! in percentage
        double precision :: mass_range_tolerance ! If the predicted mass lies outside the mass range within the tolerance,
        ! the mass will be rounded to the closest value within the grid.
        double precision, allocatable :: masses(:), xmin(:), xmax(:), ymin(:), ymax(:), xstep(:), ystep(:)
        integer, allocatable :: xsize(:), ysize(:)

!- Setup the global likelihood data object

        do i = 1, nllhs
            llhdata(i)%n2Dslices = n2Dslices(i)
!   allocate(llhdata(i)%D2llhdata(n2Dslices(i)))
        enddo

!- Try reading the data in from binary file
        open (file_id_common4, file=trim(adjustl(pathname))// &
              '/Expt_tables/tautau_llh.binary', form='unformatted')
        read (file_id_common4, iostat=ios) llhdata

        if (ios .ne. 0) then

!------Setup first likelihood -->  CMS-PAS-HIG 14-029

            analysisID = 14029
            pub = .False.
            label = 'CMS-HIG-PAS 14-029'
            citekey = 'CMS:2015mca'
            stem = 'Expt_tables/CMStables/CMS_tautau_llh_14029/'
            description = '(pp) -> phi -> tautau, using -2ln(L) reconstruction'
            expt = 'CMS'
            particle_type = Hneut
            lumi = 19.7D0
            energy = 8.0D0
            mass_sep = 0.20D0
            mass_range_tolerance = 10.0D0

            allocate (masses(n2Dslices(1)), xmin(n2Dslices(1)), xmax(n2Dslices(1)))
            allocate (ymin(n2Dslices(1)), ymax(n2Dslices(1)), xstep(n2Dslices(1)))
            allocate (ystep(n2Dslices(1)), xsize(n2Dslices(1)), ysize(n2Dslices(1)))

            masses = (/90.0D0, 100.0D0, 110.0D0, 120.0D0, 125.0D0, 130.0D0, 140.0D0, &
                       150.0D0, 160.0D0, 170.0D0, 180.0D0, 190.0D0, 200.0D0, 210.0D0, &
                       220.0D0, 230.0D0, 240.0D0, 250.0D0, 275.0D0, 300.0D0, 325.0D0, &
                       350.0D0, 375.0D0, 400.0D0, 425.0D0, 450.0D0, 475.0D0, 500.0D0, &
                       550.0D0, 600.0D0, 650.0D0, 700.0D0, 750.0D0, 800.0D0, 850.0D0, &
                       900.0D0, 950.0D0, 1000.0D0/)

!-------- Settings for the CMS-14-029 analysis ---------!
            xmin = (/0.0875000, 0.0750000, 0.0500000, 0.0150000, 0.0125000, 0.0087500, 0.0050000, &
                     0.0050000, 0.0025000, 0.0025000, 0.0015000, 0.0015000, 0.0010000, 0.0010000, &
                     0.0010000, 0.0010000, 0.0010000, 0.0003750, 0.0003750, 0.0003000, 0.0003000, &
                     0.0003000, 0.0003000, 0.0003000, 0.0003000, 0.0001750, 0.0001750, 0.0001750, &
                     0.0001750, 0.0001000, 0.0001000, 0.0000750, 0.0000750, 0.0000625, 0.0000625, &
                     0.0000375, 0.0000375, 0.0000300/)
            xmax = (/34.9125000, 29.9250000, 19.9500000, 5.9850000, 4.9875000, 3.4912500, 1.9950000, &
                     1.9950000, 0.9975000, 0.9975000, 0.5985000, 0.5985000, 0.3990000, 0.3990000, &
                     0.3990000, 0.3990000, 0.3990000, 0.1496250, 0.1496250, 0.1197000, 0.1197000, &
                     0.1197000, 0.1197000, 0.1197000, 0.1197000, 0.0698250, 0.0698250, 0.0698250, &
                     0.0698250, 0.0399000, 0.0399000, 0.0299250, 0.0299250, 0.0249375, 0.0249375, &
                     0.0149625, 0.0149625, 0.0119700/)
            ymin = (/0.0375000, 0.0250000, 0.0200000, 0.0100000, 0.0075000, 0.0050000, 0.0025000, &
                     0.0025000, 0.0015000, 0.0015000, 0.0010000, 0.0010000, 0.0007500, 0.0007500, &
                     0.0007500, 0.0007500, 0.0007500, 0.0005000, 0.0005000, 0.0003000, 0.0003000, &
                     0.0003000, 0.0003000, 0.0003000, 0.0003000, 0.0001750, 0.0001750, 0.0001000, &
                     0.0001000, 0.0000750, 0.0000750, 0.0000625, 0.0000625, 0.0000500, 0.0000500, &
                     0.0000375, 0.0000375, 0.0000300/)
            ymax = (/14.9625000, 9.9750000, 7.9800000, 3.9900000, 2.9925000, 1.9950000, 0.9975000, &
                     0.9975000, 0.5985000, 0.5985000, 0.3990000, 0.3990000, 0.2992500, 0.2992500, &
                     0.2992500, 0.2992500, 0.2992500, 0.1995000, 0.1995000, 0.1197000, 0.1197000, &
                     0.1197000, 0.1197000, 0.1197000, 0.1197000, 0.0698250, 0.0698250, 0.0399000, &
                     0.0399000, 0.0299250, 0.0299250, 0.0249375, 0.0249375, 0.0199500, 0.0199500, &
                     0.0149625, 0.0149625, 0.0119700/)
            xstep = (/0.1750000, 0.1500000, 0.1000000, 0.0300000, 0.0250000, 0.0175000, 0.0100000, &
                      0.0100000, 0.0050000, 0.0050000, 0.0030000, 0.0030000, 0.0020000, 0.0020000, &
                      0.0020000, 0.0020000, 0.0020000, 0.0007500, 0.0007500, 0.0006000, 0.0006000, &
                      0.0006000, 0.0006000, 0.0006000, 0.0006000, 0.0003500, 0.0003500, 0.0003500, &
                      0.0003500, 0.0002000, 0.0002000, 0.0001500, 0.0001500, 0.0001250, 0.0001250, &
                      0.0000750, 0.0000750, 0.0000600/)
            ystep = (/0.0750000, 0.0500000, 0.0400000, 0.0200000, 0.0150000, 0.0100000, 0.0050000, &
                      0.0050000, 0.0030000, 0.0030000, 0.0020000, 0.0020000, 0.0015000, 0.0015000, &
                      0.0015000, 0.0015000, 0.0015000, 0.0010000, 0.0010000, 0.0006000, 0.0006000, &
                      0.0006000, 0.0006000, 0.0006000, 0.0006000, 0.0003500, 0.0003500, 0.0002000, &
                      0.0002000, 0.0001500, 0.0001500, 0.0001250, 0.0001250, 0.0001000, 0.0001000, &
                      0.0000750, 0.0000750, 0.0000600/)
!------------------
            xsize = (/200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200/)
            ysize = (/200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200/)
!------------------
            do j = 1, n2Dslices(1)
                call write_metadata(analysisID, pub, label, citekey, description, expt, lumi, energy, particle_type, &
                                    mass_sep, mass_range_tolerance, masses(j), xmin(j), xmax(j), xstep(j), xsize(j), ymin(j), &
                                    ymax(j), ystep(j), ysize(j), llhdata(1)%D2llhdata(j))

                prefix = 'L_data_b_'
                call read2Ddata(masses(j), llhdata(1)%D2llhdata(j), stem, prefix, 'obs', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
                prefix = 'L_asimov_b_'
                call read2Ddata(masses(j), llhdata(1)%D2llhdata(j), stem, prefix, 'pred', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
            enddo

            deallocate (masses, xmin, xmax)
            deallocate (ymin, ymax, xstep)
            deallocate (ystep, xsize, ysize)

!------ likelihood CMS-PAS-HIG 16-037 is superseded by CMS-PAS-HIG 17-020 (see next)
!------Setup second likelihood -->  CMS-PAS-HIG 17-020

            analysisID = 17020
            pub = .False.
            label = '[hep-ex] 1803.06553 (CMS)' ! 'CMS-HIG-PAS 17-020'
            citekey = 'Sirunyan:2018zut'
            stem = 'Expt_tables/CMStables/CMS_tautau_llh_17020/'
            description = '(pp) -> phi -> tautau, using -2ln(L) reconstruction'
            expt = 'CMS'
            particle_type = Hneut
            lumi = 35.9D0
            energy = 13.0D0
            mass_sep = 0.20D0
            mass_range_tolerance = 10.0D0

            allocate (masses(n2Dslices(2)), xmin(n2Dslices(2)), xmax(n2Dslices(2)))
            allocate (ymin(n2Dslices(2)), ymax(n2Dslices(2)), xstep(n2Dslices(2)))
            allocate (ystep(n2Dslices(2)), xsize(n2Dslices(2)), ysize(n2Dslices(2)))

!-------- Settings for the CMS-17-020 analysis ---------!
            masses = (/90.0D0, 100.0D0, 110.0D0, 120.0D0, 125.0D0, 130.0D0, 140.0D0, &
                       160.0D0, 180.0D0, 200.0D0, 250.0D0, 350.0D0, 400.0D0, 450.0D0, &
                       500.0D0, 600.0D0, 700.0D0, 800.0D0, 900.0D0, 1000.0D0, 1200.0D0, &
                       1400.0D0, 1600.0D0, 1800.0D0, 2000.0D0, 2300.0D0, 2600.0D0, 2900.0D0, 3200.0D0/)

! --- no SM Higgs in BG ---!

            xmin = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
            xmax = (/26.865D0, 44.775D0, 24.875D0, 17.91D0, 14.925D0, 8.955D0, 5.97D0, 2.985D0, 2.0895D0, &
                     1.194D0, 0.796D0, 0.14925D0, 0.10945D0, 0.0796D0, 0.0796D0, 0.064675D0, 0.054725D0, &
                     0.034825D0, 0.024875D0, 0.018905D0, 0.012935D0, 0.0114425D0, 0.006965D0, 0.0064675D0, &
                     0.00597D0, 0.00597D0, 0.0054725D0, 0.0054725D0, 0.0054725D0/)
            xstep = (/0.135D0, 0.225D0, 0.125D0, 0.09D0, 0.075D0, 0.045D0, 0.03D0, 0.015D0, 0.0105D0, &
                      0.006D0, 0.004D0, 0.00075D0, 0.00055D0, 0.0004D0, 0.0004D0, 0.000325D0, 0.000275D0, &
                      0.000175D0, 0.000125D0, 9.5D-05, 6.5D-05, 5.75D-05, 3.5D-05, 3.25D-05, 3D-05, &
                      3D-05, 2.75D-05, 2.75D-05, 2.75D-05/)
            ymin = (/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
            ymax = (/24.875D0, 14.925D0, 14.925D0, 9.95D0, 7.96D0, 5.97D0, 5.4725D0, 2.4875D0, 1.791D0, &
                     1.194D0, 0.597D0, 0.14925D0, 0.08955D0, 0.0597D0, 0.04975D0, 0.024875D0, 0.01791D0, &
                     0.01194D0, 0.00995D0, 0.008955D0, 0.0064675D0, 0.0054725D0, 0.004975D0, 0.004975D0, &
                     0.0044775D0, 0.00398D0, 0.0042785D0, 0.004179D0, 0.0042785D0/)
            ystep = (/0.125D0, 0.075D0, 0.075D0, 0.05D0, 0.04D0, 0.03D0, 0.0275D0, 0.0125D0, 0.009D0, &
                      0.006D0, 0.003D0, 0.00075D0, 0.00045D0, 0.0003D0, 0.00025D0, 0.000125D0, 9D-05, &
                      6D-05, 5D-05, 4.5D-05, 3.25D-05, 2.75D-05, 2.5D-05, 2.5D-05, 2.25D-05, 2D-05, &
                      2.15D-05, 2.1D-05, 2.15D-05/)

! --- with SM Higgs in BG ---!

! xmin = (/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,&
! &        0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,&
! &        0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/)
! xmax = (/22.885D0,35.82D0,20.895D0,13.93D0,10.945D0,8.4575D0,4.975D0,2.786D0,1.791D0,&
! &         1.194D0,0.796D0,0.1592D0,0.1194D0,0.08955D0,0.084575D0,0.06965D0,0.04975D0,&
! &         0.034825D0,0.02587D0,0.020895D0,0.012935D0,0.00995D0,0.00796D0,0.006965D0,&
! &         0.00597D0,0.00597D0,0.00597D0,0.0054725D0,0.0054725D0/)
! xstep = (/0.115D0,0.18D0,0.105D0,0.07D0,0.055D0,0.0425D0,0.025D0,0.014D0,0.009D0,0.006D0,&
! &         0.004D0,0.0008D0,0.0006D0,0.00045D0,0.000425D0,0.00035D0,0.00025D0,0.000175D0,&
! &         0.00013D0,0.000105D0,6.5D-05,5D-05,4D-05,3.5D-05,3D-05,3D-05,3D-05,2.75D-05,2.75D-05/)
! ymin = (/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,&
! &        0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,&
! &        0.0D0,0.0D0,0.0D0/)
! ymax = (/23.88D0,12.935D0,12.935D0,8.955D0,7.96D0,5.4725D0,4.975D0,2.4875D0,1.791D0,&
! &        1.2935D0,0.64675D0,0.1194D0,0.0995D0,0.0597D0,0.04975D0,0.02587D0,0.01592D0,&
! &        0.01194D0,0.00995D0,0.0084575D0,0.0064675D0,0.00597D0,0.004975D0,0.004975D0,&
! &        0.0044775D0,0.0044775D0,0.0044775D0,0.0044775D0,0.0044775D0/)
! ystep = (/0.12D0,0.065D0,0.065D0,0.045D0,0.04D0,0.0275D0,0.025D0,0.0125D0,0.009D0,&
! &         0.0065D0,0.00325D0,0.0006D0,0.0005D0,0.0003D0,0.00025D0,0.00013D0,8D-05,&
! &         6D-05,5D-05,4.25D-05,3.25D-05,3D-05,2.5D-05,2.5D-05,2.25D-05,2.25D-05,&
! &         2.25D-05,2.25D-05,2.25D-05/)

!------------------
! n.b.: Maybe need to change to 199!
            xsize = (/200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, 200/)
            ysize = (/200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, &
                      200, 200, 200, 200, 200, 200, 200, 200/)
!------------------
            do j = 1, n2Dslices(2)
                call write_metadata(analysisID, pub, label, citekey, description, expt, lumi, energy, particle_type, &
                                    mass_sep, mass_range_tolerance, masses(j), xmin(j), xmax(j), xstep(j), xsize(j), ymin(j), &
                                    ymax(j), ystep(j), ysize(j), llhdata(2)%D2llhdata(j))

!    prefix='2D_scan_noSMHinBG_'
                prefix = '2D_scan_noSMHinBG_'
                call read2Ddata(masses(j), llhdata(2)%D2llhdata(j), stem, prefix, 'obs', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
!    prefix='2D_scan_noSMHinBG_Asimov_'
                prefix = '2D_scan_noSMHinBG_Asimov_'
                call read2Ddata(masses(j), llhdata(2)%D2llhdata(j), stem, prefix, 'pred', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
            enddo

            deallocate (masses, xmin, xmax)
            deallocate (ymin, ymax, xstep)
            deallocate (ystep, xsize, ysize)

!------Setup third likelihood -->  ATLAS 1709.07242

            analysisID = 170907242
            pub = .False.
            label = 'arXiv:1709.07242 [hep-ex] (ATLAS)'
            citekey = 'Aaboud:2017sjh'
            stem = 'Expt_tables/ATLtables/ATL_tautau_llh_1709.07242/'
            description = '(pp) -> phi -> tautau, using -2ln(L) reconstruction'
            expt = 'ATL'
            particle_type = Hneut
            lumi = 36.1D0
            energy = 13.0D0
            mass_sep = 0.20D0
            mass_range_tolerance = 10.0D0

            allocate (masses(n2Dslices(3)), xmin(n2Dslices(3)), xmax(n2Dslices(3)))
            allocate (ymin(n2Dslices(3)), ymax(n2Dslices(3)), xstep(n2Dslices(3)))
            allocate (ystep(n2Dslices(3)), xsize(n2Dslices(3)), ysize(n2Dslices(3)))

            masses = (/200.0D0, 250.0D0, 300.0D0, 350.0D0, 400.0D0, 500.0D0, 600.0D0, &
                       700.0D0, 800.0D0, 1000.0D0, 1200.0D0, 1500.0D0, 1750.0D0, &
                       2000.0D0, 2250.0D0/)

!-------- Settings for the ATLAS 1709.07242 ---------!
            xmin = (/0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
                     0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
                     0.0000000/)
            xmax = (/3.0000000, 1.8000000, 1.0000000, 0.7000000, 0.6000000, 0.2200000, 0.1400000, &
                     0.1000000, 0.0700000, 0.0500000, 0.0300000, 0.0300000, 0.0300000, 0.0300000, &
                     0.0300000/)
            xstep = (/3./99., 1.8/99., 1./99., 0.7/99., 0.6/99., 0.22/99., 0.14/99., &
                      0.1/99., 0.07/99., 0.05/99., 0.03/99., 0.03/99., 0.03/99., 0.03/99., &
                      0.03/99./)
            ymin = (/0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
                     0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
                     0.0000000/)
            ymax = (/3.0000000, 1.8000000, 1.0000000, 0.7000000, 0.6000000, 0.2200000, 0.1400000, &
                     0.1000000, 0.0700000, 0.0500000, 0.0300000, 0.0300000, 0.0300000, 0.0300000, &
                     0.0300000/)
            ystep = (/3./99., 1.8/99., 1./99., 0.7/99., 0.6/99., 0.22/99., 0.14/99., &
                      0.1/99., 0.07/99., 0.05/99., 0.03/99., 0.03/99., 0.03/99., 0.03/99., &
                      0.03/99./)
!------------------
! n.b.: Maybe need to change to 199!
            xsize = (/100, 100, 100, 100, 100, 100, 100, &
                      100, 100, 100, 100, 100, 100, 100, &
                      100/)
            ysize = (/100, 100, 100, 100, 100, 100, 100, &
                      100, 100, 100, 100, 100, 100, 100, &
                      100/)
!------------------
            do j = 1, n2Dslices(3)
                call write_metadata(analysisID, pub, label, citekey, description, expt, lumi, energy, particle_type, &
                                    mass_sep, mass_range_tolerance, masses(j), xmin(j), xmax(j), xstep(j), xsize(j), ymin(j), &
                                    ymax(j), ystep(j), ysize(j), llhdata(3)%D2llhdata(j))

                prefix = 'L_data_b_'
                call read2Ddata(masses(j), llhdata(3)%D2llhdata(j), stem, prefix, 'obs', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
                prefix = 'L_asimov_b_'
                call read2Ddata(masses(j), llhdata(3)%D2llhdata(j), stem, prefix, 'pred', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
            enddo

            deallocate (masses, xmin, xmax)
            deallocate (ymin, ymax, xstep)
            deallocate (ystep, xsize, ysize)

!------Setup forth likelihood -->  ATLAS 2002.12223

            analysisID = 200212223
            pub = .True.
            label = 'arXiv:2002.12223 [hep-ex] (ATLAS)'
            citekey = 'Aad:2020zxo'
            stem = 'Expt_tables/ATLtables/ATL_tautau_llh_2002.12223/'
            description = '(pp) -> phi -> tautau, using -2ln(L) reconstruction'
            expt = 'ATL'
            particle_type = Hneut
            lumi = 139.0D0
            energy = 13.0D0
            mass_sep = 0.20D0
            mass_range_tolerance = 10.0D0

            allocate (masses(n2Dslices(4)), xmin(n2Dslices(4)), xmax(n2Dslices(4)))
            allocate (ymin(n2Dslices(4)), ymax(n2Dslices(4)), xstep(n2Dslices(4)))
            allocate (ystep(n2Dslices(4)), xsize(n2Dslices(4)), ysize(n2Dslices(4)))

            masses = (/200.0D0, 250.0D0, 300.0D0, 350.0D0, 400.0D0, 500.0D0, 600.0D0, &
                       700.0D0, 800.0D0, 1000.0D0, 1200.0D0, 1500.0D0, 2000.0D0, &
                       2500.0D0/)

            !--------TODO:  Settings for the 200212223 analysis ---------!
            xmin = (/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
            xmax = (/800., 500., 500., 300., 100., 100., 30., 20., 10., 6., 5., 3., 2., 2./)
            xsize = (/100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100/)
            do i = lbound(xmin, dim=1), ubound(xmax, dim=1)
                xstep(i) = (xmax(i) - xmin(i))/(xsize(i) - 1)
            enddo
            ymin = (/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
            ymax = (/800., 500., 500., 300., 100., 100., 30., 20., 10., 6., 5., 3., 2., 2./)
            ysize = (/100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100/)
            do i = lbound(ymin, dim=1), ubound(ymax, dim=1)
                ystep(i) = (ymax(i) - ymin(i))/(ysize(i) - 1)
            enddo
            do j = 1, n2Dslices(4)
                call write_metadata(analysisID, pub, label, citekey, description, expt, &
                                    lumi, energy, particle_type, mass_sep, mass_range_tolerance, &
                                    masses(j), xmin(j), xmax(j), xstep(j), xsize(j), ymin(j), &
                                    ymax(j), ystep(j), ysize(j), llhdata(4)%D2llhdata(j))
                if (preventOverexclusion) then
                    prefix = 'L_obs_'
                else
                    prefix = 'L_obs_orig_'
                endif
                call read2Ddata(masses(j), llhdata(4)%D2llhdata(j), stem, prefix, 'obs', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
                prefix = 'L_exp_'
                call read2Ddata(masses(j), llhdata(4)%D2llhdata(j), stem, prefix, 'pred', status)
                if (status .ne. 0) stop 'Error in setup_likelihoods: Data files not found!'
            enddo

            deallocate (masses, xmin, xmax)
            deallocate (ymin, ymax, xstep)
            deallocate (ystep, xsize, ysize)

!--< NEW

!  write(*,*) 'Creating binaries...'
#ifndef WEBVERSION
            rewind (file_id_common4)
            write (file_id_common4) llhdata
#endif
        endif

        close (file_id_common4)
!--------------------------------------------------------------
    end subroutine setup_likelihoods
!--------------------------------------------------------------
    subroutine fudgelimit_CMS14029(fudge)
        implicit none
        logical, intent(in) :: fudge

        fudgelimit = fudge

    end subroutine fudgelimit_CMS14029
!--------------------------------------------------------------
    subroutine write_metadata(analysisID, pub, label, citekey, description, expt, lumi, energy, particle_type, &
                              mass_sep, mass_range_tolerance, mass, xmin, xmax, xstep, xsize, ymin, &
                              ymax, ystep, ysize, llh2D)
!--------------------------------------------------------------
        implicit none

        integer, intent(in) :: analysisID, particle_type
        logical, intent(in) :: pub
        character(LEN=*), intent(in) :: label, citekey, description, expt
        double precision, intent(in) :: mass_sep, mass, xmin, xmax, ymin, ymax, xstep, ystep, lumi, &
            energy, mass_range_tolerance
        integer, intent(in) :: xsize, ysize
        type(likelihood2D), intent(inout) :: llh2D

        llh2D%analysisID = analysisID
        llh2D%pub = pub
        llh2D%label = trim(adjustl(label))
        llh2D%citekey = trim(adjustl(citekey))
        llh2D%description = trim(adjustl(description))
        llh2D%expt = expt
        llh2D%energy = energy
        llh2D%lumi = lumi
        llh2D%particle_x = particle_type
        llh2D%mass_sep = mass_sep
        llh2D%mass_range_tolerance = mass_range_tolerance
        llh2D%mass = mass
        llh2D%xmin = xmin
        llh2D%xmax = xmax
        llh2D%xstep = xstep
        llh2D%xsize = xsize
        llh2D%ymin = ymin
        llh2D%ymax = ymax
        llh2D%ystep = ystep
        llh2D%ysize = ysize

!  if(abs((xmax-xmin)/xstep+1.-xsize).gt.small) then
!   write(*,*) 'Warning: xsize does not match: ',(xmax-xmin)/xstep+1,' vs. ', xsize
!  endif
!  if(abs((ymax-ymin)/ystep+1.-ysize).gt.small) then
!   write(*,*) 'Warning: ysize does not match: ',(ymax-ymin)/ystep+1,' vs. ', ysize
!  endif
!--------------------------------------------------------------
    end subroutine write_metadata
!--------------------------------------------------------------
    subroutine read2Ddata(mass, llh2D, stem, prefix, obspred, status)
!--------------------------------------------------------------
        use usefulbits, only: file_id_common2, vvsmall
        use install_data, only: pathname
        implicit none

        double precision, intent(in) :: mass
        type(likelihood2D), intent(inout) :: llh2D
        character(LEN=*), intent(in) :: stem, obspred, prefix
        integer, intent(out) :: status

        integer :: size, i, k, remove, posx, posy
        double precision, dimension(1:3) :: values
        type(threetuple), allocatable :: dummylikelihood(:), likelihood(:)
        double precision, dimension(200, 200) :: llhgrid
        double precision, dimension(1:3) :: bestfit
        character(LEN=4) :: intstring
        character(LEN=100) ::  filename

        bestfit = 0.0D0
        write (intstring, "(I4)") int(mass)
!   if(trim(adjustl(obspred)).eq.'obs') then
!   filename = 'L_data_SMHb_'//trim(adjustl(intstring))//'.out'
        filename = trim(adjustl(prefix))//trim(adjustl(intstring))//'.out'
!   else if(trim(adjustl(obspred)).eq.'pred') then
!   filename = 'L_asimovSMH_SMHb_'//trim(adjustl(intstring))//'.out'
!    filename = trim(adjustl(prefix))//trim(adjustl(intstring))//'.out'
!   else
!    stop 'Error in read2Ddata: obspred unknown.'
!   endif

        llh2D%mass = mass

        call get_length_of_file(trim(adjustl(pathname))//trim(adjustl(stem))// &
                                trim(adjustl(filename)), size, status)
        if (status .ne. 0) return

        open (file_id_common2, file=trim(adjustl(pathname))//trim(adjustl(stem))// &
              trim(adjustl(filename)))

        allocate (dummylikelihood(size))

        k = 0
        remove = 0
        do i = 1, size
            read (file_id_common2, *) values
            if (abs(values(3)) .lt. vvsmall) then
                bestfit = values
                remove = remove + 1
            else
                k = k + 1
                dummylikelihood(k)%val = values
            endif
        enddo

        close (file_id_common2)

        llh2D%size = k
!   write(*,*) 'size = ',llh2D%size

        allocate (likelihood(k))

        do i = 1, k
            likelihood(i)%val = dummylikelihood(i)%val
        enddo

        deallocate (dummylikelihood)

        ! write(*,*) "Reading ", obspred, "grid of ",llh2D%analysisID," with mass = ", mass
        llhgrid = 0.0D0

        do i = lbound(likelihood, dim=1), ubound(likelihood, dim=1)

            posx = nint((likelihood(i)%val(1) - llh2D%xmin)/llh2D%xstep) + 1
            posy = nint((likelihood(i)%val(2) - llh2D%ymin)/llh2D%ystep) + 1
            ! write(*,*) i, likelihood(i)%val(1), llh2D%xmin, posx
            ! write(*,*) i, likelihood(i)%val(2), llh2D%ymin, posy
! DEBUGGING:
            if (posx .gt. llh2D%xsize) then
                write (*, *) "Error in x direction: ", mass, i, likelihood(i)%val(1), llh2D%xmin, llh2D%xstep, posx
            endif
            if (posy .gt. llh2D%ysize) then
                write (*, *) "Error in y direction: ", mass, i, likelihood(i)%val(2), llh2D%ymin, llh2D%ystep, posy
            endif
!--------
            llhgrid(posx, posy) = likelihood(i)%val(3)
        enddo

!  k=0
!  do i=1,llh2D%xsize
!   do j=1, llh2D%ysize
!    if(llhgrid(i,j).lt.small) then
!     k=k+1
!!     write(*,*) i,j,llh2D%llhgrid(i,j)
!    endif
!   enddo
!  enddo
!  write(*,*) k, 'grid points are unfilled!'

        if (trim(adjustl(obspred)) .eq. 'obs') then
            llh2D%llhgrid_obs = llhgrid
            llh2D%bestfit_obs = bestfit
        else if (trim(adjustl(obspred)) .eq. 'pred') then
            llh2D%llhgrid_pred = llhgrid
            llh2D%bestfit_pred = bestfit
        endif

        deallocate (likelihood)
!--------------------------------------------------------------
    end subroutine read2Ddata
    !--------------------------------------------------------------
    subroutine get_likelihood(analysisID, jj, t, llh, M_av, nc, cbin, obspred, cbin_in)
        !--------------------------------------------------------------
        use usefulbits, only: dataset
        implicit none

        type(dataset), intent(in) :: t
        integer, intent(in) :: analysisID, jj
        character(LEN=*), intent(in) :: obspred
        double precision, intent(out) :: llh, M_av
        integer, intent(out) :: nc, cbin
        integer, optional, intent(in) :: cbin_in
        double precision :: cfact1, cfact2
        integer :: kk
        logical :: analysisfound = .False.

        llh = 0.0D0

        do kk = 1, nllhs
            if (llhdata(kk)%D2llhdata(1)%analysisID .eq. analysisID) then
                analysisfound = .True.
                if (present(cbin_in)) then
!      write(*,*) "Calling get_likelihood with cbin_in = ", cbin_in, ", obspred = ", obspred
                    call calcfact_llh(kk, jj, t, cfact1, cfact2, M_av, nc, cbin, cbin_in)
                else
                    call calcfact_llh(kk, jj, t, cfact1, cfact2, M_av, nc, cbin)
                endif
                if (cbin .ne. 0) then
                    llh = get_llh_ggbb_phi_tautau(kk, M_av, cfact1, cfact2, obspred)
                endif
            endif
        enddo

        if (.not. analysisfound) then
            write (*, *) "WARNING: No likelihood information available for analysis ", analysisID
        endif

    end subroutine get_likelihood
! ======================================
!     calculate CLs from likelihoods
! ======================================
    !! expected CL_{S+B}
    !! calculated using eq (65) of 1007.1727 for \f$\tilde{q}_mu\f$ = q_obs = \f$\tilde{q}_A\f$ = q_exp
    function CLsb_expected(q_exp)
        implicit none
        double precision, intent(in) :: q_exp
        double precision cum, ccum
        double precision CLsb_expected
        call cumnor(sqrt(q_exp), cum, ccum)
        CLsb_expected = ccum
    end function

    !! expected CL_S
    !! using CLsb_expected() and CLb_expected = 0.5
    function CLs_expected(q_exp)
        implicit none
        double precision, intent(in) :: q_exp
        double precision, parameter :: inverse_CLb_expected = 2.D0
        double precision CLs_expected
        CLs_expected = CLsb_expected(q_exp)*inverse_CLb_expected
    end function

    !! observed CL_{S+B}
    !! calculated using eq (65) of 1007.1727 for \f$\tilde{q}_\mu\f$ = q_obs and \f$(\mu/\sigma)^2\f$ = q_exp
    function CLsb_observed(q_exp, q_obs)
        implicit none
        double precision, intent(in) :: q_exp, q_obs
        double precision cum, ccum
        double precision CLsb_observed
        if (q_obs <= q_exp) then
            call cumnor(sqrt(q_obs), cum, ccum)
        else
            call cumnor((q_obs + q_exp)/(2*sqrt(q_exp)), cum, ccum)
        endif
        CLsb_observed = ccum
    end function

    !! observed CL_B
    !! calculated using eq (64) of 1007.1727 for \f$\mu' = 0\f$, \f$\tilde{q}_\mu\f$ = q_obs and \f$(\mu/\sigma)^2\f$ = q_exp
    function CLb_observed(q_exp, q_obs)
        implicit none
        double precision, intent(in) :: q_exp, q_obs
        double precision cum, ccum
        double precision CLb_observed
        if (q_obs <= q_exp) then
            call cumnor(sqrt(q_obs) - sqrt(q_exp), cum, ccum)
        else
            call cumnor((q_obs - q_exp)/(2*sqrt(q_exp)), cum, ccum)
        endif
        CLb_observed = ccum
    end function

    !! observed CL_S
    !! calculated from CLsb_observed and CLb_observed
    function CLs_observed(q_exp, q_obs)
        use usefulbits, only: div, vsmall
        implicit none
        double precision, intent(in) :: q_exp, q_obs
        double precision CLs_observed
        if (q_exp < 0.1) then ! possibly small q_exp limit
            CLs_observed = div(CLsb_observed(q_exp, q_obs), CLb_observed(q_exp, q_obs), 1.D0, 1.D0/vsmall)
        else
            CLs_observed = CLsb_observed(q_exp, q_obs)/(CLb_observed(q_exp, q_obs) + vsmall)
        endif

    end function

    !--------------------------------------------------------------
    subroutine calcpredratio_llh(c, jj, t, M_av, fact, nc, predratio)
        !--------------------------------------------------------------
        use usefulbits, only: dataset, small
        implicit none

        !prsep(h,n)%tlist,prsep(h,n)%findj,t,axis_i(n),ncomb(n),predratio(n)
        type(dataset), intent(in) :: t
        integer, intent(in) :: c, jj
        double precision, intent(out) :: M_av, predratio, fact
        integer, intent(out) ::  nc

        double precision :: cfact1, cfact2
        double precision :: sf, expllh, expllh0
        integer :: kk, cbin

        double precision sf_low, sf_high, llh_low, llh_high, sf_max

        call calcfact_llh(c, jj, t, cfact1, cfact2, M_av, nc, cbin)
        fact = cfact1 + cfact2

        expllh = get_llh_ggbb_phi_tautau(c, M_av, cfact1, cfact2, 'pred')
        expllh0 = expllh
!  write(*,*) 'c = ',c,' jj = ',jj,'nc = ',nc,' M_av = ', M_av, 'ggH = ',cfact1, 'bbH = ',cfact2

        if (cfact1 .lt. small .and. cfact2 .lt. small) then
!   write(*,*) 'Very small cross sections. CMS tautau is not competitive, expected llh =',expllh
            predratio = -1.0D0
            return
        endif

!  write(*,*) "Find predicted ratio..."
        if (expllh .lt. 0.0D0) then
            predratio = -1.0D0
            return
        endif

! OS implementation of bisection
! Maximum scale factor corresponds to min_predratio=1/sf_max
        sf_max = 1.d0
        llh_high = get_llh_ggbb_phi_tautau(c, M_av, sf_max*cfact1, sf_max*cfact2, 'pred')
        kk = 0
        do while (CLs_expected(llh_high) .gt. 5D-2 .and. kk .lt. 100)
            sf_max = sf_max*10.d0
            llh_high = get_llh_ggbb_phi_tautau(c, M_av, sf_max*cfact1, sf_max*cfact2, 'pred')
            kk = kk + 1
        enddo

        sf_low = 0.d0
        sf_high = sf_max
        kk = 0

        sf = (sf_high + sf_low)/2.

        llh_low = get_llh_ggbb_phi_tautau(c, M_av, sf_low*cfact1, sf_low*cfact2, 'pred')
        llh_high = get_llh_ggbb_phi_tautau(c, M_av, sf_high*cfact1, sf_high*cfact2, 'pred')

        ! Check that root exists in starting interval
        if (CLs_expected(llh_low) .gt. 5D-2 .and. CLs_expected(llh_high) .lt. 5D-2) then
            expllh = get_llh_ggbb_phi_tautau(c, M_av, sf*cfact1, sf*cfact2, 'pred')

!   write(*,*) kk, sf_low, sf_high, sf, expllh
            do while (abs(CLs_expected(expllh) - 5D-2) .gt. 1D-4 .and. kk .lt. 100)

                ! Five digit precision on scale factor in case of oscillating solutions
                if ((sf_high - sf_low) .LT. 1D-5) exit

                kk = kk + 1
                if (CLs_expected(expllh) .gt. 5D-2) sf_low = sf
                if (CLs_expected(expllh) .lt. 5D-2) sf_high = sf
                sf = (sf_high + sf_low)/2.

                expllh = get_llh_ggbb_phi_tautau(c, M_av, sf*cfact1, sf*cfact2, 'pred')

!    write(*,*) kk, sf_low, sf_high, sf, expllh, CLs_expected(expllh)
            enddo
        else
            sf = sf_max
        endif

        if (kk .ge. 100) write (*, *) 'Warning: Maximum number of iterations reached with no convergence:', &
            'predratio might be unreliable.'
!      write(*,*) kk, sf_low, sf_high, sf, expllh
        if (sf .le. 0.0D0) sf = 1.0D-10

        predratio = 1./sf

    end subroutine calcpredratio_llh
    !--------------------------------------------------------------
    subroutine check_against_bound_llh(c, jj, t, M_av, nc, obsratio)
        !--------------------------------------------------------------
        use usefulbits, only: dataset, vsmall
        implicit none

        type(dataset), intent(in) :: t
        integer, intent(in) :: c, jj, nc
        double precision, intent(in) :: M_av
        double precision, intent(out) :: obsratio

        double precision :: cfact1, cfact2
        double precision :: sf, obsllh, expllh, M_av2, interval
        double precision :: diff(2)
        integer :: kk, nc2, cbin

        double precision sf_low, sf_high, llh_low, llh_high, llh_exp_low, llh_exp_high, sf_max

        call calcfact_llh(c, jj, t, cfact1, cfact2, M_av2, nc2, cbin)

        if (abs(M_av2 - M_av) .gt. vsmall .or. nc .ne. nc2) then
            write (*, *) M_av2, M_av, nc, nc2
            stop 'Error in subroutine check_against_bound_llh !'

        endif

        obsllh = get_llh_ggbb_phi_tautau(c, M_av, cfact1, cfact2, 'obs')
!  write(*,*) "Find observed ratio..."
        if (obsllh .lt. 0.0D0) then
            obsratio = -1.0D0
            return
        endif

        sf = 1.0D0
        kk = 0
        interval = 0.50D0
        diff = 0.0D0

!! OS implementation of bisection
        sf_max = 1.d0
        llh_high = get_llh_ggbb_phi_tautau(c, M_av, sf_max*cfact1, sf_max*cfact2, 'obs')
        llh_exp_high = get_llh_ggbb_phi_tautau(c, M_av, sf_max*cfact1, sf_max*cfact2, 'pred')
        kk = 0
        do while (CLs_observed(llh_exp_high, llh_high) .gt. 5D-2 .and. kk .lt. 100)
            sf_max = sf_max*10.d0
            llh_high = get_llh_ggbb_phi_tautau(c, M_av, sf_max*cfact1, sf_max*cfact2, 'obs')
            llh_exp_high = get_llh_ggbb_phi_tautau(c, M_av, sf_max*cfact1, sf_max*cfact2, 'pred')
            kk = kk + 1
        enddo

        sf_low = 0.d0
        sf_high = sf_max
        kk = 0

        sf = (sf_high + sf_low)/2.

        llh_low = get_llh_ggbb_phi_tautau(c, M_av, sf_low*cfact1, sf_low*cfact2, 'obs')
        llh_high = get_llh_ggbb_phi_tautau(c, M_av, sf_high*cfact1, sf_high*cfact2, 'obs')
        llh_exp_low = get_llh_ggbb_phi_tautau(c, M_av, sf_low*cfact1, sf_low*cfact2, 'pred')
        llh_exp_high = get_llh_ggbb_phi_tautau(c, M_av, sf_high*cfact1, sf_high*cfact2, 'pred')
        ! Check that root exists in starting interval
        if (CLs_observed(llh_exp_low, llh_low) .gt. 5D-2 .and. CLs_observed(llh_exp_high, llh_high) .lt. 5D-2) then
            obsllh = get_llh_ggbb_phi_tautau(c, M_av, sf*cfact1, sf*cfact2, 'obs')
            expllh = get_llh_ggbb_phi_tautau(c, M_av, sf*cfact1, sf*cfact2, 'pred')

            do while (abs(CLs_observed(expllh, obsllh) - 5D-2) .gt. 1D-4 .and. kk .lt. 100)

                ! Five digit precision on scale factor in case of oscillating solutions
                if ((sf_high - sf_low) .LT. 1D-5) exit

                kk = kk + 1
                if (CLs_observed(expllh, obsllh) .gt. 5D-2) sf_low = sf
                if (CLs_observed(expllh, obsllh) .lt. 5D-2) sf_high = sf
                sf = (sf_high + sf_low)/2.

                obsllh = get_llh_ggbb_phi_tautau(c, M_av, sf*cfact1, sf*cfact2, 'obs')
                expllh = get_llh_ggbb_phi_tautau(c, M_av, sf*cfact1, sf*cfact2, 'pred')
            enddo
        else
            sf = sf_max
        endif

        if (kk .ge. 100) write (*, *) 'Warning: Maximum number of iterations reached with no convergence:', &
            'obsratio might be unreliable.'

        if (sf .le. 0.0D0) sf = 1.0D-10
        obsratio = 1./sf

    end subroutine check_against_bound_llh
!--------------------------------------------------------------
    subroutine calcfact_llh(c, jj, t, cfact1, cfact2, M_av, nc, cbin, cbin_in)
        !--------------------------------------------------------------
        ! This routine calculates the predicted rates for each process
        ! c is the number of the likelihood, following the order in the setup
        use usefulbits, only: dataset, np, div, vvsmall, small
        use theory_BRfunctions
        use theory_XS_SM_functions

        implicit none
        !--------------------------------------input
        type(dataset), intent(in) :: t
        integer, intent(in) :: c, jj
        integer, optional, intent(in) :: cbin_in
        !-----------------------------------output
        double precision, intent(out) :: cfact1, cfact2, M_av
        integer, intent(out) :: nc, cbin
        !-----------------------------------internal
        double precision, allocatable :: mass(:), fact1(:), fact2(:)
        integer :: npart, f, j !number of particles
        double precision :: numerator, denominator
        logical, allocatable :: Havail(:)

        cfact1 = 0.0D0
        cfact2 = 0.0D0

!   if(c.eq.1) then ! this is the CMS tautau likelihood
        npart = np(llhdata(c)%D2llhdata(1)%particle_x)
        allocate (mass(npart), fact1(npart), fact2(npart))

        mass(:) = t%particle(llhdata(c)%D2llhdata(1)%particle_x)%M(:)
        fact1 = 0.0D0
        fact2 = 0.0D0

        allocate (Havail(npart))

        if (present(cbin_in)) then
            call available_Higgses(Havail, npart, cbin_in)
        else
            Havail = .True.
        endif

!   write(*,*) 'calcfact_llh: mass = ', mass
        cbin = 0

        !write(*,*) "Calling calcfact_llh..."

        do j = 1, npart
            if (Havail(j) .and. Havail(jj)) then
!   if(abs(mass(jj)-mass(j)).le.CMS_llhdata(1)%mass_sep*mass(jj))then
!   if((abs(mass(jj)-mass(j)).le.CMS_llhdata(1)%mass_sep*100.0D0)) then
                if (abs(mass(jj) - mass(j)) .le. llhdata(c)%D2llhdata(1)%mass_sep*max(mass(jj), mass(j))) then
                    if (abs(llhdata(c)%D2llhdata(1)%energy - 8.0D0) .le. small) then
                        fact1(j) = t%lhc8%channelrates(j, 6, 4)*t%lhc8%XS_gg_H_SM(j)
                        fact2(j) = t%lhc8%channelrates(j, 7, 4)*t%lhc8%XS_bb_H_SM(j)

                        if (fudgelimit) then
!-- Using an NLO acceptance factor (-50%) for CMS 14-029 in the gg->phi cross section for mphi < 200 GeV.
                            if ((llhdata(c)%D2llhdata(1)%analysisID .eq. 14029) .and. (mass(j) .le. 200.0D0)) then
                                fact1(j) = t%lhc8%channelrates(j, 6, 4)*t%lhc8%XS_gg_H_SM(j)*0.5D0
                            endif
!-- Using an NLO acceptance factor (-40%) for CMS 14-029 in the bb->phi cross section.
                            if ((llhdata(c)%D2llhdata(1)%analysisID .eq. 14029)) then
                                fact2(j) = t%lhc8%channelrates(j, 7, 4)*t%lhc8%XS_bb_H_SM(j)*0.6D0
                            endif
                        endif
!--
!        fact1(j)=t%lhc8%XS_gg_hj_ratio(j)*XS_lhc8_gg_H_SM(mass(j))*t%BR_hjtautau(j)
!        fact2(j)=t%lhc8%XS_bb_hj_ratio(j)*XS_lhc8_bb_H_SM(mass(j))*t%BR_hjtautau(j)
                    else if (abs(llhdata(c)%D2llhdata(1)%energy - 13.0D0) .le. small) then
!        fact1(j)=t%lhc13%XS_gg_hj_ratio(j)*XS_lhc13_gg_H_SM(mass(j))*t%BR_hjtautau(j)
!        fact2(j)=t%lhc13%XS_bb_hj_ratio(j)*XS_lhc13_bb_H_SM(mass(j))*t%BR_hjtautau(j)
                        if (llhdata(c)%D2llhdata(1)%analysisID == 200212223) then
                            ! Convert from pb to fb
                            fact1(j) = 1000.0D0*t%lhc13%channelrates(j, 6, 4)*t%lhc13%XS_gg_H_SM(j)
                            fact2(j) = 1000.0D0*t%lhc13%channelrates(j, 7, 4)*t%lhc13%XS_bb_H_SM(j)
                        else
                            fact1(j) = t%lhc13%channelrates(j, 6, 4)*t%lhc13%XS_gg_H_SM(j)
                            fact2(j) = t%lhc13%channelrates(j, 7, 4)*t%lhc13%XS_bb_H_SM(j)
                        endif
                    else
                        write (*, *) 'WARNING: No likelihood data available at this energy!'
                    endif
                    cbin = cbin + 2**(j - 1)
!     write(*,*) cbin
                endif
            endif
        enddo

!   write(*,*) 'calcfact_llh: fact1 = ', fact1
!   write(*,*) 'calcfact_llh: fact2 = ', fact2

        if (fact1(jj) .le. vvsmall .and. fact2(jj) .le. vvsmall) then
!Higgs jj doesn't contribute - wait until another call of this subroutine before
!looking at nearby masses
            M_av = mass(jj)
            nc = 0
            cfact1 = 0.0D0
            cfact2 = 0.0D0
        else
!find M_av weighted by the rates (only using higgs which have non-zero fact):
            f = 0
            numerator = 0.0D0
            denominator = 0.0D0
            do j = 1, npart
                if (fact1(j) .gt. vvsmall .or. fact2(j) .gt. vvsmall) then
                    f = f + 1
                    numerator = numerator + mass(j)*(fact1(j) + fact2(j))
                    denominator = denominator + (fact1(j) + fact2(j))
                endif
            enddo

            if (abs(denominator) .gt. 0.0D0) then
                M_av = numerator/denominator
            else
                M_av = 0.0D0
                write (*, *) 'Warning: could not find average mass for tautau analysis', llhdata(c)%D2llhdata(1)%analysisID
            endif

            nc = f !f will always be > 0 because we've already made sure that fact(jj)>0.0D0
            cfact1 = sum(fact1)
            cfact2 = sum(fact2)
        endif

        deallocate (mass, fact1, fact2)
        deallocate (Havail)

!   endif
    end subroutine calcfact_llh
!--------------------------------------------------------------
    function get_llh_from_interpolation(llh2D, xval0, yval0, obspred)
        ! This routine returns the likelihood from grid interpolation
!--------------------------------------------------------------
        use usefulbits, only: small
        type(likelihood2D), intent(in) :: llh2D
        double precision, intent(in) :: xval0, yval0
        character(LEN=*), intent(in) :: obspred
        double precision :: get_llh_from_interpolation
        integer :: posx1, posx2, posy1, posy2, posxlow, posxup
        double precision :: xval, yval, xval1, xval2, yval1, yval2, fvaly1, fvaly2, fval
        integer :: i, max_expansion
        logical :: lowerxfound, upperxfound
        double precision, dimension(200, 200) :: llhgrid
! Josh Sayre extension --
        double precision :: truex, truey
! --

        if (trim(adjustl(obspred)) .eq. 'obs') then
            llhgrid = llh2D%llhgrid_obs
        else if (trim(adjustl(obspred)) .eq. 'pred') then
            llhgrid = llh2D%llhgrid_pred
        else
            stop 'Error in get_llh_from_interpolation: obspred unknown.'
        endif

        ! This number is the maximal number of tries, that the algorithm steps further up or
        ! downwards the grid in case of non-existent grid points.
        max_expansion = 3

        if (isnan(xval0)) then
            xval = 0.0D0
        else
            xval = xval0
        endif

        if (isnan(yval0)) then
            yval = 0.0D0
        else
            yval = yval0
        endif

! write(*,*) 'input x,y vals = ', xval0, xval, yval0,  yval

! Josh Sayre extension --
        truex = xval
        truey = yval
! --

        if (xval .lt. llh2D%xmin) then
!  write(*,*) 'Warning: ggH rate outside grid -- value too small.'
            xval = llh2D%xmin
        else if (xval .gt. llh2D%xmax) then
!  write(*,*) 'Warning: ggH rate outside grid -- value too large.'
            xval = llh2D%xmax - small
        endif
        if (yval .lt. llh2D%ymin) then
!  write(*,*) 'Warning: bbH rate outside grid -- value too small.'
            yval = llh2D%ymin
        else if (yval .gt. llh2D%ymax) then
!  write(*,*) 'Warning: bbH rate outside grid -- value too large.'
            yval = llh2D%ymax - small
        endif

! Get coordinates of neighboring points
        posx1 = int((xval - llh2D%xmin)/llh2D%xstep) + 1
        if (posx1 .lt. llh2D%xsize) then
            posx2 = int((xval - llh2D%xmin)/llh2D%xstep) + 2
        else
!   write(*,*) "Warning: on the edge of the grid in x-position!"
            posx1 = int((xval - llh2D%xmin)/llh2D%xstep)
            posx2 = int((xval - llh2D%xmin)/llh2D%xstep) + 1
        endif

        posy1 = int((yval - llh2D%ymin)/llh2D%ystep) + 1
        if (posy1 .lt. llh2D%ysize) then
            posy2 = int((yval - llh2D%ymin)/llh2D%ystep) + 2
        else
!   write(*,*) "Warning: on the edge of the grid in y-position!"
            posy1 = int((yval - llh2D%ymin)/llh2D%ystep)
            posy2 = int((yval - llh2D%ymin)/llh2D%ystep) + 1
        endif

        ! Get x,y values at these points
        xval1 = (posx1 - 1)*llh2D%xstep + llh2D%xmin
        xval2 = (posx2 - 1)*llh2D%xstep + llh2D%xmin
        yval1 = (posy1 - 1)*llh2D%ystep + llh2D%ymin
        yval2 = (posy2 - 1)*llh2D%ystep + llh2D%ymin

!  write(*,*) "llh2D steps, min: ", llh2D%xstep, llh2D%ystep, llh2D%xmin, llh2D%ymin
!  write(*,*) "Coordinates: ", xval1, xval2, yval1, yval2
!  write(*,*)  "Indices: ", posx1, posx2, posy1, posy2
!  write(*,*) "Values (1,1), (2,1): ", llhgrid(posx1,posy1), llhgrid(posx2,posy1)
!  write(*,*) "Values (1,2), (2,2): ", llhgrid(posx1,posy2), llhgrid(posx2,posy2)

        ! Do bilinear interpolation and retrieve likelihood value
        if (abs(llhgrid(posx1, posy1)) .gt. 0.0D0 .and. abs(llhgrid(posx2, posy1)) .gt. 0.0D0) then
! Josh Sayre extension --
            fvaly1 = (xval2 - truex)/(xval2 - xval1)*llhgrid(posx1, posy1) + &
                     (truex - xval1)/(xval2 - xval1)*llhgrid(posx2, posy1)
!   write(*,*) "fvaly1 = ", fvaly1
! --
!   fvaly1 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posx1,posy1) + &
!   &        (xval - xval1) / (xval2 - xval1) * llhgrid(posx2,posy1)
        else if (abs(llhgrid(posx1, posy1)) .gt. 0.0D0) then
            fvaly1 = llhgrid(posx1, posy1)
        else if (abs(llhgrid(posx2, posy1)) .gt. 0.0D0) then
            fvaly1 = llhgrid(posx2, posy1)
        else
            write (*, *) 'Warning in interpolating grid: Both gridpoints in x-direction missing!'
            write (*, *) '(', xval1, yval1, '), (', xval2, yval1, ')'
            lowerxfound = .False.
            upperxfound = .False.
            do i = 1, max_expansion
                if (posx1 - i .gt. 0) then
                    if (abs(llhgrid(posx1 - i, posy1)) .gt. 0.0D0) then
                        xval1 = (posx1 - i - 1)*llh2D%xstep + llh2D%xmin
                        lowerxfound = .True.
                        posxlow = posx1 - i
                    endif
                endif
                if (posx2 + i .le. ubound(llhgrid, dim=1)) then
                    if (abs(llhgrid(posx2 + i, posy1)) .gt. 0.0D0) then
                        xval2 = (posx2 + i - 1)*llh2D%xstep + llh2D%xmin
                        upperxfound = .True.
                        posxup = posx2 + i
                    endif
                endif
                if (lowerxfound .or. upperxfound) exit
            enddo
            if (lowerxfound .and. upperxfound) then
! Josh Sayre extension --
                fvaly1 = (xval2 - truex)/(xval2 - xval1)*llhgrid(posxlow, posy1) + &
                         (truex - xval1)/(xval2 - xval1)*llhgrid(posxup, posy1)
! --
!     fvaly1 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posxlow,posy1) + &
!     &        (xval - xval1) / (xval2 - xval1) * llhgrid(posxup,posy1)
            else if (lowerxfound) then
                fvaly1 = llhgrid(posxlow, posy1)
            else if (upperxfound) then
                fvaly1 = llhgrid(posxup, posy1)
            else
                write (*, *) 'Found no surrounding points!'
                fvaly1 = 99999.0D0
            endif
        endif

        if (abs(llhgrid(posx1, posy2)) .gt. 0.0D0 .and. abs(llhgrid(posx2, posy2)) .gt. 0.0D0) then
! Josh Sayre extension --
            fvaly2 = (xval2 - truex)/(xval2 - xval1)*llhgrid(posx1, posy2) + &
                     (truex - xval1)/(xval2 - xval1)*llhgrid(posx2, posy2)
!   write(*,*) "fvaly2 = ", fvaly2

! --
!   fvaly2 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posx1,posy2) + &
!   &        (xval - xval1) / (xval2 - xval1) * llhgrid(posx2,posy2)
        else if (abs(llhgrid(posx1, posy2)) .gt. 0.0D0) then
            fvaly2 = llhgrid(posx1, posy2)
        else if (abs(llhgrid(posx2, posy2)) .gt. 0.0D0) then
            fvaly2 = llhgrid(posx2, posy2)
        else
            write (*, *) 'Warning in interpolating grid: Both gridpoints in x-direction missing!'
            write (*, *) '(', xval1, yval2, '), (', xval2, yval2, ')'
            lowerxfound = .False.
            upperxfound = .False.
            do i = 1, max_expansion
                if (posx1 - i .gt. 0) then
                    if (abs(llhgrid(posx1 - i, posy2)) .gt. 0.0D0) then
                        xval1 = (posx1 - i - 1)*llh2D%xstep + llh2D%xmin
                        lowerxfound = .True.
                        posxlow = posx1 - i
                    endif
                endif
                if (posx2 + i .le. ubound(llhgrid, dim=1)) then
                    if (abs(llhgrid(posx2 + i, posy2)) .gt. 0.0D0) then
                        xval2 = (posx2 + i - 1)*llh2D%xstep + llh2D%xmin
                        upperxfound = .True.
                        posxup = posx2 + i
                    endif
                endif
!    write(*,*) i, lowerxfound, upperxfound
                if (lowerxfound .or. upperxfound) exit
            enddo
            if (lowerxfound .and. upperxfound) then
! Josh Sayre extension --
                fvaly2 = (xval2 - truex)/(xval2 - xval1)*llhgrid(posxlow, posy2) + &
                         (truex - xval1)/(xval2 - xval1)*llhgrid(posxup, posy2)
! --
!     fvaly2 = (xval2 - xval) / (xval2 - xval1) * llhgrid(posxlow,posy2) + &
!     &        (xval - xval1) / (xval2 - xval1) * llhgrid(posxup,posy2)
            else if (lowerxfound) then
                fvaly2 = llhgrid(posxlow, posy2)
            else if (upperxfound) then
                fvaly2 = llhgrid(posxup, posy2)
            else
!    write(*,*) 'Found no surrounding points!'
                fvaly2 = 99999.0D0
            endif
        endif

        fval = (yval2 - yval)/(yval2 - yval1)*fvaly1 + &
               (yval - yval1)/(yval2 - yval1)*fvaly2

        ! we want to obtain - 2 * Delta log(likelihood) \simeq chi^2
        !  get_llh_from_interpolation = 2.0D0*fval

        ! Change(20/11/2017,TS): Just retrieve the grid value here, modify later in get_llh_ggbb_phi_tautau if needed.
        get_llh_from_interpolation = fval
!   write(*,*) "get_llh_from_interpolation = ", get_llh_from_interpolation

        return

    end function get_llh_from_interpolation
!--------------------------------------------------------------
    function get_llh_ggbb_phi_tautau(c, mpred, rate1, rate2, obspred)
!--------------------------------------------------------------
        use theory_XS_SM_functions, only: XS_lhc8_gg_H_SM, XS_lhc8_bb_H_SM, XS_lhc13_gg_H_SM, XS_lhc13_bb_H_SM
        use theory_BRfunctions, only: BRSM_Htautau
        use usefulbits, only: small

        implicit none

        double precision, intent(in) :: mpred, rate1, rate2
        character(LEN=*), intent(in) :: obspred
        integer, intent(in) :: c

        double precision :: get_llh_ggbb_phi_tautau
        double precision :: ggH1, bbH1, llh1, ggH2, bbH2, llh2, mpred_tmp
        integer :: i

        get_llh_ggbb_phi_tautau = 0.0D0
        mpred_tmp = mpred

        if (mpred_tmp .lt. (llhdata(c)%D2llhdata(1)%mass - llhdata(c)%D2llhdata(1)%mass_range_tolerance)) then
!   write(*,*) 'Warning: predicted mass below lowest mass value of CMS analysis.'
            get_llh_ggbb_phi_tautau = 0.0D0
        else if (mpred_tmp .gt. (llhdata(c)%D2llhdata(llhdata(c)%n2Dslices)%mass &
                                 + llhdata(c)%D2llhdata(1)%mass_range_tolerance)) then
!   write(*,*) 'Warning: predicted mass above highest mass value of CMS analysis.'
            get_llh_ggbb_phi_tautau = 0.0D0
        else
            do i = 1, llhdata(c)%n2Dslices - 1
                if (mpred_tmp .lt. llhdata(c)%D2llhdata(1)%mass) then
                    mpred_tmp = llhdata(c)%D2llhdata(1)%mass
                else if (mpred_tmp .gt. llhdata(c)%D2llhdata(llhdata(c)%n2Dslices)%mass) then
                    mpred_tmp = llhdata(c)%D2llhdata(llhdata(c)%n2Dslices)%mass
                endif
                if (((mpred_tmp - llhdata(c)%D2llhdata(i)%mass) .ge. 0.0D0) &
                    .and. (mpred_tmp - llhdata(c)%D2llhdata(i + 1)%mass .lt. 0.0D0)) then
                    ! Rescale to grid mass values by SM rate ratios
                    if (rescale) then
                        if ((llhdata(c)%D2llhdata(1)%energy - 8.0D0) .le. small) then
                            ggH1 = rate1*XS_lhc8_gg_H_SM(llhdata(c)%D2llhdata(i)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i)%mass) &
                                   /(XS_lhc8_gg_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                            bbH1 = rate2*XS_lhc8_bb_H_SM(llhdata(c)%D2llhdata(i)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i)%mass) &
                                   /(XS_lhc8_bb_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                        else if ((llhdata(c)%D2llhdata(1)%energy - 13.0D0) .le. small) then
                            ggH1 = rate1*XS_lhc13_gg_H_SM(llhdata(c)%D2llhdata(i)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i)%mass) &
                                   /(XS_lhc13_gg_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                            bbH1 = rate2*XS_lhc13_bb_H_SM(llhdata(c)%D2llhdata(i)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i)%mass) &
                                   /(XS_lhc13_bb_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                        else
                            write (*, *) "WARNING: Unknown energy for likelihood data!"
                        endif
                    else
                        ggH1 = rate1
                        bbH1 = rate2
                    endif
                    llh1 = get_llh_from_interpolation(llhdata(c)%D2llhdata(i), ggH1, bbH1, obspred)

                    if (rescale) then
                        if ((llhdata(c)%D2llhdata(1)%energy - 8.0D0) .le. small) then
                            ggH2 = rate1*XS_lhc8_gg_H_SM(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   /(XS_lhc8_gg_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                            bbH2 = rate2*XS_lhc8_bb_H_SM(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   /(XS_lhc8_bb_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                        else if ((llhdata(c)%D2llhdata(1)%energy - 13.0D0) .le. small) then
                            ggH2 = rate1*XS_lhc13_gg_H_SM(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   /(XS_lhc13_gg_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                            bbH2 = rate2*XS_lhc13_bb_H_SM(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   *BRSM_Htautau(llhdata(c)%D2llhdata(i + 1)%mass) &
                                   /(XS_lhc13_bb_H_SM(mpred_tmp)*BRSM_Htautau(mpred_tmp))
                        else
                            write (*, *) "WARNING: Unknown energy for likelihood data!"
                        endif
                    else
                        ggH2 = rate1
                        bbH2 = rate2
                    endif

                    llh2 = get_llh_from_interpolation(llhdata(c)%D2llhdata(i + 1), ggH2, bbH2, obspred)

!      write(*,*) "Mass values = ", mpred_tmp, llhdata(c)%D2llhdata(i)%mass, llhdata(c)%D2llhdata(i+1)%mass
!      write(*,*) "Rates in Interpolation = ", ggH1, bbH1,ggH2, bbH2, llh1, llh2, i

                    get_llh_ggbb_phi_tautau = llh1 + (mpred_tmp - llhdata(c)%D2llhdata(i)%mass)/ &
                                              (llhdata(c)%D2llhdata(i + 1)%mass - llhdata(c)%D2llhdata(i)%mass)* &
                                              (llh2 - llh1)
                    exit
                else if (mpred_tmp - llhdata(c)%D2llhdata(i + 1)%mass .le. small) then
                    get_llh_ggbb_phi_tautau = get_llh_from_interpolation(llhdata(c)%D2llhdata(i + 1), &
                                                                         rate1, rate2, obspred)
!      write(*,*) "Rates = ", rate1, rate2, obspred, i+1, get_llh_ggbb_phi_tautau
                endif
            enddo
        endif

!-- CMS presents chi^2/2 values in data!
        if (llhdata(c)%D2llhdata(1)%analysisID .eq. 14029 .or. &
            llhdata(c)%D2llhdata(1)%analysisID .eq. 16037 .or. &
            llhdata(c)%D2llhdata(1)%analysisID .eq. 17020) then
            get_llh_ggbb_phi_tautau = get_llh_ggbb_phi_tautau*2.0D0
        endif
!--

!   write(*,*)  "Likelihood (channel, obspred) = ",get_llh_ggbb_phi_tautau, " (",c,",",obspred,")"
        return
!--------------------------------------------------------------
    end function get_llh_ggbb_phi_tautau
!--------------------------------------------------------------
    subroutine outputproc_llh(proc, file_id_common, descrip)
!--------------------------------------------------------------
        use usefulbits, only: listprocesses

        type(listprocesses), intent(in) :: proc
        integer, intent(in) :: file_id_common
        character(LEN=200), intent(out) :: descrip

        character(LEN=45) :: label
        integer :: i, j
        character(LEN=1) :: jj

        if (file_id_common .ne. 0) continue ! silence unused variable warning

        i = proc%findi
        j = proc%findj

        write (jj, '(I1)') j

        label = '('//trim(llhdata(proc%tlist)%D2llhdata(1)%label)//')'

        descrip = '(pp)->h'//jj//'->tautau, using -2ln(L) reconstruction  '//label

    end subroutine outputproc_llh
!--------------------------------------------------------------
    subroutine get_length_of_file(filename, n, status)
!--------------------------------------------------------------
        use usefulbits, only: file_id_common3
        implicit none

        character(LEN=*), intent(in) :: filename
        integer, intent(out) :: n, status
        integer :: error

        open (file_id_common3, file=trim(adjustl(filename)), action='read', status='old', iostat=status)
        if (status .ne. 0) then
            write (*, *) 'Bad status', status, 'with the following file:'
            write (*, *) trim(adjustl(filename))
        endif
        n = 0
        do
            read (file_id_common3, *, iostat=error)
            if (error == -1) exit
            n = n + 1
        enddo
        close (file_id_common3)
    end subroutine get_length_of_file
!--------------------------------------------------------------
    subroutine available_Higgses(Havail, nH, cbin_in)
        implicit none

        integer, intent(in) :: nH, cbin_in
        logical, intent(inout) :: Havail(:)

        integer :: c, i, r

        if (size(Havail, dim=1) .ne. nH) then
            write (*, *) "Warning in subroutine available_Higgses: sizes do not match!"
            Havail = .True.
        else
            c = cbin_in
            do i = 1, nH
                r = c/(2**(nH - i))
                if (r .eq. 1) then
                    Havail(nH - i + 1) = .False.
                elseif (r .eq. 0) then
                    Havail(nH - i + 1) = .True.
                else
                    stop "Error in subroutine available_Higgses: not valid binary!"
                endif
                c = mod(c, 2**(nH - i))
            enddo
!  write(*,*) "cbin = ",cbin_in, " Havail = ",Havail
        endif

    end subroutine available_Higgses

    subroutine cumnor(arg, cum, ccum)

        !*****************************************************************************80
        !
        !! CUMNOR computes the cumulative normal distribution.
        !
        !  Discussion:
        !
        !    This function evaluates the normal distribution function:
        !
        !                              / x
        !                     1       |       -t*t/2
        !          P(x) = ----------- |      e       dt
        !                 sqrt(2 pi)  |
        !                             /-oo
        !
        !    This transportable program uses rational functions that
        !    theoretically approximate the normal distribution function to
        !    at least 18 significant decimal digits.  The accuracy achieved
        !    depends on the arithmetic system, the compiler, the intrinsic
        !    functions, and proper selection of the machine dependent
        !    constants.
        !
        !  Author:
        !
        !    William Cody
        !    Mathematics and Computer Science Division
        !    Argonne National Laboratory
        !    Argonne, IL 60439
        !
        !  Reference:
        !
        !    William Cody,
        !    Rational Chebyshev approximations for the error function,
        !    Mathematics of Computation,
        !    1969, pages 631-637.
        !
        !    William Cody,
        !    Algorithm 715:
        !    SPECFUN - A Portable FORTRAN Package of Special Function Routines
        !    and Test Drivers,
        !    ACM Transactions on Mathematical Software,
        !    Volume 19, Number 1, 1993, pages 22-32.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) ARG, the upper limit of integration.
        !
        !    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
        !    complementary CDF.
        !
        !  Local Parameters:
        !
        !    Local, real ( kind = 8 ) EPS, the argument below which anorm(x)
        !    may be represented by 0.5 and above which  x*x  will not underflow.
        !    A conservative value is the largest machine number X
        !    such that   1.0D+00 + X = 1.0D+00   to machine precision.
        !
        implicit none

        double precision, parameter, dimension(5) :: a = (/2.2352520354606839287D+00, &
                                                           1.6102823106855587881D+02, &
                                                           1.0676894854603709582D+03, &
                                                           1.8154981253343561249D+04, &
                                                           6.5682337918207449113D-02/)
        double precision, parameter, dimension(4) :: b = (/4.7202581904688241870D+01, &
                                                           9.7609855173777669322D+02, &
                                                           1.0260932208618978205D+04, &
                                                           4.5507789335026729956D+04/)
        double precision, parameter, dimension(9) :: c = (/3.9894151208813466764D-01, &
                                                           8.8831497943883759412D+00, &
                                                           9.3506656132177855979D+01, &
                                                           5.9727027639480026226D+02, &
                                                           2.4945375852903726711D+03, &
                                                           6.8481904505362823326D+03, &
                                                           1.1602651437647350124D+04, &
                                                           9.8427148383839780218D+03, &
                                                           1.0765576773720192317D-08/)
        double precision, parameter, dimension(8) :: d = (/2.2266688044328115691D+01, &
                                                           2.3538790178262499861D+02, &
                                                           1.5193775994075548050D+03, &
                                                           6.4855582982667607550D+03, &
                                                           1.8615571640885098091D+04, &
                                                           3.4900952721145977266D+04, &
                                                           3.8912003286093271411D+04, &
                                                           1.9685429676859990727D+04/)
        double precision, parameter, dimension(6) :: p = (/2.1589853405795699D-01, &
                                                           1.274011611602473639D-01, &
                                                           2.2235277870649807D-02, &
                                                           1.421619193227893466D-03, &
                                                           2.9112874951168792D-05, &
                                                           2.307344176494017303D-02/)
        double precision, parameter, dimension(5) :: q = (/1.28426009614491121D+00, &
                                                           4.68238212480865118D-01, &
                                                           6.59881378689285515D-02, &
                                                           3.78239633202758244D-03, &
                                                           7.29751555083966205D-05/)
        double precision, parameter :: root32 = 5.656854248D+00
        double precision, parameter :: sixten = 16.0D+00
        double precision, parameter :: sqrpi = 3.9894228040143267794D-01
        double precision, parameter :: thrsh = 0.66291D+00
        double precision, intent(in) :: arg
        double precision, intent(out) :: ccum
        double precision, intent(out) :: cum
        double precision del
        double precision eps
        integer i
        double precision temp
        double precision x
        double precision xden
        double precision xnum
        double precision y
        double precision xsq
        !
        !  Machine dependent constants
        !
        eps = epsilon(1.0D+00)*0.5D+00

        x = arg
        y = abs(x)

        if (y <= thrsh) then
            !
            !  Evaluate  anorm  for  |X| <= 0.66291
            !
            if (eps < y) then
                xsq = x*x
            else
                xsq = 0.0D+00
            end if

            xnum = a(5)*xsq
            xden = xsq
            do i = 1, 3
                xnum = (xnum + a(i))*xsq
                xden = (xden + b(i))*xsq
            end do
            cum = x*(xnum + a(4))/(xden + b(4))
            temp = cum
            cum = 0.5D+00 + temp
            ccum = 0.5D+00 - temp
            !
            !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
            !
        else if (y <= root32) then

            xnum = c(9)*y
            xden = y
            do i = 1, 7
                xnum = (xnum + c(i))*y
                xden = (xden + d(i))*y
            end do
            cum = (xnum + c(8))/(xden + d(8))
            xsq = aint(y*sixten)/sixten
            del = (y - xsq)*(y + xsq)
            cum = exp(-xsq*xsq*0.5D+00)*exp(-del*0.5D+00)*cum
            ccum = 1.0D+00 - cum

            if (0.0D+00 < x) then
                temp = cum
                cum = ccum
                ccum = temp
            end if
            !
            !  Evaluate ANORM for sqrt(32) < |X|.
            !
        else

            cum = 0.0D+00
            xsq = 1.0D+00/(x*x)
            xnum = p(6)*xsq
            xden = xsq
            do i = 1, 4
                xnum = (xnum + p(i))*xsq
                xden = (xden + q(i))*xsq
            end do

            cum = xsq*(xnum + p(5))/(xden + q(5))
            cum = (sqrpi - cum)/y
            xsq = aint(x*sixten)/sixten
            del = (x - xsq)*(x + xsq)
            cum = exp(-xsq*xsq*0.5D+00) &
                  *exp(-del*0.5D+00)*cum
            ccum = 1.0D+00 - cum

            if (0.0D+00 < x) then
                temp = cum
                cum = ccum
                ccum = temp
            end if
        end if

        if (cum < tiny(cum)) then
            cum = 0.0D+00
        end if

        if (ccum < tiny(ccum)) then
            ccum = 0.0D+00
        end if

        return
    end subroutine
!******************************************************************
end module likelihoods
!******************************************************************
