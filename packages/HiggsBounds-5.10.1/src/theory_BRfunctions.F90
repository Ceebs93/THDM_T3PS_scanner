! This file is part of HiggsBounds
!  -KW
!******************************************************************
module theory_BRfunctions
! Note that these are the Standard Model branching ratios used internally by HiggsBounds
! NOT the Standard Model branching ratio functions available to the user
! (see access_SM.f90 for those)
!******************************************************************

    use interpolate
    use S95tables_type1

    implicit none

    !table type 1-----------------------------
    type(table1), allocatable :: BRSM(:)
    !------------------------------------------

    double precision :: BRSMt1Mhmax, BRSMt1Mhmin

contains

    !******************************************************
    subroutine setup_BRSM
        ! reads in the Standard Model Branching ratios and total decay width from files
        !******************************************************
        use install_data, only: pathname, pathname_length
        use usefulbits, only: file_id_common2
        implicit none
        !------------------------------------internal
        integer :: x, xbeg, xend
        character(len=100), allocatable :: filename(:)
        character(LEN=pathname_length + 150) :: fullfilename
        integer :: col, ios
        !--------------------------------------------

!   allocate(BRSM(2))
        allocate (BRSM(4))

        xbeg = lbound(BRSM, dim=1)
        xend = ubound(BRSM, dim=1)

        allocate (filename(xbeg:xend))
        x = xbeg - 1

        x = x + 1
        BRSM(x)%xmin = 0.0D0
        BRSM(x)%xmax = 1000.0D0
        BRSM(x)%sep = 0.1D0
!  filename(x)='br.sm1_HDecay'
!   filename(x)='br.sm1_LHCHiggsXSWG_HDecay'
        filename(x) = 'br.bsm1_YR4.dat'

        x = x + 1
        BRSM(x)%xmin = 0.0D0
        BRSM(x)%xmax = 1000.0D0
        BRSM(x)%sep = 0.1D0
!  filename(x)='br.sm2_HDecay'
!   filename(x)='br.sm2_LHCHiggsXSWG_HDecay'
        filename(x) = 'br.bsm2_YR4.dat'

        x = x + 1
        BRSM(x)%xmin = 120.0D0
        BRSM(x)%xmax = 130.0D0
        BRSM(x)%sep = 0.1D0
        filename(x) = 'br.sm1_YR4.dat'

        x = x + 1
        BRSM(x)%xmin = 120.0D0
        BRSM(x)%xmax = 130.0D0
        BRSM(x)%sep = 0.1D0
        filename(x) = 'br.sm2_YR4.dat'

        ! checks we've filled the whole array
        if (x .ne. xend) then
            stop 'error in setup_BRSM (a)'
        endif

        col = 7
        ! do loop to read in S95 tables
        do x = xbeg, xend
            BRSM(x)%nx = nint((BRSM(x)%xmax - BRSM(x)%xmin)/BRSM(x)%sep) + 1
            allocate (BRSM(x)%dat(BRSM(x)%nx, col - 1))
        enddo

        open (file_id_common2, file=trim(adjustl(pathname))//'Theory_tables/'// &
              'BRSM.binary', form='unformatted')

        read (file_id_common2, iostat=ios) BRSM(xbeg)%dat

        if (ios .eq. 0) then

            do x = xbeg + 1, xend
                read (file_id_common2) BRSM(x)%dat
            enddo

        else
            rewind (file_id_common2)
            do x = xbeg, xend
                fullfilename = trim(adjustl(pathname))//'Theory_tables/YR4/' &
                               //trim(filename(x))

                call read_tabletype1(BRSM(x), 3, col, fullfilename)
#ifndef WEBVERSION
                write (file_id_common2) BRSM(x)%dat
#endif
            enddo
        endif
        close (file_id_common2)

        ! we want the (smallest) range of these functions
        ! initial (extreme) values
        BRSMt1Mhmin = 1.0D6
        BRSMt1Mhmax = 3.0D3

        do x = xbeg, xend
            if (BRSM(x)%xmax .gt. BRSMt1Mhmax) BRSMt1Mhmax = BRSM(x)%xmax
            if (BRSM(x)%xmin .lt. BRSMt1Mhmin) BRSMt1Mhmin = BRSM(x)%xmin
        enddo

        deallocate (filename)

    end subroutine setup_BRSM

    !******************************************************
    subroutine out_of_BRSM_range(funcname, res, strict)
        !******************************************************
        implicit none
        character(LEN=*), intent(in) :: funcname
        double precision, intent(out) :: res
        logical, optional :: strict

        !write(*,*)'The SM Higgs mass given as the'
        !write(*,*)'argument to function '//funcname//' should be'
        !write(*,*)'between',BRSM(1)%xmin,'and',BRSM(1)%xmax
        !write(*,*)'(in units of GeV).'

        ! set res to an initial value
        if (trim(adjustl(funcname)) .eq. 'BRSM_GammaTot') then
            res = 1.0D10
        else
            res = 0.0D0
        endif

        if (present(strict)) then
            if (strict) then
                write (*, *) 'The SM Higgs mass given as the'
                write (*, *) 'argument to function '//funcname//' is'
                write (*, *) 'out of range.'

                res = -1.0D0
            endif
        endif

    end subroutine out_of_BRSM_range
    !******************************************************
    function BRSM_HWW(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_HWW
        double precision :: gamma
        character(LEN=20):: functionname = 'BRSM_HWW'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! WW is in br.sm2_HDecay, column 5 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(2), 5 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(4), 5 - 1, interpol)
                endif
            endif

            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_HWW(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif

        BRSM_HWW = interpol

    end function BRSM_HWW

    !******************************************************
    function BRSM_HZZ(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_HZZ
        double precision :: gamma
        character(LEN=20):: functionname = 'BRSM_HZZ'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! ZZ is in br.sm2_HDecay, column 6 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(2), 6 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(4), 6 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_HZZ(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_HZZ = interpol

    end function BRSM_HZZ
    !******************************************************
    function BRSM_Hbb(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: gamma
        double precision :: BRSM_Hbb
        character(LEN=20):: functionname = 'BRSM_Hbb'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! bb is in br.sm1_HDecay, column 2 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(1), 2 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(3), 2 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_Hbb(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Hbb = interpol

    end function BRSM_Hbb
    !******************************************************
    function BRSM_Htautau(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: gamma
        double precision :: BRSM_Htautau
        character(LEN=20):: functionname = 'BRSM_Htautau'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! tautau is in br.sm1_HDecay, column 3 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(1), 3 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(3), 3 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_Htautau(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Htautau = interpol

    end function BRSM_Htautau
    !******************************************************
    function BRSM_Hgaga(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: gamma
        double precision :: BRSM_Hgaga
        character(LEN=20):: functionname = 'BRSM_Hgaga'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! photon photon is in br.sm2_HDecay, column 3 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(2), 3 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(4), 3 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_Hgaga(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Hgaga = interpol

    end function BRSM_Hgaga
    !******************************************************
    function BRSM_Hgg(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: gamma
        double precision :: BRSM_Hgg
        character(LEN=20):: functionname = 'BRSM_Hgg'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! gluon gluon is in br.sm2_HDecay, column 2 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(2), 2 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(4), 2 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_Hgg(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Hgg = interpol

    end function BRSM_Hgg
    !******************************************************
    function BRSM_Htoptop(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: gamma
        double precision :: BRSM_Htoptop
        character(LEN=20):: functionname = 'BRSM_Htoptop'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! top top is in br.sm1_HDecay, column 7 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(1), 7 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(3), 7 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_Gamma_Htoptop(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma/BRSM_extrapol_GammaTot(x)
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Htoptop = interpol

    end function BRSM_Htoptop

    !******************************************************
    function BRSM_Hcc(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_Hcc
        character(LEN=20):: functionname = 'BRSM_Hcc'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! charm charm is in br.sm1_HDecay, column 6 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(1), 6 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(3), 6 - 1, interpol)
                endif
            endif

            if (interpol .lt. -1.0D-7) then
                ! Test if we are still in "viable" mass region:
                if (BRSM_extrapol_GammaTot(x) .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    ! Set BR to zero (no extrapolation beyond 1 TeV available)
                    interpol = 0.0D0
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Hcc = interpol

    end function BRSM_Hcc
    !******************************************************
    function BRSM_Hss(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_Hss
        character(LEN=20):: functionname = 'BRSM_Hss'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! strange strange is in br.sm1_HDecay, column 5 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(1), 5 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(3), 5 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                ! Test if we are still in "viable" mass region:
                if (BRSM_extrapol_GammaTot(x) .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    ! Set BR to zero (no extrapolation beyond 1 TeV available)
                    interpol = 0.0D0
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Hss = interpol

    end function BRSM_Hss
    !******************************************************
    function BRSM_Hmumu(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_Hmumu
        character(LEN=20):: functionname = 'BRSM_Hmumu'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! muon muon is in br.sm1_HDecay, column 4 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(1), 4 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(3), 4 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                ! Test if we are still in "viable" mass region:
                if (BRSM_extrapol_GammaTot(x) .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    ! Set BR to zero (no extrapolation beyond 1 TeV available)
                    interpol = 0.0D0
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_Hmumu = interpol

    end function BRSM_Hmumu
    !******************************************************
    function BRSM_HZga(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_HZga
        character(LEN=20):: functionname = 'BRSM_HZga'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! Z gamma is in br.sm2_HDecay, column 4 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(2), 4 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(4), 4 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                ! Test if we are still in "viable" mass region:
                if (BRSM_extrapol_GammaTot(x) .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    ! Set BR to zero (no extrapolation beyond 1 TeV available)
                    interpol = 0.0D0
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_HZga = interpol

    end function BRSM_HZga

    !******************************************************
    function BRSM_Gamma_tWpb(mt, strict)
        !******************************************************
        use usefulbits, only: mbmb, MW, GF, pi, alphas
        ! numbers and equation read from http://pdg.lbl.gov/ 22.10.2009
        ! BRSM_Gamma_tWpb is in GeV/c^2
        implicit none
        double precision, intent(in) :: mt
        double precision :: BRSM_Gamma_tWpb
        double precision :: mbpole
        double precision :: MW2bymt2
        character(LEN=20):: functionname = 'BRSM_Gamma_tWpb'
        logical, optional :: strict

        !eq 17, review of quark masses:
        mbpole = mbmb*(1.0D0 + 0.09D0 + 0.05D0 + 0.03D0)

        MW2bymt2 = MW**2.0D0/mt**2.0D0

        if (mt .gt. (MW + mbpole)) then
            !eq 1, review of top quark
            BRSM_Gamma_tWpb = &
                (GF*mt**3.0D0)/(8.0D0*pi*sqrt(2.0D0)) &
                *(1.0D0 - MW2bymt2)**2.0D0 &
                *(1.0D0 + 2.0D0*MW2bymt2) &
                *(1.0D0 - (2.0D0*alphas)/(3.0D0*pi) &
                  *((2.0D0*pi**2.0D0/3.0D0) - (5.0D0/2.0D0)))
        else
            call out_of_BRSM_range(functionname, BRSM_Gamma_tWpb, strict)
        endif

    end function BRSM_Gamma_tWpb
    !******************************************************
    function BRSM_GammaTot(x, strict, EWcorr)
        !******************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: interpol
        double precision :: BRSM_GammaTot
        double precision :: gamma
        character(LEN=20):: functionname = 'BRSM_GammaTot'
        logical, optional :: strict
        logical, optional :: EWcorr

        if (x .ge. 0.0D0) then
            ! GammaTot is in br.sm2_HDecay, column 7 (first column is Higgs mass)
            call interpolate_tabletype1(x, BRSM(2), 7 - 1, interpol)
            if (present(EWcorr)) then
                if (EWcorr .and. (abs(x - 125.0D0) .le. 5.0D0)) then
                    call interpolate_tabletype1(x, BRSM(4), 7 - 1, interpol)
                endif
            endif
            if (interpol .lt. -1.0D-7) then
                gamma = BRSM_extrapol_GammaTot(x)
                if (gamma .lt. -1.0D-7) then
                    call out_of_BRSM_range(functionname, interpol, strict)
                else
                    interpol = gamma
                endif
            endif
        else
            call out_of_BRSM_range(functionname, interpol, strict)
        endif
        BRSM_GammaTot = interpol

    end function BRSM_GammaTot
    !************************************************************
    function BRSM_extrapol_Gamma_Hbb(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_Hbb

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_Hbb = 0.00599731512D0 + 7.7064197D-09*x**2
        else
            BRSM_extrapol_Gamma_Hbb = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_Hbb
    !************************************************************
    function BRSM_extrapol_Gamma_Htoptop(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_Htoptop

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_Htoptop = 0.03269275D0*x + 9.19076988D-06*x**2
        else
            BRSM_extrapol_Gamma_Htoptop = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_Htoptop
    !************************************************************
    function BRSM_extrapol_Gamma_Htautau(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_Htautau

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_Htautau = 0.000763876605D0 + 1.3666917658D-09*x**2
        else
            BRSM_extrapol_Gamma_Htautau = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_Htautau
    !************************************************************
    function BRSM_extrapol_Gamma_Hgg(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_Hgg

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_Hgg = 0.051098502D0 + 1.3378912984D-08*x**2
        else
            BRSM_extrapol_Gamma_Hgg = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_Hgg
    !************************************************************
    function BRSM_extrapol_Gamma_Hgaga(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_Hgaga

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_Hgaga = 0.46982516D-03 - 1.6225307D-06*x + 1.41882959D-09*x**2
        else
            BRSM_extrapol_Gamma_Hgaga = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_Hgaga
    !************************************************************
    function BRSM_extrapol_Gamma_HWW(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_HWW

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_HWW = -0.068683437327D0*x + 0.23825627195D-03*x**2 + 1.433138585D-10*x**4
        else
            BRSM_extrapol_Gamma_HWW = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_HWW
    !************************************************************
    function BRSM_extrapol_Gamma_HZZ(x)
        ! This function can be used to extrapolate the partial width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_Gamma_HZZ

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_Gamma_HZZ = -0.03436449766D0*x + 0.115788832D-03*x**2 + 7.3407115787D-11*x**4
        else
            BRSM_extrapol_Gamma_HZZ = -1.0D0
        endif
    end function BRSM_extrapol_Gamma_HZZ
    !************************************************************
    function BRSM_extrapol_GammaTot(x)
        ! This function can be used to extrapolate the total width
        ! in the mass region [1:3] TeV. It was fitted in the large mass
        ! region to the YR4 predictions.
        ! (TS 27/05/2016)
        !************************************************************
        implicit none
        double precision, intent(in) :: x
        double precision :: BRSM_extrapol_GammaTot

        if (x .ge. 1000D0 .and. x .le. 3000D0) then
            BRSM_extrapol_GammaTot = BRSM_extrapol_Gamma_Hbb(x) + BRSM_extrapol_Gamma_Htoptop(x) + &
                                     BRSM_extrapol_Gamma_Htautau(x) + BRSM_extrapol_Gamma_Hgg(x) + &
                                     BRSM_extrapol_Gamma_Hgaga(x) + BRSM_extrapol_Gamma_HWW(x) + &
                                     BRSM_extrapol_Gamma_HZZ(x)
        else
            BRSM_extrapol_GammaTot = -1.0D0
        endif
    end function BRSM_extrapol_GammaTot
    !************************************************************
    subroutine deallocate_BRSM
        !************************************************************
        implicit none
        !-----------------------------------internal
        integer x
        !-------------------------------------------
        do x = lbound(BRSM, dim=1), ubound(BRSM, dim=1)
            deallocate (BRSM(x)%dat)
        enddo

        deallocate (BRSM)

    end subroutine deallocate_BRSM
    !************************************************************
end module theory_BRfunctions
!*******************************************************************
