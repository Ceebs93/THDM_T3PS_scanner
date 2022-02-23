!> @file
!> Higgsbounds user subroutines.
!!
!! These are the main routines needed to use HiggsBounds.
!! They cover both the input of model predictions in different formats and the
!! actual running of HiggsBounds. See the subroutine [documentation](doc/subroutine.md) for a
!! detailed description of how to use this interface.
!!
!! If you want to use HiggsBounds with C or C++ please take a look at @ref HiggsBounds_subroutines.h

!> Initializes HiggsBounds.
!! It calls subroutines to read in the Standard Model data,
!! data from LEP, Tevatron and LHC experiments, and
!! sets up lists of processes which should be tested.
!! @param nHiggsneut Number of neutral Higgs bosons in the model
!! @param nHiggsplus Number of charged Higgs bosons in the model
!! @param whichanalyses_in selection of experimental data to be used
!! | value | experimental data                                    |
!! |-------|------------------------------------------------------|
!! | onlyL | only LEP data                                        |
!! | onlyH | only Tevatron and LHC data                           |
!! | LandH | LEP, Tevatron, LHC data                              |
!! | onlyP | only published results (i.e. those with an arXiv-id) |
!! | list  | only data defined in usefulbits::analysislist        |
subroutine initialize_HiggsBounds(nHiggsneut, nHiggsplus, whichanalyses_in)
    use usefulbits, only: np, Hneut, Hplus, Chineut, Chiplus, debug, inputmethod, &
                          theo, whichanalyses, HiggsBounds_info, just_after_run, BRdirectinput, &
                          file_id_debug1, file_id_debug2, run_HB_classic, anyH
    use input, only: setup_input, check_number_of_particles, check_whichanalyses
    use S95tables, only: setup_S95tables
    use likelihoods, only: setup_likelihoods
    use theory_BRfunctions, only: setup_BRSM
    use theory_XS_SM_functions, only: setup_XSSM
    use channels, only: setup_channels
    use output, only: setup_output
#ifdef enableCHISQ
    use S95tables_type3, only: clsb_t3, fillt3needs_M2_gt_2M1
    use usefulbits, only: allocate_if_stats_required
    use S95tables, only: S95_t2
#endif

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !--------------------------------------input
    integer, intent(in) :: nHiggsneut
    integer, intent(in) :: nHiggsplus
    character(LEN=5), intent(in) :: whichanalyses_in
    !-----------------------------------internal
    logical :: messages
    !-------------------------------------------

#ifdef DEBUGGING
    debug = .True.
#else
    debug = .False.
#endif

    messages = debug .or. (inputmethod == 'datfile')

    np(Hneut) = nHiggsneut
    np(Hplus) = nHiggsplus
    np(anyH) = np(Hneut) + np(Hplus)

    np(Chineut) = 0! do not change this without contacting us first!
    np(Chiplus) = 0! do not change this without contacting us first!

    whichanalyses = whichanalyses_in

    if (inputmethod == 'subrout') then
        if (allocated(theo)) then
            stop 'subroutine initialize_HiggsBounds has already been called once'
        endif

        if (messages) write (*, *) 'doing other preliminary tasks...'; call flush (6)
        call setup_input
    endif

    if (inputmethod .ne. 'datfile') call HiggsBounds_info
    if (run_HB_classic .EQV. .True.) then
        PRINT *, "run_HB_classic=True - HiggsBounds is running in classic mode"
    endif

    if (messages) write (*, *) 'reading in Standard Model tables...'; call flush (6)
    call setup_BRSM

    call setup_XSSM

    if (messages) write (*, *) 'reading in S95tables...'; call flush (6)
    call setup_S95tables

    if (messages) write (*, *) 'reading in likelihoods...'; call flush (6)
    call setup_likelihoods

    if (messages) then
        open (file_id_debug2, file='debug_predratio.txt')
        open (file_id_debug1, file='debug_channels.txt')
    endif

    if (messages) write (*, *) 'sorting out processes to be checked...'; call flush (6)
    call setup_channels

    if (messages) write (*, *) 'preparing output arrays...'; call flush (6)
    call setup_output

#ifdef enableCHISQ
    if (allocated(allocate_if_stats_required)) then
        call fillt3needs_M2_gt_2M1(clsb_t3, S95_t2)
    endif
#endif

    just_after_run = .False.
    BRdirectinput = .False.
end subroutine initialize_HiggsBounds

!> Calls \ref higgsbounds_subroutines::initialize_higgsbounds.
!! Experimental data set is specified by an integer.
!! @param nHiggsneut Number of neutral Higgs bosons in the model
!! @param nHiggsplus Number of charged Higgs bosons in the model
!! @param flag selection of experimental data to be used
!!
!! | flag | value | experimental data                                    |
!! |------|-------|------------------------------------------------------|
!! |  1   | onlyL | only LEP data                                        |
!! |  2   | onlyH | only Tevatron and LHC data                           |
!! |  3   | LandH | LEP, Tevatron, LHC data                              |
!! |  4   | onlyP | only published results (i.e. those with an arXiv-id) |
!! |  5   | list  | only data defined in usefulbits::analysislist        |
subroutine initialize_HiggsBounds_int(nHiggsneut, nHiggsplus, flag)

    implicit none

    integer, intent(in) :: nHiggsneut, nHiggsplus, flag

    interface
        subroutine initialize_HiggsBounds(nHiggsneut, nHiggsplus, whichanalyses_in)
            integer, intent(in) :: nHiggsneut
            integer, intent(in) :: nHiggsplus
            character(LEN=5), intent(in) :: whichanalyses_in
        end subroutine initialize_HiggsBounds
    end interface

    IF (flag .EQ. 1) then
        call initialize_HiggsBounds(nHiggsneut, nHiggsplus, "onlyL")
    elseif (flag .EQ. 2) then
        call initialize_HiggsBounds(nHiggsneut, nHiggsplus, "onlyH")
    elseif (flag .EQ. 3) then
        call initialize_HiggsBounds(nHiggsneut, nHiggsplus, "LandH")
    elseif (flag .EQ. 4) then
        call initialize_HiggsBounds(nHiggsneut, nHiggsplus, "onlyP")
    elseif (flag .EQ. 5) then
        call initialize_HiggsBounds(nHiggsneut, nHiggsplus, "list ")
    else
        stop "Illegal value for flag in call to initialize_HB"
    endif
end subroutine

!> Inputs an SLHA file to HiggsBounds.
!! @param infile Name of the SLHA input file.
!! Note: the SLHA file has to contain the HiggsBounds input blocks (see manual) and Higgs boson decay tables.
subroutine HiggsBounds_input_SLHA(infile)
    use usefulbits, only: whichinput, infile1, theo, effC, just_after_run
    use extra_bits_for_SLHA
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    character(len=300), intent(in) :: infile
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    whichinput = 'SLHA'
    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif
    infile1 = infile
    call getSLHAdata(theo(n), effC(n), infile1)

    just_after_run = .False.
end subroutine HiggsBounds_input_SLHA

!> Input of neutral Higgs boson's properties (mass, total width, CP) to HiggsBounds.
!!
!! Note: When using the effective coupling input the total width can be internally derived
!! from the provided input (effective couplings and branching ratios of Higgs decays to non-SM final states)
!! under the assumption of no additional BSM contributions to the total width beyond those provided in the input.
!! This derivation is activated when providing a negative value for the total width.
!!
!! @param Mh mass values (in GeV) of the neutral Higgs bosons
!! @param GammaTotal_hj total decay width (in GeV) of the neutral Higgs bosons
!! @param CP_value CP properties of neutral Higgs bosons (+1: CP-even, -1: CP-odd, 0: CP-mixed)

subroutine HiggsBounds_neutral_input_properties(Mh, GammaTotal_hj, CP_value)
    use usefulbits, only: theo, np, Hneut, anyH, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: Mh(np(Hneut)), GammaTotal_hj(np(Hneut))
    integer, intent(in) :: CP_value(np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_mass_width should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_mass_width'
    endif

    theo(n)%particle(Hneut)%M = Mh
    theo(n)%particle(Hneut)%Mc = Mh
    theo(n)%particle(Hneut)%GammaTot = GammaTotal_hj
    theo(n)%CP_value = CP_value

    theo(n)%particle(anyH)%M(:np(Hneut)) = Mh
    theo(n)%particle(anyH)%Mc(:np(Hneut)) = Mh
    theo(n)%particle(anyH)%GammaTot(:np(Hneut)) = GammaTotal_hj

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_properties

!> Input of neutral Higgs boson's effective couplings (aka scale factors) to HiggsBounds.
!! These are used to obtain production cross sections and branching ratios by
!! approximately rescaling the corresponding predictions for a SM Higgs boson.
!! @param ghjss_s Scalar (SM-normalized) effective Higgs couplings to strange quarks
!! @param ghjss_p Pseudoscalar (SM-normalized) effective Higgs coupling to strange quarks
!! @param ghjcc_s Scalar (SM-normalized) effective Higgs coupling to charm quarks
!! @param ghjcc_p Pseudoscalar (SM-normalized) effective Higgs coupling to charm quarks
!! @param ghjbb_s Scalar (SM-normalized) effective Higgs coupling to bottom quarks
!! @param ghjbb_p Pseudoscalar (SM-normalized) effective Higgs coupling to bottom quarks
!! @param ghjtt_s Scalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghjtt_p Pseudoscalar (SM-normalized) effective Higgs coupling to top quarks
!! @param ghjmumu_s Scalar (SM-normalized) effective Higgs coupling to muons
!! @param ghjmumu_p Pseudoscalar (SM-normalized) effective Higgs coupling to muons
!! @param ghjtautau_s Scalar (SM-normalized) effective Higgs coupling to tau leptons
!! @param ghjtautau_p Pseudoscalar (SM-normalized) effective Higgs coupling to tau leptons
!! @param ghjWW (SM-normalized) effective Higgs coupling to W bosons
!! @param ghjZZ (SM-normalized) effective Higgs coupling to Z bosons
!! @param ghjZga (SM-normalized) effective Higgs coupling to a Z boson and a photon
!! @param ghjgaga (SM-normalized) effective Higgs coupling to photons
!! @param ghjgg (SM-normalized) effective Higgs coupling to gluons
!! @param ghjhiZ (SM-normalized) effective Higgs-Higgs-Z boson coupling
subroutine HiggsBounds_neutral_input_effC( &
    ghjss_s, ghjss_p, ghjcc_s, ghjcc_p, &
    ghjbb_s, ghjbb_p, ghjtt_s, ghjtt_p, &
    ghjmumu_s, ghjmumu_p, &
    ghjtautau_s, ghjtautau_p, &
    ghjWW, ghjZZ, ghjZga, &
    ghjgaga, ghjgg, ghjhiZ)

    use usefulbits, only: theo, np, Hneut, effC, whichinput, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: &
        ghjss_s(np(Hneut)), ghjss_p(np(Hneut)), &
        ghjcc_s(np(Hneut)), ghjcc_p(np(Hneut)), &
        ghjbb_s(np(Hneut)), ghjbb_p(np(Hneut)), &
        ghjtt_s(np(Hneut)), ghjtt_p(np(Hneut)), &
        ghjmumu_s(np(Hneut)), ghjmumu_p(np(Hneut)), &
        ghjtautau_s(np(Hneut)), ghjtautau_p(np(Hneut)), &
        ghjWW(np(Hneut)), ghjZZ(np(Hneut)), ghjZga(np(Hneut)), &
        ghjgaga(np(Hneut)), ghjgg(np(Hneut)), &
        ghjhiZ(np(Hneut), np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    whichinput = 'effC'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_effC should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_effC'
    endif

    effC(n)%hjss_s = ghjss_s
    effC(n)%hjss_p = ghjss_p
    effC(n)%hjcc_s = ghjcc_s
    effC(n)%hjcc_p = ghjcc_p
    effC(n)%hjbb_s = ghjbb_s
    effC(n)%hjbb_p = ghjbb_p
    effC(n)%hjtt_s = ghjtt_s
    effC(n)%hjtt_p = ghjtt_p
    effC(n)%hjmumu_s = ghjmumu_s
    effC(n)%hjmumu_p = ghjmumu_p
    effC(n)%hjtautau_s = ghjtautau_s
    effC(n)%hjtautau_p = ghjtautau_p

    effC(n)%hjWW = ghjWW
    effC(n)%hjZZ = ghjZZ
    effC(n)%hjZga = ghjZga
    effC(n)%hjgaga = ghjgaga
    effC(n)%hjgg = ghjgg

    effC(n)%hjhiZ = ghjhiZ

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_effC

!> Input of neutral Higgs boson's branching ratios for Standard Model final states to HiggsBounds.
!! @param BR_hjss Branching ratio of Higgs decay to strange quarks, \f$ h_j \to s \bar{s} \f$
!! @param BR_hjcc Branching ratio of Higgs decay to charm quarks, \f$ h_j \to c \bar{c} \f$
!! @param BR_hjbb Branching ratio of Higgs decay to bottom quarks, \f$ h_j \to b \bar{b} \f$
!! @param BR_hjtt Branching ratio of Higgs decay to top quarks, \f$ h_j \to t \bar{t} \f$
!! @param BR_hjmumu Branching ratio of Higgs decay to muons, \f$ h_j \to \mu^+ \mu^- \f$
!! @param BR_hjtautau Branching ratio of Higgs decay to tau leptons, \f$ h_j \to \tau^+ \tau^- \f$
!! @param BR_hjWW Branching ratio of Higgs decay to W bosons, \f$ h_j \to W^+W^- \f$
!! @param BR_hjZZ Branching ratio of Higgs decay to Z bosons, \f$ h_j \to ZZ \f$
!! @param BR_hjZga Branching ratio of Higgs decay to a Z boson and a photon, \f$ h_j \to Z\gamma \f$
!! @param BR_hjgaga Branching ratio of Higgs decay to photons, \f$ h_j \to \gamma\gamma \f$
!! @param BR_hjgg Branching ratio of Higgs decay to gluons, \f$ h_j \to gg \f$
!!
subroutine HiggsBounds_neutral_input_SMBR(BR_hjss, BR_hjcc, BR_hjbb, &
                                          BR_hjtt, BR_hjmumu, &
                                          BR_hjtautau, BR_hjWW, &
                                          BR_hjZZ, BR_hjZga, BR_hjgaga, &
                                          BR_hjgg)

    use usefulbits, only: theo, np, Hneut, just_after_run, BRdirectinput

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: BR_hjss(np(Hneut)), BR_hjcc(np(Hneut)), &
        BR_hjbb(np(Hneut)), BR_hjtt(np(Hneut)), &
        BR_hjmumu(np(Hneut)), BR_hjtautau(np(Hneut)), &
        BR_hjWW(np(Hneut)), BR_hjZZ(np(Hneut)), &
        BR_hjZga(np(Hneut)), BR_hjgaga(np(Hneut)), &
        BR_hjgg(np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_SMBR should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_SMBR'
    endif

    theo(n)%BR_hjss = BR_hjss
    theo(n)%BR_hjcc = BR_hjcc
    theo(n)%BR_hjbb = BR_hjbb
    theo(n)%BR_hjtt = BR_hjtt
    theo(n)%BR_hjmumu = BR_hjmumu
    theo(n)%BR_hjtautau = BR_hjtautau
    theo(n)%BR_hjWW = BR_hjWW
    theo(n)%BR_hjZZ = BR_hjZZ
    theo(n)%BR_hjZga = BR_hjZga
    theo(n)%BR_hjgaga = BR_hjgaga
    theo(n)%BR_hjgg = BR_hjgg

    just_after_run = .False.
    BRdirectinput = .True.
end subroutine HiggsBounds_neutral_input_SMBR

!> Input of neutral Higgs boson's branching ratios for non-Standard Model final states to HiggsBounds.
!! @param BR_hjinvisible Branching ratio of Higgs decay into invisible final states, \f$ h_j \to \mathrm{invisible} \f$
!! @param BR_hkhjhi Branching ratio of Higgs decay into two other Higgs bosons, \f$h_k\to h_i h_j\f$
!! @param BR_hjhiZ Branching ratio of Higgs decay into a Higgs boson and a Z boson, \f$h_j\to h_i Z\f$
!! @param BR_hjemu Branching ratio of lepton-flavor-violating Higgs decay into electron and muon, \f$h_j\to e^\pm \mu^\mp\f$
!! @param BR_hjetau Branching ratio of lepton-flavor-violating Higgs decay into electron and tau lepton, \f$h_j\to e^\pm \tau^\mp\f$
!! @param BR_hjmutau Branching ratio of lepton-flavor-violating Higgs decay into muon and tau lepton, \f$h_j\to \mu^\pm \tau^\mp\f$
!! @param BR_hjHpiW Branching ratio of Higgs decay into a charged Higgs boson and a W boson, \f$h_j\to h^\pm_i W^\mp\f$
subroutine HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi, BR_hjhiZ, &
                                             BR_hjemu, BR_hjetau, BR_hjmutau, BR_hjHpiW)
    use usefulbits, only: theo, np, Hneut, Hplus, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: BR_hjinvisible(np(Hneut)), &
        BR_hkhjhi(np(Hneut), np(Hneut), np(Hneut)), &
        BR_hjhiZ(np(Hneut), np(Hneut)), &
        BR_hjemu(np(Hneut)), &
        BR_hjetau(np(Hneut)), &
        BR_hjmutau(np(Hneut))
    double precision, intent(in) :: BR_hjHpiW(np(Hneut), np(Hplus))
    !--------------------------------------internal
    integer, parameter :: n = 1

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_nonSMBR should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_nonSMBR'
    endif

    theo(n)%BR_hjinvisible = BR_hjinvisible
    theo(n)%full_BR_inv = .false.
    theo(n)%BR_hkhjhi = BR_hkhjhi
    theo(n)%BR_hjhiZ = BR_hjhiZ
    theo(n)%BR_hjemu = BR_hjemu
    theo(n)%BR_hjetau = BR_hjetau
    theo(n)%BR_hjmutau = BR_hjmutau
    theo(n)%BR_hjHpiW = BR_hjHpiW

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_nonSMBR

!> Input of neutral Higgs boson's effective couplings (aka scale factors) to HiggsBounds.
!! These are used to obtain production cross sections.
!! @param ghjuu_s Scalar (SM-normalized) effective Higgs couplings to up quarks
!! @param ghjuu_p Pseudoscalar (SM-normalized) effective Higgs coupling to up quarks
!! @param ghjdd_s Scalar (SM-normalized) effective Higgs coupling to down quarks
!! @param ghjdd_p Pseudoscalar (SM-normalized) effective Higgs coupling to down quarks
!! @param ghjee_s Scalar (SM-normalized) effective Higgs coupling to electrons
!! @param ghjee_p Pseudoscalar (SM-normalized) effective Higgs coupling to electrons
subroutine HiggsBounds_neutral_input_effC_firstgen( &
    ghjuu_s, ghjuu_p, ghjdd_s, ghjdd_p, ghjee_s, ghjee_p)

    use usefulbits, only: theo, np, Hneut, effC, whichinput, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: &
        ghjuu_s(np(Hneut)), ghjuu_p(np(Hneut)), &
        ghjdd_s(np(Hneut)), ghjdd_p(np(Hneut)), &
        ghjee_s(np(Hneut)), ghjee_p(np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    whichinput = 'effC'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_effC_firstgen should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_effC_firstgen'
    endif

    effC(n)%hjuu_s = ghjuu_s
    effC(n)%hjuu_p = ghjuu_p
    effC(n)%hjdd_s = ghjdd_s
    effC(n)%hjdd_p = ghjdd_p
    effC(n)%hjee_s = ghjee_s
    effC(n)%hjee_p = ghjee_p

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_effC_firstgen

!> Input of neutral Higgs boson's flavor violating couplings to HiggsBounds.
!! These dimensionless couplings are absent in the SM.
!! @param ghjuc_s Scalar flavor violating Higgs coupling to up and charm quark
!! @param ghjuc_p Pseudoscalar flavor violating Higgs coupling to up and charm quark
!! @param ghjut_s Scalar flavor violating Higgs coupling to up and top quark
!! @param ghjut_p Pseudoscalar flavor violating Higgs coupling to up and top quark
!! @param ghjct_s Scalar flavor violating Higgs coupling to charm and top quark
!! @param ghjct_p Pseudoscalar flavor violating Higgs coupling to charm and top quark
!! @param ghjds_s Scalar flavor violating Higgs coupling to down and strange quark
!! @param ghjds_p Pseudoscalar flavor violating Higgs coupling to down and strange quark
!! @param ghjdb_s Scalar flavor violating Higgs coupling to down and bottom quark
!! @param ghjdb_p Pseudoscalar flavor violating Higgs coupling to down and bottom quark
!! @param ghjsb_s Scalar flavor violating Higgs coupling to strange and bottom quark
!! @param ghjsb_p Pseudoscalar flavor violating Higgs coupling to strange and bottom quark

subroutine HiggsBounds_neutral_input_effC_FV( &
    ghjuc_s, ghjuc_p, ghjut_s, ghjut_p, ghjct_s, ghjct_p, &
    ghjds_s, ghjds_p, ghjdb_s, ghjdb_p, ghjsb_s, ghjsb_p)

    use usefulbits, only: theo, np, Hneut, effC, whichinput, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: &
        ghjuc_s(np(Hneut)), ghjuc_p(np(Hneut)), &
        ghjut_s(np(Hneut)), ghjut_p(np(Hneut)), &
        ghjct_s(np(Hneut)), ghjct_p(np(Hneut)), &
        ghjds_s(np(Hneut)), ghjds_p(np(Hneut)), &
        ghjdb_s(np(Hneut)), ghjdb_p(np(Hneut)), &
        ghjsb_s(np(Hneut)), ghjsb_p(np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    whichinput = 'effC'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_effC_FV should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_effC_FV'
    endif

    effC(n)%hjuc_s = ghjuc_s
    effC(n)%hjuc_p = ghjuc_p
    effC(n)%hjut_s = ghjut_s
    effC(n)%hjut_p = ghjut_p
    effC(n)%hjct_s = ghjct_s
    effC(n)%hjct_p = ghjct_p
    effC(n)%hjds_s = ghjds_s
    effC(n)%hjds_p = ghjds_p
    effC(n)%hjdb_s = ghjdb_s
    effC(n)%hjdb_p = ghjdb_p
    effC(n)%hjsb_s = ghjsb_s
    effC(n)%hjsb_p = ghjsb_p

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_effC_FV

!> Input of neutral Higgs boson's branching ratios to first generation SM fermions.
!! @param BR_hjuu Branching ratio of Higgs decay into up-quarks, \f$ h_j \to u \bar{u} \f$
!! @param BR_hjuu Branching ratio of Higgs decay into down-quarks, \f$ h_j \to d \bar{d} \f$
!! @param BR_hjuu Branching ratio of Higgs decay into electrons, \f$ h_j \to e^+e^- \f$
subroutine HiggsBounds_neutral_input_firstgenBR(BR_hjuu, BR_hjdd, BR_hjee)
    use usefulbits, only: theo, np, Hneut, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: BR_hjuu(np(Hneut)), BR_hjdd(np(Hneut)), BR_hjee(np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_nonSMBR should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_nonSMBR'
    endif

    theo(n)%BR_hjuu = BR_hjuu
    theo(n)%BR_hjdd = BR_hjdd
    theo(n)%BR_hjee = BR_hjee

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_firstgenBR

!> Input of neutral Higgs boson's branching ratios to first generation SM fermions.
!! @param BR_hjuc Branching ratio of Higgs decay into up-quark and charm quark, \f$ h_j \to u \bar{d} + c.c. \f$
!! @param BR_hjds Branching ratio of Higgs decay into down-quark and strange quark, \f$ h_j \to d \bar{s} + c.c. \f$
!! @param BR_hjut Branching ratio of Higgs decay into up-quark and top quark, \f$ h_j \to u \bar{t} + c.c. \f$
!! @param BR_hjdb Branching ratio of Higgs decay into down-quark and bottom quark, \f$ h_j \to d \bar{b} + c.c. \f$
!! @param BR_hjct Branching ratio of Higgs decay into charm-quark and top quark, \f$ h_j \to c \bar{t} + c.c. \f$
!! @param BR_hjsb Branching ratio of Higgs decay into strange-quark and bottom quark, \f$ h_j \to s \bar{b} + c.c. \f$
subroutine HiggsBounds_neutral_input_FVBR(BR_hjuc, BR_hjds, BR_hjut, BR_hjdb, BR_hjct, BR_hjsb)
    use usefulbits, only: theo, np, Hneut, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: BR_hjuc(np(Hneut)), BR_hjds(np(Hneut)), &
        BR_hjut(np(Hneut)), BR_hjdb(np(Hneut)), &
        BR_hjct(np(Hneut)), BR_hjsb(np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_FVBR should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_FVBR'
    endif

    theo(n)%BR_hjuc = BR_hjuc
    theo(n)%BR_hjds = BR_hjds
    theo(n)%BR_hjut = BR_hjut
    theo(n)%BR_hjdb = BR_hjdb
    theo(n)%BR_hjct = BR_hjct
    theo(n)%BR_hjsb = BR_hjsb

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_FVBR

!> Input of neutral Higgs boson's cross sections for lepton collider (LEP). (See manual for details on normalization.)
!! @param XS_ee_hjZ_ratio (SM-normalized) cross section for Higgs-Strahlung, \f$ e^+e^- \to h_j Z\f$
!! @param XS_ee_bbhj_ratio (SM-normalized) cross section for \f$ b\bar{b}\f$-associated production, \f$ e^+e^- \to h_j b\bar{b}\f$
!! @param XS_ee_tautauhj_ratio (SM-normalized) cross section for \f$ \tau^+\tau^-\f$-associated production, \f$ e^+e^- \to h_j \tau^+\tau^-\f$
!! @param XS_ee_hjhi_ratio (2HDM-normalized) cross section for double Higgs production, \f$ e^+e^- \to h_i h_j\f$
subroutine HiggsBounds_neutral_input_LEP(XS_ee_hjZ_ratio, XS_ee_bbhj_ratio, &
                                         XS_ee_tautauhj_ratio, XS_ee_hjhi_ratio)

    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run
    implicit none
    !---------------------------------------------
    double precision, intent(in) :: XS_ee_hjZ_ratio(np(Hneut)), &
        XS_ee_bbhj_ratio(np(Hneut)), XS_ee_tautauhj_ratio(np(Hneut)), &
        XS_ee_hjhi_ratio(np(Hneut), np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_LEP should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_LEP'
    endif

    theo(n)%lep%XS_hjZ_ratio = XS_ee_hjZ_ratio
    theo(n)%lep%XS_bbhj_ratio = XS_ee_bbhj_ratio
    theo(n)%lep%XS_tautauhj_ratio = XS_ee_tautauhj_ratio
    theo(n)%lep%XS_hjhi_ratio = XS_ee_hjhi_ratio

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_LEP

!> Input of neutral Higgs boson's cross sections for hadron colliders.
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param CS_hj_ratio (SM-normalized) cross section for single Higgs production
!! @param CS_gg_hj_ratio (SM-normalized) cross section for gluon fusion production, \f$gg\to h_j\f$
!! @param CS_bb_hj_ratio (SM-normalized) cross section for Higgs production in association with bottom quarks, \f$b\bar{b}\to h_j\f$
!! @param CS_hjW_ratio (SM-normalized) cross section for Higgs production in association with a W boson
!! @param CS_hjZ_ratio (SM-normalized) cross section for Higgs production in association with a Z boson
!! @param CS_vbf_ratio (SM-normalized) cross section for Higgs production in vector boson fusion
!! @param CS_tthj_ratio (SM-normalized) cross section for Higgs production in association with a top quark pair
!! @param CS_thj_tchan_ratio (SM-normalized) cross section for Higgs production in association with single top quark (t-channel process)
!! @param CS_thj_schan_ratio (SM-normalized) cross section for Higgs production in association with single top quark (s-channel process)
!! @param CS_qq_hjZ_ratio (SM-normalized) cross section for qq-bar-initiated Higgs production in association with a Z boson
!! @param CS_gg_hjZ_ratio (SM-normalized) cross section for gg-initiated Higgs production in association with a Z boson
!! @param CS_tWhj_ratio (SM-normalized) cross sections fo Higgs prodoction in association with a top quark and a W
!! @param CS_hjhi Cross section (in pb) for non-resonant double Higgs production
subroutine HiggsBounds_neutral_input_hadr(collider, CS_hj_ratio, &
                                          CS_gg_hj_ratio, CS_bb_hj_ratio, &
                                          CS_hjW_ratio, CS_hjZ_ratio, &
                                          CS_vbf_ratio, CS_tthj_ratio, &
                                          CS_thj_tchan_ratio, CS_thj_schan_ratio, &
                                          CS_qq_hjZ_ratio, CS_gg_hjZ_ratio, &
                                          CS_tWhj_ratio, &
                                          CS_hjhi)

    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run, hadroncolliderdataset

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    double precision, intent(in) :: CS_hj_ratio(np(Hneut)), &
        CS_gg_hj_ratio(np(Hneut)), CS_bb_hj_ratio(np(Hneut)), &
        CS_hjW_ratio(np(Hneut)), CS_hjZ_ratio(np(Hneut)), &
        CS_vbf_ratio(np(Hneut)), CS_tthj_ratio(np(Hneut)), &
        CS_thj_tchan_ratio(np(Hneut)), CS_thj_schan_ratio(np(Hneut)), &
        CS_qq_hjZ_ratio(np(Hneut)), CS_gg_hjZ_ratio(np(Hneut)), &
        CS_tWhj_ratio(np(Hneut)), &
        CS_hjhi(np(Hneut), np(Hneut))
    integer, intent(in) :: collider
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_hadr should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_hadr'
    endif

    select case (collider)
    case (2)
        call set_input(theo(n)%tev)
    case (7)
        call set_input(theo(n)%lhc7)
    case (8)
        call set_input(theo(n)%lhc8)
    case (13)
        call set_input(theo(n)%lhc13)
    case default
        stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr'
    end select

    just_after_run = .False.

contains

    subroutine set_input(dataset)

        implicit none
        type(hadroncolliderdataset) :: dataset

        dataset%XS_hj_ratio = CS_hj_ratio
        dataset%XS_gg_hj_ratio = CS_gg_hj_ratio
        dataset%XS_bb_hj_ratio = CS_bb_hj_ratio
        dataset%XS_hjW_ratio = CS_hjW_ratio
        dataset%XS_hjZ_ratio = CS_hjZ_ratio
        dataset%XS_gg_hjZ_ratio = CS_qq_hjZ_ratio
        dataset%XS_qq_hjZ_ratio = CS_gg_hjZ_ratio
        dataset%XS_vbf_ratio = CS_vbf_ratio
        dataset%XS_tthj_ratio = CS_tthj_ratio
        dataset%XS_thj_tchan_ratio = CS_thj_tchan_ratio
        dataset%XS_thj_schan_ratio = CS_thj_schan_ratio
        dataset%XS_tWhj_ratio = CS_tWhj_ratio
        dataset%XS_hjhi = CS_hjhi

    end subroutine set_input
end subroutine HiggsBounds_neutral_input_hadr

! TODO: subroutine HiggsBounds_neutral_input_ZHprod(collider,CS_qq_hjZ_ratio,CS_gg_hjZ_ratio)

!> Input of neutral Higgs boson's channel rates.
!! Input of neutral Higgs boson's signal rates for all hadron collider channels (i.e. production x decay mode).
!! This routine can be used to circumvent the narrow width approximation (that is otherwise being employed)
!! for specific collider channels. All non-zero elements of the input matrix channelrates will be used
!! to overwrite the channel rates that are otherwise obtained from the conventional HiggsBounds input
!! in the narrow-width approximation.
!!
!! The input is given in terms of the signal rate (for production times decay), normalized to the
!! corresponding production cross section (without decay rate!) in the SM at the same Higgs mass,
!! for instance, the (normalized) channel rate \f$ \frac{\sigma(gg\to h_3 \to \tau^+\tau^- + X)}{\sigma_\text{SM}(gg \to H)}\f$
!! for Higgs boson \f$h_3\f$ is given by the matrix element \f$\text{channelrates}(3,6,4)\f$ (see below).
!!
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param channelrates Input matrix for each neutral Higgs boson.
!! First index enumerates \f$i\f$ the Higgs bosons, second index \f$i_\text{prod}\f$ enumerates the production processes,
!! third index \f$i_\text{decay}\f$ enumerates the decay processes.
!!  | Index |   Production mode                  | | Index |   Decay mode                 |
!!  |-------|------------------------------------|-|-------|------------------------------|
!!  |   1   |    single Higgs production         | |   1   | \f$h_i \to \gamma\gamma \f$  |
!!  |   2   |    vector boson fusion             | |   2   | \f$h_i \to W W \f$           |
!!  |   3   | \f$h_i W\f$ production             | |   3   | \f$h_i \to Z Z \f$           |
!!  |   4   | \f$h_i Z\f$ production             | |   4   | \f$h_i \to \tau^+\tau^- \f$  |
!!  |   5   | \f$h_i t\bar{t}\f$ production      | |   5   | \f$h_i \to b\bar{b} \f$      |
!!  |   6   | \f$gg \to h_i\f$ production        | |   6   | \f$h_i \to Z\gamma \f$       |
!!  |   7   | \f$b\bar{b} \to h_i\f$ production  | |   7   | \f$h_i \to c\bar{c} \f$      |
!!  |   8   | \f$h_i t\f$ production (t-channel) | |   8   | \f$h_i \to \mu^+\mu^- \f$    |
!!  |   9   | \f$h_i t\f$ production (s-channel) | |   9   | \f$h_i \to gg \f$            |
!!  |   10  | \f$q\bar{q} \to h_i Z\f$ production| |   10  | \f$h_i \to s\bar{s} \f$      |
!!  |   11  | \f$gg \to h_i Z\f$ production      | |   11  | \f$h_i \to t\bar{t} \f$      |
subroutine HiggsBounds_neutral_input_hadr_channelrates(collider, channelrates)
    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run, hadroncolliderdataset, &
                          Nprod, Ndecay

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    double precision, intent(in) :: channelrates(np(Hneut), Nprod, Ndecay)
    integer, intent(in) :: collider
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_hadr_channelrates should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_hadr'
    endif

    select case (collider)
    case (2)
        theo(n)%tev%channelrates_tmp = channelrates
    case (7)
        theo(n)%lhc7%channelrates_tmp = channelrates
    case (8)
        theo(n)%lhc8%channelrates_tmp = channelrates
    case (13)
        theo(n)%lhc13%channelrates_tmp = channelrates
    case default
        stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_channelrates'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_hadr_channelrates

!> Input of charged Higgs boson's properties (mass, width, LEP cross section, BRs) to HiggsBounds
!! Input of charged Higgs boson's masses, total widths, production cross section at lepton colliders, and
!! branching ratios.
!!
!! @param Mhplus mass values (in GeV) of the charged Higgs bosons
!! @param GammaTotal_Hpj total widths (in GeV) of the charged Higgs bosons
!! @param CS_ee_HpjHmj_ratio (2HDM-normalized) cross section for the LEP process \f$ e^+e^- \to h^+_j h^-_j\f$ (see manual for details on normalization)
!! @param BR_tWpb Branching ratio for the top quark decay, \f$ t\to W^+ b \f$
!! @param BR_tHpjb Branching ratio for the top quark decay to charged Higgs bosons, \f$t\to h^+_j b\f$
!! @param BR_Hpjcs Branching ratio for the charged Higgs decay \f$h^+_j \to c \bar{s}\f$
!! @param BR_Hpjcb Branching ratio for the charged Higgs decay \f$h^+_j \to c \bar{b}\f$
!! @param BR_Hpjtaunu Branching ratio for the charged Higgs decay \f$h^+_j \to \tau \nu \f$
!! @param BR_Hpjtb Branching ratio for the charged Higgs decay \f$h^+_j \to t\bar{b} \f$
!! @param BR_HpjWZ Branching ratio for the charged Higgs decay \f$h^+_j \to W^+ Z \f$
!! @param BR_HpjhiW Branching ratio for the charged Higgs decay to a neutral Higgs boson, \f$h^+_j \to h_i W^+ \f$
subroutine HiggsBounds_charged_input(Mhplus, GammaTotal_Hpj, &
                                     CS_ee_HpjHmj_ratio, &
                                     BR_tWpb, BR_tHpjb, &
                                     BR_Hpjcs, BR_Hpjcb, BR_Hpjtaunu, BR_Hpjtb, &
                                     BR_HpjWZ, BR_HpjhiW)

    use usefulbits, only: theo, np, Hplus, Hneut, anyH, just_after_run

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: Mhplus(np(Hplus)), GammaTotal_Hpj(np(Hplus)), &
        CS_ee_HpjHmj_ratio(np(Hplus)), &
        BR_tWpb, BR_tHpjb(np(Hplus)), &
        BR_Hpjcs(np(Hplus)), BR_Hpjcb(np(Hplus)), BR_Hpjtaunu(np(Hplus)), &
        BR_Hpjtb(np(Hplus)), BR_HpjWZ(np(Hplus))
    double precision, intent(in) :: BR_HpjhiW(np(Hplus), np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hplus) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_charged_input should'
        write (*, *) 'only be called if np(Hplus)>0'
        stop 'error in subroutine HiggsBounds_charged_input'
    endif

    theo(n)%particle(Hplus)%M = Mhplus
    theo(n)%particle(Hplus)%Mc = Mhplus
    theo(n)%particle(Hplus)%GammaTot = GammaTotal_Hpj

    theo(n)%particle(anyH)%M(np(Hneut) + 1:) = Mhplus
    theo(n)%particle(anyH)%Mc(np(Hneut) + 1:) = Mhplus
    theo(n)%particle(anyH)%GammaTot(np(Hneut) + 1:) = GammaTotal_Hpj

    theo(n)%lep%XS_HpjHmj_ratio = CS_ee_HpjHmj_ratio

    theo(n)%BR_tWpb = BR_tWpb
    theo(n)%BR_tHpjb = BR_tHpjb

    theo(n)%BR_Hpjcs = BR_Hpjcs
    theo(n)%BR_Hpjcb = BR_Hpjcb
    theo(n)%BR_Hpjtaunu = BR_Hpjtaunu
    theo(n)%BR_Hpjtb = BR_Hpjtb
    theo(n)%BR_HpjWZ = BR_HpjWZ
    theo(n)%BR_HpjhiW = BR_HpjhiW

    just_after_run = .False.
end subroutine HiggsBounds_charged_input

!> Input of charged scalar boson's exotic branching ratios
!! Branching ratios of scalar decays to final states involving first-generation quarks or light leptons (e, mu)
!!
!! @param BR_Hpjud Branching ratio for the charged Higgs decay \f$h^+_j \to u \bar{d}\f$
!! @param BR_Hpjus Branching ratio for the charged Higgs decay \f$h^+_j \to u \bar{s}\f$
!! @param BR_Hpjcd Branching ratio for the charged Higgs decay \f$h^+_j \to c \bar{d}\f$
!! @param BR_Hpjub Branching ratio for the charged Higgs decay \f$h^+_j \to u \bar{b}\f$

subroutine HiggsBounds_charged_input_exoticBR(BR_Hpjud, BR_Hpjus, BR_Hpjcd, BR_Hpjub, BR_Hpjenu, BR_Hpjmunu)

    use usefulbits, only: theo, np, Hplus, just_after_run

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: BR_Hpjud(np(Hplus)), BR_Hpjus(np(Hplus)), &
        BR_Hpjcd(np(Hplus)), BR_Hpjub(np(Hplus)), &
        BR_Hpjenu(np(Hplus)), BR_Hpjmunu(np(Hplus))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hplus) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_charged_input should'
        write (*, *) 'only be called if np(Hplus)>0'
        stop 'error in subroutine HiggsBounds_charged_input'
    endif

    theo(n)%BR_Hpjud = BR_Hpjud
    theo(n)%BR_Hpjus = BR_Hpjus
    theo(n)%BR_Hpjcd = BR_Hpjcd
    theo(n)%BR_Hpjub = BR_Hpjub
    theo(n)%BR_Hpjenu = BR_Hpjenu
    theo(n)%BR_Hpjmunu = BR_Hpjmunu

    just_after_run = .False.
end subroutine HiggsBounds_charged_input_exoticBR

!> Input of charged Higgs boson's production cross sections at hadron colliders.
!! Input of charged Higgs boson's production cross sections at hadron colliders
!! to HiggsBounds. The cross section is given in pb and corresponds to the sum of
!! the quoted process and its charge-conjugate.
!!
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param CS_Hpmjtb Cross section for \f$h_j^+ \bar{t} b\f$ (+ c.c.) production
!! @param CS_Hpmjcb Cross section for \f$h_j^+ \bar{c} b\f$ (+ c.c.) production
!! @param CS_Hpmjbjet Cross section for \f$h_j^+ b + j\f$ (+ c.c.) production (with light-flavor jet \f$j\f$)
!! @param CS_Hpmjcjet Cross section for \f$h_j^+ c + j\f$ (+ c.c.) production (with light-flavor jet \f$j\f$)
!! @param CS_qq_Hpm Cross section for \f$s\f$-channel single \f$h_j^\pm\f$ production
!! @param CS_HpmjW Cross section for \f$h_j^+ W^-\f$ (+ c.c.) production
!! @param CS_HpmjZ Cross section for \f$h_j^+ Z\f$ (+ c.c.) production
!! @param CS_vbf_Hpmj Cross section for \f$h_j^+\f$ (+ c.c.) production in vector boson fusion (VBF)
!! @param CS_HpjHmj Cross section for \f$h_j^+ h_j^-\f$ production
!! @param CS_Hpmjhi Cross section for \f$h_j^+ h_i\f$ (+c.c.) production
subroutine HiggsBounds_charged_input_hadr(collider, CS_Hpmjtb, CS_Hpmjcb, &
                                          CS_Hpmjbjet, CS_Hpmjcjet, CS_qq_Hpmj, CS_HpmjW, &
                                          CS_HpmjZ, CS_vbf_Hpmj, CS_HpjHmj, CS_Hpmjhi)

    use usefulbits, only: theo, np, Hplus, Hneut, just_after_run, hadroncolliderdataset
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: CS_Hpmjtb(np(Hplus)), CS_Hpmjcb(np(Hplus)), &
        CS_Hpmjbjet(np(Hplus)), CS_Hpmjcjet(np(Hplus)), &
        CS_qq_Hpmj(np(Hplus)), &
        CS_HpmjW(np(Hplus)), CS_HpmjZ(np(Hplus)), &
        CS_vbf_Hpmj(np(Hplus)), CS_HpjHmj(np(Hplus))
    integer, intent(in) :: collider
    double precision, intent(in) :: CS_Hpmjhi(np(Hplus), np(Hneut))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hplus) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_charged_input should'
        write (*, *) 'only be called if np(Hplus)>0'
        stop 'error in subroutine HiggsBounds_charged_input'
    endif

    select case (collider)
    case (2)
        call set_input(theo(n)%tev)
    case (7)
        call set_input(theo(n)%lhc7)
    case (8)
        call set_input(theo(n)%lhc8)
    case (13)
        call set_input(theo(n)%lhc13)
    case default
        stop 'wrong input for collider to subroutine HiggsBounds_charged_input_hadr'
    end select

    just_after_run = .False.
contains

    subroutine set_input(dataset)
        type(hadroncolliderdataset) :: dataset
        dataset%XS_Hpmjtb = CS_Hpmjtb
        dataset%XS_Hpmjcb = CS_Hpmjcb
        dataset%XS_Hpmjbjet = CS_Hpmjbjet
        dataset%XS_Hpmjcjet = CS_Hpmjcjet
        dataset%XS_qq_Hpmj = CS_qq_Hpmj
        dataset%XS_vbf_Hpmj = CS_vbf_Hpmj
        dataset%XS_HpmjW = CS_HpmjW
        dataset%XS_HpmjZ = CS_HpmjZ
        dataset%XS_HpjHmj = CS_HpjHmj
        dataset%XS_Hpmjhi = CS_Hpmjhi
    end subroutine set_input
end subroutine HiggsBounds_charged_input_hadr


!> Couplings of the charged Higgs bosons to fermions.
!! The couplings are defines as in eq. 2 of the beyond Higgs paper.
subroutine HiggsBounds_charged_input_effC_fermions( &
    hcjud_L, hcjud_R, hcjcs_L, hcjcs_R, hcjtb_L, hcjtb_R, &
    hcjus_L, hcjus_R, hcjub_L, hcjub_R, hcjcd_L, hcjcd_R, &
    hcjcb_L, hcjcb_R, hcjtd_L, hcjtd_R, hcjts_L, hcjts_R)

    use usefulbits, only: theo, np, Hplus, effC_Hc, whichinput, just_after_run
#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------input
    double precision, intent(in) :: &
        hcjud_L(np(Hplus)), hcjud_R(np(Hplus)), &
        hcjcs_L(np(Hplus)), hcjcs_R(np(Hplus)), &
        hcjtb_L(np(Hplus)), hcjtb_R(np(Hplus)), &
        hcjus_L(np(Hplus)), hcjus_R(np(Hplus)), &
        hcjub_L(np(Hplus)), hcjub_R(np(Hplus)), &
        hcjcd_L(np(Hplus)), hcjcd_R(np(Hplus)), &
        hcjcb_L(np(Hplus)), hcjcb_R(np(Hplus)), &
        hcjtd_L(np(Hplus)), hcjtd_R(np(Hplus)), &
        hcjts_L(np(Hplus)), hcjts_R(np(Hplus))
    !--------------------------------------internal
    integer, parameter :: n = 1
    !----------------------------------------------

    whichinput = 'effC'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hplus) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_charged_input_effC_fermions should'
        write (*, *) 'only be called if np(Hplus)>0'
        stop 'error in subroutine HiggsBounds_charged_input_effC_fermions'
    endif

    effC_Hc(n)%hcjud_L = hcjud_L
    effC_Hc(n)%hcjud_R = hcjud_R
    effC_Hc(n)%hcjcs_L = hcjcs_L
    effC_Hc(n)%hcjcs_R = hcjcs_R
    effC_Hc(n)%hcjtb_L = hcjtb_L
    effC_Hc(n)%hcjtb_R = hcjtb_R
    effC_Hc(n)%hcjus_L = hcjus_L
    effC_Hc(n)%hcjus_R = hcjus_R
    effC_Hc(n)%hcjub_L = hcjub_L
    effC_Hc(n)%hcjub_R = hcjub_R
    effC_Hc(n)%hcjcd_L = hcjcd_L
    effC_Hc(n)%hcjcd_R = hcjcd_R
    effC_Hc(n)%hcjcb_L = hcjcb_L
    effC_Hc(n)%hcjcb_R = hcjcb_R
    effC_Hc(n)%hcjtd_L = hcjtd_L
    effC_Hc(n)%hcjtd_R = hcjtd_R
    effC_Hc(n)%hcjts_L = hcjts_L
    effC_Hc(n)%hcjts_R = hcjts_R

    just_after_run = .False.
end subroutine HiggsBounds_charged_input_effC_fermions

!> Get neutral Higgs boson's (SM-normalized) hadronic cross section.
!! This subroutines returns the SM-normalized hadronic cross sections for the neutral Higgs
!! boson \f$h_i\f$ and the specified hadron collider from the internal data in HiggsBounds.
!! @param i Index of the neutral Higgs boson
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param singleH (SM-normalized) cross section for single Higgs production
!! @param ggH (SM-normalized) cross section for Higgs production in gluon fusion, \f$ gg\to h_i\f$
!! @param bbH (SM-normalized) cross section for Higgs production in association with bottom quarks, \f$ b\bar{b}\to h_i\f$
!! @param VBF (SM-normalized) cross section for Higgs production in vector boson fusion
!! @param WH (SM-normalized) cross section for Higgs production in association with a W boson
!! @param ZH (SM-normalized) cross section for Higgs production in association with a Z boson
!! @param ttH (SM-normalized) cross section for Higgs production in association with a top quark pair
!! @param tH_tchan (SM-normalized) cross section for Higgs production in association with a top quark (t-channel)
!! @param tH_schan (SM-normalized) cross section for Higgs production in association with a top quark (s-channel)
!! @param qqZH (SM-normalized) cross section for qq-bar-initiated Higgs production in association with a Z boson
!! @param ggZH (SM-normalized) cross section for gg-initiated Higgs production in association with a Z boson
subroutine HiggsBounds_get_neutral_hadr_CS(i, collider, &
                                           singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH, ggZH)

    use usefulbits, only: theo, np, Hneut, hadroncolliderdataset
    implicit none
    integer, intent(in) :: i, collider
    double precision, intent(out) :: singleH, ggH, bbH, VBF, WH, ZH, ttH, tH_tchan, tH_schan, qqZH, ggZH

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (i .gt. np(Hneut)) then
        write (*, "(A,I2,A)") 'WARNING: Requested neutral Higgs h', i, ' not part of the model!'
    else
        select case (collider)
        case (2)
            call get_cross_section(theo(1)%tev)
        case (7)
            call get_cross_section(theo(1)%lhc7)
        case (8)
            call get_cross_section(theo(1)%lhc8)
        case (13)
            call get_cross_section(theo(1)%lhc13)
        case default
            stop 'wrong input for collider to subroutine HiggsBounds_get_neutral_SMnormalizedCS'
        end select
    endif

contains

    subroutine get_cross_section(dataset)
        type(hadroncolliderdataset) :: dataset

        singleH = dataset%XS_hj_ratio(i)
        ggH = dataset%XS_gg_hj_ratio(i)
        bbH = dataset%XS_bb_hj_ratio(i)
        VBF = dataset%XS_vbf_ratio(i)
        WH = dataset%XS_hjW_ratio(i)
        ZH = dataset%XS_hjZ_ratio(i)
        ttH = dataset%XS_tthj_ratio(i)
        tH_tchan = dataset%XS_thj_tchan_ratio(i)
        tH_schan = dataset%XS_thj_schan_ratio(i)
        qqZH = dataset%XS_qq_hjZ_ratio(i)
        ggZH = dataset%XS_gg_hjZ_ratio(i)

    end subroutine get_cross_section
end subroutine HiggsBounds_get_neutral_hadr_CS

!> Get neutral Higgs boson's branching ratios (for Higgs decays into SM final states).
!! This subroutine returns the branching ratios of the neutral Higgs boson \f$h_i\f$ decays
!! into SM final states from the internal data in HiggsBounds.
!! @param i Index of the neutral Higgs boson
!! @param BR_hjss Branching ratio of Higgs decay to strange quarks, \f$ h_i \to s \bar{s} \f$
!! @param BR_hjcc Branching ratio of Higgs decay to charm quarks, \f$ h_i \to c \bar{c} \f$
!! @param BR_hjbb Branching ratio of Higgs decay to bottom quarks, \f$ h_i \to b \bar{b} \f$
!! @param BR_hjtt Branching ratio of Higgs decay to top quarks, \f$ h_i \to t \bar{t} \f$
!! @param BR_hjmumu Branching ratio of Higgs decay to muons, \f$ h_i \to \mu^+ \mu^- \f$
!! @param BR_hjtautau Branching ratio of Higgs decay to tau leptons, \f$ h_i \to \tau^+ \tau^- \f$
!! @param BR_hjWW Branching ratio of Higgs decay to W bosons, \f$ h_i \to W^+W^- \f$
!! @param BR_hjZZ Branching ratio of Higgs decay to Z bosons, \f$ h_i \to ZZ \f$
!! @param BR_hjZga Branching ratio of Higgs decay to a Z boson and a photon, \f$ h_i \to Z\gamma \f$
!! @param BR_hjgaga Branching ratio of Higgs decay to photons, \f$ h_i \to \gamma\gamma \f$
!! @param BR_hjgg Branching ratio of Higgs decay to gluons, \f$ h_i \to gg \f$
subroutine HiggsBounds_get_neutral_BR(i, BR_hjss, BR_hjcc, BR_hjbb, &
                                      BR_hjtt, BR_hjmumu, BR_hjtautau, BR_hjWW, BR_hjZZ, BR_hjZga, &
                                      BR_hjgaga, BR_hjgg)

    use usefulbits, only: theo, np, Hneut
    implicit none

    integer, intent(in) :: i
    double precision, intent(out) :: BR_hjss, BR_hjcc, BR_hjbb, &
        BR_hjtt, BR_hjmumu, BR_hjtautau, BR_hjWW, BR_hjZZ, BR_hjZga, &
        BR_hjgaga, BR_hjgg

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (i .gt. np(Hneut)) then
        write (*, "(A,I2,A)") 'WARNING: Requested neutral Higgs h', i, ' not part of the model!'
    else
        BR_hjss = theo(1)%BR_hjss(i)
        BR_hjcc = theo(1)%BR_hjcc(i)
        BR_hjbb = theo(1)%BR_hjbb(i)
        BR_hjtt = theo(1)%BR_hjtt(i)
        BR_hjmumu = theo(1)%BR_hjmumu(i)
        BR_hjtautau = theo(1)%BR_hjtautau(i)
        BR_hjWW = theo(1)%BR_hjWW(i)
        BR_hjZZ = theo(1)%BR_hjZZ(i)
        BR_hjZga = theo(1)%BR_hjZga(i)
        BR_hjgaga = theo(1)%BR_hjgaga(i)
        BR_hjgg = theo(1)%BR_hjgg(i)
    endif
end subroutine HiggsBounds_get_neutral_BR

!> Get neutral Higgs boson's total width.
!! This subroutine returns the total width of the neutral Higgs boson \f$h_i\f$.
subroutine HiggsBounds_get_neutral_total_width(i, GammaTot, GammaTot_SM)

    use usefulbits, only: theo, np, Hneut
    implicit none

    integer, intent(in) :: i
    double precision, intent(out) :: GammaTot, GammaTot_SM

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (i .gt. np(Hneut)) then
        write (*, "(A,I2,A)") 'WARNING: Requested neutral Higgs h', i, ' not part of the model!'
    else
        GammaTot = theo(1)%particle(Hneut)%GammaTot(i)
        GammaTot_SM = theo(1)%particle(Hneut)%GammaTot_SM(i)
    endif
end subroutine HiggsBounds_get_neutral_total_width

!> Input of theoretical Higgs mass uncertainties to HiggsBounds.
!! @param dMhneut Theoretical mass uncertainties (in GeV) for the neutral Higgs bosons
!! @param dMhch Theoretical mass uncertainties (in GeV) for the charged Higgs bosons
subroutine HiggsBounds_set_mass_uncertainties(dMhneut, dMhch)
    use usefulbits, only: theo, np, Hneut, Hplus

    implicit none

    double precision, intent(in) :: dMhneut(np(Hneut))
    double precision, intent(in) :: dMhch(np(Hplus))

    theo(1)%particle(Hneut)%dMh = dMhneut
    theo(1)%particle(Hplus)%dMh = dMhch
end subroutine HiggsBounds_set_mass_uncertainties

!> Internal subroutine
!! @private
subroutine get_mass_variation_param(n)
    use usefulbits, only: theo, np, Hneut, Hplus, diffMhneut, diffMhch, ndmh, dmhsteps, small_mh
    implicit none

    integer, intent(in) :: n

    double precision :: dMhneut(np(Hneut))
    double precision :: dMhch(np(Hplus))
    integer :: km(np(Hneut) + np(Hplus))
    integer :: dm(dmhsteps**(np(Hneut) + np(Hplus)), np(Hneut) + np(Hplus))

    integer i, j, k, kp

    if (np(Hneut) .gt. 0) dMhneut = theo(n)%particle(Hneut)%dMh
    if (np(Hplus) .gt. 0) dMhch = theo(n)%particle(Hplus)%dMh

    if (modulo(dmhsteps, 2) .NE. 1) then
        stop 'Wrong number of steps in set_mass_uncertainty: must be odd (>=3)'
    endif

    ndmh = 0
    do i = 1, np(Hneut)

        IF (dMhneut(i) .GT. small_mh) THEN
            ndmh = ndmh + 1
        ENDIF
        km(i) = -(dmhsteps - 1) / 2
    enddo

    do i = 1, np(Hplus)
        IF (dMhch(i) .GT. small_mh) ndmh = ndmh + 1
        km(i + np(Hneut)) = -(dmhsteps - 1) / 2
    enddo

    IF (ndmh .EQ. 0) THEN
        RETURN
    ENDIF

    if (allocated(diffMhneut)) deallocate (diffMhneut)
    if (allocated(diffMhch)) deallocate (diffMhch)

    allocate (diffMhneut(dmhsteps**(np(Hneut) + np(Hplus)), np(Hneut)))
    allocate (diffMhch(dmhsteps**(np(Hneut) + np(Hplus)), np(Hplus)))

    k = 1
    do i = 1, dmhsteps**ndmh

        do j = 1, ndmh
            dm(i, j) = km(j)
        enddo

        km(k) = km(k) + 1

        do j = 2, ndmh
            IF (modulo(i, dmhsteps**(j - 1)) .EQ. 0) THEN
                km(j) = km(j) + 1
                km(j - 1) = -1
            ENDIF
        ENDDO

    enddo

    do i = 1, dmhsteps**ndmh
        k = 1

        do j = 1, np(Hneut)
            IF (dMhneut(j) .GT. small_mh) THEN
                diffMhneut(i, j) = theo(n)%particle(Hneut)%M(j) + dm(i, k) * dMhneut(k) / ((dmhsteps - 1) / 2)
                k = k + 1
            ELSE
                diffMhneut(i, j) = theo(n)%particle(Hneut)%M(j)
            ENDIF
        enddo

        kp = k
        do j = 1, np(Hplus)
            IF (dMhch(j) .GT. small_mh) THEN
                diffMhch(i, j) = theo(n)%particle(Hplus)%M(j) + dm(i, k) * dMhch(k - (kp - 1)) / ((dmhsteps - 1) / 2)
                k = k + 1
            ELSE
                diffMhch(i, j) = theo(n)%particle(Hplus)%M(j)
            ENDIF
        enddo
    enddo
end subroutine get_mass_variation_param

!> Runs HiggsBounds.
!! Runs HiggsBounds on the previously specified model point.
!! Per default, it calls \ref higgsbounds_subroutines::run_higgsbounds_full, i.e. for each Higgs boson in the model,
!! HiggsBounds first finds the most sensitive experimental analysis (out of the dataset specified earlier
!! by whichanalyses in #initialize_higgsbounds), and then
!! tests the predicted rate against the observed limit in this analysis. The output is then connected
!! by a logical "OR", i.e. if one Higgs boson of the model is excluded, then the model is regarded excluded.
!!
!! If the logical parameter usefulbits::run_HB_classic is set True, however, the old (version 3) method
!! of determining the output only from the experimental analysis that is most sensitive among all Higgs bosons
!! is being used.
!!
!! Note, that if many data points are tested at the same time (as for inputmethod==datfiles),
!! this subroutine only returns the results of the last datapoint. The full results are saved in fullHBres.
!!
!! @param HBresult Main binary HiggsBounds result,
!!  |   HBresult   |   description      |
!!  |--------------|--------------------|
!!  |      0       |   excluded (95%CL) |
!!  |      1       |   allowed (95%CL)  |
!!  |      -1      |   invalid point    |
!! @param chan Number of the channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!! @param obsratio Ratio of the predicted rate over the observed limit for this channel
!! @param ncombined Number of Higgs bosons that have been combined in order to calculate the obsratio
subroutine run_HiggsBounds(HBresult, chan, obsratio, ncombined)
    use usefulbits, only: np, Hneut, Hplus, run_HB_classic

    implicit none
    integer HBresult, chan, ncombined
    double precision obsratio

    integer hbres(0:np(Hneut) + np(Hplus)), hbchan(0:np(Hneut) + np(Hplus)), hbcomb(0:np(Hneut) + np(Hplus))
    double precision hbobs(0:np(Hneut) + np(Hplus))

    if (run_HB_classic .EQV. .True.) then
        call run_HiggsBounds_classic(HBresult, chan, obsratio, ncombined)
        return
    endif

    call run_HiggsBounds_full(hbres, hbchan, hbobs, hbcomb)

    HBresult = hbres(0)
    chan = hbchan(0)
    obsratio = hbobs(0)
    ncombined = hbcomb(0)
end subroutine run_HiggsBounds

!> Runs HiggsBounds and returns results for a single Higgs boson.
!! Runs HiggsBounds by internally calling \ref higgsbounds_subroutines::run_higgsbounds
!! and extracting the result for a single Higgs boson specified by the index h.
!! @param h Index of the Higgs boson for which the results should be given
!! @param HBresult Main binary HiggsBounds result,
!!  |   HBresult   |   description      |
!!  |--------------|--------------------|
!!  |      0       |   excluded (95%CL) |
!!  |      1       |   allowed (95%CL)  |
!!  |      -1      |   invalid point    |
!! @param chan Number of the channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!! @param obsratio Ratio of the predicted rate over the observed limit for this channel
!! @param ncombined Number of Higgs bosons that have been combined in order to calculate the obsratio
subroutine run_HiggsBounds_single(h, HBresult, chan, obsratio, ncombined)
    use usefulbits, only: np, Hneut, Hplus

    implicit none
    integer, intent(in) :: h
    integer, intent(out) :: HBresult, chan, ncombined
    double precision, intent(out) :: obsratio

    integer hbres(0:np(Hneut) + np(Hplus)), hbchan(0:np(Hneut) + np(Hplus)), hbcomb(0:np(Hneut) + np(Hplus))
    double precision hbobs(0:np(Hneut) + np(Hplus))

    IF (h .LT. 0) stop "Illegal number of Higgs boson: h < 0"
    if (h .GT. np(Hneut) + np(Hplus)) stop "Illegal number of Higgs boson"

    call run_HiggsBounds_full(hbres, hbchan, hbobs, hbcomb)

    HBresult = hbres(h)
    chan = hbchan(h)
    obsratio = hbobs(h)
    ncombined = hbcomb(h)
end subroutine run_HiggsBounds_single

!> Runs HiggsBounds.
!! Runs HiggsBounds on the previously specified model point. For each Higgs boson in the model,
!! HiggsBounds first finds the most sensitive experimental analysis (out of the dataset specified earlier
!! by usefulbits::whichanalyses in \ref higgsbounds_subroutines::initialize_higgsbounds), and then
!! tests the predicted rate against the observed limit in this analysis.
!!
!! The output arrays are of (length number of neutral and charged Higgs bosons) + 1.
!! The zeroth entry represents the global result (i.e. a logical *AND* combination of the
!! outcomes of the individual Higgs bosons), the remaining entries enumerate the outcome for
!! the neutral and charged Higgs bosons.
!! @param HBresult Main binary HiggsBounds result,
!!  |   HBresult   |   description      |
!!  |--------------|--------------------|
!!  |      0       |   excluded (95%CL) |
!!  |      1       |   allowed (95%CL)  |
!!  |      -1      |   invalid point    |
!! @param chan Number of the channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!! @param obsratio Ratio of the predicted rate over the observed limit for this channel
!! @param ncombined Number of Higgs bosons that have been combined in order to calculate the obsratio
subroutine run_HiggsBounds_full(HBresult, chan, &
                                obsratio, ncombined)

    use usefulbits, only: theo, res, just_after_run, ndmh, numres, &
                          np, Hneut, Hplus, dmhsteps, ndat, fullHBres, small_mh, &
                          HBresult_all, ncombined_all, chan_all, obsratio_all, predratio_all
    use channels, only: check_channels
    use theo_manip, only: HB5_complete_theo, HB5_recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------output
    integer, intent(out)::  HBresult(0:np(Hneut) + np(Hplus))
    integer, intent(out)::  chan(0:np(Hneut) + np(Hplus))
    integer, intent(out)::  ncombined(0:np(Hneut) + np(Hplus))
    double precision, intent(out) :: obsratio(0:np(Hneut) + np(Hplus))

    double precision :: Mhneut(np(Hneut))
    double precision :: Mhch(np(Hplus))
    !-------------------------------------internal
    integer :: n, i, j, ind, part, k
    !---------------------------------------------

    if (lbound(HBresult, dim=1) .NE. 0) stop "run_HiggsBounds_full: Array HBresult must begin with element 0"
    if (ubound(HBresult, dim=1) .NE. (np(Hneut) + np(Hplus))) then
        stop "run_HiggsBounds_full: Upper limit of array HBresult must be equal to number of Higgses"
    endif
    if (lbound(chan, dim=1) .NE. 0) stop "run_HiggsBounds_full: Array chan must begin with element 0"
    if (ubound(chan, dim=1) .NE. (np(Hneut) + np(Hplus))) then
        stop "run_HiggsBounds_full: Upper limit of array chan must be equal to number of Higgses"
    endif
    if (lbound(obsratio, dim=1) .NE. 0) stop "run_HiggsBounds_full: Array obsratio must begin with element 0"
    if (ubound(obsratio, dim=1) .NE. (np(Hneut) + np(Hplus))) then
        stop "run_HiggsBounds_full: Upper limit of array obsratio must be equal to number of Higgses"
    endif
    if (lbound(ncombined, dim=1) .NE. 0) stop "run_HiggsBounds_full: Array ncombined must begin with element 0"
    if (ubound(ncombined, dim=1) .NE. (np(Hneut) + np(Hplus))) then
        stop "run_HiggsBounds_full: Upper limit of array ncombined must be equal to number of Higgses"
    endif

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (.not. allocated(HBresult_all)) allocate (HBresult_all(0:np(Hneut) + np(Hplus), numres))
    if (.not. allocated(chan_all)) allocate (chan_all(0:np(Hneut) + np(Hplus), numres))
    if (.not. allocated(ncombined_all)) allocate (ncombined_all(0:np(Hneut) + np(Hplus), numres))
    if (.not. allocated(obsratio_all)) allocate (obsratio_all(0:np(Hneut) + np(Hplus), numres))
    if (.not. allocated(predratio_all)) allocate (predratio_all(0:np(Hneut) + np(Hplus), numres))

    call HB5_complete_theo

    do n = 1, ndat

        theo(n)%particle(Hneut)%Mc = theo(n)%particle(Hneut)%M
        theo(n)%particle(Hplus)%Mc = theo(n)%particle(Hplus)%M

        do i = 0, ubound(Hbresult, dim=1)
            obsratio_all(i, :) = -999d0
            predratio_all(i, :) = -999d0
            HBresult_all(i, :) = 1
            chan_all(i, :) = -999
            ncombined_all(i, :) = -999
            obsratio(i) = -999d0
            HBresult(i) = 1
            chan(i) = -999
            ncombined(i) = -999
        enddo

        do i = 1, np(Hneut)
            if (theo(n)%particle(Hneut)%dMh(i) .GT. small_mh) then
                ndmh = ndmh + 1
            endif
        enddo
        do i = 1, np(Hplus)
            if (theo(n)%particle(Hplus)%dMh(i) .GT. small_mh) then
                ndmh = ndmh + 1
            endif
        enddo

!  Do we have mass uncertainties to take care off
        IF (ndmh .GT. 0) THEN

            if (np(Hneut) .ne. 0) Mhneut = theo(n)%particle(Hneut)%M
            if (np(Hplus) .ne. 0) Mhch = theo(n)%particle(Hplus)%M

!   Loop over all Higgses
            do i = 1, np(Hneut) + np(Hplus)
                obsratio_all(i, :) = 1.D23

                IF (i .LE. np(Hneut)) THEN
                    ind = i
                    part = Hneut
                ELSE
                    ind = i - np(Hneut)
                    part = Hplus
                ENDIF

!     Check for mass steps for this particular Higgs boson
                IF (theo(n)%particle(part)%dMh(ind) .GT. small_mh) THEN
                    theo(n)%particle(part)%M(ind) = theo(n)%particle(part)%M(ind) &
                                                    - theo(n)%particle(part)%dMh(ind)
                    do j = 1, dmhsteps
!                        write(*,*) "Mass variation: ", theo(n)%particle(part)%M(ind), j
                        call HB5_recalculate_theo_for_datapoint(n)

                        call check_channels(theo(n), res(n), i)

                        do k = 1, size(res(n)%obsratio)
                            IF (res(n)%obsratio(k) .LT. obsratio_all(i, k)) THEN
                                HBresult_all(i, k) = res(n)%allowed95(k)
                                chan_all(i, k) = res(n)%chan(k)
                                obsratio_all(i, k) = res(n)%obsratio(k)
                                predratio_all(i, k) = res(n)%predratio(k)
                                ncombined_all(i, k) = res(n)%ncombined(k)
                            ENDIF
                        enddo
                        theo(n)%particle(part)%M(ind) = theo(n)%particle(part)%M(ind) &
                                                        + theo(n)%particle(part)%dMh(ind) / (dmhsteps - 1) * 2

                    enddo
                else
                    call HB5_recalculate_theo_for_datapoint(n)
                    call check_channels(theo(n), res(n), i)

                    do k = 1, size(res(n)%obsratio)
                        HBresult_all(i, k) = res(n)%allowed95(k)
                        chan_all(i, k) = res(n)%chan(k)
                        obsratio_all(i, k) = res(n)%obsratio(k)
                        predratio_all(i, k) = res(n)%predratio(k)
                        ncombined_all(i, k) = res(n)%ncombined(k)

                    enddo
                endif

                HBresult(i) = HBresult_all(i, 1)
                chan(i) = chan_all(i, 1)
                obsratio(i) = obsratio_all(i, 1)
                ncombined(i) = ncombined_all(i, 1)

!         Logical OR between exclusions (one Higgs excluded = combined exclusion)
                HBresult(0) = HBresult(0) * HBresult(i)

!    Save the data for the Higgs that has the highest ratio of theory/obs
                IF (obsratio(i) .GT. obsratio(0)) THEN
                    chan(0) = chan(i)
                    obsratio(0) = obsratio(i)
                    ncombined(0) = ncombined(i)
                ENDIF

                theo(n)%particle(Hneut)%M = Mhneut
                theo(n)%particle(Hplus)%M = Mhch

            enddo

        ELSE

            call HB5_recalculate_theo_for_datapoint(n)

            do i = 1, np(Hneut) + np(Hplus)
                call check_channels(theo(n), res(n), i)

                do k = 1, size(res(n)%obsratio)
                    HBresult_all(i, k) = res(n)%allowed95(k)
                    chan_all(i, k) = res(n)%chan(k)
                    obsratio_all(i, k) = res(n)%obsratio(k)
                    predratio_all(i, k) = res(n)%predratio(k)
                    ncombined_all(i, k) = res(n)%ncombined(k)
                enddo

                HBresult(i) = HBresult_all(i, 1)
                chan(i) = chan_all(i, 1)
                obsratio(i) = obsratio_all(i, 1)
                ncombined(i) = ncombined_all(i, 1)

                HBresult(0) = HBresult(0) * res(n)%allowed95(1)

                IF (obsratio(i) .GT. obsratio(0)) THEN
                    chan(0) = res(n)%chan(1)
                    obsratio(0) = res(n)%obsratio(1)
                    ncombined(0) = res(n)%ncombined(1)
                ENDIF
            enddo
        ENDIF

        fullHBres(n)%allowed95 = HBresult(0)
        fullHBres(n)%chan = chan(0)
        fullHBres(n)%obsratio = obsratio(0)
        fullHBres(n)%ncombined = ncombined(0)
    enddo

    just_after_run = .True.
end subroutine run_HiggsBounds_full

!> This routine provides information about the most sensitive channel(s) for a specific Higgs boson.
!! After the HiggsBounds run this subroutine provides information about the most sensitive channels
!! (with the position in the ranking specified by the user, up to maximal rank \ref usefulbits::numres)
!! for a given Higgs boson. The routine returns both the ratio of predicted rate over the expected limit (predratio),
!! and over the observed limit (obsratio).
!! @param nH Index of the Higgs boson, enumerating first over the neutral Higgs bosons then the charged Higgs bosons
!! @param pos Rank in sensitivity of the channel (with 1 being the most sensitive channel)
!! @param HBresult The HiggsBounds result (0,1 or -1) from this specific channel and Higgs boson
!! @param chan Number of the channel, as defined in Key.dat
!! @param obsratio Ratio of the predicted rate over the observed limit for this channel
!! @param predratio Ratio of the predicted rate over the expected limit for this channel
!! @param ncombined Number of Higgs bosons that have been combined in order to calculate the predicted rate
subroutine HiggsBounds_get_most_sensitive_channels_per_Higgs(nH, pos, HBresult, chan, obsratio, predratio, ncombined)
    use usefulbits, only: HBresult_all, obsratio_all, chan_all, ncombined_all, predratio_all, &
                          just_after_run, np, Hneut, Hplus, numres

    implicit none

    integer, intent(in) :: nH, pos
    integer, intent(out) :: HBresult, chan, ncombined
    double precision, intent(out) :: obsratio, predratio

    HBresult = 0
    chan = 0
    obsratio = 0
    predratio = 0
    ncombined = 0

    if (just_after_run .and. allocated(HBresult_all)) then
        if (nH .le. np(Hneut) + np(Hplus)) then
            if (pos .le. numres) then
                HBresult = HBresult_all(nH, pos)
                chan = chan_all(nH, pos)
                obsratio = obsratio_all(nH, pos)
                predratio = predratio_all(nH, pos)
                ncombined = ncombined_all(nH, pos)
            else
                write (*, *) 'WARNING: request exceeds the number of stored most sensitive channels (', numres, ')'
            endif
        else
            write (*, *) 'WARNING: requested Higgs boson is invalid (choose between 1 and ', np(Hneut) + np(Hplus), '!)'
        endif
    else
        write (*, *) 'WARNING: Please call run_HiggsBounds or run_HiggsBounds_full before calling', &
            ' HiggsBounds_get_most_sensitive_channels!'
    endif
end subroutine HiggsBounds_get_most_sensitive_channels_per_Higgs

!> This routine provides information about the overall most sensitive channel(s).
!! After the HiggsBounds run this provides information about the most sensitive channels
!! (with the position in the ranking specified by the user, up to maximal rank #usefulbits::numres).
!! The routine returns both the ratio of predicted rate over the expected limit (predratio),
!! and over the observed limit (obsratio).
!! @param pos Rank in sensitivity of the channel (with 1 being the most sensitive channel)
!! @param HBresult The HiggsBounds result (0,1 or -1) from this specific channel and Higgs boson
!! @param chan Number of the channel, as defined in Key.dat
!! @param obsratio Ratio of the predicted rate over the observed limit for this channel
!! @param predratio Ratio of the predicted rate over the expected limit for this channel
!! @param ncombined Number of Higgs bosons that have been combined in order to calculate the predicted rate
subroutine HiggsBounds_get_most_sensitive_channels(pos, HBresult, chan, obsratio, predratio, ncombined)
    use usefulbits, only: HBresult_all, obsratio_all, predratio_all, chan_all, ncombined_all, &
                          just_after_run, np, Hneut, Hplus, numres

    implicit none

    integer, intent(in) :: pos
    integer, intent(out) :: HBresult, chan, ncombined
    double precision, intent(out) :: obsratio, predratio
    integer :: i, j, k, count

    integer, allocatable :: nH_rank(:), pos_rank(:), posflat(:)
    double precision, allocatable :: predratio_tmp(:)

    allocate (nH_rank(numres), pos_rank(numres), posflat(numres), predratio_tmp(numres * (np(Hneut) + np(Hplus))))

    HBresult = 0
    chan = 0
    obsratio = 0
    ncombined = 0

    predratio_tmp = 0

    count = 0
    if (just_after_run .and. allocated(HBresult_all)) then
        if (pos .le. numres) then
            do j = 1, np(Hneut) + np(Hplus)
                do i = 1, numres
                    count = count + 1
                    predratio_tmp(count) = predratio_all(j, i)
                enddo
            enddo

            do i = 1, numres
                posflat(i) = maxloc(predratio_tmp, 1)
                predratio_tmp(posflat(i)) = -1.0D0
            enddo

            count = 0

            do j = 1, np(Hneut) + np(Hplus)
                do i = 1, numres
                    count = count + 1
                    do k = 1, numres
                        if (count .eq. posflat(k)) then
                            nH_rank(k) = j
                            pos_rank(k) = i
                        endif
                    enddo
                enddo
            enddo

            HBresult = HBresult_all(nH_rank(pos), pos_rank(pos))
            chan = chan_all(nH_rank(pos), pos_rank(pos))
            obsratio = obsratio_all(nH_rank(pos), pos_rank(pos))
            predratio = predratio_all(nH_rank(pos), pos_rank(pos))
            ncombined = ncombined_all(nH_rank(pos), pos_rank(pos))

        else
            write (*, *) 'WARNING: request exceeds the number of stored most sensitive channels (', numres, ')'
        endif
    else
        write (*, *) 'WARNING: Please call run_HiggsBounds or run_HiggsBounds_full before calling', &
            ' HiggsBounds_get_most_sensitive_channels!'
    endif

    deallocate (nH_rank, pos_rank, posflat, predratio_tmp)
end subroutine HiggsBounds_get_most_sensitive_channels

!> Run HiggsBounds in classic mode.
!! Uses the old (version 3) method of determining the output only from
!! the experimental analysis that is most sensitive among all Higgs bosons.
!!
!! **This should only be used for comparison purposes
!! or as a prerequisite to the LEP Chisq extension (see higgsbounds_get_lepchisq()).**
!!
!! @param HBresult Main binary HiggsBounds result,
!!  |   HBresult   |   description      |
!!  |--------------|--------------------|
!!  |      0       |   excluded (95%CL) |
!!  |      1       |   allowed (95%CL)  |
!!  |      -1      |   invalid point    |
!! @param chan Number of the channel predicted to have the highest statistical sensitivity, as defined in Key.dat
!! @param obsratio Ratio of the predicted rate over the observed limit for this channel
!! @param ncombined Number of Higgs bosons that have been combined in order to calculate the obsratio
subroutine run_HiggsBounds_classic(HBresult, chan, obsratio, ncombined)
    use usefulbits, only: theo, res, debug, just_after_run, ndmh, diffmhneut, diffmhch, &
                          np, Hneut, Hplus, full_dmth_variation, dmhsteps, ndat, fullHBres
    use channels, only: check_channels
    use theo_manip, only: HB5_complete_theo, HB5_recalculate_theo_for_datapoint

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    !----------------------------------------output
    integer, intent(out)::     HBresult, chan, ncombined
    double precision, intent(out) :: obsratio

    double precision :: Mhneut(np(Hneut))
    double precision :: Mhch(np(Hplus))
    !-------------------------------------internal
    integer :: n, i
    integer :: HBresult_tmp, chan_tmp, ncombined_tmp
    double precision :: obsratio_tmp
    !---------------------------------------------

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    HBresult_tmp = 0 ! take care of incorrect maybe-uninitialized warnings
    chan_tmp = 0 ! take care of incorrect maybe-uninitialized warnings

    call HB5_complete_theo

    do n = 1, ndat

        theo(n)%particle(Hneut)%Mc = theo(n)%particle(Hneut)%M

        call get_mass_variation_param(n)

        IF (ndmh .GT. 0) THEN

            if (np(Hneut) .ne. 0) Mhneut = theo(n)%particle(Hneut)%M
            if (np(Hplus) .ne. 0) Mhch = theo(n)%particle(Hplus)%M

            obsratio_tmp = 10.0E6 ! Set to very large initial value
            do i = 1, dmhsteps**ndmh
                theo(n)%particle(Hneut)%M = diffMhneut(i, :)
                theo(n)%particle(Hplus)%M = diffMhch(i, :)

                if (debug) write (*, *) 'manipulating input...'; call flush (6)
                call HB5_recalculate_theo_for_datapoint(n)

                if (debug) write (*, *) 'compare each data point to the experimental bounds...'; call flush (6)
                call check_channels(theo(n), res(n), 0)

                HBresult = res(n)%allowed95(1)
                chan = res(n)%chan(1)
                obsratio = res(n)%obsratio(1)
                ncombined = res(n)%ncombined(1)

                IF (.NOT. full_dmth_variation) THEN
                    IF (HBresult .EQ. 1) THEN
                        just_after_run = .True.
                        exit
                    ENDIF
                ELSE
                    IF (obsratio .lt. obsratio_tmp) THEN
                        HBresult_tmp = HBresult
                        chan_tmp = chan
                        obsratio_tmp = obsratio
                        ncombined_tmp = ncombined
                    ENDIF
                ENDIF
            enddo

            IF (full_dmth_variation) THEN
                HBresult = HBresult_tmp
                chan = chan_tmp
                obsratio = obsratio_tmp
                ncombined = ncombined
                just_after_run = .True.
            ENDIF

            theo(n)%particle(Hneut)%M = Mhneut
            theo(n)%particle(Hplus)%M = Mhch
            call HB5_recalculate_theo_for_datapoint(n)
            call check_channels(theo(n), res(n), 0)

        ELSE

            if (debug) write (*, *) 'manipulating input...'; call flush (6)
            call HB5_recalculate_theo_for_datapoint(n)

            if (debug) write (*, *) 'compare each data point to the experimental bounds...'; call flush (6)

            call check_channels(theo(n), res(n), 0)

            HBresult = res(n)%allowed95(1)
            chan = res(n)%chan(1)
            obsratio = res(n)%obsratio(1)
            ncombined = res(n)%ncombined(1)

            just_after_run = .True.
        ENDIF

        fullHBres(n)%allowed95 = HBresult
        fullHBres(n)%chan = chan
        fullHBres(n)%obsratio = obsratio
        fullHBres(n)%ncombined = ncombined

    enddo

    just_after_run = .True.
end subroutine run_HiggsBounds_classic

!> Get exclusion likelihoods.
!! Finds the combination of Higgs bosons (cluster) that gives the maximal likelihood
!! and returns the likelihood for this combination. See 1507.06706 for details.
!! @param analysisID likelihood for which analysis
!! @param Hindex index of the Higgs bosons that originated the dominant cluster
!! @param nc number of Higgs bosons in the dominant cluster
!! @param cbin binary code denoting which Higgs bosons are contained in the cluster.
!!             Calculated as \f$ \sum_\mathrm{Higgs} 2^{i-1} \f$ where the sum
!!             goes over all Higgs bosons in the clusters and \f$i\f$ is their index.
!!             For example, for the indexing \f$h1=h, h2=H, h3=A\f$, the combination
!!             H+A would give `cbin = 6` and a cluster formed only by the h gives `cbin = 1`.
!! @param M averaged mass value
!! @param llh value of the likelihood
!! @param obspred which likelihood to return, 'obs' = observed, 'pred' = expected/predicted
subroutine HiggsBounds_get_likelihood(analysisID, Hindex, nc, cbin, M, llh, obspred)
    use usefulbits, only: theo, np, Hneut
    use theo_manip, only: HB5_complete_theo
    use likelihoods, only: get_likelihood, calcpredratio_llh
    implicit none

    integer, intent(in) :: analysisID
    integer, intent(out) :: Hindex, nc, cbin
    double precision, intent(out) :: llh, M
    character(LEN=*), intent(in) :: obspred

    integer :: i
    double precision, allocatable :: expllh(:)

    double precision, allocatable :: mass(:)
    integer, allocatable :: nclist(:)

    allocate (expllh(np(Hneut)), mass(np(Hneut)), nclist(np(Hneut)))
    expllh = 0.0D0

    call HB5_complete_theo

! Determine most sensitive combination
    do i = 1, np(Hneut)
        call get_likelihood(analysisID, i, theo(1), expllh(i), mass(i), nclist(i), cbin, 'pred')
    enddo

    Hindex = maxloc(expllh, dim=1)

    call get_likelihood(analysisID, Hindex, theo(1), llh, M, nc, cbin, obspred)

    deallocate (mass, nclist, expllh)
end subroutine HiggsBounds_get_likelihood

!> Get exclusion likelihoods involving the specified Higgs boson.
!! Similar to higgsbounds_get_likelihood_for_comb() but only consider clusters
!! involving the Higgs specified by Hindex.
!! @param analysisID likelihood for which analysis
!! @param cbin_in binary code indicating which Higgs bosons are **not** to
!!                be included in the clustering,
!!                see higgsbounds_get_likelihood() for details.
!! @param Hindex index of the Higgs bosons to always include
!! @param nc number of Higgs bosons in the dominant cluster
!! @param cbin binary code denoting which Higgs bosons are contained in the cluster,
!!             see higgsbounds_get_likelihood() for details.
!! @param M averaged mass value
!! @param llh value of the likelihood
!! @param obspred which likelihood to return, 'obs' = observed, 'pred' = expected/predicted
subroutine HiggsBounds_get_likelihood_for_Higgs(analysisID, cbin_in, Hindex, nc, cbin, M, llh, obspred)
    use usefulbits, only: theo
    use theo_manip, only: HB5_complete_theo
    use likelihoods, only: get_likelihood, calcpredratio_llh
    implicit none

    integer, intent(in) :: analysisID, Hindex
    integer, intent(out) ::  nc, cbin
    double precision, intent(out) :: llh, M
    integer, intent(in) :: cbin_in
    character(LEN=*), intent(in) :: obspred

    call HB5_complete_theo

    call get_likelihood(analysisID, Hindex, theo(1), llh, M, nc, cbin, obspred, cbin_in)
end subroutine HiggsBounds_get_likelihood_for_Higgs

!> Get exclusion likelihoods with some Higgs bosons excluded from the clustering.
!! Finds the combination of Higgs bosons (cluster) that gives the best exclusion
!! without including any of the Higgs bosons indivated by cbin_in
!! and returns the likelihood for this combination. See 1507.06706 for details.
!! @param analysisID likelihood for which analysis
!! @param cbin_in binary code indicating which Higgs bosons are **not** to
!!                be included in the clustering,
!!                see higgsbounds_get_likelihood() for details.
!! @param Hindex index of the Higgs bosons that originated the dominant cluster
!! @param nc number of Higgs bosons in the dominant cluster
!! @param cbin binary code denoting which Higgs bosons are contained in the cluster,
!!             see higgsbounds_get_likelihood() for details.
!! @param M averaged mass value
!! @param llh value of the likelihood
!! @param obspred which likelihood to return, 'obs' = observed, 'pred' = expected/predicted
subroutine HiggsBounds_get_likelihood_for_comb(analysisID, cbin_in, Hindex, nc, cbin, M, llh, obspred)
    use usefulbits, only: theo, np, Hneut
    use theo_manip, only: HB5_complete_theo
    use likelihoods, only: get_likelihood, calcpredratio_llh
    implicit none

    integer, intent(in) :: analysisID, cbin_in
    integer, intent(out) :: Hindex, nc, cbin
    double precision, intent(out) :: llh, M
    character(LEN=*), intent(in) :: obspred

    integer :: i
    double precision, allocatable :: obsllh(:)
    double precision, allocatable :: mass(:)
    integer, allocatable :: nclist(:), cbinlist(:)

    allocate (obsllh(np(Hneut)), mass(np(Hneut)), nclist(np(Hneut)), cbinlist(np(Hneut)))
    obsllh = 0.0D0

    call HB5_complete_theo

! Determine most sensitive combination
    do i = 1, np(Hneut)
        call get_likelihood(analysisID, i, theo(1), obsllh(i), mass(i), nclist(i), cbinlist(i), obspred, cbin_in)
    enddo

    Hindex = maxloc(obsllh, dim=1)
    llh = obsllh(Hindex)
    M = mass(Hindex)
    nc = nclist(Hindex)
    cbin = cbinlist(Hindex)

    deallocate (mass, nclist, obsllh, cbinlist)
end subroutine HiggsBounds_get_likelihood_for_comb

!> Outputs the HiggsBounds results to the SLHA file.
!! Can only be used if higgsbounds_input_slha() was used beforehand.
!! Writes to the file specified there.
subroutine HiggsBounds_SLHA_output
    use usefulbits, only: whichinput, just_after_run
    use output, only: do_output
    implicit none

    if (.not. just_after_run) then
        stop 'subroutine run_HiggsBounds should be called before subroutine HiggsBounds_SLHA_output'
    endif

    select case (whichinput)
    case ('SLHA')
        call do_output
    case default
        stop 'The subroutine HiggsBounds_SLHA_output should only be used when whichinput=SLHA'
    end select
end subroutine HiggsBounds_SLHA_output

#ifdef enableCHISQ

!> Initialize the LEP Chisq extension.
!! If the extension is used this needs to be called exactly once and
!! before the call to #initialize_higgsbounds
subroutine initialize_HiggsBounds_chisqtables
    use S95tables_type3
    use usefulbits, only: allocate_if_stats_required, theo
    implicit none

    if (allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds_chisqtables should be called before subroutine initialize_HiggsBounds'
    elseif (allocated(clsb_t3)) then
        stop 'subroutine initialize_HiggsBounds_chisqtables has already been called once'
    endif

    allocate (clsb_t3(ntable3))

    call initializetables_type3_blank(clsb_t3)
    call initializetables3(clsb_t3)

    call readclsbfiles_binary

    if (allocated(allocate_if_stats_required)) then
        stop 'error in subroutine initialize_HiggsBounds_chisqtables'
    else
        allocate (allocate_if_stats_required(1))
    endif
end subroutine initialize_HiggsBounds_chisqtables

!> Deallocate data structures used by the LEP Chisq extension.
!! This needs to be called exactly once when you are done using HiggsBounds.
subroutine finish_HiggsBounds_chisqtables
    use S95tables_type3
    use usefulbits, only: allocate_if_stats_required
    implicit none
    integer :: x

    if (.not. allocated(clsb_t3)) then
        stop 'initialize_HiggsBounds_chisqtables should be called first'
    endif

    do x = lbound(clsb_t3, dim=1), ubound(clsb_t3, dim=1)
        deallocate (clsb_t3(x)%dat)
    enddo
    deallocate (filename)
    deallocate (clsb_t3)
    deallocate (allocate_if_stats_required)
end subroutine finish_HiggsBounds_chisqtables

!> Calculate the LEP chisq value.
!! Evaluates the \f$\chi^2\f$ of the tabulated LEP results.
!! Only run this after calling run_higgsbounds_classic().
!! @param theory_uncertainty_1s 1 sigma mass uncertainty to use in the chisq calculation.
!!                              This is a separate value from the one set by higgsbounds_set_mass_uncertainties().
!! @param chisq_withouttheory \f$\chi^2\f$-value ignoring the theory_uncertainty_1s
!! @param chisq_withtheory \f$\chi^2\f$-value including the theory_uncertainty_1s
!! @param channel code indicating which analysis was used to derive the \f$\chi^2\f$-value
subroutine HiggsBounds_get_LEPChisq(theory_uncertainty_1s, chisq_withouttheory, chisq_withtheory, channel)
    use usefulbits, only: res, theo, pr, just_after_run, vsmall
    use interpolate
    use S95tables_type1
    use S95tables_type3
    use S95tables
    use extra_bits_for_chisquared
    implicit none

    integer, intent(out) :: channel
    integer :: x, c, z, y
    integer :: id
    double precision, intent(in) :: theory_uncertainty_1s
    double precision, intent(out) :: chisq_withouttheory, chisq_withtheory
    double precision :: low_chisq, sigma

    x = 1
    low_chisq = 1.0D-2

    if (.not. allocated(theo)) then
        stop 'subroutine HiggsBounds_initialize must be called first'
    elseif (.not. allocated(clsb_t3)) then
        stop 'subroutine initialize_HiggsBounds_chisqtables must be called first'
    elseif (.not. just_after_run) then
        stop 'subroutine run_HiggsBounds must be called first'
    endif

    sigma = theory_uncertainty_1s
    if (sigma .lt. vsmall) then
        write (*, *) 'Warning: will not calculate chi^2 with theory uncertainty'
    endif

    chisq_withtheory = -2.0D0
    chisq_withouttheory = -2.0D0

    z = 2;
    c = res(x)%chan(z)
    channel = c

    if (res(x)%allowed95(z) .eq. -1) then! labels an unphysical parameter point
        chisq_withtheory = -1.0D0
        chisq_withouttheory = -1.0D0
    elseif (c .gt. 0) then ! labels a physical parameter point and a real channel

        id = S95_t1_or_S95_t2_idfromelementnumber(pr(c)%ttype, pr(c)%tlist)
        y = clsb_t3elementnumber_from_S95table(pr(c)%ttype, id)

        if (y .gt. 0) then
            call get_chisq(sigma, res(x)%axis_i(z), res(x)%axis_j(z), res(x)%sfactor(z), &
                           y, chisq_withouttheory, chisq_withtheory)
        else
            write (*, *) 'hello y=', y
            stop 'problem here with y'
        endif

    else
        chisq_withtheory = 0.0D0
        chisq_withouttheory = 0.0D0
    endif
end subroutine HiggsBounds_get_LEPChisq
#endif

!> Deallocate data structures and close files.
!! This needs to be called exactly once when you are done using HiggsBounds.
subroutine finish_HiggsBounds
    use usefulbits, only: deallocate_usefulbits, debug, theo, debug, &
                          file_id_debug1, file_id_debug2
    use S95tables, only: deallocate_S95tables
    use theory_BRfunctions, only: deallocate_BRSM
    use theory_XS_SM_functions, only: deallocate_XSSM

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    if (debug) then
        close (file_id_debug2)
        close (file_id_debug1)
    endif

    if (.not. allocated(theo)) then
        stop 'initialize_HiggsBounds  should be called first'
    endif

    if (debug) write (*, *) 'finishing off...'; call flush (6)
    call deallocate_BRSM
    call deallocate_XSSM
    call deallocate_S95tables
    call deallocate_usefulbits

    if (debug) write (*, *) 'finished'; call flush (6)
end subroutine finish_HiggsBounds

!> Set a single effective coupling to the given value(s).
!! Only sets the specified couling leaving the other ones unchanged.
!! See also higgsbounds_neutral_input_effc() and higgsbounds_neutral_input_effc_double().
!! @param quantity which quantity to set, valid values for the different coupligns are
!! identical to the argument names of higgsbounds_neutral_input_effc()
!! @param val values to set
subroutine HiggsBounds_neutral_input_effC_single(quantity, val)
    use usefulbits, only: theo, np, Hneut, effC, whichinput, just_after_run
    implicit none
    !---------------------------------------------
    character(LEN=*), intent(in) :: quantity
    double precision, intent(in) :: val(np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'effC'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_effC_single should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_effC_single'
    endif

    select case (trim(adjustl(quantity)))
    case ("ghjcc_s")
        effC(n)%hjcc_s = val
    case ("ghjcc_p")
        effC(n)%hjcc_p = val
    case ("ghjss_s")
        effC(n)%hjss_s = val
    case ("ghjss_p")
        effC(n)%hjss_p = val
    case ("ghjbb_s")
        effC(n)%hjbb_s = val
    case ("ghjbb_p")
        effC(n)%hjbb_p = val
    case ("ghjtt_s")
        effC(n)%hjtt_s = val
    case ("ghjtt_p")
        effC(n)%hjtt_p = val
    case ("ghjmumu_s")
        effC(n)%hjmumu_s = val
    case ("ghjmumu_p")
        effC(n)%hjmumu_p = val
    case ("ghjtautau_s")
        effC(n)%hjtautau_s = val
    case ("ghjtautau_p")
        effC(n)%hjtautau_p = val
    case ("ghjWW")
        effC(n)%hjWW = val
    case ("ghjZZ")
        effC(n)%hjZZ = val
    case ("ghjZga")
        effC(n)%hjZga = val
    case ("ghjgaga")
        effC(n)%hjgaga = val
    case ("ghjgg")
        effC(n)%hjgg = val

    case default
        stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_effC_single'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_effC_single

!> Set a single 2d coupling matrix of the effective coupling input.
!! Only sets the specified couling leaving the other ones unchanged.
!! See also higgsbounds_neutral_input_effc() and higgsbounds_neutral_input_effc_double().
!! @param quantity which coupling to set, currently only 'ghjhiZ' is supported
!! @param val value to set
subroutine HiggsBounds_neutral_input_effC_double(quantity, val)
    use usefulbits, only: theo, np, Hneut, effC, whichinput, just_after_run
    implicit none
    !---------------------------------------------
    character(LEN=*), intent(in) :: quantity
    double precision, intent(in) :: val(np(Hneut), np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'effC'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_effC_double should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_effC_double'
    endif

    select case (trim(adjustl(quantity)))
    case ("ghjhiZ")
        effC(n)%hjhiZ = val
    case default
        stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_effC_double'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_effC_double

!> Set a single LEP cross section to the given value(s).
!! Only sets the specified cxn leaving the other ones unchanged.
!! See also higgsbounds_neutral_input_lep() and higgsbounds_neutral_input_lep_double().
!! @param quantity which quantity to set, valid values for the different cross sections are
!! identical to the argument names of higgsbounds_neutral_input_lep()
!! @param val values to set
subroutine HiggsBounds_neutral_input_LEP_single(quantity, val)
    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run
    implicit none
    !---------------------------------------------
    character(LEN=*), intent(in) :: quantity
    double precision, intent(in) :: val(np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_LEP_single should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_LEP_single'
    endif

    select case (trim(adjustl(quantity)))
    case ("XS_hjZ_ratio")
        theo(n)%lep%XS_hjZ_ratio = val
    case ("XS_bbhj_ratio")
        theo(n)%lep%XS_bbhj_ratio = val
    case ("XS_tautauhj_ratio")
        theo(n)%lep%XS_tautauhj_ratio = val
    case default
        stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_LEP_single'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_LEP_single

!> Set a single 2d LEP cross section to the given value(s).
!! Only sets the specified cxn leaving the other ones unchanged.
!! See also higgsbounds_neutral_input_lep() and higgsbounds_neutral_input_lep_single().
!! @param quantity which quantity to set, currently only 'XS_ee_hjhi_ratio' is supported
!! @param val values to set
subroutine HiggsBounds_neutral_input_LEP_double(quantity, val)
    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run
    implicit none
    !---------------------------------------------
    character(LEN=*), intent(in) :: quantity
    double precision, intent(in) :: val(np(Hneut), np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_LEP_double should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_LEP_double'
    endif

    select case (trim(adjustl(quantity)))
    case ("XS_hjhi_ratio")
        theo(n)%lep%XS_hjhi_ratio = val
    case default
        stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_LEP_double'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_LEP_double

!> Set a single hadron collider cross section to the given value(s).
!! Only sets the specified cxn leaving the other ones unchanged.
!! See also higgsbounds_neutral_input_hadr() and higgsbounds_neutral_input_hadr_double().
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param quantity which quantity to set, valid values for the different cross sections are
!! identical to the argument names of higgsbounds_neutral_input_hadr()
!! @param val values to set
subroutine HiggsBounds_neutral_input_hadr_single(collider, quantity, val)
    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run, &
                          hadroncolliderdataset
    implicit none
    !---------------------------------------------
    integer, intent(in) :: collider
    character(LEN=*), intent(in) :: quantity
    double precision, intent(in) :: val(np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_hadr_single should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_hadr_single'
    endif

    select case (collider)
    case (2)
        call set_input(theo(n)%tev)
    case (7)
        call set_input(theo(n)%lhc7)
    case (8)
        call set_input(theo(n)%lhc8)
    case (13)
        call set_input(theo(n)%lhc13)
    case default
        stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_single'
    end select

    just_after_run = .False.

contains

    subroutine set_input(dataset)
        type(hadroncolliderdataset) :: dataset
        select case (trim(adjustl(quantity)))
        case ("XS_hj_ratio")
            dataset%XS_hj_ratio = val
        case ("XS_gg_hj_ratio")
            dataset%XS_gg_hj_ratio = val
        case ("XS_bb_hj_ratio")
            dataset%XS_bb_hj_ratio = val
            dataset%XS_hjb_ratio = val
        case ("XS_vbf_ratio")
            dataset%XS_vbf_ratio = val
        case ("XS_hjZ_ratio")
            dataset%XS_hjZ_ratio = val
        case ("XS_gg_hjZ_ratio")
            dataset%XS_gg_hjZ_ratio = val
        case ("XS_qq_hjZ_ratio")
            dataset%XS_qq_hjZ_ratio = val
        case ("XS_hjW_ratio")
            dataset%XS_hjW_ratio = val
        case ("XS_tthj_ratio")
            dataset%XS_tthj_ratio = val
        case ("XS_thj_tchan_ratio")
            dataset%XS_thj_tchan_ratio = val
        case ("XS_thj_schan_ratio")
            dataset%XS_thj_schan_ratio = val
        case ("XS_tWhj_ratio")
            dataset%XS_tWhj_ratio = val
        case default
            stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_hadr_single'
        end select
    end subroutine set_input
end subroutine HiggsBounds_neutral_input_hadr_single

!> Set a single 2d hadron collider cross section to the given value(s).
!! Only sets the specified cxn leaving the other ones unchanged.
!! See also higgsbounds_neutral_input_hadr() and higgsbounds_neutral_input_hadr_single().
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param quantity which quantity to set, currently only 'XS_hjhi' is supported
!! @param val values to set
subroutine HiggsBounds_neutral_input_hadr_double(collider, quantity, val)
    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run, &
                          hadroncolliderdataset
    implicit none
    !---------------------------------------------
    integer, intent(in) :: collider
    character(LEN=*), intent(in) :: quantity
    double precision, intent(in) :: val(np(Hneut), np(Hneut))
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'

    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (np(Hneut) .eq. 0) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_hadr_double should'
        write (*, *) 'only be called if np(Hneut)>0'
        stop 'error in subroutine HiggsBounds_neutral_input_hadr_double'
    endif

    select case (trim(adjustl(quantity)))
    case ("XS_hjhi")
        select case (collider)
        case (2)
            theo(n)%tev%XS_hjhi = val
        case (7)
            theo(n)%lhc7%XS_hjhi = val
        case (8)
            theo(n)%lhc8%XS_hjhi = val
        case (13)
            theo(n)%lhc13%XS_hjhi = val
        case default
            stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_double'
        end select
    case default
        stop 'wrong input for quantity to subroutine HiggsBounds_neutral_input_hadr_double'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_hadr_double

!> Sets a single channelrate.
!! Input an element of the channelrate matrix, see #higgsbounds_neutral_input_hadr_channelrates()
!! Elements of the matrix channelrates with values < 0 will be overwritten by XS times BR
!! using the narrow width approximation.
!! @param collider selects the collider experiment for which the input is given
!!  | collider | experiment    |
!!  |----------|---------------|
!!  |    2     | Tevatron      |
!!  |    7     | LHC at 7 TeV  |
!!  |    8     | LHC at 8 TeV  |
!!  |    13     | LHC at 13 TeV|
!! @param nHiggs Index of the Higgs boson
!! @param p Index of the production mode
!! @param d Index of the decay mode
!! @param val value to set
subroutine HiggsBounds_neutral_input_hadr_channelrates_single(collider, nHiggs, p, d, val)
    use usefulbits, only: theo, np, Hneut, whichinput, just_after_run, hadroncolliderdataset

#if defined(NAGf90Fortran)
    use F90_UNIX_IO, only: flush
#endif

    implicit none
    double precision, intent(in) :: val
    integer, intent(in) :: collider, p, d, nHiggs
    !-------------------------------------internal
    integer, parameter :: n = 1
    !---------------------------------------------
    whichinput = 'hadr'
    if (.not. allocated(theo)) then
        stop 'subroutine initialize_HiggsBounds must be called first'
    endif

    if (nHiggs .gt. np(Hneut)) then
        write (*, *) 'subroutine HiggsBounds_neutral_input_hadr_channelrates_single should'
        write (*, *) 'only be called with nHiggs <= np(Hneut)'
        stop 'error in subroutine HiggsBounds_neutral_input_hadr_channelrates_single'
    endif

    select case (collider)
    case (2)
        theo(n)%tev%channelrates_tmp(nHiggs, p, d) = val
    case (7)
        theo(n)%lhc7%channelrates_tmp(nHiggs, p, d) = val
    case (8)
        theo(n)%lhc8%channelrates_tmp(nHiggs, p, d) = val
    case (13)
        theo(n)%lhc13%channelrates_tmp(nHiggs, p, d) = val
    case default
        stop 'wrong input for collider to subroutine HiggsBounds_neutral_input_hadr_channelrates_single'
    end select

    just_after_run = .False.
end subroutine HiggsBounds_neutral_input_hadr_channelrates_single

!> Resets the channelrates.
!! Use to undo values set by HiggsBounds_neutral_input_hadr_channelrates_single()
!! and HiggsBounds_neutral_input_hadr_channelrates()
subroutine HiggsBounds_neutral_input_hadr_channelrates_clean
    use theo_manip, only: clean_channelrates
    implicit none
    call clean_channelrates
end subroutine HiggsBounds_neutral_input_hadr_channelrates_clean
