program HSLimit
    use usefulbits, only: inputmethod, np, Hneut, Hplus, whichanalyses, vsmall, theo
    use input, only: do_input
    use usefulbits_hs, only: output_level, HiggsSignals_info, Exptdir, Nparam
    use io, only: setup_input_for_hs, do_output_for_hs

    implicit none
    double precision :: sHat, sUp, CLs, Chisq_mu, CHisq_mh, Chisq, Pvalue
    double precision :: mh, BRinv, kappa, kappa0
    integer :: i, j, nobs
    integer :: outfile = 11111

    mh = 125.09
    kappa = 1.0D0

    open (outfile, file="results/HSLimitEffC.dat")
    write (outfile, "(6(A16,10X))") "kappa", "BRinv", "sHat", "sUp", "CLs", "chisq"

    call initialize_HiggsSignals(1, 0, "latestresults")

    do i = 0, 20, 2
        BRinv = i*1D-2
        do j = 0, 10, 1
            kappa = 0.9 + j*1D-2
            call HiggsBounds_neutral_input_properties(mh, -1D0, 1)
            call HiggsBounds_neutral_input_effC( &
                kappa, 0D0, kappa, 0D0, &
                kappa, 0D0, kappa, 0D0, &
                kappa, 0D0, &
                kappa, 0D0, &
                kappa, kappa, kappa, &
                kappa, kappa, 0D0)
            call HiggsBounds_neutral_input_nonSMBR(BRinv, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0)

            call set_HiggsSignals_limit(sHat, sUp, CLs)
            call run_HiggsSignals_full(Chisq_mu, CHisq_mh, Chisq, nobs, Pvalue)
            write (outfile, *) kappa, BRinv, sHat, sUp, CLs, Chisq
            print *, kappa, BRinv
        enddo
    enddo
    call finish_HiggsSignals
end program HSLimit
