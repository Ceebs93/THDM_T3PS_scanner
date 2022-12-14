PROGRAM AllAnalyses
    use S95tables
    use likelihoods
    use theory_BRfunctions, only: setup_BRSM, deallocate_BRSM
    use usefulbits, only: listprocesses

    integer I
    character(LEN=200) :: descrip
    type(listprocesses) :: proc

    call setup_BRSM
    call setup_S95tables
    call setup_likelihoods

    PRINT *, ""
    PRINT *, " LISTING ALL HIGGSBOUNDS ANALYSES"
    PRINT *, ""
    PRINT *, " Tables of Type I", ubound(S95_t1, dim=1), "total"
    PRINT *, ""
    PRINT *, "      ID   expt     E   lumi   part  label                             " &
        //"           cite key                               Mmin    Mmax"
    PRINT *, " -----------------------------------------------------------------------------------" &
        //"-----------------------------------------------------"
    proc%ttype = 1
    DO I = lbound(S95_t1, dim=1), ubound(S95_t1, dim=1)
        proc%tlist = I
        call outputproc(proc, 21, descrip, 0)
        WRITE (*, '(I10,A,A,A,F6.2,F6.1,I6,A,A,A,2F8.1,A,A)') &
            S95_t1(I)%id, "  ", adjustl(S95_t1(I)%expt), "  ", &
            S95_t1(I)%energy, S95_t1(I)%lumi, &
            S95_t1(I)%particle_x, "   ", ADJUSTL(S95_t1(I)%label), &
            ADJUSTL(S95_t1(I)%citekey), S95_t1(I)%xmin, S95_t1(I)%xmax, "   # ", trim(adjustl(descrip))
    ENDDO
    PRINT *, " -----------------------------------------------------------------------------------" &
        //"-----------------------------------------------------"

    PRINT *, ""
    PRINT *, " Tables of Type II", ubound(S95_t2, dim=1), "total"
    PRINT *, ""
    PRINT *, "      ID   expt     E   lumi   part  label                                     " &
        //"   cite key                               Mmin1   Mmax1   Mmin2   Mmax2"
    PRINT *, " -----------------------------------------------------------------------------------" &
        //"---------------------------------------------------------------------"
    proc%ttype = 2
    DO I = lbound(S95_t2, dim=1), ubound(S95_t2, dim=1)
        proc%tlist = I
        call outputproc(proc, 21, descrip, 0)
        WRITE (*, '(I10,A,A,A,F6.2,F6.1,I4,I3,A,A,A,4F8.1,A,A)') &
            S95_t2(I)%id, "  ", adjustl(S95_t2(I)%expt), "  ", S95_t2(I)%energy, &
            S95_t2(I)%lumi, S95_t2(I)%particle_x1, S95_t2(I)%particle_x2, &
            "  ", ADJUSTL(S95_t2(I)%label), ADJUSTL(S95_t2(I)%citekey), &
            S95_t2(I)%xmin1, S95_t2(I)%xmax1, S95_t2(I)%xmin2, S95_t2(I)%xmax2, &
            "   # ", trim(adjustl(descrip))
    ENDDO
    PRINT *, " -----------------------------------------------------------------------------------" &
        //"---------------------------------------------------------------------"
    PRINT *, ""
    PRINT *, " Likelihood Tables", nllhs, "total"
    PRINT *, ""
    PRINT *, "      ID   expt     E   lumi   part  label                                " &
        //"        citekey                                Mmin    Mmax"
    PRINT *, " -----------------------------------------------------------------------------" &
        //"-----------------------------------------------------------"
    DO I = 1, nllhs
        WRITE (*, '(I10,A,A,A,F6.2,F6.1,I6,A,A,A,2F8.1,A,A)') &
            llhdata(I)%D2llhdata(1)%analysisID, "  ", adjustl(llhdata(I)%D2llhdata(1)%expt), "  ", &
            llhdata(I)%D2llhdata(1)%energy, llhdata(I)%D2llhdata(1)%lumi, &
            llhdata(I)%D2llhdata(1)%particle_x, "   ", adjustl(llhdata(I)%D2llhdata(1)%label), &
            adjustl(llhdata(I)%D2llhdata(1)%citekey), llhdata(I)%D2llhdata(1)%mass, &
            llhdata(I)%D2llhdata(n2Dslices(I))%mass, "   # ", adjustl(llhdata(I)%D2llhdata(1)%description)
    ENDDO
    PRINT *, " ------------------------------------------------------------------------------" &
        //"-----------------------------------------------------------"
    print *, ""
END
