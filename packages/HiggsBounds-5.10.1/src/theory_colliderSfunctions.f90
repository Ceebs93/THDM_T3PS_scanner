! This file is part of HiggsBounds

module theory_colliderSfunctions

contains

    real*8 function lhc7_rHVBF_WW(x)
        !******************************************************
        !* Used HAWK with settings of input parameters of arXiv:1101.0593 [hep-ph].
        !* (PDF: MSTW 2008 NNLO)
        !* Calculated Born CS with only t-channel graphs (*),
        !* (a) with only WW-fusion graphs and
        !* (b) with all WW/ZZ fusion graphs, and calculated
        !* ratio (a)/(b).
        !* (*): s-channel graphs are neglected because they are
        !*      highly suppressed after VBF cuts.
        !******************************************************
        !* x : Higgs mass in GeV
        !* fit: strictly valid in range [90:1000], deviations from data below 0.3%
        !*      slight extrapolation to range [70:1020] is allowed.
        !* 27/4/2011, Oliver Brein (revised: 29/4/2011)
        !******************************
        implicit none
        real*8 a, b, c, d, e, f, g, h, i, j, k, x
        real*8 a0, b0, c0, d0

        a0 = 0.764629881486091
        b0 = -0.000128901014699432
        c0 = 2.248155381437e-07
        d0 = -2.38464077139572e-10

        i = 353411.391484508
        j = -13532.7423394166
        k = 216.248242966815
        a = -1.13447368036044
        b = 0.0100380489651014
        c = -3.45357638530139e-05
        d = 7.69000469151772e-08
        e = -1.10305464638084e-10
        f = 9.79724660472145e-14
        g = -4.88552475279702e-17
        h = 1.04278241770435e-20

        lhc7_rHVBF_WW = 0d0

        if (x .lt. 70d0) then
            x = 70d0
            lhc7_rHVBF_WW = a0 + b0*x + c0*x**2 + d0*x**3
        endif
        if ((x .ge. 70d0) .and. (x .lt. 200d0)) then
            lhc7_rHVBF_WW = a0 + b0*x + c0*x**2 + d0*x**3
        endif
        if ((x .ge. 200d0) .and. (x .le. 1020d0)) then
            lhc7_rHVBF_WW = i/x**3 + j/x**2 + k/x + a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6 + h*x**7
        endif
        if (x .gt. 1020d0) then
        endif

    end function

    real*8 function lhc7_rHVBF_ZZ(x)
        !******************************************************
        !* Used HAWK with settings of input parameters of arXiv:1101.0593 [hep-ph].
        !* (PDF: MSTW 2008 NNLO)
        !* Calculated Born CS with only t-channel graphs (*),
        !* (a) with only ZZ-fusion graphs and
        !* (b) with all WW/ZZ fusion graphs, and calculated
        !* ratio (a)/(b).
        !* (*): s-channel graphs are neglected because they are
        !*      highly suppressed after VBF cuts.
        !******************************************************
        !* x : Higgs mass in GeV
        !* fit: strictly valid in range [90:1000], deviations from data below 0.2%
        !*      slight extrapolation to range [70:1020] is allowed.
        !* 27/4/2011, Oliver Brein (revised: 29/4/2011)
        !******************************
        implicit none
        real*8 a, b, c, d, e, f, g, h, i, j, k, x
        real*8 a0, b0, c0, d0

        a0 = 0.237215317453696
        b0 = 0.000102401893913729
        c0 = -8.00259372973129e-08
        d0 = -1.1874470382624e-11

        i = 5170.13979049644
        j = 169.605166589112
        k = -7.27222629415314
        a = 0.32245038001123
        b = -0.000343674727914338
        c = 9.99626245300363e-07
        d = -7.7551434795945e-10
        e = -1.73547041303278e-12
        f = 4.41364218643659e-15
        g = -3.67002024110552e-18
        h = 1.08578419805009e-21

        lhc7_rHVBF_ZZ = 0d0

        if (x .lt. 70d0) then
            x = 70d0
            lhc7_rHVBF_ZZ = a0 + b0*x + c0*x**2 + d0*x**3
        endif
        if ((x .ge. 70d0) .and. (x .lt. 90d0)) then
            lhc7_rHVBF_ZZ = a0 + b0*x + c0*x**2 + d0*x**3
        endif
        if ((x .ge. 90d0) .and. (x .le. 1020d0)) then
            lhc7_rHVBF_ZZ = i/x**3 + j/x**2 + k/x + a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6 + h*x**7
        endif
    end function

    real*8 function lhc8_rHVBF_WW(x)
        !******************************************************
        !* Used HAWK-2.0.1.
        !* (PDF: NNPDF23_nlo_as_0118_qed)
        !* Calculated Born CS
        !* (a) with only SM HWW coupling,
        !* (b) with only SM HZZ coupling,
        !* (c) with both SM HWW and HZZ coupling.
        !* Calculated ratios (a)/((a)+(b)) and (b)/((a)+(b)), and
        !* performed simultaneous fit, such that the sum is 1.
        !* Cross checked that (a)+(b) ~= (c), therefore interference
        !* is negligible.
        !******************************************************
        !* x : Higgs mass in GeV
        !* fit: strictly valid in range [10:1000].
        !* 01/06/2016, TS
        !******************************
        implicit none
        real*8 am1, a0, a1, a2, x

        am1 = 9.99844541408024e-06
        a0 = 0.250107980787387
        a1 = 6.52202639729497e-05
        a2 = -3.20376357582592e-08

        lhc8_rHVBF_WW = 0d0

        if (x .ge. 0d0) then
            lhc8_rHVBF_WW = 1.0d0 - (am1/x + a0 + a1*x + a2*x**2.0d0)
        endif

        if (x .gt. 1050d0) then
            write (*, *) 'function lhc8_rHVBF_WW might not be a good fit (m_H > 1050 GeV)'
        endif

    end function

    real*8 function lhc8_rHVBF_ZZ(x)
        !******************************************************
        !* Used HAWK-2.0.1.
        !* (PDF: NNPDF23_nlo_as_0118_qed)
        !* Calculated Born CS
        !* (a) with only SM HWW coupling,
        !* (b) with only SM HZZ coupling,
        !* (c) with both SM HWW and HZZ coupling.
        !* Calculated ratios (a)/((a)+(b)) and (b)/((a)+(b)), and
        !* performed simultaneous fit, such that the sum is 1.
        !* Cross checked that (a)+(b) ~= (c), therefore interference
        !* is negligible.
        !******************************************************
        !* x : Higgs mass in GeV
        !* fit: strictly valid in range [10:1000].
        !* 01/06/2016, TS
        !******************************
        implicit none
        real*8 am1, a0, a1, a2, x

        am1 = 9.99844541408024e-06
        a0 = 0.250107980787387
        a1 = 6.52202639729497e-05
        a2 = -3.20376357582592e-08

        lhc8_rHVBF_ZZ = 0d0

        if (x .ge. 0d0) then
            lhc8_rHVBF_ZZ = am1/x + a0 + a1*x + a2*x**2.0d0
        endif

        if (x .gt. 1050d0) then
            write (*, *) 'function lhc8_rHVBF_ZZ might not be a good fit (m_H > 1050 GeV)'
        endif

    end function

    real*8 function lhc13_rHVBF_WW(x)
        !******************************************************
        !* Used HAWK-2.0.1.
        !* (PDF: NNPDF23_nlo_as_0118_qed)
        !* Calculated Born CS
        !* (a) with only SM HWW coupling,
        !* (b) with only SM HZZ coupling,
        !* (c) with both SM HWW and HZZ coupling.
        !* Calculated ratios (a)/((a)+(b)) and (b)/((a)+(b)), and
        !* performed simultaneous fit, such that the sum is 1.
        !* Cross checked that (a)+(b) ~= (c), therefore interference
        !* is negligible.
        !******************************************************
        !* x : Higgs mass in GeV
        !* fit: strictly valid in range [10:3000].
        !* 01/06/2016, TS
        !******************************
        implicit none
        real*8 a0, a1, a2, x

        a0 = 0.263427271363775
        a1 = 2.72728035702793e-05
        a2 = -3.89549044847535e-09

        lhc13_rHVBF_WW = 0d0

        if (x .ge. 0d0) then
            lhc13_rHVBF_WW = 1.0d0 - (a0 + a1*x + a2*x**2.0d0)
        endif

        if (x .gt. 3050d0) then
            write (*, *) 'function lhc13_rHVBF_WW might not be a good fit (m_H > 3050 GeV)'
        endif

    end function

    real*8 function lhc13_rHVBF_ZZ(x)
        !******************************************************
        !* Used HAWK-2.0.1.
        !* (PDF: NNPDF23_nlo_as_0118_qed)
        !* Calculated Born CS
        !* (a) with only SM HWW coupling,
        !* (b) with only SM HZZ coupling,
        !* (c) with both SM HWW and HZZ coupling.
        !* Calculated ratios (a)/((a)+(b)) and (b)/((a)+(b)), and
        !* performed simultaneous fit, such that the sum is 1.
        !* Cross checked that (a)+(b) ~= (c), therefore interference
        !* is negligible.
        !******************************************************
        !* x : Higgs mass in GeV
        !* fit: strictly valid in range [10:3000].
        !* 01/06/2016, TS
        !******************************
        implicit none
        real*8 a0, a1, a2, x

        a0 = 0.263427271363775
        a1 = 2.72728035702793e-05
        a2 = -3.89549044847535e-09

        lhc13_rHVBF_ZZ = 0d0

        if (x .ge. 0d0) then
            lhc13_rHVBF_ZZ = a0 + a1*x + a2*x**2.0d0
        endif

        if (x .gt. 3050d0) then
            write (*, *) 'function lhc13_rHVBF_ZZ might not be a good fit (m_H > 3050 GeV)'
        endif

    end function

end module theory_colliderSfunctions
