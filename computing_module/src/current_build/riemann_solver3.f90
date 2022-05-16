module riemann_solver3
    use kind_parameters
contains
    subroutine guess_p(PM, p1, u1, rho1, p4, u4, rho4, gamma)
        real(dkind), intent (out) :: PM
        real(dkind), intent (in) :: p1, rho1, u1, p4, rho4, u4, gamma
        real(dkind), allocatable :: tol, quser, cup, ppv, pmin, pmax, qmax
        real(dkind), allocatable :: cspd1, cspd4, G1, G2, G3, G4, G5, G6, G7, G8

        tol = 1e-6;
        cspd1 = sqrt(gamma*p1/rho1)
        cspd4 = sqrt(gamma*p4/rho4)
        G1 = 0.5e0*(gamma - 1.e0)/gamma
        G2 = 0.5e0*(gamma + 1.e0)/gamma
        G3 = 1.e0/G1
        G4 = 1.e0/(G1*gamma)
        G5 = 1.e0/(G2*gamma)
        G6 = G1/G2
        G7 = G1*gamma
        G8 = gamma - 1.e0

		PM = -1e20
        ! Define user-chosen pressure ratio, below which we consider
		! use of the pressure estimate from the Primitive Variable Riemann Solver
		quser = 2.e0

        ! 1) Obtain initial guess using the primitive variable Riemann solver
		cup  = 0.25e0*(rho1 + rho4)*(cspd1 + cspd4)
        ppv  = 0.5e0*(p1 + p4) + 0.5e0*(u1 - u4)*cup
        ppv  = max(0.e0, ppv)
        pmin = min(p1, p4)
        pmax = max(p1, p4)
        qmax = pmax/pmin

        ! 2) Decide whether to use the PVRS value, the Two-Rarefaction value,
		!    or the Two-Shock value
		if( (qmax <= quser).and.(pmin <= ppv).and.(ppv <= pmax) ) then
            ! Use the PVRS pressure value
			PM = ppv
        else
            if(ppv < pmin) then
                ! Use the Two-Rarefaction Riemann solver
				pq = (p1/p4) ** G1
                um = (pq*u1/cspd1 + u4/cspd4 + G4*(pq - 1.e0))/(pq/cspd1 + 1.e0/cspd4)
                ptL = 1.e0 + G7*(u1 - um)/cspd1
                ptR = 1.e0 + G7*(um - u4)/cspd4
                PM  = 0.5e0*(p1*(ptL ** G3) + p4*(ptR ** G3))
            else
                ! Use Two-Shock Riemann solver with PVRS as estimate
                ge1 = sqrt((G5/rho1)/(G6*p1 + ppv))
                ge4 = sqrt((G5/rho4)/(G6*p4 + ppv))
                PM  = (ge1*p1 + ge4*p4 - u4 + u1)/(ge1 + ge4)
            end if
        end if
    ! 3) Return the pressure estimate in the star region
    end subroutine guess_p

    subroutine prefun(F, dF, P, densK, presK, cspdK, gamma)
        real(dkind), intent (out) :: F, dF
        real(dkind), intent (in) :: P, densK, presK, cspdK, gamma
        real(dkind), allocatable :: prat, aK, bK, G1, G2, G4, G5, G6

        G1 = 0.5e0*(gamma - 1.e0)/gamma
        G2 = 0.5e0*(gamma + 1.e0)/gamma
        G4 = 1.e0/(G1*gamma)
        G5 = 1.e0/(G2*gamma)
        G6 = G1/G2
        ! Outputs
		F  = -1e20
        dF = -1e20

        ! Rarefaction wave
        if(P <= presK) then
            prat = P/presK
            F    = G4*cspdK*(prat ** G1 - 1.e0)
            dF   = (1.e0/(densK*cspdK)) * prat ** (-G2)

            ! Shock wave
		else
            aK  = G5/densK
            bK  = G6*presK
            qrt = sqrt(aK/(bK + P))
            F   = (P - presK)*qrt
            dF  = (1.e0 - 0.5e0*(P - presK)/(bK + P))*qrt
    end if
    end subroutine prefun

    subroutine solve_riemann_problem3(RP, RU, RRHO, p1, u1, rho1, p4, u4, rho4, t0, gamma)
        real(dkind), intent (out) :: RP, RU, RRHO
        real(dkind), intent (in) :: p1, rho1, u1, p4, rho4, u4, t0, gamma
        real(dkind) :: TOL, G0, G1, G2, G3, G4, G5, G6, G7, G8, cspd1, cspd4
        real(dkind) :: pstar, pstar0, ps, ustar, us, MAX_ITER, TOL_PRES
        real(dkind) :: dOut, uOut, pOut, S, fL, fR, dfL, dfR, pold, cmr
        cspd1 = sqrt(gamma*p1/rho1)
        cspd4 = sqrt(gamma*p4/rho4)
        G1 = 0.5e0*(gamma - 1.e0)/gamma
        G2 = 0.5e0*(gamma + 1.e0)/gamma
        G3 = 1.e0/G1
        G4 = 1.e0/(G1*gamma)
        G5 = 1.e0/(G2*gamma)
        G6 = G1/G2
        G7 = G1*gamma
        G8 = gamma - 1.e0
        ! call starpu -> samplePt -> Rp, RU, RHO
        pstar = -1e20
        ustar = -1e20
        ! Define solver criteria
        MAX_ITER = 100
        TOL_PRES = 1e-8
        ! 1) Compute initial pressure guess and prepare for Newton-Raphson
        call guess_p(pstar0, p1, u1, rho1, p4, u4, rho4, gamma)
        pold   = pstar0
        pcur   = pold
        deltau = u4 - u1
        do it = 1, MAX_ITER
        ! 2) Perform Newton-Raphson iteration to determine pressure in
        !    the star region
            call prefun(fL, dfL, pold, rho1, p1, cspd1, gamma)
            call prefun(fR, dfR, pold, rho4, p4, cspd4, gamma)
            pcur = pold - (fL + fR + deltau)/(dfL + dfR)
            dp = 2.e0*abs((pcur - pold)/(pcur + pold))
            if(dp <= TOL_PRES) then
                exit
            end if

            if(pcur < 0.e0) then
                pcur = TOL_PRES
            end if
            pold = pcur
        end do
        if(dp > TOL_PRES) then
            print *, 'Newton-Raphson did not converge'
        end if
        pstar = pcur
        ustar = 0.5e0*(u1 + u4 + fR - fL)
        pS = pstar
        uS = ustar



        dOut = -1e20
        uOut = -1e20
        pOut = -1e20
        S = 0
        ! Left of contact discontinuity
        if(S <= uS) then
            ! Left rarefaction
            if(pS <= p1) then
                shl = u1 - cspd1
                ! Left initial state
                if(S <= shl) then
                    dOut = rho1
                    pOut = p1
                    uOut = u1
                else
                    cml = cspd1*((pS/p1) ** G1)
                    stl = uS - cml
                    ! Left star state
                    if(S > stl) then
                        dOut = rho1*((pS/p1) ** (1.0/gamma))
                        uOut = uS
                        pOut = pS
                        ! Left fan
                    else
                        uOut = G5*(cspd1 + G7*u1 + S)
                        c    = G5*(cspd1 + G7*(u1 - S))
                        dOut = rho1*((c/cspd1) ** G4)
                        pOut = p1*((c/cspd1) ** G3)
                    end if
                end if
                ! Left shock
            else
                pml = pS/p1
                sl  = u1 - cspd1*sqrt(G2*pml + G1)
                ! Left initial state
                if(S <= sl) then
                    dOut = rho1
                    pOut = p1
                    uOut = u1
                    ! Left star state
                else
                    dOut = rho1*(pml + G6)/(G6*pml + 1.e0)
                    pOut = pS
                    uOut = uS
                end if
            end if
            ! Right of contact discontinuity
        else
            ! Right shock
            if(pS > p4) then
                pmr = pS/p4
                sr  = u4 + cspd4*sqrt(G2*pmr + G1)
                ! Right initial state
                if(S >= sr) then
                    dOut = rho4
                    pOut = p4
                    uOut = u4
                    ! Right star state
                else
                    dOut = rho4*(pmr + G6)/(G6*pmr + 1.e0)
                    pOut = pS
                    uOut = uS
                end if
                ! Right rarefaction
            else
                shr = u4 + cspd4
                ! Right initial state
                if(S >= shr) then
                    dOut = rho4
                    pOut = p4
                    uOut = u4
                else
                    cmr = cspd4*((pS/p4) ** G1)
                    stR = uS + cmr
                    ! Right star state
                    if(S <= stR) then
                        dOut = rho4*((pS/p4) ** (1.e0/gamma))
                        pOut = pS
                        uOut = uS
                        ! Right fan
                    else
                        uOut = G5*(-cspd4 + G7*u4 + S)
                        c    = G5*(cspd4 - G7*(u4 - S))
                        dOut = rho4*((c/cspd4) ** G4)
                        pOut = p4*((c/cspd4) ** G3)
                    end if
                end if
            end if
        end if
        RP = pOut
        RU = uOut
        RRhO = dOut
    end subroutine solve_riemann_problem3
end module riemann_solver3
