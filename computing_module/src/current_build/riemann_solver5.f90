module riemann_solver5
    use kind_parameters
contains
    subroutine find_p_star_newtonraphson(pStarOut, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R)
        real(dkind), intent (out) :: pStarOut
        real(dkind), intent (in) :: rho_L, u_L, p_L, rho_R, u_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R
        real(dkind) :: p_star, p_star_next, TOL, tp, tpd
        TOL = 1e-6
        ! First we set the initial guess for p_star using a simple mean-value approximation
        p_star_next = 0.5*(p_L+p_R)
        ! Now use the Newton-Raphson algorithm
        do it  = 1, 100
            p_star = p_star_next
            call total_pressure_function(tp, p_star,rho_L,u_L,p_L,rho_R,u_R,p_R, gamma_L, pinf_L, gamma_R, pinf_R)
            call total_pressure_function_deriv(tpd, p_star,rho_L,p_L,rho_R,p_R, gamma_L, pinf_L, gamma_R, pinf_R)
            p_star_next = p_star - tp / tpd
            p_star_next = max(p_star_next, TOL)
            if (abs(p_star_next - p_star)/(0.5*(p_star+p_star_next)) > TOL) then
                exit
            end if

            if (it == 100) then
                p_star_next = 0.5*(p_L+p_R)
            end if
        end do

        pStarOut = p_star_next
    end subroutine find_p_star_newtonraphson

    subroutine total_pressure_function(pOut, p_star, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R)
        real(dkind), intent (out) :: pOut
        real(dkind), intent (in) :: p_star, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R
        real(dkind) :: ff1, ff2
        call f(ff1, p_star, rho_L, p_L, gamma_L, pinf_L)
        call f(ff2, p_star, rho_R, p_R, gamma_R, pinf_R)
        pOut = ff1 + ff2 + u_R - u_L
    end subroutine total_pressure_function

    subroutine total_pressure_function_deriv(pOut, p_star, rho_L, p_L, rho_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R)
        real(dkind), intent (out) :: pOut
        real(dkind), intent (in) :: p_star, rho_L, p_L, rho_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R
        real(dkind) :: fd1, fd2
        call f_deriv(fd1, p_star, rho_L, p_L, gamma_L, pinf_L)
        call f_deriv(fd2, p_star, rho_R, p_R, gamma_R, pinf_R)
        pOut =	fd1 + fd2
    end subroutine total_pressure_function_deriv

    subroutine f(pOut, p_star, rho, p, gamma, pinf)
        real(dkind), intent (out) :: pOut
        real(dkind), intent (in) :: p_star, rho, p, gamma, pinf
        real(dkind) :: qk, aa
        if (p_star > p) then
            call Q_K(qk, p_star, rho, p, gamma, pinf)
            pOut = (p_star - p)/qk
        else
            call aaa(aa, rho,p,gamma,pinf)
            pOut = (2.0*aa/(gamma-1.0))*((((p_star + pinf)/(p + pinf)) ** ((gamma-1.0)/(2.0*gamma))) - 1.0)
        end if
    end subroutine f

    subroutine f_deriv(pOut, p_star, rho, p, gamma, pinf)
        real(dkind), intent (out) :: pOut
        real(dkind), intent (in) :: p_star, rho, p, gamma, pinf
        real(dkind) :: A, B, aa
        A = 2.0/((gamma+1.0)*rho)
        B = (p+pinf)*(gamma-1.0)/(gamma+1.0)

        if (p_star > p) then
            pOut = sqrt(A/(B+p_star+pinf))*(1.0 - ((p_star-p)/(2.0*(B+p_star+pinf))))
        else
            call aaa(aa, rho,p,gamma,pinf)
            pOut = (1.0/(rho*aa))*(((p_star+pinf)/(p+pinf)) ** (-(gamma+1.0)/(2.0*gamma)))
        end if
    end subroutine f_deriv

    subroutine set_left_rarefaction_fan_state(rhoOut, pOut, uOut, rho_L, p_L, u_L, gamma_L, pinf_L, S)
        real(dkind), intent (out) :: rhoOut, pOut, uOut
        real(dkind), intent (in) :: rho_L, u_L, p_L, gamma_L, pinf_L, S
        real(dkind) :: a_L
        call aaa(a_L, rho_L,p_L,gamma_L,pinf_L)
        rhoOut = rho_L*(((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(u_L - S)) ** (2.0/(gamma_L - 1.0)))
        uOut = (2.0/(gamma_L+1.0))*(a_L + S + ((gamma_L-1.0)/2.0)*u_L)
        pOut = (p_L + pinf_L)*(((2.0/(gamma_L+1.0)) + ((gamma_L-1.0)/(a_L*(gamma_L+1.0)))*(u_L - S)) ** ((2.0*gamma_L)/(gamma_L-1.0))) - pinf_L
    end subroutine set_left_rarefaction_fan_state

    subroutine set_right_rarefaction_fan_state(rhoOut, pOut, uOut, rho_R, p_R, u_R, gamma_R, pinf_R, S)
        real(dkind), intent (out) :: rhoOut, pOut, uOut
        real(dkind), intent (in) :: rho_R, u_R, p_R, gamma_R, pinf_R, S
        real(dkind) :: a_R
        call aaa(a_R, rho_R,p_R,gamma_R,pinf_R)
        rhoOut = rho_R*(((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(u_R - S)) ** (2.0/(gamma_R - 1.0)))
        uOut = (2.0/(gamma_R+1.0))*(- a_R + S + ((gamma_R-1.0)/2.0)*u_R)
        pOut = (p_R + pinf_R)*(((2.0/(gamma_R+1.0)) - ((gamma_R-1.0)/(a_R*(gamma_R+1.0)))*(u_R - S)) ** ((2.0*gamma_R)/(gamma_R-1.0))) - pinf_R
    end subroutine set_right_rarefaction_fan_state

    subroutine Q_K(QKOut, p_star, rho, p, gamma, pinf)
        real(dkind), intent (out) :: QKOut
        real(dkind), intent (in) :: p_star, rho, p, gamma, pinf
        real(dkind) :: A, B
        A = 2.0/((gamma+1.0)*rho)
        B = (p+pinf)*(gamma-1.0)/(gamma+1.0)
        QKOut = sqrt((p_star+pinf+B)/A)
    end subroutine Q_K

    subroutine aaa(aOut, rho, p, gamma, pinf)
        real(dkind), intent (out) :: aOut
        real(dkind), intent (in) :: rho, p, gamma, pinf
        aOut = sqrt(gamma*((p+pinf)/rho))
    end subroutine aaa

    subroutine solve_riemann_problem5(pOut, uOut, rhoOut, p_L, u_L, rho_L, p_R, u_R, rho_R, t0, gamma_L, gamma_R)
        real(dkind), intent (out) :: rhoOut, pOut, uOut
        real(dkind), intent (in) :: p_L, rho_L, u_L, p_R, rho_R, u_R, t0, gamma_L, gamma_R
        real(dkind) :: S, S_L = 0, S_HL = 0, S_TL = 0, S_HR = 0, S_TR = 0, pinf_L, pinf_R
        real(dkind) :: P_STAR, ss1, ss2, rho_star_L, rho_star_R, qk, a_L, a_R, a_star_L, a_star_R
        pinf_L = 0.0
        pinf_R = 0.0
        ! Calculate p_star
        call find_p_star_newtonraphson(P_STAR, rho_L, u_L, p_L, rho_R, u_R, p_R, gamma_L, pinf_L, gamma_R, pinf_R)
        ! Calculate u_star
        call f(ss1, P_STAR,rho_R,p_R,gamma_R,pinf_R)
        call f(ss2, P_STAR,rho_L,p_L,gamma_L,pinf_L)
        S_STAR = 0.5*(u_L+u_R) + 0.5*(ss1 - ss2)
        ! Solution now depends on character of 1st and 3rd waves
        if (P_STAR > p_L) then
            ! Left shock
            rho_star_L = rho_L*((2.0*gamma_L*pinf_L + (gamma_L+1.0)*P_STAR + (gamma_L-1.0)*p_L)/(2.0*(p_L + gamma_L*pinf_L) + (gamma_L-1.0)*P_STAR + (gamma_L-1.0)*p_L))
            call Q_K(qk, P_STAR,rho_L,p_L,gamma_L,pinf_L)
            S_L = u_L - (qk / rho_L)

        else
            ! Left rarefaction
            rho_star_L = rho_L*(((P_STAR + pinf_L)/(p_L + pinf_L)) ** (1.0/gamma_L))

            call aaa(a_L, rho_L, p_L, gamma_L, pinf_L)
            a_star_L = a_L*(((P_STAR + pinf_L)/(p_L + pinf_L)) ** ((gamma_L-1.0)/(2.0*gamma_L)))

            S_HL = u_L - a_L
            S_TL = S_STAR - a_star_L
        end if
        if (P_STAR > p_R) then
            ! Right shock
            rho_star_R = rho_R*((2.0*gamma_R*pinf_R + (gamma_R+1.0)*P_STAR + (gamma_R-1.0)*p_R)/(2.0*(p_R + gamma_R*pinf_R) + (gamma_R-1.0)*P_STAR + (gamma_R-1.0)*p_R))
            call Q_K(qk, P_STAR,rho_R,p_R,gamma_R,pinf_R)
            S_R = u_R + (qk/rho_R)
        else
            ! Right rarefaction
            rho_star_R = rho_R*(((P_STAR + pinf_R)/(p_R + pinf_R)) ** (1.0/gamma_R))
            call aaa(a_R, rho_R,p_R,gamma_R, pinf_R)
            a_star_R = a_R*(((P_STAR + pinf_R)/(p_R + pinf_R)) ** ((gamma_R-1.0)/(2.0*gamma_R)))
            S_HR = u_R + a_R
            S_TR = S_STAR + a_star_R
        end if

        S = 0

        if (S < S_STAR) then
            ! To the left of the contact
            if (P_STAR > p_L) then
                ! Left shock
                if (S < S_L) then
                    pOut = p_L
                    uOut = u_L
                    rhoOut = rho_L
                else
                    rhoOut = rho_star_L
                    uOut = S_STAR
                    pOut = P_STAR
                end if
            else
                ! Left rarefaction
                if (S < S_HL) then
                    pOut = p_L
                    uOut = u_L
                    rhoOut = rho_L
                else
                    if (S > S_TL) then
                        rhoOut = rho_star_L
                        uOut = S_STAR
                        pOut = P_STAR
                    else
                        call set_left_rarefaction_fan_state(rhoOut, pOut, uOut, rho_L, p_L, u_L, gamma_L, pinf_L, S)
                    end if
                end if
            end if
        else
            ! To the right of the contact
            if (P_STAR > p_R) then
                ! Right shock
                if (S > S_R) then
                    pOut = p_R
                    uOut = u_R
                    rhoOut = rho_R
                else
                    rhoOut = rho_star_R
                    uOut = S_STAR
                    pOut = P_STAR
                end if
            else
                ! Right rarefaction
                if (S > S_HR) then
                    pOut = p_R
                    uOut = u_R
                    rhoOut = rho_R
                else
                    if (S < S_TR) then
                        rhoOut = rho_star_R
                        uOut = S_STAR
                        pOut = P_STAR
                    else
                        call set_left_rarefaction_fan_state(rhoOut, pOut, uOut, rho_R, p_R, u_R, gamma_R, pinf_R, S)
                    end if
                end if
            end if
        end if
    end subroutine solve_riemann_problem5
end module riemann_solver5
