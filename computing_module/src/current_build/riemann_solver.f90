module riemann_solver
    use kind_parameters
contains
    subroutine solve_riemann_problem(RP, RU, RRHO, p1, u1, rho1, p4, u4, rho4, t0)
        real(dkind), intent (out) :: RP, RU, RRHO
        real(dkind), intent (in) :: p1, rho1, u1, p4, rho4, u4, t0
        real(dkind), allocatable :: tol, gamma, delta, alpha, zeta, ws, wd, a1, a2, a3, a4, p2, p3, u2, u3, ustar, u
        real(dkind), allocatable :: M1, M4, V1, dp2, dp3, x1, x2, x3, x4, x5, xUS, xUD, A, t, x
        character(len = 3) :: tipo
        integer, allocatable :: iter
        tol = 0.01
        if(abs(p1-p4)<tol.and.abs(u1-u4)<tol.and.abs(rho1-rho4)<tol) then
            RP = p1;RU = u1;RRHO = rho1;
            return
        end if
        t = t0
        x = 0
        dp3 = 0
        dp2 = 0
        iter = 0
        gamma = 1.4
        delta = (gamma - 1) / 2
        alpha = 2 * gamma / (gamma - 1)
        ws = 0
        wd = 0
        a1 = sqrt(gamma * p1 / rho1)
        a4 = sqrt(gamma * p4 / rho4)
        zeta = (p1 / p4)**(1 / alpha) * (a1 / a4)
        if (u4 - u1>a4 / delta + a1 / delta) then
            print *, 'No solution for:'
            print *, 'rho1 =', rho1, 'p1 =', p1, 'u1=', u1, &
                    'rho4 =', rho4, 'p4 =', p4, 'u4=', u4, &
                    'rho =', rrho, 'p =', rp, 'u=', ru
        else
            ustar = (zeta * (a1 + delta * u1) - a4 + delta * u4) / ((1 + zeta) * delta)
            u = ustar
            p2 = 1
            p3 = 0.5

            do while (abs(p2 - p3)>tol)
                ! Checking wave on first family
                if(u>u1) then
                    a2 = a1 + delta * (u1 - u)
                    p2 = p1 * (a2 / a1)**alpha
                    dp2 = -gamma * p2 / a2
                else if(u<u1) then
                    V1 = (gamma + 1) / 4 * (u1 - u) + sqrt(((gamma + 1) / 4)**2 * (u1 - u)**2 + a1**2)
                    M1 = V1 / a1
                    p2 = p1 * (1 + (2 * gamma) / (gamma + 1) * (M1**2 - 1))
                    a2 = a1 * sqrt((gamma + 1 + (gamma - 1) * (p2 / p1)) / (gamma + 1 + (gamma - 1) * (p1 / p2)))
                    ws = u1 - V1
                    dp2 = -2 * gamma * (p1 / a1) * ((abs(M1))**3 / (1 + M1**2))
                end if
                !Checking wave on the second family
                if (u>u4) then !S
                    V4 = (gamma + 1) / 4 * (u - u4) + sqrt(((gamma + 1) / 4)**2 * (u - u4)**2 + a4**2)
                    M4 = V4 / a4
                    p3 = p4 * (1 + (2 * gamma) / (gamma + 1) * (M4**2 - 1))
                    a3 = a4 * sqrt((gamma + 1 + (gamma - 1) * (p3 / p4)) / (gamma + 1 + (gamma - 1) * (p4 / p3)))
                    wd = u4 - V4
                    dp3 = 2 * gamma * (p4 / a4) * ((abs(M4))**3 / (1 + M4**2))
                else if (u<u4) then !%R
                    a3 = a4 + delta * (u - u4)
                    p3 = p4 * (a3 / a4)**alpha
                    dp3 = gamma * p3 / a3
                end if
                !New velocity and update of the iteration
                u = u - (p2 - p3) / (dp2 - dp3)
                iter = iter + 1
                ! Checking number of iterations
                if (iter>45000) then
                    print *, 'The problem can not be solved using this method!'
                    exit
                end if
            end do
        end if
        rho2 = gamma * p2 / a2**2
        rho3 = gamma * p3 / a3**2
        !Type of the solution
        if ((u>u1).and.(abs(p2 - p1)>tol).and.(abs(p2 - p4)>tol)) then
            if (u<u4) then
                tipo = 'RCR'
            else
                tipo = 'RCS'
            end if
        else if ((u<u1).and.(abs(p2 - p1)>tol).and.(abs(p2 - p4)>tol)) then
            if (u<u4) then
                tipo = 'SCR'
            else
                tipo = 'SCS'
            end if
        else if ((abs(u - u1)<tol).and.(abs(p2 - p1)<tol).and.(abs(u - u4)>tol)) then
            if (u<u4) then
                tipo = 'NCR'
                u = u1;p2 = p1;p3 = p1;a2 = a1;rho2 = rho1;
            else
                tipo = 'NCS'
                u = u1;p2 = p1;p3 = p1;a2 = a1;rho2 = rho1;
            end if
        else if ((abs(u - u1)>tol).and.(abs(p2 - p4)<tol).and.(abs(u - u4)<tol)) then
            if (u>u1) then
                tipo = 'RCN'
                u = u4;p2 = p4;p3 = p4;a3 = a4;rho3 = rho4;
            else
                tipo = 'SCN'
                u = u4;p2 = p4;p3 = p4;a3 = a4;rho3 = rho4;
            end if
        else
            tipo = 'NCN'
            u = u1;p2 = p1; p3 = p1;a2 = a1; a3 = a1;rho2 = rho1; rho3 = rho1;
        end if
        x1 = (u1 - a1) * t
        x2 = (u - a2) * t
        x3 = u * t
        x4 = (u + a3) * t
        x5 = (u4 + a4) * t
        xUS = -abs(ws) * t
        xUD = abs(wd) * t
        if (tipo=='RCR') then
            if(x<=x1) then
                RU = u1;RP = p1;A = a1;RRHO = rho1;
            else if((x<=x2).and.(x>=x1)) then
                RU = 1 / (1 + delta) * (a1+delta*u1+x/t)
                A = 1/(1+delta)*(a1)+delta/(delta+1)*(u1)-delta/(1+delta)*x/t
                RP = p1 * (A / a1)**(alpha)
                RRHO = rho1 * (A / a1)**(alpha / gamma);
            else if ((x>=x2).and.(x<=x3)) then
                RU = u;
                A = a2;
                RP = p2;
                RRHO = rho2;
            else if ((x>=x3).and.(x<=x4)) then
                RU = u;A = a3;RP = p3;RRHO = rho3;
            else if ((x>=x4).and.(x<=x5)) then
                RU = 1/(1+delta)*(x/t-a4+delta*u4);
                A = 1/(1+delta)*(a4)-delta/(1+delta)*(u4)+delta/(1+delta)*x/t;
                RP = p4 * (A / a4)**(alpha);
                RRHO = rho4 * (A / a4)**(alpha / gamma);
            else if (x>=x5) then
                RU = u4;RP = p4;A = a4;RRHO = rho4;
            end if
        else if (tipo=='RCS') then
            if (x <= x1) then
                RU = u1;RP = p1;A = a1;RRHO = rho1;
            else if (x<=x2.and.x >=x1) then
                RU = 1 / (1 + delta) * (a1 + delta * u1 + x / t);
                A = 1 / (1 + delta) * (a1) + delta / (delta + 1) * (u1) - delta / (1 + delta) * x / t;
                RP = p1 * (A / a1)**(alpha);
                RRHO = rho1 * (A / a1)**(alpha / gamma);
            else if (x>=x2.and.x <=x3) then
                RU = u;A = a2;RP = p2;RRHO = rho2;
            else if (x>=x3.and.x<=xUD) then
                RU = u;RP = p3;A = a3;RRHO = rho3;
            else if (x>=xUD) then
                RU = u4;RP = p4;A = a4;RRHO = rho4;
            end if
        else if (tipo=='SCR') then
            if (x<=xUS) then
                RU = u1;RP = p1;A = a1;RRHO = rho1;
            else if (x>=xUS.and.x<=x3) then
                RU = u;RP = p2;A = a2;RRHO = rho2;
            else if (x>=x3.and.x<=x4) then
                RU = u;A = a3;RP = p3;RRHO = rho3;
            else if (x>=x4.and.x<=x5) then
                RU = 1 / (1 + delta) * (x / t - a4 + delta * u4);
                A = 1 / (1 + delta) * (a4) - delta / (delta + 1) * (u4) + delta / (1 + delta) * x / t;
                RP = p4 * (A / a4)**(alpha);
                RRHO = rho4 * (A / a4)**(alpha / gamma);
            else if (x>=x5) then
                RU = u4;RP = p4;A = a4;RRHO = rho4;
            end if
        else if (tipo=='SCS') then
            if (x<=xUS) then
                RU = u1;RP = p1;A = a1;RRHO = rho1;
            else if (x>=xUS.and.x<=x3) then
                RU = u;RP = p2;A = a2;RRHO = rho2;
            else if(x>=x3.and.x<=xUD) then
                RU = u;RP = p3;A = a3;RRHO = rho3;
            else if (x>=xUD) then
                RU = u4;RP = p4;A = a4;RRHO = rho4;
            end if
        else if (tipo=='NCR') then
            if (x<=x3) then
                RU = u1;A = a1;RP = p1;RRHO = rho1;
            else if (x>=x3.and.x<=x4) then
                RU = u;A = a3;RP = p3;RRHO = rho3;
            else if (x>=x4.and.x<=x5) then
                RU = 1 / (1 + delta) * (x / t - a4 + delta * u4);
                A = 1 / (1 + delta) * (a4) - delta / (delta + 1) * (u4) + delta / (1 + delta) * x / t;
                RP = p4 * (A / a4)**(alpha);
                RRHO = rho4 * (A / a4)**(alpha / gamma);
            else if (x>=x5) then
                RU = u4;RP = p4;A = a4;RRHO = rho4;
            end if
        else if (tipo=='NCS') then
            if (x<=x3) then
                RU = u1;A = a1;RP = p1;RRHO = rho1;
            else if (x>=x3.and.x<=xUD) then
                RU = u;RP = p3;A = a3;RRHO = rho3;
            else if (x>=xUD) then
                RU = u4;RP = p4;A = a4;RRHO = rho4;
            end if
        else if (tipo=='RCN') then
            if (x <=x1) then
                RU = u1; RP = p1; A = a1;RRHO = rho1;
            else if (x>=x1.and.x<=x2) then
                RU = 1 / (1 + delta) * (a1 + delta * u1 + x / t);
                A = 1 / (1 + delta) * (a1) + delta / (delta + 1) * (u1) - delta / (1 + delta) * x / t;
                RP = p1 * (A / a1)**(alpha);
                RRHO = rho1 * (A / a1)**(alpha / gamma);
            else if (x>=x2.and.x<=x3) then
                RU = u;A = a2;RP = p2;RRHO = rho2;
            else if (x>=x3) then
                RU = u4;A = a4;RP = p4;RRHO = rho4;
            end if
        else if (tipo=='SCN') then
            if (x<=xUS) then
                RU = u1;RP = p1;A = a1;RRHO = rho1;
            else if (x>=xUS.and.x<=x3) then
                RU = u;RP = p2;A = a2;RRHO = rho2;
            else if (x>=x3) then
                RU = u4;A = a4;RP = p4;RRHO = rho4;
            end if
        else if (tipo=='NCN') then
            RU = u1;RP = p1;A = a1;RRHO = rho1;
        end if
    end subroutine solve_riemann_problem
end module riemann_solver
