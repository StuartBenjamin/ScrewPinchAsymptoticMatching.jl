
function test_ideal_stability(Bp, Bt, dpdr, R0, r0, rb; m1ncap = 7, m0ncap=3, rs0=1e-5, nmax=8, del=1e-5, integrator_reltol=10^(-20), integrator_reltol_no_rs=1e-5, verbose=true, rs_verbose=false, verify_sols=false, case=0, ignore_Suydam=false, return_psi_small=false, kwargs...)
    length(find_zeros(Bt,0.0,rb))>0 && error("MHD analysis not set up for Bt switching sign.")
    
    q = q_(Bt,Bp,R0)

    if !test_Suydam(Bt, q, dpdr, rb) && !ignore_Suydam
        verbose && print("Equilibrium ideal MHD unstable, failed Suydam criteria.\n")
        return false
    end

    m1set = -m1ncap:m1ncap
    m0set = 1:m0ncap

    m=1
    n=0

    for n in m1set
        rs = find_rs(q,m,n,rb;verbose=rs_verbose)
        rs==-2 && @warn("Multiple rational surfaces detected.")
        if rs==-1 
            verbose && print("rs outside q range. ")
            k = k_(n,R0)
            solin = integrate_psi_in(Bp, Bt, dpdr, k, m, r0, rb; integrator_reltol=integrator_reltol_no_rs, g=g_)
            solout = integrate_psi_out(Bp, Bt, dpdr, k, m, rb, r0; integrator_reltol=integrator_reltol_no_rs, g=g_)
    
            if (!ideal_stability_solin(solin, 0.0, r0, rb) || !ideal_stability_solin(solout, 0.0, r0, rb*(1-1e-10)))
                verbose && print("Equilibrium unstable to ideal mode n=$(n),m=$(m).\n")
                verbose && display(display(plot(plot(solin.t,[i[2] for i in solin.u]),plot(solout.t,[i[2] for i in solout.u]))))
                return false
            end
        else
            k = k_(n,R0)
            f__= r -> f_(k, m, Bt, Bp)(r)
            g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
            ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
            solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=false)

            if !ideal_stability_solin(solin, del, r0, rs) 
                verbose && print("Equilibrium unstable to ideal mode n=$(n),m=$(m).\n")
                verbose && Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=true)
                return false
            end

            if !ideal_stability_out(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=false, return_psi_small=return_psi_small)
                verbose && print("Equilibrium unstable to ideal mode n=$(n),m=$(m).\n")
                verbose && ideal_stability_out(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=true)
                return false
            end
        end
        verbose && print("Equilibrium stable to ideal mode n=$(n),m=$(m).\n")
    end

    m=0
    for n in m0set
        k = k_(n,R0)
        solin = integrate_psi_in(Bp, Bt, dpdr, k, m, r0, rb; integrator_reltol=integrator_reltol_no_rs, g=g_)
        solout = integrate_psi_out(Bp, Bt, dpdr, k, m, rb, r0; integrator_reltol=integrator_reltol_no_rs, g=g_)

        if (!ideal_stability_solin(solin, 0.0, r0, rb) || !ideal_stability_solin(solout, 0.0, r0, rb*(1-1e-10)))
            verbose && print("Equilibrium unstable to ideal mode n=$(n),m=$(m).\n")
            verbose && display(display(plot(plot(solin.t,[i[2] for i in solin.u]),plot(solout.t,[i[2] for i in solout.u]))))
            return false
        end

        verbose && print("Equilibrium stable to ideal mode n=$(n),m=$(m).\n")
    end

    return true
end

function test_ideal_stability(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del, q; integrator_reltol=10^(-20), verify_sols=false, case=0, return_psi_small=false, kwargs...)
    length(find_zeros(Bt,0.0,rb))>0 && error("MHD analysis not set up for Bt switching sign.")

    if !verify_sols && case == 0
        if !test_Suydam(Bt, q, dpdr, rb) 
            return false
        end

        if m!=0
            f__= r -> f_(k, m, Bt, Bp)(r)
            g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
            ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
            solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=false)

            if !ideal_stability_solin(solin, del, r0, rs) 
                return false
            end

            if !ideal_stability_out(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=false, return_psi_small=return_psi_small)
                return false
            end
        else
            k = k_(n,R0)
            solin = integrate_psi_in(Bp, Bt, dpdr, k, m, r0, rb; integrator_reltol=integrator_reltol, g=g_)
            solout = integrate_psi_out(Bp, Bt, dpdr, k, m, rb, r0; integrator_reltol=integrator_reltol, g=g_)

            if (!ideal_stability_solin(solin, 0.0, r0, rb) || !ideal_stability_solin(solout, 0.0, r0, rb*(1-1e-10)))
                return false
            end
        end
    else verify_sols
        if (!test_Suydam(Bt, q, dpdr, rb) && !(case>1)) || case == 1 
            return test_Suydam(Bt, q, dpdr, rb;plotresults=true)
        end
        
        if m!=0
            f__= r -> f_(k, m, Bt, Bp)(r)
            g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
            ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
            solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=false)

            if (!ideal_stability_solin(solin, del, r0, rs) && !(case>2)) || case==2
                Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=true)
                return ideal_stability_solin(solin, del, r0, rs)
            end

            if case==3 || !ideal_stability_out(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=false)
                return ideal_stability_out(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=true, return_psi_small=return_psi_small)
            end
        else
            k = k_(n,R0)
            solin = integrate_psi_in(Bp, Bt, dpdr, k, m, r0, rb; integrator_reltol=integrator_reltol, g=g_)
            solout = integrate_psi_out(Bp, Bt, dpdr, k, m, rb, r0; integrator_reltol=integrator_reltol, g=g_)

            if ((!ideal_stability_solin(solin, 0.0, r0, rb) || !ideal_stability_solin(solout, 0.0, r0, rb*(1-1e-10))) && !(case>2)) || case==2
                display(plot(plot(solin.t,[i[2] for i in solin.u]),plot(solout.t,[i[2] for i in solout.u])))
                return ideal_stability_solin(solin, 0.0, r0, rb)
            end
        end
    end

    return true
end

function test_Suydam(Bt::Function, q::Function, dpdr::Function, rb; plotresults=false, kwargs...)
    plotresults && return test_Suydam_(Bt, q, dpdr, rb; plotresults=true, kwargs...)

    qprime(r)=ForwardDiff.derivative(q,r)
    local_Suyd(r) = r*Bt(r)^2*(qprime(r)/q(r))^2+8*dpdr(r)

    if length(find_zeros(local_Suyd,0.0,rb))>0
        return false
    end
    return true
end

function test_Suydam_(Bt::Function, q::Function, dpdr::Function, rb; plotresults=true, rvec = range(1e-5,rb,200),ylims=(-1,5))
    qprime(r)=ForwardDiff.derivative(q,r)
    local_Suyd(r) = r*Bt(r)^2*(qprime(r)/q(r))^2+8*dpdr(r)

    Ds(r) = -2*dpdr(r)*q(r)^2/(r*Bt(r)^2*qprime(r)^2)
    local_Suyd2(r) = (1/4)*(1-4*Ds(r))

    plo=plot(rvec,[local_Suyd.(rvec),local_Suyd2.(rvec)],ylims=ylims, label=["bounded Suyd" "Reference Suyd"])
    display(plo)

    if length(find_zeros(local_Suyd,0.0,rb))>0
        return false,plo,local_Suyd,local_Suyd2
    end    
    return true,plo,local_Suyd,local_Suyd2
end

function ideal_stability_solin(solin, del, r0, rs)
    fsolin(r) = solin(r)[2]

    if length(find_zeros(fsolin,r0,rs-del))>0
        return false
    end

    return true
end

function ideal_stability_out(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=1e-20, plot_solution=false, return_psi_small=false)
    del, ξ_plus, ξ_minus, ξ_plus_prime, ξ_minus_prime, solin, solout, inScale, outScale, Δl, Δr, ξl, ξr, cl, cr, psi_l, psi_r = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; bigoutput=true, integrator_reltol=integrator_reltol, plot_solution=false, plot_soln=false)
    return ideal_stability_out(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus, ξ_minus, ξ_plus_prime, ξ_minus_prime, solin, solout, inScale, outScale, Δl, Δr, ξl, ξr, cl, cr, psi_l, psi_r; integrator_reltol=integrator_reltol, plot_solution=plot_solution, return_psi_small=return_psi_small)
end

function ideal_stability_out(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus, ξ_minus, ξ_plus_prime, ξ_minus_prime, solin, solout, inScale, outScale, Δl, Δr, ξl, ξr, cl, cr, psi_l, psi_r; integrator_reltol=1e-20, plot_solution=false, return_psi_small=false)
    psi_r_minus(r) = cr*ξ_minus(r)*F(k, m, Bt, Bp)(r)
    psi_r_plus(r) = cr*Δr*ξ_plus(r)*F(k, m, Bt, Bp)(r)

    psi_r_minus_deriv(r) = ForwardDiff.derivative(psi_r_minus,r)
    psi_r_plus_deriv(r) = ForwardDiff.derivative(psi_r_plus,r)

    #sum of the two above should equal psi_r(r), and solout(i)[2]*outScale
    psi_ode_problem_plus = SecondOrderODEProblem(make_psi_ODE_RUTH(Bp, Bt, dpdr, k, m; g=g_),[psi_r_plus_deriv(rs+del)],[psi_r_plus(rs+del)],(rs+del,rb))
    psi_plus_solout = solve(psi_ode_problem_plus,Tsit5(),reltol=integrator_reltol)
    
    psi_ode_problem_minus = SecondOrderODEProblem(make_psi_ODE_RUTH(Bp, Bt, dpdr, k, m; g=g_),[psi_r_minus_deriv(rs+del)],[psi_r_minus(rs+del)],(rs+del,rb))
    psi_minus_solout = solve(psi_ode_problem_minus,Tsit5(),reltol=integrator_reltol)

    err = abs(psi_plus_solout(rb)[2]+psi_minus_solout(rb)[2]-solout(rb)[2]*outScale)/solout(rs+del)[2]*outScale

    if plot_solution
        psiplot = plot([psi_plus_solout.t,psi_plus_solout.t],[[psi_plus_solout(i)[2]+psi_minus_solout(i)[2] for i in psi_plus_solout.t],[solout(i)[2]*outScale for i in psi_plus_solout.t]],label=["Ψ shoot out" "Ψ shooting in"])
        plusplot =  plot([psi_plus_solout.t],[[psi_plus_solout(i)[2] for i in psi_plus_solout.t]],label="Ψ small (zero-crossing = unstable)")
        display(plot(psiplot,plusplot))

        print("Ψ(rb)/Ψ(rs+del)=$(round(err;sigdigits=3)).\n")
        #print("Ψsmall(rs+del)=$(round(psi_plus_solout(rs+del)[2];sigdigits=3)).\n")
        #print("Ψsmall(rb)=$(round(psi_plus_solout(rb)[2];sigdigits=3)).\n")
    end
    if return_psi_small
        return psi_plus_solout
    end

    if length(find_zeros(r -> psi_plus_solout(r)[2],rs+del,rb)) > 0 #No need to divide by F (and use ξ_plus_solout) since F won't cross zero over this interval, and is just a scale factor
        return false
    end

    return true 
end

function ideal_nowall_external_stability(Bp, Bt, dpdr, k, m, rs, rb, del, solout;rb_diff=1e-10) #Zero since we enforce ξ_out(rb) = 0 
    k0=k0_(k,m)
    
    ξ_out(r)=solout(r)/F(k, m, Bt, Bp)(r)
    ξ_out_prime(r)=ForwardDiff.derivative(ξ_out,r)

    #Wall terms
    #K(r)=SpecialFunctions.besselk(m, ξ_out(r))
    #Kprime(r)=ForwardDiff.derivative(K,r)

    no_wall(r)=(F(k, m, Bt, Bp)(r)^2/(k0^2))*r*ξ_out_prime(r)/ξ_out(r)+F(k, m, Bt, Bp)(r)*Fdag(k, m, Bt, Bp)(r)/(k0^2)

    if (no_wall(r)*ξ_out(rb)^2)<0.0
        return false
    end

    return true
end

if false #deprecated functions
    function get_ξ_plus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=10^(-20))
        psi_plus(r) =  ξ_plus(r)*F(k, m, Bt, Bp)(r)
        psi_plus_prime(r) = ForwardDiff.derivative(psi_plus,r)

        psi_ode_problem = SecondOrderODEProblem(make_psi_ODE_RUTH(Bp, Bt, dpdr, k, m; g=g_),[psi_plus_prime(rs+del)],[psi_plus(rs+del)],(rs+del,rb))
        psi_plus_solout = solve(psi_ode_problem,Tsit5(),reltol=integrator_reltol)

        ξ_plus_solout(r) = psi_plus_solout(r)[2]/F(k, m, Bt, Bp)(r)

        return ξ_plus_solout, psi_plus_solout
    end

    function get_ξ_minus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_minus; integrator_reltol=10^(-20))
        psi_minus(r) =  ξ_minus(r)*F(k, m, Bt, Bp)(r)
        psi_minus_prime(r) = ForwardDiff.derivative(psi_minus,r)

        psi_ode_problem = SecondOrderODEProblem(make_psi_ODE_RUTH(Bp, Bt, dpdr, k, m; g=g_),[psi_minus_prime(rs+del)],[psi_minus(rs+del)],(rs+del,rb))
        psi_minus_solout = solve(psi_ode_problem,Tsit5(),reltol=integrator_reltol)

        ξ_minus_solout(r) = psi_minus_solout(r)[2]/F(k, m, Bt, Bp)(r)

        return ξ_minus_solout, psi_minus_solout
    end

    function ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=1e-20, plot_solution=false, p5=nothing)
        ξ_plus_solout, psi_plus_solout = get_ξ_plus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol)

        if plot_solution
            if p5 isa Nothing   
                plot(psi_plus_solout.t,[i[2] for i in psi_plus_solout.u])
            else
                p6=plot(psi_plus_solout.t,[i[2] for i in psi_plus_solout.u])
                display(plot(p5,p6))
            end
        end

        if length(find_zeros(ξ_plus_solout,rs+del,rb))>0
            return false
        end

        return true
    end

    function plot_ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus, ξ_minus, outScale, solout;integrator_reltol=1e-20)
        ξ_plus_solout, psi_plus_solout = get_ξ_plus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol)
        ξ_minus_solout, psi_minus_solout = get_ξ_minus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_minus; integrator_reltol=integrator_reltol)
        Δl,Δr,del2 = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del)

        ξ_solout(r) = ξ_minus_solout(r)+Δr*ξ_plus_solout(r)
        psi_solout_raw(r) = ξ_solout(r)*F(k, m, Bt, Bp)(r)
        psi_solout(r) = outScale*psi_solout_raw(r)*(solout(rs+del)[2]/psi_solout_raw(rs+del))

        #p5 = plot([psi_plus_solout.t,psi_plus_solout.t,psi_plus_solout.t],[[psi_solout(i) for i in psi_plus_solout.t],
        #[psi_plus_solout(i)[2] for i in psi_plus_solout.t],
        #[psi_minus_solout(i)[2] for i in psi_plus_solout.t]],labels=["Ψ, where Ψ(rb)/Ψ(rs+del)=$(round(psi_solout(rb)/psi_solout(rs+del);sigdigits=3))" "Ψ small" "Ψ big"])

        ppsi = plot([psi_plus_solout.t,psi_plus_solout.t],[[psi_solout(i) for i in psi_plus_solout.t],[outScale*solout(i)[2] for i in psi_plus_solout.t]],label=["Ψ, Ψ(rb)/Ψ(rs+del)=$(round(psi_solout(rb)/psi_solout(rs+del);sigdigits=3))" false],legend=:topright)

        psmall = plot(psi_plus_solout.t,[psi_plus_solout(i)[2] for i in psi_plus_solout.t],label="Ψ small (zero-crossing = unstable)")

        pbig = plot(psi_plus_solout.t,[psi_minus_solout(i)[2] for i in psi_plus_solout.t],label="Ψ big")

        p5 = plot(ppsi,psmall,pbig)
        display(plot(ppsi,psmall,pbig))

        return p5
    end
end
