
function test_ideal_stability(Bp, Bt, dpdr, R0, r0, rs, rb; m1ncap = 10, m0ncap=3, rs0=1e-5, nmax=8, del=1e-5, integrator_reltol=10^(-20), verbose=true, verify_sols=false, case=0, kwargs...)
    q = q_(Bt,Bp,R0)

    if !test_Suydam(Bt, q, dpdr, rb) 
        verbose && print("Equilibrium ideal MHD unstable, failed Suydam criteria.\n")
        return false
    end

    m1set = -m1ncap:m1ncap
    m0set = 0:m0ncap

    m=1
    n=0

    for n in m1set
        k = k_(n,R0)
        f__= r -> f_(k, m, Bt, Bp)(r)
        g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
        ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=plot_solution)

        if !ideal_stability_solin(solin, del, r0, rs) 
            verbose && print("Equilibrium ideal MHD unstable, failed Newcomb criteria for mode n=$(n),m=$(m).\n")
            return false
        end

        if !ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol, plot_solution=false)
            verbose && print("Equilibrium ideal MHD unstable, failed Newcomb criteria for mode n=$(n),m=$(m).\n")
            return false
        end
    end

    m=0
    for n in m0set
        k = k_(n,R0)
        f__= r -> f_(k, m, Bt, Bp)(r)
        g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
        ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=plot_solution)

        if !ideal_stability_solin(solin, del, r0, rs) 
            verbose && print("Equilibrium ideal MHD unstable, failed Newcomb criteria for mode n=$(n),m=$(m).\n")
            return false
        end

        if !ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol, plot_solution=false)
            verbose && print("Equilibrium ideal MHD unstable, failed Newcomb criteria for mode n=$(n),m=$(m).\n")
            return false
        end
    end

    return true
end

function test_ideal_stability(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del, q; integrator_reltol=10^(-20), verify_sols=false, case=0, kwargs...)
    if !verify_sols && case == 0
        if !test_Suydam(Bt, q, dpdr, rb) 
            return false
        end

        f__= r -> f_(k, m, Bt, Bp)(r)
        g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
        ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=false)

        if !ideal_stability_solin(solin, del, r0, rs) 
            return false
        end

        if !ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol, plot_solution=false)
            return false
        end

    else verify_sols
        if (!test_Suydam(Bt, q, dpdr, rb) && !(case>1)) || case == 1 
            return test_Suydam(Bt, q, dpdr, rb;plotresults=true)
        end
        
        f__= r -> f_(k, m, Bt, Bp)(r)
        g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
        ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=false)

        if (!ideal_stability_solin(solin, del, r0, rs) && !(case>1)) || case==2
            Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=true)
            return ideal_stability_solin(solin, del, r0, rs)
        end

        if !ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol, plot_solution=false) || case==3
            ξ_plus_solout, psi_plus_solout = get_ξ_plus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol)
            ξ_minus_solout, psi_minus_solout = get_ξ_minus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_minus; integrator_reltol=integrator_reltol)
            Δl,Δr,del2 = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del)
    
            ξ_solout(r) = ξ_minus_solout(r)+Δr*ξ_plus_solout(r)
            psi_solout(r) = ξ_solout(r)*F(k, m, Bt, Bp)(r)
    
            display(plot([psi_plus_solout.t,psi_plus_solout.t,psi_plus_solout.t],[[psi_solout(i) for i in psi_plus_solout.t],
            [psi_plus_solout(i)[2] for i in psi_plus_solout.t],
            [psi_minus_solout(i)[2] for i in psi_plus_solout.t]],labels=["Ψ, where Ψ(rb)/Ψ(rs+del)=$(round(psi_solout(rb)/psi_solout(rs+del);sigdigits=3))" "Ψ small" "Ψ big"]))

            return ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol, plot_solution=false)
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

function ideal_stability_solin(solin, del, r0, rs)
    fsolin(r) = solin(r)[2]

    if length(find_zeros(fsolin,r0,rs-del))>0
        return false
    end

    return true
end

function ideal_stability_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=1e-20, plot_solution=false)
    ξ_plus_solout, psi_plus_solout = get_ξ_plus_solout(Bp, Bt, dpdr, k, m, rs, rb, del, ξ_plus; integrator_reltol=integrator_reltol)

    plot_solution && plot(psi_plus_solout.t,[i[2] for i in psi_plus_solout.u])

    if length(find_zeros(ξ_plus_solout,rs+del,rb))>0
        return false
    end

    return true
end

function ideal_nowall_external_stability(Bp, Bt, dpdr, k, m, rs, rb, del, solout;rb_diff=1e-10)
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


#TESTING

