using TaylorSeries

#cd("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia")
include("ChandraScrewPinchEquilibrium.jl")
include("OuterIntegrator.jl")

##########################################################################################################################################################
#Ingredient Functions
f(k, m, Bt, Bp, k0) = r -> r*F(k, m, Bt, Bp)(r)^2/(k0(k,m)(r)^2)
f_(k, m, Bt, Bp) = r -> f(k, m, Bt, Bp, k0_)(r)

shift_taylor(rs,nmax::Int) = rs + Taylor1(typeof(rs),nmax)
taylor_at_rs(f_, nmax::Int, rs) = f_(shift_taylor(rs,nmax))
nth_coeff(taylor_at_rs,n) = getcoeff(taylor_at_rs, n)
denom(σ,n,f2,g0) = f2*(σ + n)*(σ + n + 1) - g0

f2_(f,rs) = getcoeff(taylor_at_rs(f, 2, rs),2)
g0_(g,rs) = getcoeff(taylor_at_rs(g, 0, rs),0)

σ_plus_(f2,g0) = -1/2 + sqrt(1/4 + g0/f2)
σ_minus_(f2,g0) = -1/2 - sqrt(1/4 + g0/f2)

f_tay = f__ -> taylor_at_rs(f__, nmax+3, rs)
g_tay = g__ -> taylor_at_rs(g__, nmax+1, rs)

function ξseries(σ, f2, g0, f_taylor, g_taylor, rs, nmax; ξ0 = 1.0)
    ξs = Number[]
    push!(ξs,ξ0)

    for n=1:(nmax+1)
        internal_sum = 0.0
        for m=0:(n-1)
            internal_sum += ξs[m+1]*(nth_coeff(f_taylor,n-m+2)*(σ+m)*(σ+n+1)-nth_coeff(g_taylor,n-m))
        end
        push!(ξs,-(internal_sum/denom(σ,n,f2,g0)))
    end

    return ξs
end

function ξseries_plus_minus(f__::Function, g__::Function, rs, nmax; ξ0 = 1.0)
    f_taylor = taylor_at_rs(f__, nmax+3, rs)
    g_taylor = taylor_at_rs(g__, nmax+1, rs)

    f2 = nth_coeff(f_taylor,2)
    g0 = nth_coeff(g_taylor,0)
    σ_plus = σ_plus_(f2,g0)
    σ_minus = σ_minus_(f2,g0)

    ξs_plus = ξseries(σ_plus, f2, g0, f_taylor, g_taylor, rs, nmax; ξ0=ξ0)
    ξs_minus = ξseries(σ_minus, f2, g0, f_taylor, g_taylor, rs, nmax; ξ0=ξ0)

    return ξs_plus, ξs_minus, σ_plus, σ_minus, f2, g0
end

function ξ_full(ξs,σ,rs,r,rs0) 
    ξ_ser = 0.0
    
    x = r-rs
    
    for i in 1:length(ξs)
        ξ_ser += x^(i-1)*ξs[i]
    end
        
    ξ_ser *= abs(x)^σ

    return ξ_ser
end

ξ(ξs,σ,rs,rs0) = r -> ξ_full(ξs,σ,rs,r,rs0) 

function ξ_plus_minus(f::Function, g::Function, rs, rs0, nmax; kwargs...)
    ξs_plus, ξs_minus, σ_plus, σ_minus, f2, g0 = ξseries_plus_minus(f, g, rs, nmax; kwargs...)

    ξ_plus = ξ(ξs_plus,σ_plus,rs,rs0)
    ξ_minus = ξ(ξs_minus,σ_minus,rs,rs0)

    return ξ_plus, ξ_minus, σ_plus, σ_minus
end

function κ_(Bp, Bt, dpdr, k, m, rs; g1=gzero,test=false)
    κ_Furth(g1) = r -> (g1(k, m, dpdr, Bt, Bp)(r)/F(k, m, Bt, Bp)(r)+(dH(k, m)(r)*dF(k, m, Bt, Bp)(r)+H(k, m)(r)*d2F(k, m, Bt, Bp)(r)))/(H(k, m)(r)*dF(k, m, Bt, Bp)(r))
    κ_Nish(g1) = r -> (g1(k, m, dpdr, Bt, Bp)(r)/(F(k, m, Bt, Bp)(r)^2)+(dH(k, m)(r)*dF(k, m, Bt, Bp)(r)+H(k, m)(r)*d2F(k, m, Bt, Bp)(r))/F(k, m, Bt, Bp)(r))/H(k, m)(r)
    κ_Glass(g1) = r -> (r-rs)*κ_Nish(g1)(r)

    if test
        rvec1 = range(rs-0.01,rs-0.00000000001,length=100)
        rvec2 = range(rs+0.01,rs+0.00000000001,length=100)
        p11 = plot([rvec1,rvec2],[κ_Furth(g1).(rvec1),κ_Furth(g1).(rvec2)],color=:black,label=["κ Furth"  false])
        p11 = plot!([rvec1,rvec2],[κ_Glass(g1).(rvec1),κ_Glass(g1).(rvec2)],ylims = (minimum(vcat(κ_Furth(g1).(rvec1),κ_Furth(g1).(rvec2))),maximum(vcat(κ_Furth(g1).(rvec1),κ_Furth(g1).(rvec2)))),color=:blue,label=["κ Glass"  false])
        p11 = vline!([rs],label="rs",linestyle=:dash)
        display(p11)
    end

    if !isnan(κ_Furth(g1)(rs))
        return κ_Furth(g1)(rs)
    else
        @warn("κ was nan, using an approximation (rel error = $((κ_Furth(g1)(rs+1e-10)-κ_Glass(g1)(rs+1e-10))/κ_Furth(g1)(rs+1e-10)))")
            #
        return κ_Furth(g1)(rs+1e-10)
    end
end

function κ_Nish(Bp, Bt, dpdr, k, m, rs; g=gzero)
    f__= r -> f_(k, m, Bt, Bp)(r)
    g__ = r -> gzero(k, m, dpdr, Bt, Bp)(r)

    g = g__(rs)
    F_ = Fdag(k, m, Bt, Bp)(rs)

    Fdag_deriv(r) = ForwardDiff.derivative(Fdag(k, m, Bt, Bp),r)

    H(k, m)
    #Fdag(k, m, Bt, Bp)(r) #could be Fdag
    #H(k, m)
    κ = (g/f_+ForwardDiff.derivate(x -> H(k, m)(x)*Fdag_deriv(x),rs))/(H(k, m)(rs)*Fdag_deriv(rs))
end

function ψ_plus_minus_zeroPressure(Bp, Bt, dpdr, k, m, rs; test_κ=false)
    κ = κ_(Bp, Bt, dpdr, k, m, rs; g1=gzero, test=test_κ)

    ψ_s = r -> 1+κ*(r-rs)*log(abs(r-rs))+(1/2)*κ^2*(r-rs)^2*log(abs(r-rs))-(3/4)*κ^2*(r-rs)^2
    ψ_b = r -> (r-rs)+(1/2)*κ*(r-rs)^2+(1/12)*κ^2*(r-rs)^3

    #ψ_s_Nish = r -> ψ_s(r)
    #ψ_b_Nish = r -> ψ_b(r)

    #ψ_s_Glass = r -> ψ_s(r) 
    #ψ_b_Glass = r -> (r-rs)+(1/2)*κ*(r-rs)^2+(1/12)*κ^2*(r-rs)^2 #probs wrong

    ψ_s_prime = r -> ForwardDiff.derivative(ψ_s,r)
    ψ_b_prime = r -> ForwardDiff.derivative(ψ_b,r)

    return ψ_s, ψ_b, ψ_s_prime, ψ_b_prime
end

function zl_zr_calculator(solin, solout, inScale, outScale, del, k, m, Bt, Bp, rs; del2 = 0.0, use_ψ=false)
    if del2 == 0.0 || (del2 != 0.0 && del2 < del)
        del2 = del
        outScale2 = outScale
    else
        outScale2 = inScale*solin(rs-del2)[2]/(solout(rs+del2)[2])
    end

    if use_ψ
        zl = solin(rs-del2)[2]/solin(rs-del2)[1]
        zr = solout(rs+del2)[2]/solout(rs+del2)[1]
    else
        ξ_num_in(r) = solin(r)[2]/F(k, m, Bt, Bp)(r)
        ξ_num_out(r) = solout(r)[2]/F(k, m, Bt, Bp)(r)

        ξprime_num_in(r) =  (solin(r)[1]*F(k, m, Bt, Bp)(r)-solin(r)[2]*dF(k, m, Bt, Bp)(r))/(F(k, m, Bt, Bp)(r)^2)
        ξprime_num_out(r) =  (solout(r)[1]*F(k, m, Bt, Bp)(r)-solout(r)[2]*dF(k, m, Bt, Bp)(r))/(F(k, m, Bt, Bp)(r)^2)

        zl = ξ_num_in(rs-del2)/ξprime_num_in(rs-del2)
        zr = ξ_num_out(rs+del2)/ξprime_num_out(rs+del2)
    end

    return zl,zr,del2,outScale2
end

function plot_full_Psis_w_frobenius(solin, solout, inScale, outScale, ξl, ξr, del, nmax, rs; plotwidth=10, plot_soln_1=false, res = 100)

    closevec_in = range(rs-plotwidth*del,rs-del,length=res)
    closevec_out = range(rs+plotwidth*del,rs+del,length=res)

    p1 = plot([closevec_in  closevec_out],[[solin(i)[2]*inScale for i in closevec_in]  [solout(i)[2]*outScale for i in closevec_out]], label = "Psi")
    !plot_soln_1 && (p1 = plot!([closevec_in  closevec_out],[[ξl(i) for i in closevec_in]  [ξr(i) for i in closevec_out]], xlims = (rs-plotwidth*del,rs+plotwidth*del), label="Psi expansion O($(nmax))"))

    display(p1)

    return p1
end

#I converted psi (integrated) to ξ
    #Did matching with ξ
    #Convert the matched ξ's back to Psi's
    #All occurs within zl_zr_calculator, Δl_Δr_calculator.

function Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=10, plot_soln = true, plot_soln_1=false, kwargs...)
    f__= r -> f_(k, m, Bt, Bp)(r)
    g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
    ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rs, rs0, nmax; kwargs...)
    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=integrator_reltol, plot_solution=plot_solution)

    ξ_plus_prime(r) = ForwardDiff.derivative(ξ_plus,r)
    ξ_minus_prime(r) = ForwardDiff.derivative(ξ_minus,r)

    zl, zr, del2, outScale2 = zl_zr_calculator(solin, solout, inScale, outScale, del, k, m, Bt, Bp, rs; del2 = del2)

    Δl = (zl*ξ_minus_prime(rs-del2)-ξ_minus(rs-del2))/(ξ_plus(rs-del2)-zl*ξ_plus_prime(rs-del2))
    Δr = (zr*ξ_minus_prime(rs+del2)-ξ_minus(rs+del2))/(ξ_plus(rs+del2)-zr*ξ_plus_prime(rs+del2))

    ############################################################################################################################################################
    #cl = solin(rs-del2)[2]*inScale/(ξ_minus(rs-del2)+Δl*ξ_plus(rs-del2))
    #cr = solout(rs+del2)[2]*outScale2/(ξ_minus(rs+del2)+Δr*ξ_plus(rs+del2))
    #ξl(r) = cl*(ξ_minus(r)+Δl*ξ_plus(r))
    #ξr(r) = cr*(ξ_minus(r)+Δr*ξ_plus(r))

    #plot_soln && plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, ξl, ξr, del2, nmax, rs; plotwidth=plotwidth, plot_soln_1=plot_soln_1)
    ############################################################################################################################################################

    ############################################################################################################################################################
    ξl(r) = ξ_minus(r)+Δl*ξ_plus(r)
    ξr(r) = ξ_minus(r)+Δr*ξ_plus(r)

    psi_l_raw(r) = ξl(r)*F(k, m, Bt, Bp)(r)
    psi_r_raw(r) = ξr(r)*F(k, m, Bt, Bp)(r)

    cl = solin(rs-del2)[2]*inScale/psi_l_raw(rs-del2)
    cr = solout(rs+del2)[2]*outScale2/psi_r_raw(rs+del2)

    psi_l(r) = cl*psi_l_raw(r)
    psi_r(r) = cr*psi_r_raw(r)

    plot_soln && plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=plotwidth, plot_soln_1=plot_soln_1)
    ############################################################################################################################################################

    return Δl,Δr,del2
end

function Δl_Δr_calculator(Bp, Bt, dpdr, k, m, rs, ξ_plus, ξ_minus, solin, solout, inScale, outScale, del, nmax, del2; plotwidth=10, plot_soln = true, plot_soln_1=false, kwargs...)
    if del2 < del
        @error("del2 must be greater than or equal to del")
    end

    ξ_plus_prime(r) = ForwardDiff.derivative(ξ_plus,r)
    ξ_minus_prime(r) = ForwardDiff.derivative(ξ_minus,r)

    zl, zr, del2, outScale2 = zl_zr_calculator(solin, solout, inScale, outScale, del, k, m, Bt, Bp, rs; del2 = del2)

    Δl = (zl*ξ_minus_prime(rs-del2)-ξ_minus(rs-del2))/(ξ_plus(rs-del2)-zl*ξ_plus_prime(rs-del2))
    Δr = (zr*ξ_minus_prime(rs+del2)-ξ_minus(rs+del2))/(ξ_plus(rs+del2)-zr*ξ_plus_prime(rs+del2))

    ############################################################################################################################################################
    ξl(r) = ξ_minus(r)+Δl*ξ_plus(r)
    ξr(r) = ξ_minus(r)+Δr*ξ_plus(r)

    psi_l_raw(r) = ξl(r)*F(k, m, Bt, Bp)(r)
    psi_r_raw(r) = ξr(r)*F(k, m, Bt, Bp)(r)

    cl = solin(rs-del2)[2]*inScale/psi_l_raw(rs-del2)
    cr = solout(rs+del2)[2]*outScale2/psi_r_raw(rs+del2)

    psi_l(r) = cl*psi_l_raw(r)
    psi_r(r) = cr*psi_r_raw(r)

    plot_soln && plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=plotwidth, plot_soln_1=plot_soln_1)
    ############################################################################################################################################################

    return Δl,Δr,del2
end

#Must use g_int = gzero
function Δl_Δr_calculator_zeroPressure(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, del; test_κ=false, integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=10, plot_soln = true, plot_soln_1=false, kwargs...)
    ψ_plus, ψ_minus, ψ_plus_prime, ψ_minus_prime = ψ_plus_minus_zeroPressure(Bp, Bt, dpdr, k, m, rs; test_κ=test_κ)
    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=gzero, integrator_reltol=integrator_reltol, plot_solution=plot_solution)

    zl, zr, del2, outScale2 = zl_zr_calculator(solin, solout, inScale, outScale, del, k, m, Bt, Bp, rs; del2 = del2, use_ψ=true)

    Δl = (zl*ψ_minus_prime(rs-del2)-ψ_minus(rs-del2))/(ψ_plus(rs-del2)-zl*ψ_plus_prime(rs-del2))
    Δr = (zr*ψ_minus_prime(rs+del2)-ψ_minus(rs+del2))/(ψ_plus(rs+del2)-zr*ψ_plus_prime(rs+del2))

    ############################################################################################################################################################
    #cl = solin(rs-del2)[2]*inScale/(ξ_minus(rs-del2)+Δl*ξ_plus(rs-del2))
    #cr = solout(rs+del2)[2]*outScale2/(ξ_minus(rs+del2)+Δr*ξ_plus(rs+del2))
    #ξl(r) = cl*(ξ_minus(r)+Δl*ξ_plus(r))
    #ξr(r) = cr*(ξ_minus(r)+Δr*ξ_plus(r))

    #plot_soln && plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, ξl, ξr, del2, nmax, rs; plotwidth=plotwidth, plot_soln_1=plot_soln_1)
    ############################################################################################################################################################

    ############################################################################################################################################################
    psi_l_raw(r) = ψ_minus(r)+Δl*ψ_plus(r)
    psi_r_raw(r) = ψ_minus(r)+Δr*ψ_plus(r)

    cl = solin(rs-del2)[2]*inScale/psi_l_raw(rs-del2)
    cr = solout(rs+del2)[2]*outScale2/psi_r_raw(rs+del2)

    psi_l(r) = cl*psi_l_raw(r)
    psi_r(r) = cr*psi_r_raw(r)

    nmax = 2
    plot_soln && plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=plotwidth, plot_soln_1=plot_soln_1)
    ############################################################################################################################################################

    return Δl,Δr,del2
end

function Δl_Δr_calculator_zeroPressure(rs, ψ_plus, ψ_minus, ψ_plus_prime, ψ_minus_prime, solin, solout, inScale, outScale, del, del2; test_κ=false, plot_solution=true, plotwidth=10, plot_soln = true, plot_soln_1=false, kwargs...)
    zl, zr, del2, outScale2 = zl_zr_calculator(solin, solout, inScale, outScale, del, 0.0, 0.0, x->nothing, x->nothing, rs; del2 = del2, use_ψ=true)

    Δl = (zl*ψ_minus_prime(rs-del2)-ψ_minus(rs-del2))/(ψ_plus(rs-del2)-zl*ψ_plus_prime(rs-del2))
    Δr = (zr*ψ_minus_prime(rs+del2)-ψ_minus(rs+del2))/(ψ_plus(rs+del2)-zr*ψ_plus_prime(rs+del2))

    ############################################################################################################################################################
    psi_l_raw(r) = ψ_minus(r)+Δl*ψ_plus(r)
    psi_r_raw(r) = ψ_minus(r)+Δr*ψ_plus(r)

    cl = solin(rs-del2)[2]*inScale/psi_l_raw(rs-del2)
    cr = solout(rs+del2)[2]*outScale2/psi_r_raw(rs+del2)

    psi_l(r) = cl*psi_l_raw(r)
    psi_r(r) = cr*psi_r_raw(r)

    nmax = 2
    plot_soln && plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=plotwidth, plot_soln_1=plot_soln_1)
    ############################################################################################################################################################

    return Δl,Δr,del2
end

#delvec = 10 .^ range(-7,-0.5, length = 400)
function Δl_Δr_xmin(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec; big_n=20, xmin_tol = 1e-5, integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=10, plot_soln = false, verbose = false, kwargs...)
    f__= r -> f_(k, m, Bt, Bp)(r)
    g__ = r -> g_(k, m, dpdr, Bt, Bp)(r)
    ξ_plus0, ξ_minus0, σ_plus0, σ_minus0 = ξ_plus_minus(f__, g__, rs, rs0, 0; kwargs...)
    ξ_plus20, ξ_minus20, σ_plus20, σ_minus20 = ξ_plus_minus(f__, g__, rs, rs0, big_n; kwargs...)

    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, minimum(delvec); g=g_, integrator_reltol=integrator_reltol, plot_solution=plot_solution)

    Δl_Δr_calc(ξ_plus,ξ_minus) = del2 -> Δl_Δr_calculator(Bp, Bt, dpdr, k, m, rs, ξ_plus, ξ_minus, solin, solout, inScale, outScale, del, nmax, del2; plotwidth=plotwidth, plot_soln = plot_soln, plot_soln_1=false)

    Δvec0 = []
    Δvec20 = []
    for i in 1:length(delvec)
        Δl0,Δr0,del2 = Δl_Δr_calc(ξ_plus0,ξ_minus0)(delvec[i])
        Δl20,Δr20,del2 = Δl_Δr_calc(ξ_plus20,ξ_minus20)(delvec[i])

        push!(Δvec0,(Δl0+Δr0))
        push!(Δvec20,(Δl20+Δr20))
    end

    Δdiffvec = Δvec20 .- Δvec0
    min_Arg = 0
    for j = 1:length(delvec)
        if abs(Δdiffvec[end+1-j]) < xmin_tol
            min_Arg = length(delvec)+1-j
            break
        end
        verbose && print(abs(Δdiffvec[end+1-j]),"\n")
    end
    if min_Arg == 0
        @warn("xmin_tol too small, outer expansions didn't converge")
        min_Arg = 1
        plot_soln = true
    end

    output = Δl_Δr_calc(ξ_plus0,ξ_minus0)(delvec[min_Arg])

    pdel = Plots.plot(delvec, Δvec0, xaxis=:log, legend=:bottomleft, label = "O(0)", title = "Δl + Δr as a function of matching width", xlabel="del")
    pdel = plot!(delvec, Δvec20, xaxis=:log, legend=:bottomleft, label = "O($(big_n))", title = "Δl + Δr as a function of matching width", xlabel="del",  ylabel = "Δl + Δr")
    pdel = vline!([delvec[min_Arg]],label="xmin",linestyle=:dash)
    plot_soln && display(pdel)

    return delvec[min_Arg], output, Δdiffvec[min_Arg]
end

#delvec = 10 .^ range(-7,-0.5, length = 400)
#Integrates separately for each del in delvec. #Δl_Δr_xmin (new vers above) integrates once into the smallest del, then uses del2 functionality to match at successively larger dels without having to re-integrate
function Δl_Δr_xmin_deprecated(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec; xmin_tol = 1e-5, integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=10, plot_soln = true, verbose = false, kwargs...)
    Δl_Δr_calc1(nmax) = del -> Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false)

    Δvec0 = []
    Δvec20 = []
    for i in 1:length(delvec)
        Δl0,Δr0,del2 = Δl_Δr_calc1(0)(delvec[i])
        Δl20,Δr20,del2 = Δl_Δr_calc1(20)(delvec[i])

        push!(Δvec0,(Δl0+Δr0))
        push!(Δvec20,(Δl20+Δr20))
    end

    Δdiffvec = Δvec20 .- Δvec0
    min_Arg = 0
    for j = 1:length(delvec)
        if abs(Δdiffvec[end+1-j]) < xmin_tol
            min_Arg = length(delvec)+1-j
            break
        end
        verbose && print(abs(Δdiffvec[end+1-j]),"\n")
    end
    if min_Arg == 0
        @warn("xmin_tol too small, outer expansions didn't converge")
        min_Arg = 1
        plot_soln = true
    end

    pdel = plot(delvec, Δvec0, xaxis=:log, legend=:bottomleft, label = "O($(0))", title = "Δl + Δr as a function of matching width", xlabel="del")
    pdel = plot!(delvec, Δvec20, xaxis=:log, legend=:bottomleft, label = "O($(20))", title = "Δl + Δr as a function of matching width", xlabel="del",  ylabel = "Δl + Δr")
    pdel = vline!([delvec[min_Arg]],label="xmin",linestyle=:dash)
    plot_soln && display(pdel)

    return delvec[min_Arg], Δl_Δr_calc1(0)(Δdiffvec[min_Arg]), Δdiffvec[min_Arg]
end

#delvec = 10 .^ range(-7,-0.5, length = 400)
function Δl_Δr_xmin_visualiser(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, delvec; finite_pressure=true, test_κ=false, big_n=20, xmin_tol = 1e-5, integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=10, plot_soln = false, verbose = false, kwargs...)
    if finite_pressure
        nmax = -1 #deprecated variable 
        return Δl_Δr_xmin(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec; big_n=big_n, xmin_tol = xmin_tol, integrator_reltol=integrator_reltol, plot_solution=true, del2 = del2, plotwidth=plotwidth, plot_soln = true, verbose = verbose, kwargs...)
    else
        ψ_plus, ψ_minus, ψ_plus_prime, ψ_minus_prime = ψ_plus_minus_zeroPressure(Bp, Bt, dpdr, k, m, rs; test_κ=test_κ)
        
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, minimum(delvec); g=gzero, integrator_reltol=integrator_reltol, plot_solution=plot_solution)
    
        Δl_Δr_calc = del2 -> Δl_Δr_calculator_zeroPressure(rs, ψ_plus, ψ_minus, ψ_plus_prime, ψ_minus_prime, solin, solout, inScale, outScale, del, del2; test_κ=false, integrator_reltol=integrator_reltol, plot_solution=plot_solution, plot_soln = plot_soln, plot_soln_1=false)
    
        Δvec = []
        for i in 1:length(delvec)
            Δl,Δr,del2 = Δl_Δr_calc(delvec[i])
            push!(Δvec,(Δl+Δr))
        end
    
        pdel = plot(delvec, Δvec, xaxis=:log, legend=:bottomleft, label = "O($(0))", title = "Δl + Δr as a function of matching width", xlabel="del")
        display(pdel)
    
        return delvec, Δvec
    end
end
