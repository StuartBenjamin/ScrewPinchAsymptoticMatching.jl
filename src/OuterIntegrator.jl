using LinearAlgebra
using Statistics
using Interpolations
using StaticArrays
using DifferentialEquations
using OrdinaryDiffEq
using ForwardDiff
using Plots
using Optim

#cd("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia")
include("ChandraScrewPinchEquilibrium.jl")

##########################################################################################################################################################
#Inputs
k_(n,R0) = n/R0;
#m = 2;
#n = 1;
#qtest = m/n;
#c0 = 1;
#r0 = 0.001;

##########################################################################################################################################################
#Setting up cylindrical Newcomb equation
F(k, m, Bt, Bp) = r -> k*Bt(r) - (m/r)*Bp(r)
Fdag(k, m, Bt, Bp) = r -> k*Bt(r) + (m/r)*Bp(r)
H(k, m) = r -> r^3/(k^2*r^2 + m^2)

k0_(k,m) = r -> sqrt(k^2+m^2/r^2)
g(k, m, dpdr, Bt, Bp, k0) = r -> 2*dpdr(r)*k^2/(k0(k,m)(r)^2) + (r^2*k0(k,m)(r)^2 - 1)*r*F(k, m, Bt, Bp)(r)^2/(r^2*k0(k,m)(r)^2) + 2*(k^2/(r*k0(k,m)(r)^4))*F(k, m, Bt, Bp)(r)*Fdag(k, m, Bt, Bp)(r) #From ideal MHD
g_(k, m, dpdr, Bt, Bp) = r -> g(k, m, dpdr, Bt, Bp, k0_)(r)  #From ideal MHD
gzero(k, m, dpdr, Bt, Bp) = r -> 2*0.0*k^2/(k0_(k,m)(r)^2) + (r^2*k0_(k,m)(r)^2 - 1)*r*F(k, m, Bt, Bp)(r)^2/(r^2*k0_(k,m)(r)^2) + 2*(k^2/(r*k0_(k,m)(r)^4))*F(k, m, Bt, Bp)(r)*Fdag(k, m, Bt, Bp)(r) 

dH(k, m) = r -> -2*k^2*r^4/((m^2 + k^2*r^2)^2) + 3*r^2/(m^2 + k^2*r^2)
dF(k, m, Bt, Bp) = r -> ForwardDiff.derivative(F(k, m, Bt, Bp), r)
d2F(k, m, Bt, Bp) = r -> ForwardDiff.derivative(dF(k, m, Bt, Bp), r)

##########################################################################################################################################################
#Solving for rs
function find_rs_Optim(q,m,n,rb,q0,ν,rs0;verbose=true)
    f1= x -> abs(q(first(x))-m/n)
    qtest = m/n
    res = optimize(f1,[q_Furth_find_rs(q0,ν,rs0)(qtest)],LBFGS())
    #res = optimize(ftest,xb/2,rb)


    if Optim.converged(res) #&& !(abs(res.minimizer-a)<0.1*a || abs(res.minimizer)<0.1)
        rs = res.minimizer
        p1 = plot(0:(rb/300):rb,[q(i) for i in 0:(rb/300):rb],title="q profile",label=false)
        p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="q = $(m)/$(n)",line=:dash)
        p1 = vline!([res.minimizer],label="resonant surface location",line=:dash)
    else    
        @warn "Your chosen rbtional surface lies outside your minor rbdius!"
        p1 = plot(0:(rb/300):rb,[q(i) for i in 0:(rb/300):rb],title="q profile",label=false)
        p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="proposed rational q-value",line=:dash)
        p1 = vline!([res.minimizer],label="proposed resonant surface location",line=:dash)
    end

    verbose && display(p1)

    return res.minimizer[1], p1
end

function find_rs(q,m,n,rb,q0,ν,rs0; verbose=true, useOptim=false)
    useOptim && (return find_rs_Optim(q,m,n,rb,q0,ν,rs0;verbose=verbose))

    qtest = m/n
    f1 = x -> abs(q(first(x))-qtest)

    zeros = find_zeros(f1,0.0,rb)

    if length(zeros)>1
        print("Multiple rs values detected: \n") 
        print("rs = {$(zeros[1])")
        for i in 2:length(zeros)
            print(",$(zeros[i])")
        end
        print("}\n")
    elseif length(zeros)==0
        print("No zeros detected.\n")
        @warn "Your chosen rational surface lies outside your minor radius!"
        return find_rs_Optim(q,m,n,rb,q0,ν,rs0;verbose=verbose)
    end

    rs = zeros[1]

    p1 = plot(0:(rb/300):rb,[q(i) for i in 0:(rb/300):rb],title="q profile",label=false)
    p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="q = $(m)/$(n)",line=:dash, xlabel="r (m)")
    p1 = vline!([rs],label="resonant surface location",line=:dash)

    verbose && display(p1)

    return rs, p1
end

##########################################################################################################################################################
#Integrator functions

@inline function psi_ODE_RUTH(psi, psi_deriv, r, Bp::Function, Bt::Function, dpdr::Function, k, m::Int; g::Function = g_)
    return psi[1]*(g(k, m, dpdr, Bt, Bp)(r)/(F(k, m, Bt, Bp)(r)^2)+(dH(k, m)(r)*dF(k, m, Bt, Bp)(r)+H(k, m)(r)*d2F(k, m, Bt, Bp)(r))/F(k, m, Bt, Bp)(r))/H(k, m)(r) - psi_deriv[1]*dH(k, m)(r)/H(k, m)(r)
end

function make_psi_ODE_RUTH(Bp::Function, Bt::Function, dpdr::Function, k, m::Int; g::Function = g_)
    psi_ode = function f(du,u,p,t)
        ddu = [psi_ODE_RUTH(u[1], du[1], t, Bp, Bt, dpdr, k, m; g=g)]
    end
    return psi_ode
end

function integrate_psi_in(Bp::Function, Bt::Function, dpdr::Function, k, m::Int, r0, r1; integrator_reltol=10^(-20), g::Function = g_)
    psi_ode_problem = SecondOrderODEProblem(make_psi_ODE_RUTH(Bp, Bt, dpdr, k, m; g=g),[(m - 1)*r0^(m - 2)],[r0^(m - 1)],(r0,r1))
    sol = solve(psi_ode_problem,Tsit5(),reltol=integrator_reltol)
    return sol
end

function integrate_psi_out(Bp::Function, Bt::Function, dpdr::Function, k, m::Int, rb, r3; integrator_reltol=10^(-20), g::Function = g_)
    psi_ode_problem = SecondOrderODEProblem(make_psi_ODE_RUTH(Bp, Bt, dpdr, k, m; g=g),[-1.0],[0.0],(rb,r3))
    sol = solve(psi_ode_problem,Tsit5(),reltol=integrator_reltol)
    return sol
end

##########################################################################################################################################################
#Double integration w matching (constant psi)
maxdel(sol) = maximum([i[2] for i in sol.u])

function Psi_w_scales(Bp::Function, Bt::Function, dpdr::Function, k, m::Int, r0, rs, rb, del; integrator_reltol=10^(-20), g::Function = g_, plot_solution=true)
    solin = integrate_psi_in(Bp, Bt, dpdr, k, m, r0, rs-del; integrator_reltol=integrator_reltol, g=g)
    solout = integrate_psi_out(Bp, Bt, dpdr, k, m, rb, rs+del; integrator_reltol=integrator_reltol, g=g)

    maxin = maxdel(solin)
    outScale = solin.u[end][2]/(maxin*solout.u[end][2])

    if plot_solution
        display(plot_full_Psis(solin, solout, (1/maxin), outScale))
    end

    return solin, solout, (1/maxin), outScale, del
end

##########################################################################################################################################################
#Plotting (constant psi)

function plot_full_Psis(solin, solout, inScale, outScale)
    p1 = plot([solin.t;reverse(solout.t)],[[i[2]*inScale for i in solin.u];reverse([i[2]*outScale for i in solout.u])])
    return p1
end

##########################################################################################################################################################
#Raw delta prime  (constant psi)
function raw_delta_prime(Bp::Function, Bt::Function, dpdr::Function, k, m::Int, r0, rs, rb, del;del2=0.0, integrator_reltol=10^(-20), g::Function = g_, plot_solution=true)
    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; integrator_reltol=integrator_reltol, g=g)

    if del2 == 0.0 || (del2 != 0.0 && del2 < del)
        del2 = del
        outScale2 = outScale
    else
        outScale2 = inScale*solin(rs-del2)[2]/(solout(rs+del2)[2])
    end

    if plot_solution
        display(plot_full_Psis(solin, solout, inScale, outScale))
    end

    delta_prime = (solout(rs+del2)[1]*outScale2-solin(rs-del2)[1]*inScale)/(0.5*(solin(rs-del2)[2]*inScale+solout(rs+del2)[2]*outScale2))

    return (delta_prime, del2), solin, solout, inScale, outScale, del 
end

function raw_delta_prime(solin, solout, inScale, outScale, rs, del; del2=0.0)
    if del2 == 0.0 || (del2 != 0.0 && del2 < del)
        del2 = del
        outScale2 = outScale
    else
        outScale2 = inScale*solin(rs-del2)[2]/(solout(rs+del2)[2])
    end

    delta_prime = (solout(rs+del2)[1]*outScale2-solin(rs-del2)[1]*inScale)/(0.5*(solin(rs-del2)[2]*inScale+solout(rs+del2)[2]*outScale2))

    return (delta_prime, del2)
end
