using LinearAlgebra
using Statistics
using Interpolations
using StaticArrays
using DifferentialEquations
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
function find_rs(q,m,n,rb,q0,ν,rs0;verbose=true)
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

    return res.minimizer, p1
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

##########################################################################################################################################################
#Tests:
if false #Checking functions are working properly
    β = 0.0000001 
    rs0 = 1.0 
    R0 = 20 
    ν = 1.0  
    xb = 2.0 
    Bp0 = 0.0 
    q0 = 1.0

    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

    m = 2;
    n = 1;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.001;

    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]

    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, 0.0001; integrator_reltol=10^(-20), plot_solution=true)
    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, 0.0001; g=gzero, integrator_reltol=10^(-20), plot_solution=true)
    delta_prime_raw, solin, solout, inScale, outScale, del =raw_delta_prime(Bp, Bt, dpdr, n/R0, m, r0, rs, rb, 0.01; del2=0.0, integrator_reltol=10^(-20));
    delta_prime_raw, solin, solout, inScale, outScale, del =raw_delta_prime(Bp, Bt, dpdr, n/R0, m, r0, rs, rb, 0.01; g=gzero, del2=0.0, integrator_reltol=10^(-20));
    raw_delta_prime(solin, solout, inScale, outScale, del; del2=0.02)
end
if false #Sorting out different ways to write g
    β = 0.1 
    rs0 = 2.0 
    R0 = 20 
    ν = 1.0  
    xb = 1.0 
    Bp0 = 1.0 
    q0 = 1.1

    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

    m = 2;
    n = 1;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.001;
    rvec = range(0.000001,rb,200)

    gCop(k, n, m, dpdr, Bt, Bp) = r -> (2*n^2*k^2*r^2/(n^2*k^2*r^2+m^2))*dpdr(r) + (1/r)*(n*k*r*Bt(r)-m*Bp(r))^2*(n^2*k^2*r^2+m^2-1)/(n^2*k^2*r^2+m^2)+(2*n^2*k^2*r)/((n^2*k^2*r^2+m^2)^2)*(n^2*k^2*r^2*Bt(r)^2-m^2*Bp(r)^2)  #From Coppi '66
    g_old(k, m, dpdr, Bt, Bp) = r -> ((m^2 - 1)*r*F(k, m, Bt, Bp)(r)^2 + k^2*r^2*(2*dpdr(r) + r*F(k, m, Bt, Bp)(r)^2 + F(k, m, Bt, Bp)(r)*2*(k*r*Bt(r) - m*Bp(r))/(k^2*r^2 + m^2)))/(k^2*r^2 + m^2)  #From Furth '73

    #gCop is identical to g_ from Ideal MHD. g_old (from Furth) is very similar (small difference on order of <1%).

    p2b = plot([rvec,rvec,rvec],[g_(k, m, dpdr, Bt, Bp).(rvec),g_old(k, m, dpdr, Bt, Bp).(rvec),gCop(k, n, m, dpdr, Bt, Bp).(rvec)];label=["g_"   "g_old" "gCop"])
    plot(p1b,p2b)
    p3b = plot([rvec,rvec],[(g_(k, m, dpdr, Bt, Bp).(rvec)-gCop(k, n, m, dpdr, Bt, Bp).(rvec))/g_(k, m, dpdr, Bt, Bp)(2.0),(g_(k, m, dpdr, Bt, Bp).(rvec)-g_old(k, m, dpdr, Bt, Bp).(rvec))/g_(k, m, dpdr, Bt, Bp)(2.0)];label=["g_ - gCop (normalised)"   "g_ - g_old (normalised)"])

    ##########################################################################################################################################################
    #This below (k0_big,g_big) is wrong, from an old error I made
    k0_big(k,m) = r -> k^2+m^2/r^2 #wrong, an old error I made
    g_big(k, m, dpdr, Bt, Bp) = r -> g(k, m, dpdr, Bt, Bp, k0_big)(r) #wrong, see above
    p1b = plot([rvec,rvec,rvec,rvec],[g_(k, m, dpdr, Bt, Bp).(rvec),g_old(k, m, dpdr, Bt, Bp).(rvec),gCop(k, n, m, dpdr, Bt, Bp).(rvec),g_big(k, m, dpdr, Bt, Bp).(rvec)];label=["g_"   "g_old" "gCop"  "g_big"])
    ##########################################################################################################################################################

    plot(p1b,p2b,p3b)
end