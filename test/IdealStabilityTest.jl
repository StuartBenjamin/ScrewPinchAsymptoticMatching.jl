
using TaylorSeries
using MatrixEquations
using ForwardDiff
using Distances
using Plots
using Statistics
using NaNStatistics
using LinearAlgebra
using RootsAndPoles
using SpecialFunctions
using JLD2
using CSV
using DataFrames
using Optim
using Roots
using Random
using OrdinaryDiffEq

using Pkg
Pkg.add(PackageSpec(url="https://github.com/StuartBenjamin/ScrewPinchAsymptoticMatching.jl"))
using ScrewPinchAsymptoticMatching
ScrewPinchAsymptoticMatching.mod_path=chop(pathof(ScrewPinchAsymptoticMatching);tail=31)


##########################################################################################################################################################
#GEN EQUIL
##########################################################################################################################################################

    mu0 = ScrewPinchAsymptoticMatching.mu0


    if true
        q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
        rs0 = 1 #CESAR BE 2 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt) 
        R0 = 20 #Controls ratio of Bt to Bp
        ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
        xb = 2 #CESAR BE 1 #Widen device to increase q at edge (peaks current more stongly)
        Bp0 = 0.5 #Size of Bp (ratio of Bp to Bt unchanged)

        β = 1e-2#14; #Gives a true 0-pressure Chandra Equil

        m = 0;
        n = 2;
        k = ScrewPinchAsymptoticMatching.k_(n,R0)
        c0 = 1;
        r0 = 0.000001;
        del=1e-4;
        nmax=3;
    end

    Furth_Equil_qvar(plot_equil) = q0 -> Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=plot_equil, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))
    Chandra_Equil_qvar(plot_equil) = q0 -> Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=plot_equil, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

    Furth_Equil_qvar(true)(1.1)
    Chandra_Equil_qvar(true)(1.1)

    η_diff = 1.371e-5 #Checked w Cesar
    ρ = 1e20*2*1.6605e-27 #Deuterium!!

    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil_qvar(true)(1.2)
    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]


    Bpc,Btc,qc,dpdrc,pc,Jtc,Jpc,rb,outerp6 = Chandra_Equil_qvar(true)(1.2)
    rsc, rs_plot = find_rs(qc,m,n,rb,q0,ν,rs0)
    rsc = rsc[1]


    Jt_50 = ScrewPinchAsymptoticMatching.CubicSpline(range(0.0,2.0;step=2.0/50),Jt.(range(0.00000000001,2.0;length=51)))
    p_50 = ScrewPinchAsymptoticMatching.CubicSpline(range(0.0,2.0;step=2.0/50),p.(range(0.00000000001,2.0;length=51)))
    p_50 = ScrewPinchAsymptoticMatching.CubicSpline(range(0.0,2.0;step=2.0/7),vcat([pc(0.0),0.98*pc(0.0),0.9*pc(0.0)],pc.(range((2.0/7)*3,2.0;length=5))))

    Bp_s50,Bt_s50,q_s50,dpdr_s50,p_s50,Jt_s50,Jp_s50,rb_s50,outerp7,Jp2 = Spline_Equil(Jt_50,p_50,Bt(0.0),R0)

    #Bp_s50,Bt_s50,q_s50,dpdr_s50,p_s50,Jt_s50,Jp_s50,rb_s50,outerp7,Jp2 = Spline_Equil(Jt_50,p,Bt(0.0),R0;dpdr=dpdr)

    rs_s50, rs_plot = ScrewPinchAsymptoticMatching.find_rs(q_s50,m,n,rb)
    rs_s50 = rs_s50[1]

##########################################################################################################################################################
#Test stability:
##########################################################################################################################################################

F = ScrewPinchAsymptoticMatching.F
make_psi_ODE_RUTH = ScrewPinchAsymptoticMatching.make_psi_ODE_RUTH
integrate_psi_in = ScrewPinchAsymptoticMatching.integrate_psi_in

test_Suydam(Bt, q, dpdr, rb;plotresults=false)
test_Suydam(Bt, q, dpdr, rb;plotresults=true)

test_Suydam(Bt_s50, q_s50, dpdr_s50, rb;plotresults=false)
test_Suydam(Bt_s50, q_s50, dpdr_s50, rb;plotresults=true)


test_ideal_stability(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del, q)
test_ideal_stability(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del, q;case=1)
test_ideal_stability(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del, q;case=2)
test_ideal_stability(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del, q;case=3)

test_ideal_stability(Bpc, Btc, dpdrc, k, m, r0, rsc, rb, rs0, nmax, del, qc)
test_ideal_stability(Bpc, Btc, dpdrc, k, m, r0, rsc, rb, rs0, nmax, del, qc;case=1)
test_ideal_stability(Bpc, Btc, dpdrc, k, m, r0, rsc, rb, rs0, nmax, del, qc;case=2)
test_ideal_stability(Bpc, Btc, dpdrc, k, m, r0, rsc, rb, rs0, nmax, del, qc;case=3)

test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50)
test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50;case=1)
test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50;case=2)
test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50;case=3)

Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del)
Δl_Δr_calculator(Bpc, Btc, dpdrc, k, m, r0, rsc, rb, rs0, nmax, del)
Δl_Δr_calculator(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del)

test_ideal_stability(Bp, Bt, dpdr, R0, r0, rb;verbose=true, ignore_Suydam=true)
test_ideal_stability(Bpc, Btc, dpdrc, R0, r0, rb;verbose=true, ignore_Suydam=false)
test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, R0, r0, rb;verbose=true, ignore_Suydam=true)

f__= r -> f_(k, m, Btc, Bpc)(r)
g__ = r -> g_(k, m, dpdrc, Btc, Bpc)(r)
ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f__, g__, rsc, rs0, nmax)

plot_ideal_stability_solout(Bpc, Btc, dpdrc, k, m, rsc, rb, del, ξ_plus)

#test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50;case=3)
#test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50;case=3)
#test_ideal_stability(Bp_s50, Bt_s50, dpdr_s50, k, m, r0, rs_s50, rb, rs0, nmax, del, q_s50;case=3)