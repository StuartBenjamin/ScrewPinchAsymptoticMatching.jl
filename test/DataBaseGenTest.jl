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
using QuadGK
using Distributions
using StatsBase

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

    m = 2;
    n = 1;
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
rs = find_rs(q,m,n,rb,q0,ν,rs0)


Bpc,Btc,qc,dpdrc,pc,Jtc,Jpc,rb,outerp6 = Chandra_Equil_qvar(true)(1.2)
rsc = find_rs(qc,m,n,rb,q0,ν,rs0)


Jt_50 = ScrewPinchAsymptoticMatching.CubicSpline(range(0.0,2.0;step=2.0/50),Jt.(range(0.00000000001,2.0;length=51)))
p_50 = ScrewPinchAsymptoticMatching.CubicSpline(range(0.0,2.0;step=2.0/50),p.(range(0.00000000001,2.0;length=51)))
p_50 = ScrewPinchAsymptoticMatching.CubicSpline(range(0.0,2.0;step=2.0/7),vcat([pc(0.0),0.98*pc(0.0),0.9*pc(0.0)],pc.(range((2.0/7)*3,2.0;length=5))))

Bp_s50,Bt_s50,q_s50,dpdr_s50,p_s50,Jt_s50,Jp_s50,rb_s50,outerp7,Jp2 = Spline_Equil(Jt_50,p_50,Bt(0.0),R0)

#Bp_s50,Bt_s50,q_s50,dpdr_s50,p_s50,Jt_s50,Jp_s50,rb_s50,outerp7,Jp2 = Spline_Equil(Jt_50,p,Bt(0.0),R0;dpdr=dpdr)

rs_s50 = find_rs(q_s50,m,n,rb)

Jtotmax = 1.1*total_plasma_current(Jt,rb)
Jtotmin = 0.9*total_plasma_current(Jt,rb)
Jtotrange=[Jtotmax, Jtotmin]

P0bounds = [1.2*p_s50(0.0),0.8*p_s50(0.0)]

##########################################################################################################################################################
#Make Jts:
##########################################################################################################################################################

Jts = gen_n_clean_Jts(Jtotmax,rb,5000,80,generate_cleaners(rb ; Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=false,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange)); 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)])
plot_profiles(Jts, rb;ylims = (0.0,1.4*Jt(0.0)))

Jts2 = gen_n_clean_Jts(Jtotmax,rb,5000,80,generate_cleaners(rb ; Jtot_range = Jtotrange,maxgrad_width = rb/70,use_fine=true,use_coarse=true, monotonic=true, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange)); 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.6*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)])
plot_profiles(Jts2, rb;ylims = (0.0,1.4*Jt(0.0)))


##########################################################################################################################################################
#Make Pts: 
##########################################################################################################################################################

P0bounds = [1.2*p_s50(0.0),0.8*p_s50(0.0)]
P0bounds2 = 1e20

ps = gen_n_clean_pressure_profiles(P0bounds2,rb,500,8,generate_cleaners(rb ; maxgrad_width = rb/1000, minval = 0.0, use_fine=true,use_coarse=true, monotonic=true, wobbles=10, max_ref=P0bounds2); maxbatches=100, verbose=true,knotmax=6,knotmin=5)
plot_profiles(ps, rb;ylims = (0.0,1.4*P0bounds2))

##########################################################################################################################################################
#Make Equilibria with Jts:
##########################################################################################################################################################

equila1= gen_equilibria_Jts(Jts,p; Bt0=10.0, dpdr=dpdrc, R0=R0)
plot_equil(equila1;case=3)

equila2= gen_equilibria_Jts(Jts,pc; dpdr=dpdrc, R0=R0)
equila3= gen_equilibria_Jts(Jts,p_s50; Bt0=1.0, R0=R0)

plot_equil(equila1[1])
Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil_qvar(true)(1.2)

plot_equil(equila3[1])

plot_equil_short(equila2[1])

plot_equil(equila1;case=3)
plot_equil(equila2)
plot_equil(equila3)


equilas1 = gen_clean_equilibria(Jts[1:3],p; dpdr=dpdr, Bt0=10, R0=R0, ideal_verbose=true, qtest=1.9);
equilas2 = gen_clean_equilibria(Jts[1:5],p; dpdr=dpdr, Bt0=10, R0=3, m1ncap=20, ideal_verbose=true) 
equilas3 = gen_clean_equilibria(Jts[1:3],p_s50; Bt0=10, R0=R0, ideal_verbose=true, qtest=3.0);
plot_Suydam(equila1[3])

gen_clean_equilibria(Jts[1:3],p_s50; Bt0=10, R0=3, ideal_verbose=true,maxq=0.01) 
plot_equil(equila1[3])

plot_equil(equilas1)
plot_equil(equilas3)

plot_profiles(Jts, rb;ylims = (0.0,1.4*Jt(0.0)))


equilas11=gen_n_clean_equilibria(Jtotmax,p,rb,5000,50; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=3.0, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);


plot_equil(equilas11)
outequils1, inds = run_Δl_Δr_calculator(equilas11, 3, 1, 1e-4, 8, 1e-5)

outequils1_w0pressure, inds = run_Δl_Δr_calculator(equilas11, 3, 1, 1e-4, 8, 1e-5;run_zero_pressure=true)

i=3

print(outequils1_w0pressure[i].Δprimes.Δprimezero,outequils1_w0pressure[i].Δprimes.Δprime)

run_Δl_Δr_calculator(equilas11[inds[1:3]], 3, 1, 1e-4, 8, 1e-5;plot_soln_equil=true)


run_Δl_Δr_calculator(equilas11[inds[1:3]], 3, 1, 1e-4, 1, 1e-5;plot_soln_equil=true)


plot_equil(equilas11[inds[argmax(deltaprimes)]])
plot_equil(equilas11[inds[argmin(deltaprimes)]])

test_Suydam(equilas11[inds[argmin(deltaprimes)]].Bt, equilas11[inds[argmin(deltaprimes)]].q, equilas11[inds[argmin(deltaprimes)]].dpdr, equilas11[inds[argmin(deltaprimes)]].rb;plotresults=true) 

run_Δl_Δr_calculator([equilas11[inds[argmax(deltaprimes)]]], 3, 1, 1e-5, 1, 1e-2;plot_soln_equil=true)
run_Δl_Δr_calculator([equilas11[8]], 3, 1, 1e-5, 1, 1e-2;plot_soln_equil=true)

##########################################################################################################################################################
#Mass Generating Equilibria:
##########################################################################################################################################################
m=3
n=1
nmax=1 #important for stability

"""
m=3
n=1
nmax=1 #important for stability

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,1000000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,1000000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure2.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,1000000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure3.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,1000000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure4.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,1000000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure5.jld2" resistive_equils_store equils_and_inds
end
"""

##########################################################################################################################################################
#Mass Generating Equilibria:
##########################################################################################################################################################
m=3
n=1
nmax=1 #important for stability
"""

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 80

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,250; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure6.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 80

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,250; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure7.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 80

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,250; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure8.jld2" resistive_equils_store equils_and_inds
end

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 80

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,250; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure9.jld2" resistive_equils_store equils_and_inds
end
"""

##########################################################################################################################################################
#Mass Generating Equilibria:
##########################################################################################################################################################
m=3
n=1
nmax=1 #important for stability

"""
cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure6.jld2" resistive_equils_store equils_and_inds
end
cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure7.jld2" resistive_equils_store equils_and_inds
end
cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_store = []
equils_and_inds = []

loop_big_batch = 20

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,100000,1000; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(resistive_equils_store, resistive_equils)
    push!(equils_and_inds, (equils_loop, inds))

    @save "n3_study_fin_pressure8.jld2" resistive_equils_store equils_and_inds
end
"""


##########################################################################################################################################################
#Analysis of n=3:
##########################################################################################################################################################


cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils = ResistiveEquilibrium[]
@load "n3_study_fin_pressure.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils,i)
end
@load "n3_study_fin_pressure3.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils,i)
end
@load "n3_study_fin_pressure5.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils,i)
end
@load "n3_study_fin_pressure6.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils,i)
end

deltaprimes = [i.Δprimes.Δprime for i in resistive_equils]
Jts = [i.equilibrium.Jt for i in resistive_equils]
Jps = [i.equilibrium.Jp for i in resistive_equils]
rss = [i.Δprimes.rs for i in resistive_equils]
Jtots = [total_plasma_current(i.equilibrium.Jt,i.equilibrium.rb) for i in resistive_equils]
Jtgrad_at_rs = [ScrewPinchAsymptoticMatching.gradient(Jts[i],rss[i]) for i in 1:length(Jts)]
Jpgrad_at_rs = [ForwardDiff.derivative(Jps[i],rss[i]) for i in 1:length(Jps)]
Jt_at_rs = [Jts[i](rss[i]) for i in 1:length(Jts)]
shear_at_rs = [ForwardDiff.derivative(resistive_equils[i].equilibrium.q,rss[i]) for i in 1:length(Jts)]
p_at_rs = [mu0*resistive_equils[i].equilibrium.p(rss[i]) for i in 1:length(resistive_equils)]
dpdr_at_rs = [resistive_equils[i].equilibrium.dpdr(rss[i]) for i in 1:length(resistive_equils)]



plotrvec = range(0.000001,rb,200)
grad_J_in = [(minimum(Jts[i].(plotrvec))-Jts[i](resistive_equils[i].equilibrium.rb))/(plotrvec[argmin(Jts[i].(plotrvec))]-resistive_equils[i].equilibrium.rb)/(minimum(Jts[i].(plotrvec))) for i in 1:length(Jts)]
grad_J_out = [(Jts[i](resistive_equils[i].equilibrium.rb/2)-minimum(Jts[i].(plotrvec)))/(resistive_equils[i].equilibrium.rb/2-plotrvec[argmin(Jts[i].(plotrvec))])/(Jts[i](resistive_equils[i].equilibrium.rb/2)) for i in 1:length(Jts)]

filtinds = filter(i->!isnan(grad_J_in[i]), 1:length(grad_J_in))

#Optax
#Autodiff to find min points
    #Throw to an optimiser

histogram(deltaprimes; xlabel="Δ'", label=false)
histogram(rss; xlabel="3/1 rational surface location", label=false)
histogram(Jtots; xlabel="total toroidal current (A)", label=false)

histogram2d(rss,deltaprimes; xlabel="3/1 rational surface location (m)",ylabel = "Δ'")
cor(rss,deltaprimes)
histogram2d(Jtgrad_at_rs,deltaprimes; xlabel="∇Jt (A/m^3) at rs",ylabel = "Δ'")
cor(Jtgrad_at_rs,deltaprimes)

histogram2d(Jpgrad_at_rs,deltaprimes; xlabel="∇Jp (A/m^3) at rs",ylabel = "Δ'")
cor(Jpgrad_at_rs,deltaprimes)

histogram2d(shear_at_rs,deltaprimes; xlabel="∇q at rs",ylabel = "Δ'")
cor(shear_at_rs,deltaprimes)
histogram2d(Jtots,deltaprimes; xlabel="Total current (A)",ylabel = "Δ'")
cor(Jtots,deltaprimes)
histogram2d(Jt_at_rs,deltaprimes; xlabel="Jt at rs (A/m^2)",ylabel = "Δ'")
cor(Jt_at_rs,deltaprimes)
histogram2d(p_at_rs,deltaprimes; xlabel="Pressure at rs (sqrt(T))",ylabel = "Δ'")
cor(p_at_rs,deltaprimes) 
histogram2d(dpdr_at_rs,deltaprimes; xlabel="Pressure gradient at rs (sqrt(T)/m)",ylabel = "Δ'")
cor(dpdr_at_rs,deltaprimes)

histogram2d(grad_J_in[filtinds],deltaprimes[filtinds]; xlabel="∇J/J in",ylabel = "Δ'", title="∇J/Jin correlation with Δ'")
cor(grad_J_in[filtinds],deltaprimes[filtinds])
histogram2d(grad_J_out,deltaprimes; xlabel="∇J/J out",ylabel = "Δ'", title="∇J/Jout correlation with Δ'")
cor(grad_J_out,deltaprimes)


sorted_inds = sortperm(deltaprimes)
deltaprimes[sorted_inds[1:10]]
deltaprimes[sorted_inds[end-100:end]]

plotrvec = range(0.000001,rb,200)
for (io,i) in enumerate(sorted_inds[1:10])
    if io==1 
        plot(plotrvec,Jts[i].(plotrvec);label=false)
        display(vline!([rss[i]];label=false))
    else
        plot!(plotrvec,Jts[i].(plotrvec);title = "10 most stable current profiles",label=false, xlabel="r (m)", ylabel="Jt current density (A/m^2)")
        display(vline!([rss[i]];label=false))
    end
end

deltaprimes[sorted_inds[1:10]]


plotrvec = range(0.000001,rb,200)
for (io,i) in enumerate(sorted_inds[end-10:end])
    if io==1 
        plot(plotrvec,Jts[i].(plotrvec);label=false)
        display(vline!([rss[i]];label=false))
    else
        plot!(plotrvec,Jts[i].(plotrvec);title = "10 most unstable current profiles",label=false, xlabel="r (m)", ylabel="Jt current density (A/m^2)")
        display(vline!([rss[i]];label=false))
    end
end
deltaprimes[sorted_inds[end-10:end]]

plotrvec = range(0.000001,rb,200)
for (io,i) in enumerate(sorted_inds[14000:14200])
    if io==1 
        plot(plotrvec,Jts[i].(plotrvec);label=false)
        display(vline!([rss[i]];label=false))
    else
        plot!(plotrvec,Jts[i].(plotrvec);title = "100 current profiles with Δ' ~ -1",label=false, xlabel="r (m)", ylabel="Jt current density (A/m^2)")
        display(vline!([rss[i]];label=false))
    end
end
deltaprimes[sorted_inds[14000:14200]]

plotrvec = range(0.000001,rb,200)
for (io,i) in enumerate(sorted_inds[25200:25300])
    if io==1 
        plot(plotrvec,Jts[i].(plotrvec);label=false)
        display(vline!([rss[i]];label=false))
    else
        plot!(plotrvec,Jts[i].(plotrvec);title = "100 current profiles with Δ' ~ 1",label=false, xlabel="r (m)", ylabel="Jt current density (A/m^2)")
        display(vline!([rss[i]];label=false))
    end
end
deltaprimes[sorted_inds[25200:25300]]

plotrvec = range(0.000001,rb,200)
for (io,i) in enumerate(sorted_inds[end-2000:end-1900])
    if io==1 
        plot(plotrvec,Jts[i].(plotrvec);label=false)
        display(vline!([rss[i]]))
    else
        plot!(plotrvec,Jts[i].(plotrvec);label=false)
        display(vline!([rss[i]];label=false))
    end
end

deltaprimes[sorted_inds[end-2000:end-1900]]

##########################################################################################################################################################
#Mass generating of m=2,3 for n=1: Let's get rockin
##########################################################################################################################################################

mu0 = ScrewPinchAsymptoticMatching.mu0

if true #These are a mix of ARC and SPARC parameters 
    q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
    rs0 = 0.7 #Controls qedge to be around 3.3 (SPARC is 3.4)
    R0 = 3.0 #Major radius
    ν = 1.0  #Current peaking is fine, it'll change anyway 
                #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
    xb = 1/rs0 #Sets minor radius to 1.0
    Bp0 = 2.5 #Sets Bt to around 12T (~ARC/SPARC) 
    β = 1e-2 #This is a SPARC number #14; #Gives a true 0-pressure Chandra Equil

    m = 2;
    ms = [2,3,4,5]
    n = 1;
    k = ScrewPinchAsymptoticMatching.k_(n,R0)
    #c0 = 1;
end

Chandra_Equil_qvar(plot_equil) = q0 -> Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=plot_equil, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))
#Furth_Equil_qvar(plot_equil) = q0 -> Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=plot_equil, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))
#Furth_Equil_qvar(true)(1.1)

Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Chandra_Equil_qvar(true)(1.1)
rs = find_rs(q,m,n,rb,q0,ν,rs0)

total_plasma_current(Jt,1.0)

########################### Current generation inputs
Jtotmax = 1.1*total_plasma_current(Jt,rb)
Jtotmin = 0.9*total_plasma_current(Jt,rb)
Jtotrange=[Jtotmax, Jtotmin]
J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)]
Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]
Jt_upper_bound=maximum(J0bounds)
knotmax=10
knotmin=5
dirichlet_alpha=10
maxgrad_width = rb/2
max_ref=Jt(0.0)
wobbles=3


r0 = 1e-2;
del=1e-5;
nmax=1;

f_gen_equilibria(batch_size, num_clean) = gen_n_clean_equilibria(Jt_upper_bound,p,rb,batch_size,num_clean; dpdr=dpdr,Bt0=Bt(0.0), R0=R0, ideal_verbose=false, qtest=m/n, 
Jtot_range = Jtotrange, use_fine=true, use_coarse=true, monotonic=false, wobbles=wobbles, max_ref=max_ref,maxgrad_width = maxgrad_width, maxbatches=100, 
J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha, knotmax=knotmax, knotmin=knotmin);

f_gen_equilibria_6_knots(batch_size, num_clean) = gen_n_clean_equilibria(Jt_upper_bound,p,rb,batch_size,num_clean; dpdr=dpdr,Bt0=Bt(0.0), R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange, use_fine=true, use_coarse=true, monotonic=false, wobbles=wobbles, max_ref=max_ref,maxgrad_width = maxgrad_width, maxbatches=100, 
J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha, knotmax=knotmax, knotmin=6);

f_gen_equilibria_7_knots(batch_size, num_clean) = gen_n_clean_equilibria(Jt_upper_bound,p,rb,batch_size,num_clean; dpdr=dpdr,Bt0=Bt(0.0), R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange, use_fine=true, use_coarse=true, monotonic=false, wobbles=wobbles, max_ref=max_ref,maxgrad_width = maxgrad_width, maxbatches=100, 
J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha, knotmax=knotmax, knotmin=7);

f_gen_equilibria_8_knots(batch_size, num_clean) = gen_n_clean_equilibria(Jt_upper_bound,p,rb,batch_size,num_clean; dpdr=dpdr,Bt0=Bt(0.0), R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange, use_fine=true, use_coarse=true, monotonic=false, wobbles=wobbles, max_ref=max_ref,maxgrad_width = maxgrad_width, maxbatches=100, 
J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha, knotmax=knotmax, knotmin=8);

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
path="/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase"

batch_size=100000
batch_size6=4*batch_size
batch_size7=100*batch_size
batch_size8=1000*batch_size
num_clean=1000
loop_big_batch=5

gen_equil_run_Δl_Δr(f_gen_equilibria_6_knots, batch_size6, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_6knot_1.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_6_knots, batch_size6, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_6knot_2.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_7_knots, batch_size7, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_7knot_1.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_7_knots, batch_size7, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_7knot_2.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_8_knots, batch_size8, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_8knot_1.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_6_knots, batch_size6, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_6knot_3.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_7_knots, batch_size7, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_7knot_3.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_8_knots, batch_size8, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_8knot_2.jld2", loop_big_batch=loop_big_batch, path=path)
gen_equil_run_Δl_Δr(f_gen_equilibria_8_knots, batch_size8, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_8knot_3.jld2", loop_big_batch=loop_big_batch, path=path)

#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure1.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure2.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure3.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure4.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure5.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure6.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure7.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure8.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure9.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure10.jld2", loop_big_batch=loop_big_batch, path=path)

########################### Analysis
cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
resistive_equils_n2 = ResistiveEquilibrium[]
@load "n2_study_fin_pressure1.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure2.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end

########################### Removing scenarios without n=3 rs (just the tail end of Jtot distribution)
###########################
knotnums = [length(i.equilibrium.Jt.xs) for i in resistive_equils_n2]
deltaprimes_n2 = [i.Δprimes[1].Δprime for i in resistive_equils_n2]
#deltaprimes_n2 = filter(i -> resistive_equils_n2[i].Δprimes[1].m==2, 1:length(deltaprimes_n2))
deltaprimes_n3_inds = filter(i -> (length(resistive_equils_n2[i].Δprimes)>1 && resistive_equils_n2[i].Δprimes[2].m==3), 1:length(resistive_equils_n2))
deltaprimes_no_n3_inds = filter(i -> (length(resistive_equils_n2[i].Δprimes)==1), 1:length(resistive_equils_n2))
deltaprimes_n3 = [i.Δprimes[2].Δprime for i in resistive_equils_n2[deltaprimes_n3_inds]]

Jtots = [total_plasma_current(i.equilibrium.Jt,i.equilibrium.rb) for i in resistive_equils_n2]
mean(Jtots[deltaprimes_n3_inds])
mean(Jtots[deltaprimes_no_n3_inds])
plot_current_profiles(resistive_equils_n2[deltaprimes_no_n3_inds])
plot_current_profiles(resistive_equils_n2[deltaprimes_n3_inds][1:300])

resistive_equils_n2_n3 = filter(i -> (length(i.Δprimes)>1), resistive_equils_n2)
deltaprimes_n2 = [i.Δprimes[1].Δprime for i in resistive_equils_n2_n3]
deltaprimes_n3 = [i.Δprimes[2].Δprime for i in resistive_equils_n2_n3]
deltaprime_diff = deltaprimes_n2-deltaprimes_n3

########################### Removing n2 delta' outliers (they're obvious)
sorted_inds_n2 = sortperm(deltaprimes_n2)
sorted_inds_n3 = sortperm(deltaprimes_n3)

deltaprimes_n2[sorted_inds_n2[20:40]]
deltaprimes_n2[sorted_inds_n2[20:28]]
deltaprimes_n2[sorted_inds_n2[28]]

resistive_equils_n2_n3_clean = filter(i -> (i.Δprimes[1].Δprime>deltaprimes_n2[sorted_inds_n2[28]]), resistive_equils_n2_n3)
deltaprimes_n2 = [i.Δprimes[1].Δprime for i in resistive_equils_n2_n3_clean]
deltaprimes_n3 = [i.Δprimes[2].Δprime for i in resistive_equils_n2_n3_clean]
deltaprime_diff = deltaprimes_n2-deltaprimes_n3

########################### n2, n3 histograms
histogram(deltaprimes_n2)
histogram(filter(x->x>-200,deltaprimes_n3))
n_3_inds = filter(i -> (deltaprimes_n3[i]>-200), 1:length(deltaprimes_n3))

Jts = [i.equilibrium.Jt for i in resistive_equils_n2_n3_clean]
Jps = [i.equilibrium.Jp for i in resistive_equils_n2_n3_clean]
rss_n2 = [i.Δprimes[1].rs for i in resistive_equils_n2_n3_clean]
rss_n3 = [i.Δprimes[2].rs for i in resistive_equils_n2_n3_clean]
Jtots = [total_plasma_current(i.equilibrium.Jt,i.equilibrium.rb) for i in resistive_equils_n2_n3_clean]
Jtgrad_at_rs_n2 = [ScrewPinchAsymptoticMatching.gradient(Jts[i],rss_n2[i]) for i in 1:length(Jts)]
Jtgrad_at_rs_n3 = [ScrewPinchAsymptoticMatching.gradient(Jts[i],rss_n3[i]) for i in 1:length(Jts)]
Jpgrad_at_rs_n2 = [ForwardDiff.derivative(Jps[i],rss_n2[i]) for i in 1:length(Jps)]
Jpgrad_at_rs_n3 = [ForwardDiff.derivative(Jps[i],rss_n3[i]) for i in 1:length(Jps)]
Jt_at_rs_n2 = [Jts[i](rss_n2[i]) for i in 1:length(Jts)]
Jt_at_rs_n3 = [Jts[i](rss_n3[i]) for i in 1:length(Jts)]
shear_at_rs_n2 = [ForwardDiff.derivative(resistive_equils_n2_n3_clean[i].equilibrium.q,rss_n2[i]) for i in 1:length(Jts)]
shear_at_rs_n3 = [ForwardDiff.derivative(resistive_equils_n2_n3_clean[i].equilibrium.q,rss_n3[i]) for i in 1:length(Jts)]
p_at_rs_n2 = [mu0*resistive_equils_n2_n3_clean[i].equilibrium.p(rss_n2[i]) for i in 1:length(resistive_equils_n2_n3_clean)]
p_at_rs_n3 = [mu0*resistive_equils_n2_n3_clean[i].equilibrium.p(rss_n3[i]) for i in 1:length(resistive_equils_n2_n3_clean)]
dpdr_at_rs_n2 = [resistive_equils_n2_n3_clean[i].equilibrium.dpdr(rss_n2[i]) for i in 1:length(resistive_equils_n2_n3_clean)]
dpdr_at_rs_n3 = [resistive_equils_n2_n3_clean[i].equilibrium.dpdr(rss_n3[i]) for i in 1:length(resistive_equils_n2_n3_clean)]

histogram(rss_n2; xlabel="2/1 rational surface location", label=false, xlims=(0.0,1.0))
histogram(rss_n3; xlabel="3/1 rational surface location", label=false, xlims=(0.0,1.0))
histogram2d(rss_n2[negm2posm3_inds],rss_n3[negm2posm3_inds],xlabel="2/1 rational surface location (m)",ylabel = "3/1 rational surface location (m)")

histogram(Jtots; xlabel="total toroidal current (A)", label=false)

histogram2d(rss_n2,deltaprimes_n2; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
cor(rss_n2,deltaprimes_n2)
histogram2d(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds]; xlabel="3/1 rational surface location (m)",ylabel = "Δ' 3/1")
cor(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds])
histogram2d(Jtgrad_at_rs,deltaprimes; xlabel="∇Jt (A/m^3) at rs",ylabel = "Δ'")
cor(Jtgrad_at_rs,deltaprimes)

histogram2d(Jpgrad_at_rs,deltaprimes; xlabel="∇Jp (A/m^3) at rs",ylabel = "Δ'")
cor(Jpgrad_at_rs,deltaprimes)

histogram2d(shear_at_rs,deltaprimes; xlabel="∇q at rs",ylabel = "Δ'")
cor(shear_at_rs,deltaprimes)
histogram2d(Jtots,deltaprimes; xlabel="Total current (A)",ylabel = "Δ'")
cor(Jtots,deltaprimes)
histogram2d(Jt_at_rs,deltaprimes; xlabel="Jt at rs (A/m^2)",ylabel = "Δ'")
cor(Jt_at_rs,deltaprimes)
histogram2d(p_at_rs,deltaprimes; xlabel="Pressure at rs (sqrt(T))",ylabel = "Δ'")
cor(p_at_rs,deltaprimes) 
histogram2d(dpdr_at_rs,deltaprimes; xlabel="Pressure gradient at rs (sqrt(T)/m)",ylabel = "Δ'")
cor(dpdr_at_rs,deltaprimes)

histogram2d(grad_J_in[filtinds],deltaprimes[filtinds]; xlabel="∇J/J in",ylabel = "Δ'", title="∇J/Jin correlation with Δ'")
cor(grad_J_in[filtinds],deltaprimes[filtinds])
histogram2d(grad_J_out,deltaprimes; xlabel="∇J/J out",ylabel = "Δ'", title="∇J/Jout correlation with Δ'")
cor(grad_J_out,deltaprimes)

########################### n2, n3 extreme Cases
sorted_inds_n2 = sortperm(deltaprimes_n2)
sorted_inds_n3 = sortperm(deltaprimes_n3)

deltaprimes_n2[sorted_inds_n2[1:10]]
deltaprimes_n3[sorted_inds_n3[1:10]]

deltaprimes_n2[sorted_inds_n2[end-10:end]]
deltaprimes_n3[sorted_inds_n3[end-10:end]]

plot_current_profiles(resistive_equils_n2_n3_clean[sorted_inds_n2[end-10:end]];m_ind=1)
plot_equil(resistive_equils_n2_n3_clean[sorted_inds_n2[end-10:end]];case=1, Jt_ref = 1.1*Jt_upper_bound)

########################### Cases where n2 is stable and n3 is unstable (rare)
negm2posm3_inds = filter(i -> deltaprimes_n2[i]<0.0&& deltaprimes_n3[i]>0.0, 1:length(resistive_equils_n2_n3_clean))
negm2posm3_equils = filter(i -> i.Δprimes[1].Δprime<0.0 && i.Δprimes[2].Δprime>0.0, resistive_equils_n2_n3_clean)
plot_current_profiles(negm2posm3_equils)

########################### Cases where n_2 rs values aren't in super low-rs region:
histogram2d(rss_n2,deltaprimes_n2; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
cor(rss_n2,deltaprimes_n2)
histogram2d(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds]; xlabel="3/1 rational surface location (m)",ylabel = "Δ' 3/1")
cor(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds])

non_rs_inds = filter(i -> rss_n2[i]>0.3, 1:length(resistive_equils_n2_n3_clean))
non_rs_equil = filter(i -> i.Δprimes[1].rs>0.3, resistive_equils_n2_n3_clean)

########################### no_rs extreme cases:
deltaprimes_n2_nrs= [i.Δprimes[1].Δprime for i in non_rs_equil]
deltaprimes_n3_nrs = [i.Δprimes[2].Δprime for i in non_rs_equil]
sorted_inds_n2_nrs = sortperm(deltaprimes_n2_nrs)
sorted_inds_n3_nrs = sortperm(deltaprimes_n3_nrs)

Jts_nrs = [i.equilibrium.Jt for i in non_rs_equil]
Jps = [i.equilibrium.Jp for i in non_rs_equil]
rss_n2_nrs = [i.Δprimes[1].rs for i in non_rs_equil]
rss_n3_nrs = [i.Δprimes[2].rs for i in non_rs_equil]
Jtgrad_at_rs_n2_nrs = [ScrewPinchAsymptoticMatching.gradient(Jts_nrs[i],rss_n2_nrs[i]) for i in 1:length(Jts_nrs)]
Jtgrad_at_rs_n3_nrs = [ScrewPinchAsymptoticMatching.gradient(Jts_nrs[i],rss_n3_nrs[i]) for i in 1:length(Jts_nrs)]

plot_current_profiles(non_rs_equil[sorted_inds_n2_nrs[2]];m_ind=1)

histogram2d(Jtgrad_at_rs_n2[non_rs_inds],deltaprimes_n2[non_rs_inds]; xlabel="∇Jt (A/m^3) at rs",ylabel = "Δ'")
cor(Jtgrad_at_rs_n2[non_rs_inds],deltaprimes_n2[non_rs_inds])
corspearman(Jtgrad_at_rs_n2[non_rs_inds],deltaprimes_n2[non_rs_inds])

histogram(deltaprimes_n2_nrs)
histogram(deltaprimes_n3_nrs)

#What makes the current profile more stable... in a cylinder, its the instantaneous current gradient w.r.t. the bingo bongo


########################### TESTING Having a minimum:
plotrvec = range(0.000001,rb,200)
grad_J_in = [(minimum(Jts_nrs[i].(plotrvec))-Jts_nrs[i](non_rs_equil[i].equilibrium.rb))/(abs(plotrvec[argmin(Jts_nrs[i].(plotrvec))]-non_rs_equil[i].equilibrium.rb))/(abs(minimum(Jts_nrs[i].(plotrvec)))) for i in 1:length(Jts_nrs)]
grad_J_out = [(Jts_nrs[i](non_rs_equil[i].equilibrium.rb/2)-minimum(Jts_nrs[i].(plotrvec)))/(non_rs_equil[i].equilibrium.rb/2-plotrvec[argmin(Jts_nrs[i].(plotrvec))])/(Jts_nrs[i](non_rs_equil[i].equilibrium.rb/2)) for i in 1:length(Jts_nrs)]

filtinds = filter(i->!isnan(grad_J_in[i]), 1:length(grad_J_in))

histogram2d(grad_J_in[filtinds],deltaprimes_n2_nrs[filtinds]; xlabel="∇Jt (A/m^3) at rs",ylabel = "Δ'", xlims = (-10,10))
cor(grad_J_in[filtinds],deltaprimes_n2_nrs[filtinds])
corspearman(grad_J_in[filtinds],deltaprimes_n2_nrs[filtinds])

########################### TESTING being near global min (they don't get near cos q2 not near edge. Should repeat with q=3):
min_closeness = 0.15 #Norm to 1
sep_from_edge = 0.05
near_well_equils1 = filter(i -> (abs(i.Δprimes[1].rs-plotrvec[argmin(i.equilibrium.Jt.(plotrvec))])<min_closeness*i.equilibrium.rs) && (abs(plotrvec[argmin(i.equilibrium.Jt.(plotrvec))]-i.equilibrium.rs)>sep_from_edge*i.equilibrium.rs), resistive_equils_n2_n3_clean)
deltaprimes_n2_near_well1= [i.Δprimes[1].Δprime for i in near_well_equils1]
plot_current_profiles(near_well_equils1[1:50];m_ind=1)


non_rs_inds = filter(i -> rss_n2[i]>0.3, 1:length(resistive_equils_n2_n3_clean))
non_rs_equil = filter(i -> i.Δprimes[1].rs>0.3, resistive_equils_n2_n3_clean)

########################### TESTING being near local min:

#Just any well isn't important... 
    #What if rs is in a well
    min_closeness = 0.01 #Norm to 1
    sep_from_edge = 0.05 #Norm to 1

well_equils, a_to_well, well_width, well_gradient, Jt_double_derivative, well_inds = rs_near_local_min(non_rs_equil,min_closeness,sep_from_edge;return_well_info=true)
plot_current_profiles(well_equils[1:100];m_ind=1)

min_closeness = 0.03 #Norm to 1
right_in_well_inds = filter(i -> a_to_well[i]<min_closeness, 1:length(a_to_well))
histogram2d(well_gradient[right_in_well_inds],[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])
cor(well_gradient[right_in_well_inds],[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])
corspearman(float.(well_gradient[right_in_well_inds]),[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])

a_to_well[right_in_well_inds][1:10]
well_gradient[right_in_well_inds][1:10]
plot_current_profiles(well_equils[right_in_well_inds[6]];m_ind=1)

for i in right_in_well_inds[1:10]
    print("Well depth: ",round(well_gradient[i];digits=2), " Δ': ",well_equils[i].Δprimes[1].Δprime,"\n")
end

near_well_equils = filter(i -> rs_near_local_min(i,min_closeness,sep_from_edge;m_ind=1), non_rs_equil)
plot_current_profiles(near_well_equils[1:100];m_ind=1)

deltaprimes_n2_near_well = [i.Δprimes[1].Δprime for i in near_well_equils]





########################### TESTs
#randomJt(Jt_upper_bound, 100, rb; plot_profs=true, knotmax=knotmax, knotmin=knotmin, J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha)

Jts = gen_n_clean_Jts(Jt_upper_bound,rb,500,10,generate_cleaners(rb ; 
Jtot_range = Jtotrange, use_fine=true, use_coarse=true, monotonic=false, wobbles=wobbles, max_ref=max_ref,maxgrad_width = maxgrad_width); maxbatches=100, 
J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha, knotmax=knotmax, knotmin=knotmin)
plot_profiles(Jts, rb;ylims = (0.0,1.4*Jt(0.0)))

test_equils = gen_n_clean_equilibria(Jt_upper_bound,p,rb,500,10; dpdr=dpdr,Bt0=Bt(0.0), R0=R0, ideal_verbose=true, qtest=m/n, 
Jtot_range = Jtotrange, use_fine=true, use_coarse=true, monotonic=false, wobbles=wobbles, max_ref=max_ref,maxgrad_width = maxgrad_width, maxbatches=100, 
J0bounds=J0bounds, Jedgebounds=Jedgebounds, dirichlet_alpha=dirichlet_alpha, knotmax=knotmax, knotmin=knotmin);

plot_current_profiles(test_equils)
plot_equil(test_equils)

########################### Calculating Delta's
equilibria = test_equils
resistive_equils, inds = run_Δl_Δr_calculator(test_equils, ms, n, r0, nmax, del;report_err=true)
########################### Test function run
batch_size=5000
num_clean=10
gen_equil_run_Δl_Δr(f_gen_equilibria, batch_size, num_clean, ms, n, r0, nmax, del; filename = "testdog1", loop_big_batch=3, path="/Users/sbenjamin/Desktop/PHD")
cd("/Users/sbenjamin/Desktop/PHD")
@load "testdog1" resistive_equils_store equils_and_inds
plot_current_profiles(resistive_equils_store[1])
plot_equil(resistive_equils_store[1])
########################### Big function run

"""

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
@load "n3_study_test1.jld2" resistive_equils_store equils_and_inds

plot_equil(resistive_equils_store[1][1].equilibrium)
resistive_equils_store[1][1].Δprimes.Δprimezero


length(bigstore[1][1])

run_Δl_Δr_calculator([bigstore[1][1][1]], m, n, 1e-4, nmax, 1e-5; report_err=true)
run_Δl_Δr_calculator([bigstore[1][1][3]], m, n, 1e-4, nmax, 1e-5; report_err=true)

q_one_rs(bigstore[1][1][1], 3/1)

equilas1_=gen_n_clean_equilibria(Jtotmax,p,rb,5000,20;dpdr=dpdr, Bt0=0.3, Jtot_range = Jtotrange, maxgrad_width = rb/7,use_fine=false,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange)) #all fail suydm
equilas2_=gen_n_clean_equilibria(Jtotmax,pc,rb,500,20;dpdr=dpdrc, Bt0=0.3,Jtot_range = Jtotrange, maxgrad_width = rb/7,use_fine=false,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange)) #all fail suydm
equilas3_=gen_n_clean_equilibria(Jtotmax,p_s50,rb,500,20; Bt0=0.3,Jtot_range = Jtotrange, maxgrad_width = rb/7,use_fine=false,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange))

##########################################################################################################################################################
#Make Equilibria from scratch:
##########################################################################################################################################################




batch_size=1000
num_clean=3

gen_n_clean_equilibria(Jtotmax,p_s50,rb,batch_size,num_clean;max_batches=10,area_integrated_max=maximum(Jtotrange))

gen_n_clean_equilibria(Jtotmax,p,rb,batch_size,num_clean;max_ref=Jt(0.000000001), dpdr=dpdr)



#Jts = gen_n_clean_Jts(Jtotmax,rb,5000,80,generate_cleaners(rb ; Jtot_range = Jtotrange,maxgrad_width = rb/100,use_fine=false,use_coarse=true, monotonic=false, wobbles=1, area_integrated_max=maximum(Jtotrange));
#                maxbatches=10, J0bounds=[0.8*Jt(0.0),1.2*Jt(0.0)], Jedgebounds=[1.2*Jt(2.0),0.8*Jt(2.0)])
#plot_profiles(Jts, rb;ylims = (0.0,1.4*Jt(0.0))) 

#ps = gen_n_clean_pressure_profiles(1e20,rb,500,10,generate_cleaners(rb, [Jtotmax, Jtotmin]; maxgrad_width = rb/100,use_fine=false,use_coarse=true, monotonic=false, wobbles=1);
#maxbatches=10, J0bounds=[0.8*Jt(0.0),1.2*Jt(0.0)], Jedgebounds=[1.2*Jt(2.0),0.8*Jt(2.0)])

#plot_profiles(randomJt(Jtotmax, 1000, rb;J0bounds=5e5, Jedgebounds=1e5),rb)

#plot_profiles(randomJt(Jtotmax, 1000, rb;J0bounds=5e5, Jedgebounds=1e5),rb)
"""
print("cat")

