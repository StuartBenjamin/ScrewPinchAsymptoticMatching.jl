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
equilas2 = gen_clean_equilibria(Jts,p; dpdr=dpdr, Bt0=10, R0=3, m1ncap=20, ideal_verbose=true) 
equilas3 = gen_clean_equilibria(Jts[1:3],p_s50; Bt0=10, R0=R0, ideal_verbose=true, qtest=3.0);
plot_Suydam(equila1[3])

gen_clean_equilibria(Jts[1:3],p_s50; Bt0=10, R0=3, ideal_verbose=true,maxq=0.01) 
plot_equil(equila1[3])

plot_equil(equilas1)
plot_equil(equilas3)


equilas11=gen_n_clean_equilibria(Jtotmax,p,rb,5000,50; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=3.0, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);


plot_equil(equilas11)
deltaprimes, outmatrix, inds = run_Δl_Δr_calculator(equilas11, 3, 1, 1e-4, 8, 1e-5)


run_Δl_Δr_calculator(equilas11[inds[1:3]], 3, 1, 1e-4, 8, 1e-5;plot_soln_equil=true)


run_Δl_Δr_calculator(equilas11[inds[1:3]], 3, 1, 1e-4, 8, 1e-5;plot_soln_equil=true)


plot_equil(equilas11[inds[argmax(deltaprimes)]])
plot_equil(equilas11[inds[argmin(deltaprimes)]])

test_Suydam(equilas11[inds[argmin(deltaprimes)]].Bt, equilas11[inds[argmin(deltaprimes)]].q, equilas11[inds[argmin(deltaprimes)]].dpdr, equilas11[inds[argmin(deltaprimes)]].rb;plotresults=true) 

run_Δl_Δr_calculator([equilas11[inds[argmax(deltaprimes)]]], 3, 1, 1e-5, 1, 1e-2;plot_soln_equil=true)
run_Δl_Δr_calculator([equilas11[8]], 3, 1, 1e-5, 1, 1e-2;plot_soln_equil=true)

##########################################################################################################################################################
#Mass Generating Equilibria:
##########################################################################################################################################################

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")

equils_store = []
outputs_finitep_store = []

bigstore = []

loop_big_batch = 5

m=3
n=1
nmax=1 #important for stability

for i in 1:loop_big_batch
    equils_loop=gen_n_clean_equilibria(Jtotmax,p,rb,10000,10; dpdr=dpdr,Bt0=10, R0=R0, ideal_verbose=false, qtest=m/n, Jtot_range = Jtotrange,maxgrad_width = rb/7,use_fine=true,use_coarse=true, monotonic=false, wobbles=3, area_integrated_max=1.0*maximum(Jtotrange), 
maxbatches=100, J0bounds=[1.2*Jt(0.0),0.5*Jt(0.0)],Jedgebounds=[0.3*Jt(0.0),0.0*Jt(0.0)]);

    deltaprimes, outmatrix, inds = run_Δl_Δr_calculator(equils_loop, m, n, 1e-4, nmax, 1e-5)

    push!(equils_store, equils_loop)
    push!(outputs_finitep_store, (deltaprimes, outmatrix, inds))
    #@save "n3_study_1.jld2" equils_store outputs_finitep_store 

    push!(bigstore, (equils_loop, deltaprimes, outmatrix, inds))
    @save "n3_study_1.jld2" bigstore
end

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

