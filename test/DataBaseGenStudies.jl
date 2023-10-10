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
using DelimitedFiles
using Printf
using LaTeXStrings

using Pkg
Pkg.add(PackageSpec(url="https://github.com/StuartBenjamin/ScrewPinchAsymptoticMatching.jl"))
using ScrewPinchAsymptoticMatching
ScrewPinchAsymptoticMatching.mod_path=chop(pathof(ScrewPinchAsymptoticMatching);tail=31)

include("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia/ScrewPinchAsymptoticMatching.jl/src/DataBaseGen.jl")
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
ref_equil=Equilibrium(Bp,Bt,q,dpdr,p,r->Jt(r),Jp,rb,R0,rs,rs0)
plot_equil(ref_equil)

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

#gen_equil_run_Δl_Δr(f_gen_equilibria_6_knots, batch_size6, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_6knot_1.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_6_knots, batch_size6, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_6knot_2.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_7_knots, batch_size7, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_7knot_1.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_7_knots, batch_size7, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_7knot_2.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_8_knots, batch_size8, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_8knot_1.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_6_knots, batch_size6, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_6knot_3.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_7_knots, batch_size7, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_7knot_3.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_8_knots, batch_size8, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_8knot_2.jld2", loop_big_batch=loop_big_batch, path=path)
#gen_equil_run_Δl_Δr(f_gen_equilibria_8_knots, batch_size8, num_clean, ms, n, r0, nmax, del; filename = "n2_study_fin_pressure_8knot_3.jld2", loop_big_batch=loop_big_batch, path=path)

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
@load "n2_study_fin_pressure3.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure4.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure5.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure6.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure7.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure_6knot_1.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure_6knot_2.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure_7knot_1.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end
@load "n2_study_fin_pressure_7knot_2.jld2" resistive_equils_store equils_and_inds
for i in resistive_equils_store
    append!(resistive_equils_n2,i)
end

length(resistive_equils_n2)

########################### Removing scenarios without n=3 rs (just the tail end of Jtot distribution -> largest currents cause qedge to drop down)
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

########################### Removing n2 delta' outliers (they're at very small rs, integrator often broken at these rs)

rss_n2 = [i.Δprimes[1].rs for i in resistive_equils_n2_n3]
sorted_inds_n2 = sortperm(deltaprimes_n2)
sorted_inds_n3 = sortperm(deltaprimes_n3)

histogram2d(rss_n2,deltaprimes_n2; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
corspearman(rss_n2,deltaprimes_n2)
cor(rss_n2,deltaprimes_n2)


deltaprimes_n2[sorted_inds_n2[1:90]]
maximum(rss_n2[sorted_inds_n2[1:90]])

deltaprimes_n2[sorted_inds_n2[80:90]]
rss_n2[sorted_inds_n2[80:90]]


iend = 87
deltaprimes_n2[sorted_inds_n2[iend]]
maximum(rss_n2[sorted_inds_n2[1:iend]])

scatter(rss_n2[sorted_inds_n2[1:iend]],deltaprimes_n2[sorted_inds_n2[1:iend]]; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
corspearman(rss_n2[sorted_inds_n2[1:iend]],deltaprimes_n2[sorted_inds_n2[1:iend]])
cor(rss_n2[sorted_inds_n2[1:iend]],deltaprimes_n2[sorted_inds_n2[1:iend]])

resistive_equils_n2_n3_clean = filter(i -> (i.Δprimes[1].Δprime>deltaprimes_n2[sorted_inds_n2[iend]]), resistive_equils_n2_n3)
deltaprimes_n2 = [i.Δprimes[1].Δprime for i in resistive_equils_n2_n3_clean]
deltaprimes_n3 = [i.Δprimes[2].Δprime for i in resistive_equils_n2_n3_clean]
deltaprime_diff = deltaprimes_n2-deltaprimes_n3

########################### Cases where n2 is stable and n3 is unstable (rare)
negm2posm3_inds = filter(i -> deltaprimes_n2[i]<0.0&& deltaprimes_n3[i]>0.0, 1:length(resistive_equils_n2_n3_clean))
negm2posm3_equils = filter(i -> i.Δprimes[1].Δprime<0.0 && i.Δprimes[2].Δprime>0.0, resistive_equils_n2_n3_clean)
plot_current_profiles(negm2posm3_equils)

########################### general analysis, n2, n3 histograms
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
Drs_n2 = [Dr(i)[1] for i in resistive_equils_n2_n3_clean]
Drs_n3 = [Dr(i)[2] for i in resistive_equils_n2_n3_clean]

xs_s_raw = [i.equilibrium.Jt.xs for i in resistive_equils_n2_n3_clean]
knotnum_lengths = length.(xs_s_raw)

rando = rand(1:length(resistive_equils_n2_n3_clean),340)
plot_current_profiles(resistive_equils_n2_n3_clean[rando])
plot_q_profiles(resistive_equils_n2_n3_clean[rando])

histogram(knotnum_lengths; label=false, title="Distribution of radial knot number in database", thickness_scaling=1.0)
histogram(xs_s; label=false, title="Radial distribution of knots in database", thickness_scaling=1.0, bins = 200)
save()
xs_7s_raw = filter(i -> length(i)>6, xs_s_raw)
xs_s = Number[]
xs_7 = Number[]
for i in xs_s_raw
    append!(xs_s, i[2:end-1])
end
for i in xs_7s_raw
    append!(xs_7, i[2:end-1])
end
histogram(xs_s)
histogram(xs_7)

histogram(rss_n2; xlabel="radius (m)", label=false, xlims=(0.0,1.0), bins=200, title="Radial distribution of 2/1 rational surface")
histogram(rss_n3; xlabel="radius (m)", label=false, xlims=(0.0,1.0), bins=200, title="Radial distribution of 3/1 rational surface")
histogram2d(rss_n2[negm2posm3_inds],rss_n3[negm2posm3_inds],xlabel="2/1 rational surface location (m)",ylabel = "3/1 rational surface location (m)")

histogram(Jtots; xlabel="Total Jt (A)", label=false, title="Distribution of total toroidal current", bins=200)
histogram(knotnums)

########################### Important correlation: rs approaching edge of plasma

histogram2d(rss_n2,deltaprimes_n2; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
cor(rss_n2,deltaprimes_n2)
corspearman(rss_n2,deltaprimes_n2)
histogram2d(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds]; xlabel="3/1 rational surface location (m)",ylabel = "Δ' 3/1") 
cor(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds])
corspearman(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds])
#^DFINITELY A CONSEQUENCE OF ENFORCING WE HAVE AN INTERNAL MODE (perturbation = 0 at edge)

########################### n2, n3 extreme Cases
sorted_inds_n2 = sortperm(deltaprimes_n2)
sorted_inds_n3 = sortperm(deltaprimes_n3)

deltaprimes_n2[sorted_inds_n2[1:10]]
deltaprimes_n3[sorted_inds_n3[1:10]]

deltaprimes_n2[sorted_inds_n2[end-10:end]]
deltaprimes_n3[sorted_inds_n3[end-10:end]]

histogram(deltaprimes_n2)

plot_current_profiles(resistive_equils_n2_n3_clean[sorted_inds_n2[end-10:end]];m_ind=1) #THESE ARE EXTREMELY FLAT CURRENT PROFILES pushing rs_n2 towards 0 -> not interesting
plot_equil(resistive_equils_n2_n3_clean[sorted_inds_n2[end-10:end]];case=1, Jt_ref = 1.1*Jt_upper_bound)
plot_current_profiles(resistive_equils_n2_n3_clean[sorted_inds_n2[1:10]];m_ind=1) #THESE ARE the rs_n2 on the hill -> INTERESTING
plot_current_profiles(resistive_equils_n2_n3_clean[sorted_inds_n3[1:10]];m_ind=2) #THESE rs_n3 towards rb stabilising mode -> not interesting
plot_current_profiles(resistive_equils_n2_n3_clean[sorted_inds_n3[end-10:end]];m_ind=2)  #THESE ARE the rs_n3 after long flat period, unstable n3 -> INTERESTING

########################### Cases where n_2 rs values aren't in super low-rs region:
histogram2d(rss_n2,deltaprimes_n2; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
cor(rss_n2,deltaprimes_n2)
histogram2d(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds]; xlabel="3/1 rational surface location (m)",ylabel = "Δ' 3/1")
cor(rss_n3[n_3_inds],deltaprimes_n3[n_3_inds])

non_rs_inds = filter(i -> resistive_equils_n2_n3[i].Δprimes[1].rs>0.3, 1:length(resistive_equils_n2_n3_clean)) #getting rid of any correlation with position for n2
non_rs_equil = filter(i -> i.Δprimes[1].rs>0.3, resistive_equils_n2_n3_clean)

########################### First big thing: local gradient correlates with Delta' for 2/1 (evidence it might be absolute total local current gradient)
deltaprimes_n2_nrs= [i.Δprimes[1].Δprime for i in non_rs_equil]
deltaprimes_n3_nrs = [i.Δprimes[2].Δprime for i in non_rs_equil]
sorted_inds_n2_nrs = sortperm(deltaprimes_n2_nrs)
sorted_inds_n3_nrs = sortperm(deltaprimes_n3_nrs)

Jts_nrs = [i.equilibrium.Jt for i in non_rs_equil]
Jps_nrs = [i.equilibrium.Jp for i in non_rs_equil]
rss_n2_nrs = [i.Δprimes[1].rs for i in non_rs_equil]
rss_n3_nrs = [i.Δprimes[2].rs for i in non_rs_equil]
Jtots_nrs = [total_plasma_current(i.equilibrium.Jt,i.equilibrium.rb) for i in non_rs_equil]

histogram(deltaprimes_n2_nrs;xlabel="Δ' 2/1", label=false, bins=1000, title = "Δ' 2/1 distribution")

histogram2d(rss_n2_nrs,deltaprimes_n2_nrs; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
cor(rss_n2_nrs,deltaprimes_n2_nrs)
corspearman(rss_n2_nrs,deltaprimes_n2_nrs)
n_3_inds_temp = filter(i -> (deltaprimes_n3_nrs[i]>-200), 1:length(deltaprimes_n3_nrs))
histogram2d(rss_n3_nrs[n_3_inds_temp],deltaprimes_n3_nrs[n_3_inds_temp]; xlabel="3/1 rational surface location (m)",ylabel = "Δ' 3/1")
cor(rss_n3_nrs[n_3_inds_temp],deltaprimes_n3_nrs[n_3_inds_temp])
corspearman(rss_n3_nrs[n_3_inds_temp],deltaprimes_n3_nrs[n_3_inds_temp])

Jtgrad_at_rs_n2_nrs = [ScrewPinchAsymptoticMatching.gradient(Jts_nrs[i],rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
Jtgrad_at_rs_n3_nrs = [ScrewPinchAsymptoticMatching.gradient(Jts_nrs[i],rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]
Jpgrad_at_rs_n2_nrs = [ForwardDiff.derivative(Jps_nrs[i],rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
Jpgrad_at_rs_n3_nrs = [ForwardDiff.derivative(Jps_nrs[i],rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]

plot_current_profiles(non_rs_equil[sorted_inds_n2_nrs[1:10]];m_ind=1)
deltaprimes_n2_nrs[sorted_inds_n2_nrs[1:10]]

plot_current_profiles(non_rs_equil[sorted_inds_n3_nrs[1:10]];m_ind=2)
deltaprimes_n3_nrs[sorted_inds_n3_nrs[1:10]]

########################### Terrible: 
plot_current_profiles(non_rs_equil[sorted_inds_n2_nrs[end-10:end]];m_ind=1) 
deltaprimes_n2_nrs[sorted_inds_n2_nrs[end-10:end]]


###GENERAL STABLE
i=8700
plot_current_profiles(non_rs_equil[sorted_inds_n2_nrs[i-100:i]];m_ind=1) 
deltaprimes_n2_nrs[sorted_inds_n2_nrs[i-100:i]]

###GENERAL UNSTABLE
i=95000
plot_current_profiles(non_rs_equil[sorted_inds_n2_nrs[i-100:i]];m_ind=1) 
deltaprimes_n2_nrs[sorted_inds_n2_nrs[i-100:i]]


#Jt n2 big correlation
histogram2d(Jtgrad_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="∇Jt (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(Jtgrad_at_rs_n2_nrs,deltaprimes_n2_nrs)
corspearman(Jtgrad_at_rs_n2_nrs,deltaprimes_n2_nrs)

#Jp n2 smaller correlation, might come from Jt (since they're all correlated)
histogram2d(Jpgrad_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="∇Jp (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(Jpgrad_at_rs_n2_nrs,deltaprimes_n2_nrs)
corspearman(Jpgrad_at_rs_n2_nrs,deltaprimes_n2_nrs)

########################### Does local gradient correlate with Delta' for 3/1? Yes but other way (is it absolute total local current gradient?)
non_rs__n3_inds = filter(i -> resistive_equils_n2_n3[i].Δprimes[1].rs>0.3 && resistive_equils_n2_n3[i].Δprimes[2].rs<0.7, 1:length(resistive_equils_n2_n3_clean)) #getting rid of any correlation with position for n2
non_rs__n3_equil = filter(i -> i.Δprimes[1].rs>0.3 && i.Δprimes[2].rs<0.7, resistive_equils_n2_n3_clean)

deltaprimes_n2_nrs_n3= [i.Δprimes[1].Δprime for i in non_rs__n3_equil]
deltaprimes_n3_nrs_n3 = [i.Δprimes[2].Δprime for i in non_rs__n3_equil]
sorted_inds_n2_nrs_n3 = sortperm(deltaprimes_n2_nrs_n3)
sorted_inds_n3_nrs_n3 = sortperm(deltaprimes_n3_nrs_n3)

Jts_nrs_n3 = [i.equilibrium.Jt for i in non_rs__n3_equil]
Jps_nrs_n3 = [i.equilibrium.Jp for i in non_rs__n3_equil]
rss_n2_nrs_n3 = [i.Δprimes[1].rs for i in non_rs__n3_equil]
rss_n3_nrs_n3 = [i.Δprimes[2].rs for i in non_rs__n3_equil]
Jtots_nrs_n3 = [total_plasma_current(i.equilibrium.Jt,i.equilibrium.rb) for i in non_rs__n3_equil]

#There's a small number of equilibria where rs_n3 is small enough that the correlation with the edge goes away
mean(Jtots_nrs)
mean(Jtots_nrs_n3)
histogram(Jtots_nrs)
histogram!(Jtots_nrs_n3)
#However, there is strong trend in likelihood of rs being towards outer edge, and likelihood of total current being small
histogram(rss_n3_nrs_n3)
#these effects may be combining to confound the 'opposite trend in local gradient correlating with Delta'. However we can't be sure. Results are as they stand.

histogram2d(rss_n2_nrs_n3,deltaprimes_n2_nrs_n3; xlabel="2/1 rational surface location (m)",ylabel = "Δ' 2/1")
cor(rss_n2_nrs_n3,deltaprimes_n2_nrs_n3)
corspearman(rss_n2_nrs_n3,deltaprimes_n2_nrs_n3)
histogram2d(rss_n3_nrs_n3,deltaprimes_n3_nrs_n3; xlabel="3/1 rational surface location (m)",ylabel = "Δ' 3/1")
cor(rss_n3_nrs_n3,deltaprimes_n3_nrs_n3)
corspearman(rss_n3_nrs_n3,deltaprimes_n3_nrs_n3)

Jtgrad_at_rs_n2_nrs_n3 = [ScrewPinchAsymptoticMatching.gradient(Jts_nrs_n3[i],rss_n2_nrs_n3[i]) for i in 1:length(non_rs__n3_equil)]
Jtgrad_at_rs_n3_nrs_n3 = [ScrewPinchAsymptoticMatching.gradient(Jts_nrs_n3[i],rss_n3_nrs_n3[i]) for i in 1:length(non_rs__n3_equil)]
Jpgrad_at_rs_n2_nrs_n3 = [ForwardDiff.derivative(Jps_nrs_n3[i],rss_n2_nrs_n3[i]) for i in 1:length(non_rs__n3_equil)]
Jpgrad_at_rs_n3_nrs_n3 = [ForwardDiff.derivative(Jps_nrs_n3[i],rss_n3_nrs_n3[i]) for i in 1:length(non_rs__n3_equil)]

plot_current_profiles(non_rs__n3_equil[sorted_inds_n3_nrs_n3[1:10]];m_ind=2)

###Recheck n2s: Still there but weaker
histogram2d(Jtgrad_at_rs_n2_nrs_n3,deltaprimes_n2_nrs_n3; xlabel="∇Jt (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(Jtgrad_at_rs_n2_nrs_n3,deltaprimes_n2_nrs_n3)
corspearman(Jtgrad_at_rs_n2_nrs_n3,deltaprimes_n2_nrs_n3)

histogram2d(Jpgrad_at_rs_n2_nrs_n3,deltaprimes_n2_nrs_n3; xlabel="∇Jp (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(Jpgrad_at_rs_n2_nrs_n3,deltaprimes_n2_nrs_n3)
corspearman(Jpgrad_at_rs_n2_nrs_n3,deltaprimes_n2_nrs_n3)

###Recheck n3s: NOW STRONG
histogram2d(Jtgrad_at_rs_n3_nrs_n3,deltaprimes_n3_nrs_n3; xlabel="∇Jt (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(Jtgrad_at_rs_n3_nrs_n3,deltaprimes_n3_nrs_n3)
corspearman(Jtgrad_at_rs_n3_nrs_n3,deltaprimes_n3_nrs_n3)
##
histogram2d(Jpgrad_at_rs_n3_nrs_n3,deltaprimes_n3_nrs_n3; xlabel="∇Jp (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(Jpgrad_at_rs_n3_nrs_n3,deltaprimes_n3_nrs_n3)
corspearman(Jpgrad_at_rs_n3_nrs_n3,deltaprimes_n3_nrs_n3)

plot_current_profiles(non_rs__n3_equil[sorted_inds_n3_nrs_n3[1:10]];m_ind=2) #Good = on a hill
plot_current_profiles(non_rs__n3_equil[sorted_inds_n3_nrs_n3[1:10]];m_ind=2) #Terrible = in a stepp well, or on the side of a steep well

########################### Checking Absolute current gradient: 

#Jt n2 big correlation is unchanged
histogram2d(abs.(Jtgrad_at_rs_n2_nrs),deltaprimes_n2_nrs; xlabel="∇Jt (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(abs.(Jtgrad_at_rs_n2_nrs),deltaprimes_n2_nrs)
corspearman(abs.(Jtgrad_at_rs_n2_nrs),deltaprimes_n2_nrs)

#Jp n2 sligthly worsened -> it's all correlated idk
histogram2d(abs.(Jpgrad_at_rs_n2_nrs),deltaprimes_n2_nrs; xlabel="∇Jp (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(abs.(Jpgrad_at_rs_n2_nrs),deltaprimes_n2_nrs)
corspearman(abs.(Jpgrad_at_rs_n2_nrs),deltaprimes_n2_nrs)

#N3 Jt bout same
histogram2d(abs.(Jtgrad_at_rs_n3_nrs_n3),deltaprimes_n3_nrs_n3; xlabel="∇Jt (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(abs.(Jtgrad_at_rs_n3_nrs_n3),deltaprimes_n3_nrs_n3)
corspearman(abs.(Jtgrad_at_rs_n3_nrs_n3),deltaprimes_n3_nrs_n3)
#N3 Jp bout same
histogram2d(abs.(Jpgrad_at_rs_n3_nrs_n3),deltaprimes_n3_nrs_n3; xlabel="∇Jp (\$\\textrm{A/m}^{2}\$) at rs of 2/1",ylabel = "Δ' of 2/1")
cor(abs.(Jpgrad_at_rs_n3_nrs_n3),deltaprimes_n3_nrs_n3)
corspearman(abs.(Jpgrad_at_rs_n3_nrs_n3),deltaprimes_n3_nrs_n3)


########################### Second big thing: being on hill correlates with stability. For points in the center of the hill, hill steepness strong correlation with stability
min_closeness=0.0
sep_from_edge=0.05
#hill_equils, a_to_hill, hill_width, hill_gradient, Jt_double_derivative_in_hill, a_to_hill_normed, hill_inds = rs_near_local_minmax(non_rs_equil,min_closeness,sep_from_edge;well=false, return_wellhill_info=true)

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
#@save "hill_data_n1to7_6677.jld2" hill_equils a_to_hill hill_width hill_gradient Jt_double_derivative_in_hill a_to_hill_normed hill_inds
#@load "hill_data_n1to7_6677.jld2" hill_equils a_to_hill hill_width hill_gradient Jt_double_derivative_in_hill a_to_hill_normed hill_inds



#@save "hill_data_n1to4_667.jld2" hill_equils a_to_hill hill_width hill_gradient Jt_double_derivative_in_hill a_to_hill_normed hill_inds
@load "hill_data_n1to4_667.jld2" hill_equils a_to_hill hill_width hill_gradient Jt_double_derivative_in_hill a_to_hill_normed hill_inds

deltaprime_n2_hill= [i.Δprimes[1].Δprime for i in hill_equils]

hill_gradient[1:10]
plot_current_profiles(hill_equils[1:10];m_ind=1)

min_normed_closeness = 0.1 #Norm to 1
right_in_hill_inds = filter(i -> a_to_hill_normed[i]<min_normed_closeness, 1:length(a_to_hill_normed))
scatter(hill_gradient[right_in_hill_inds],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]])
cor(hill_gradient[right_in_hill_inds],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]])
corspearman(float.(hill_gradient[right_in_hill_inds]),[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]])

min_closeness = 0.008 #Stronger correlations, use this
right_in_hill_inds = filter(i -> a_to_hill[i]<min_closeness, 1:length(a_to_hill))
scatter(hill_gradient[right_in_hill_inds],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]]; xlabel="∇J/\$\\textrm{J}_{hill}\$", ylabel="Δ' of 2/1", color=:black, label=false)
cor(hill_gradient[right_in_hill_inds],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]])
corspearman(float.(hill_gradient[right_in_hill_inds]),[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]])

    ##########CHOSING M3DC1 runs
        randidxs = rand(1:length(right_in_hill_inds),5)
        scatter(hill_gradient[right_in_hill_inds],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds]])
        #scatter!(hill_gradient[right_in_hill_inds[randidxs]],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds[randidxs]]])
        is=[30,49,51,46,79]
        scatter!(hill_gradient[right_in_hill_inds[is]],[i.Δprimes[1].Δprime for i in hill_equils[right_in_hill_inds[is]]])

        #49##
        #51##
        #46##
        #79##
        #30##
        hill_m3dc1_equils = hill_equils[right_in_hill_inds[is]]
        plot_current_profiles(hill_m3dc1_equils[4:5];m_ind=1)
        deltaprime_n2_hill= [i.Δprimes[1].Δprime for i in hill_m3dc1_equils]

########################### What about being in a well? The strongest correlation... why? because well = bad but low current gradient = good, so goes from robustly stable to very unstable as well depth increases from 0. 
#well_equils, a_to_well, well_width, well_gradient, Jt_double_derivative_in_well, a_to_well_normed, well_inds = rs_near_local_minmax(non_rs_equil,min_closeness,sep_from_edge;well=true,return_wellhill_info=true)

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
#@save "well_data_n1to7_6677.jld2" well_equils a_to_well well_width well_gradient Jt_double_derivative_in_well a_to_well_normed well_inds
#@load "well_data_n1to7_6677.jld2" well_equils a_to_well well_width well_gradient Jt_double_derivative_in_well a_to_well_normed well_inds


#@save "well_data_n1to4_667.jld2" well_equils a_to_well well_width well_gradient Jt_double_derivative_in_well a_to_well_normed well_inds
@load "well_data_n1to4_667.jld2" well_equils a_to_well well_width well_gradient Jt_double_derivative_in_well a_to_well_normed well_inds

plot_current_profiles(well_equils[1:100];m_ind=1)

deltaprime_n2_well= [i.Δprimes[1].Δprime for i in well_equils]

well_gradient[1:10]
plot_current_profiles(well_equils[1:10];m_ind=1)

min_normed_closeness = 0.1 #Norm to 1
right_in_well_inds = filter(i -> a_to_well_normed[i]<min_normed_closeness, 1:length(a_to_well_normed))
scatter(well_gradient[right_in_well_inds],[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])
cor(well_gradient[right_in_well_inds],[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])
corspearman(float.(well_gradient[right_in_well_inds]),[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])

min_closeness = 0.008 #Stronger correlations, use this
right_in_well_inds = filter(i -> a_to_well[i]<min_closeness, 1:length(a_to_well))
scatter(well_gradient[right_in_well_inds],[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]]; xlabel="∇J/\$\\textrm{J}_{well}\$", ylabel="Δ' of 2/1", color=:black, label=false)
cor(well_gradient[right_in_well_inds],[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])
corspearman(float.(well_gradient[right_in_well_inds]),[i.Δprimes[1].Δprime for i in well_equils[right_in_well_inds]])


right_in_well_equils = collect(well_equils[right_in_well_inds])
right_in_well_gradient = collect(well_gradient[right_in_well_inds])

sorted_well_equils = collect(right_in_well_equils[sortperm([i.Δprimes[1].Δprime for i in right_in_well_equils])])
sorted_well_grad = collect(right_in_well_gradient[sortperm([i.Δprimes[1].Δprime for i in right_in_well_equils])])

i=1
scatter(sorted_well_grad,[i.Δprimes[1].Δprime for i in sorted_well_equils])
scatter!(sorted_well_grad[i:i],[j.Δprimes[1].Δprime for j in sorted_well_equils[i:i]]; label=false)

iswell=[4,41,60,80,100,150,230,284]
well_m3dc1_equils=collect(sorted_well_equils[iswell])
plot_current_profiles(well_m3dc1_equils;m_ind=1)
deltaprime_n2_hill= [i.Δprimes[1].Δprime for i in well_m3dc1_equils]


########################### Some general correlations:
Jts_nrs = [i.equilibrium.Jt for i in non_rs_equil]                                          #
Jps_nrs = [i.equilibrium.Jp for i in non_rs_equil]                                          #
rss_n2_nrs = [i.Δprimes[1].rs for i in non_rs_equil]                                        #
rss_n3_nrs = [i.Δprimes[2].rs for i in non_rs_equil]                                        #
Jtots_nrs = [total_plasma_current(i.equilibrium.Jt,i.equilibrium.rb) for i in non_rs_equil] #
Jt_at_rs_n2_nrs = [Jts_nrs[i](rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]               
Jt_at_rs_n3_nrs = [Jts_nrs[i](rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]
Jp_at_rs_n2_nrs = [Jps_nrs[i](rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
Jp_at_rs_n3_nrs = [Jps_nrs[i](rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]
shear_at_rs_n2_nrs = [ForwardDiff.derivative(non_rs_equil[i].equilibrium.q,rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
shear_at_rs_n3_nrs = [ForwardDiff.derivative(non_rs_equil[i].equilibrium.q,rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]
q_at_rs_n2_nrs = [non_rs_equil[i].equilibrium.q(rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
q_at_rs_n3_nrs = [non_rs_equil[i].equilibrium.q(rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]
p_at_rs_n2_nrs= [mu0*non_rs_equil[i].equilibrium.p(rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
p_at_rs_n3_nrs = [mu0*non_rs_equil[i].equilibrium.p(rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]
dpdr_at_rs_n2_nrs = [non_rs_equil[i].equilibrium.dpdr(rss_n2_nrs[i]) for i in 1:length(non_rs_equil)]
dpdr_at_rs_n3_nrs = [non_rs_equil[i].equilibrium.dpdr(rss_n3_nrs[i]) for i in 1:length(non_rs_equil)]


histogram2d(shear_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="∇q at rs 2/1",ylabel = "Δ' 2/1") #SMALL CORRELATION
cor(shear_at_rs_n2_nrs,deltaprimes_n2_nrs)
corspearman(shear_at_rs_n2_nrs,deltaprimes_n2_nrs)
histogram2d(Jtots_nrs,deltaprimes_n2_nrs; xlabel="Total current (A)",ylabel = "Δ' 2/1") #NO CORRELATION
cor(Jtots_nrs,deltaprimes_n2_nrs)
corspearman(Jtots_nrs,deltaprimes_n2_nrs)
histogram2d(Jt_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="Jt at rs 2/1 (\$\\textrm{A/m}^{2}\$)",ylabel = "Δ' 2/1") #SMALL CORRELATION
cor(Jt_at_rs_n2_nrs,deltaprimes_n2_nrs)
corspearman(Jt_at_rs_n2_nrs,deltaprimes_n2_nrs)
histogram2d(Jp_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="Jp at rs 2/1 (\$\\textrm{A/m}^{2}\$)",ylabel = "Δ' 2/1") #ALMOST NO CORRELATION
cor(Jp_at_rs_n2_nrs,deltaprimes_n2_nrs)
corspearman(Jp_at_rs_n2_nrs,deltaprimes_n2_nrs) 
histogram2d(p_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="Pressure at rs 2/1 (sqrt(T))",ylabel = "Δ' 2/1")  #SMALL CORRELATION (parabolic pressure profile its about location of rs)
cor(p_at_rs_n2_nrs,deltaprimes_n2_nrs) 
corspearman(p_at_rs_n2_nrs,deltaprimes_n2_nrs) 
histogram2d(dpdr_at_rs_n2_nrs,deltaprimes_n2_nrs; xlabel="Pressure gradient at rs 2/1 (sqrt(T)/m)",ylabel = "Δ' 2/1")   #SMALL CORRELATION LIKE p (parabolic pressure profile its about location of rs)
cor(dpdr_at_rs_n2_nrs,deltaprimes_n2_nrs)
corspearman(dpdr_at_rs_n2_nrs,deltaprimes_n2_nrs)








#What makes the current profile more stable... in a cylinder, its being on a Jt hill, not being in a Jt well.
    #Also the instantaneous current Jt gradient being small 
    
#Cases to confirm:

#CHOOSING EQUILIBRIA
# P1   1. A general stability to instability scan - randomly take a couple of equilibria from 5-10 segments of the delta' scale (for n=2, with n=3 stable for simplicity)
    
# P2   2. Test some interesting cases:
            #Some n=2 stable, n=3 unstable cases
            #A couple n=2, n=3 unstable cases
            #A couple n=2 marginal, n=3 wildly unstable cases
# P111 3. A sequence of well and hill equilibria (since this is my strongest argument)
# P2   4. Remake Cesar's stability to instability case 

hill_m3dc1_equils
well_m3dc1_equils

cd("/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase")
rvec = range(0.01,rb;step=0.01)
#print_equil_data(hill_m3dc1_equils; directoryprefactor="hill_equil", rvec=rvec)
#print_equil_data(well_m3dc1_equils; directoryprefactor="well_equil", rvec=rvec)

path="/Users/sbenjamin/Desktop/PHD/ScrewPinchDataBase"
for i in 1:length(hill_m3dc1_equils)
    mkdir("hill_equil_$(i)")
    cd("hill_equil_$(i)")

    print_equil_data_(hill_m3dc1_equils[i]; rvec=rvec)
    cd(path)
end

for i in 1:length(well_m3dc1_equils)
    mkdir("well_equil_$(i)")
    cd("well_equil_$(i)")

    print_equil_data_(well_m3dc1_equils[i]; rvec=rvec)
    cd(path)
end






#make list
#make function to gen_inputs
#download the equilibria onto julia on the cluster


#RUN PLAN:
# P2   . Do some 6-field to confirm instability 






















#USE THESE FOR TITLES:
#histogram2d(grad_J_in[filtinds],deltaprimes[filtinds]; xlabel="∇J/J in",ylabel = "Δ'", title="∇J/Jin correlation with Δ'")
#cor(grad_J_in[filtinds],deltaprimes[filtinds])
#histogram2d(grad_J_out,deltaprimes; xlabel="∇J/J out",ylabel = "Δ'", title="∇J/Jout correlation with Δ'") 
#cor(grad_J_out,deltaprimes)