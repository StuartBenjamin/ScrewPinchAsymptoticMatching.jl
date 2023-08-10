
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

#cd("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia")
#include("InnerAsymptotics.jl")
#include("InnerIntegrator.jl")
#include("InnerOuterMatching.jl")

#cd("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia")
#include("ChandraScrewPinchEquilibrium.jl")
#include("OuterIntegrator.jl")
#include("OuterAsymptotics.jl")

#include("GeneratingInnerTerms.jl")

using Pkg
Pkg.add(PackageSpec(url="https://github.com/StuartBenjamin/ScrewPinchAsymptoticMatching.jl"))
using ScrewPinchAsymptoticMatching

############################################################################################################################################################
############################################################################################################################################################

############################################################################################################################################################
############################################################################################################################################################

#To Do by 29th:
    #Do some constant-psi predictions -> make a plot to go w ya Cesar's result
    #Do some finite pressure Furth Predictions 
    #Do some zero-pressure Furth predictions (write the big and small solutions, see if you can get the plots working?)

#To Do afterwards
    #Do some finite pressure Scaffidi Predictions
        #Get some work from Cesar

##########################################################################################################################################################
##########################################################################################################################################################

mu0 = ScrewPinchAsymptoticMatching.mu0

CesarPlots = true
Furth73_Equil = false 

if Furth73_Equil
    q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
    rs0 = 1 #CESAR BE 2 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt) 
    R0 = 20 #Controls ratio of Bt to Bp
    ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
    xb = 2 #CESAR BE 1 #Widen device to increase q at edge (peaks current more stongly)
    Bp0 = 1.0 #Size of Bp (ratio of Bp to Bt unchanged)
    β = 1e-14; #Gives a true 0-pressure Chandra Equil

    m = 2;
    n = 1;
    k = ScrewPinchAsymptoticMatching.k_(n,R0)
    c0 = 1;
    r0 = 0.000001;
else #This replicates Cesar's Modular Pinch Furth mathematica script
    q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
    rs0 = 2.0 #CESAR BE 2 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt) 
    R0 = 20 #Controls ratio of Bt to Bp
    ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
    xb = 1.0 #CESAR BE 1 #Widen device to increase q at edge (peaks current more stongly)
    Bp0 = 1.0 #Size of Bp (ratio of Bp to Bt unchanged)
    β = 1e-14; #Gives a true 0-pressure Chandra Equil

    m = 2;
    n = 1;
    k = ScrewPinchAsymptoticMatching.k_(n,R0)
    c0 = 1;
    r0 = 0.1;
end

Furth_Equil_qvar(plot_equil) = q0 -> Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=plot_equil, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))
Zero_Furth_Equil_qvar(plot_equil) = q0 -> Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=plot_equil, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

Furth_Equil_qvar(true)(1.1)
Zero_Furth_Equil_qvar(true)(1.1)

η_diff = 1.371e-5 #Checked w Cesar
ρ = 1e20*2*1.6605e-27 #Deuterium!!

q0vec = []

rsvec = []
D_Qvec = []
find_xmax_funcs = []
delminvec = []
Δl_zero_Δr_zero_vec = []
Δl_vec = []
Δr_vec = []
Δl_zero_vec = []
Δr_zero_vec = []
Δl_Δr_vec = []
Δl_Δr_vec_xmin = []
inner_terms = []
scen_start = []

rsvec_Chand = []
D_Qvec_Chand = []
find_xmax_funcs_Chand = []
delminvec_Chand = []
Δl_zero_Δr_zero_vec_Chand = []
Δl_Δr_vec_Chand = []
Δl_Δr_vec_xmin_Chand = []
inner_terms_Chand = []
scen_start_Chand = []
Drvec = []

dellength = 10
q0length = 10
delvec = 10 .^ range(-7,-0.5,length=dellength)

Δ_raw_matrix = Matrix{Float64}(undef,dellength,q0length)
Δ_raw_matrix_Chand = Matrix{Float64}(undef,dellength,q0length)

Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil_qvar(false)(0.4)
rs, rs_plot = ScrewPinchAsymptoticMatching.find_rs(q,m,n,rb,q0,ν,rs0)
rs = rs[1]
capE,capF,capH,capM,capK,capG,Lr,Qr,capE_plus_capF = generateInnerTerms_ScrewPinch_Scaffidi(Bt, Bp, q, p, dpdr, rs, k, η_diff, n; γ=5/3, ρ=ρ)

nmax = 0 #Deprecated
xmax = 1e7

t_a = (rs0*xb)/(11/sqrt(mu0*ρ))
t_r=mu0*(rs0*xb)^2/η_diff

#(t_r/t_a)^(-2/5)

for i in 1:q0length
    print(i,"\n")
    Furth73_Equil ? (q0 = 0.44 + (i-1)*(1.0/q0length)) : (q0 = 1.1 + (i-1)*(0.2/q0length))
    push!(q0vec, q0)

    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil_qvar(false)(q0)
    rs, rs_plot = ScrewPinchAsymptoticMatching.find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]
    push!(rsvec,rs)

    if true
        solin0, solout0, inScale0, outScale0, del0 = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, minimum(delvec); integrator_reltol=10^(-20), plot_solution=true)
        #delta_tuple, solin, solout, inScale, outScale, del  = raw_delta_prime(Bp, Bt, dpdr, k, m, r0, rs, rb, minimum(delvec); integrator_reltol=10^(-20), g = g_, plot_solution=true)

        raw_delta_prime_del_var = del2 -> raw_delta_prime(solin0, solout0, inScale0, outScale0, rs, del0; del2=del2)[1]

        Δ_raw_matrix[:,i] =  raw_delta_prime_del_var.(delvec)
        
        Δl_zero,Δr_zero,del2 = Δl_Δr_calculator_zeroPressure(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, 1e-5; g_int=ScrewPinchAsymptoticMatching.g_, test_κ=false, integrator_reltol=10^(-20), plot_solution=false, del2 = 0.0, plotwidth=3, plot_soln = false, plot_soln_1=false)

        Δl_,Δr_,del2 = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, 1e-5; integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=100, plot_soln = true, plot_soln_1=false)

        push!(Δl_vec,Δl_)
        push!(Δr_vec,Δr_)

        push!(Δl_zero_vec,Δl_zero)
        push!(Δr_zero_vec,Δr_zero)
        push!(Δl_zero_Δr_zero_vec,Δl_zero-Δr_zero)
        push!(Δl_Δr_vec,Δl_+Δr_)
        #nmax = 20
        del = 10^-7

        #Δl,Δr,del2 = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false)
        
        
        del_min,Δl_Δr_xmin_output,abs_err = Δl_Δr_xmin(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec; big_n=20, xmin_tol = 1e-5, integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false, verbose = false)

        Δl_Δr_xmin_visualiser(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, delvec; finite_pressure=false, plot_solution=false, plotwidth=10, plot_soln = false)

        push!(delminvec,del_min)
        push!(Δl_Δr_vec_xmin,Δl_Δr_xmin_output[1]+Δl_Δr_xmin_output[2])

        #Δvec,Δl_Δr_xmin_output = Δl_Δr_xmin(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec; big_n=20, xmin_tol = 1e-3, integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false, verbose = false)

        capE,capF,capH,capM,capK,capG,Lr,Qr,capE_plus_capF = generateInnerTerms_ScrewPinch_Scaffidi(Bt, Bp, q, p, dpdr, rs, k, η_diff, n; γ=5/3, ρ=ρ)
        push!(inner_terms, generateInnerTerms_ScrewPinch_Scaffidi(Bt, Bp, q, p, dpdr, rs, k, η_diff, n; γ=5/3, ρ=ρ))

        Dr = capE+capF+capH^2
        push!(Drvec,(Dr,(pi*Dr/4)^(2/3)))
        #NEED XMAX...
        #NEED XMIN...

        push!(find_xmax_funcs,Q -> findXmax(Q,k1=3,k2=2,capE,capF,capH,capG,capK; search_range = [0,3], grain = 1000, convergence_tol = 10.0^(-7), verbose = false, truncate_terms=[2,1,1]))
        D_Q = generate_D_Q(xmax,Δr_,Δl_,capE,capF,capH,capG,capK; k1=3,k2=2, N=2000, truncate_terms=[2;1;1]) 

        #push!(delta_vec, Δl+Δr)
        push!(D_Qvec,D_Q)

        Qstart, Dr, scalediff, scen = cylinder_root_start_Glass75(Δr_,Δl_,Lr,capE,capF,capH,capG,capK; Vs=1, scalediff=0.1)
        push!(scen_start,[Qstart, Dr, scalediff, scen])
        
    end
end

if false
    q01 = 1.142299999999999
    q02 = 1.1447999999999998
    #vline([q01],color=:1,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #vline!([q02],color=:3,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #pC = plot([1.16,1.18],[0.0,0.0],color=:red,linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label="M3D-C1 (zero-pressure) stability transition")
    pC = plot([1.16,1.16],[-0.05,0.05],color=:red, linestyle=:solid, linewidth=2,ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label="M3D-C1 (zero-pressure) stability transition",dpi=500)
    #pC = plot!([q01,q01],[-0.05,0.05],color=:black, linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    pC = plot!([q02,q02],[-0.05,0.05],color=:black,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)

    #pC = plot!(q0vec,Δl_Δr_vec; color=:1, xlabel = "q0",ylabel="Δ'", label="finite-pressure Δ'")#  transition at q0 = 1.142") #O($(nmax)) expansion")
    pC = plot!(q0vec,Δl_zero_Δr_zero_vec, color=:3, label="zero-pressure Δ'",xlabel = "q0",ylabel="Δ'")#   transition at q0 = 1.145")#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))   #382-190=0.5 => 1.0 = 384.0, (384-308)/384 = 0.2 (where the zero-crossing is in Nish.) Zero crossing in Furth is around 0.2 as well

    #pC = plot!([1.18,1.18],[-0.05,0.05],color=:red,linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false,dpi=300)#, size=(750,750))
    png(pC,"pC1d")
end
if false
    q01 = 1.142299999999999
    q02 = 1.1447999999999998
    #vline([q01],color=:1,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #vline!([q02],color=:3,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    pC = plot([1.16,1.18],[0.0,0.0],color=:red,linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label="M3D-C1 (zero-pressure) stability transition")
    pC = plot!([q01,q01],[-0.2,0.2],color=:black, linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #pC = plot!([q02,q02],[-0.05,0.05],color=:black,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)

    pC = plot!(q0vec,Δl_Δr_vec; color=:1, xlabel = "q0",ylabel="Δ'", label="finite-pressure Δ'")#  transition at q0 = 1.142") #O($(nmax)) expansion")
    #pC = plot!(q0vec,Δl_zero_Δr_zero_vec, color=:3, label="zero-pressure Δ'")#   transition at q0 = 1.145")#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))   #382-190=0.5 => 1.0 = 384.0, (384-308)/384 = 0.2 (where the zero-crossing is in Nish.) Zero crossing in Furth is around 0.2 as well
    pC = plot!([1.16,1.16],[-0.2,0.2],color=:red, linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    pC = plot!([1.18,1.18],[-0.2,0.2],color=:red,linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false,dpi=300)#, size=(750,750))
    png(pC,"pC2")

    q01 = 1.142299999999999
    q02 = 1.1447999999999998
    #vline([q01],color=:1,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #vline!([q02],color=:3,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #pC = plot([1.16,1.18],[0.0,0.0],color=:red,linestyle=:solid,label="M3D-C1 (zero-pressure) stability transition")
    pC = plot([q01,q01],[-0.5,0.5],color=:black, linestyle=:dot,label=false)#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    pC = plot!([q02,q02],[-0.5,0.5],color=:black,linestyle=:dot,label=false)#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)

    pC = plot!(q0vec,Δl_Δr_vec; color=:1, xlabel = "q0",ylabel="Δ'", label="finite-pressure Δ'")#  transition at q0 = 1.142") #O($(nmax)) expansion")
    pC = plot!(q0vec,Δl_zero_Δr_zero_vec, color=:3, label="zero-pressure Δ'",dpi=300)#   transition at q0 = 1.145")#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))   #382-190=0.5 => 1.0 = 384.0, (384-308)/384 = 0.2 (where the zero-crossing is in Nish.) Zero crossing in Furth is around 0.2 as well
    #pC = plot!([1.16,1.16],[-0.2,0.2],color=:red, linestyle=:solid,label=false)
    #pC = plot!([1.18,1.18],[-0.2,0.2],color=:red,linestyle=:solid,label=false,dpi=300)#, size=(750,750))
    png(pC,"pC3")
end
if false
    q01 = 1.142299999999999
    q02 = 1.1447999999999998
    #vline([q01],color=:1,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    #vline!([q02],color=:3,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    plot([1.16,1.18],[0.0,0.0],color=:red,linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label="M3D-C1 (zero-pressure) stability transition")
    plot!([q01,q01],[-0.05,0.05],color=:black, linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    plot!([q02,q02],[-0.05,0.05],color=:black,linestyle=:dot, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)

    plot!(q0vec,Δl_Δr_vec; color=:1, xlabel = "q0",ylabel="Δ'", label="finite-pressure Δ'")#  transition at q0 = 1.142") #O($(nmax)) expansion")
    plot!(q0vec,Δl_zero_Δr_zero_vec, color=:3, label="zero-pressure Δ'")#   transition at q0 = 1.145")#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))   #382-190=0.5 => 1.0 = 384.0, (384-308)/384 = 0.2 (where the zero-crossing is in Nish.) Zero crossing in Furth is around 0.2 as well
    plot!([1.16,1.16],[-0.05,0.05],color=:red, linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    plot!([1.18,1.18],[-0.05,0.05],color=:red,linestyle=:solid, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)),label=false)
    

    #plot!(q0vec,Δl_Δr_vec_xmin)#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))
    

    #plot(rsvec,Δl_Δr_vec)#,xlims = (0.0,xb))
    #plot!(rsvec,Δl_zero_Δr_zero_vec)#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))
    #plot!(rsvec,Δl_Δr_vec_xmin)#, ylims = (minimum(Δl_Δr_vec),maximum(Δl_Δr_vec)))
    #plot(rsvec,Δl_zero_Δr_zero_vec)
end

if false  #Δe Δo orig plot
    #Checking Xmax
    qlength=10
    QvecNew = range(-5,5,length=qlength)
    xmax = 1e5
    checkXmax_Q = Q -> checkXmax(Q,xmax,3,2,capE,capF,capH,capG,capK; convergence_tol = 10.0^(-7), truncate_terms=[2,1,1])
    checkXmax_Q_bool = Q -> checkXmax_Q(Q)[1]
    sum(checkXmax_Q_bool.(QvecNew))
    qlength

    #Cookin 
    d_Q_Evaluated = []
    Δo_Evaluated = []
    Δe_Evaluated = []
    Qvec_plots = []

    xmaxtemp = 1000.0

    """
    for i in 1:length(q0vec)

        capE,capF,capH,capM,capK,capG,Lr,Qr,capE_plus_capF = inner_terms[i]
        dQ = generate_D_Q(xmaxtemp,Δl_vec[i],Δr_vec[i],capE,capF,capH,capG,capK; k1=5,k2=3, N=2000, truncate_terms=[-1;-1;-1])
        Δo,Δe = generate_inner_ratios(xmaxtemp,capE,capF,capH,capG,capK; k1=5,k2=3, N=2000, truncate_terms=[-1;-1;-1]) 

        push!(d_Q_Evaluated,dQ.(QvecNew))
        push!(Δo_Evaluated,Δo.(QvecNew))
        push!(Δe_Evaluated,Δe.(QvecNew))
        push!(Qvec_plots,QvecNew)
        print(i,"\n")
    end
    

        pdidlle=[]
        for i in 1:length(q0vec)
            if i==1
                pdidlle= plot(QvecNew,abs.(d_Q_Evaluated[i]))#,xaxis=:log)
            else
                pdidlle= plot!(QvecNew,abs.(d_Q_Evaluated[i]))
            end
        end
        display(pdidlle)

        Δe_Case = Δe_Evaluated[1]
        Δo_Case = Δo_Evaluated[1]
        Q_case = d_Q_Evaluated[1]
        QvecNew_case = QvecNew
    """

    @load "q1p1_D_Qdata.jld2" Q_case Δo_Case Δe_Case QvecNew_case

    plot(QvecNew_case,real.(Q_case))

    plot(QvecNew_case,real.(Δo_Case))
    plot(QvecNew_case,imag.(Δo_Case))

    plot(QvecNew_case,real.(Δe_Case))
    plot(QvecNew_case,imag.(Δe_Case))

    plot(QvecNew_case,real.(Δo_Case),ylims = (-10.0,10.0),title="Odd-solution asymptotic ratio Δo",ylabel="Δo",xlabel="Q",dpi=500)
    png("q0_1p1_Δo")
    plot(QvecNew_case,real.(Δo_Case),title="Odd-solution asymptotic ratio Δo",ylabel="Δo",xlabel="Q",dpi=500)
    png("q0_1p1_Δo1")

    plot(QvecNew_case,real.(Δe_Case),title="Even-solution asymptotic ratio Δe",ylabel="Δe",xlabel="Q",dpi=500)
    png("q0_1p1_Δe")
    
    sortperm(abs.(real.(Δo_Case));reverse=true)
    vline!(QvecNew_case[sortperm(abs.(real.(Δo_Case));rev=true)[1:3]])

    Qmin1 = QvecNew_case[sortperm(abs.(real.(Δo_Case));rev=true)[1]]
    Qmin2 = QvecNew_case[sortperm(abs.(real.(Δo_Case));rev=true)[3]]

    plot(QvecNew_case,real.(Q_case),ylims = (-5*1e4,5*1e4), title = "Solving dispersion relation for growth-rate Q", xlabel="Q", ylabel="Δ(Q)",label=false)
    vline!(QvecNew_case[sortperm(abs.(Q_case))[1:1]], linewidth=2, color=:black, linestyle=:dot,dpi=500,label = "Q = 1.09")
    png("q0_1p1_D_Q_crossing")
    
    #hline!([0.0], color=:black, label = false)

    Q1 = QvecNew_case[sortperm(abs.(Q_case))[1:1]]
    Q1 = QvecNew_case[sortperm(abs.(Q_case))[1:1]]

    #@save "q1p1_D_Qdata.jld2" Q_case Δo_Case Δe_Case QvecNew_case
end

if true  #Δe Δo positive only (for comparison with deltac)
    #Checking Xmax

    scan_x0=0.00001
    scan_x1=2.0
    scan_nstep = 100

    log_scan_x0=log10(scan_x0)
    qstep=(log10(scan_x1)-log_scan_x0)/scan_nstep
    qlogvec = Float64[]
    qscanvec = Float64[]
    for istep=0:scan_nstep
        qlog=log_scan_x0+istep*qstep
        q_scan=10^qlog
        push!(qlogvec,qlog)
        push!(qscanvec,q_scan)
    end

    #qlength=100
    #Qvec_test_vs_deltac= range(0,5,length=qlength)
    xmax = 3e7
    checkXmax_Q = Q -> checkXmax(Q,xmax,3,2,capE,capF,capH,capG,capK; convergence_tol = 10.0^(-7), truncate_terms=[2,1,1])
    checkXmax_Q_bool = Q -> checkXmax_Q(Q)[1]
    sum(checkXmax_Q_bool.(qscanvec))
    qlength = length(qscanvec)

    #Cookin 
    d_Q_Evaluated = []
    Δo_Evaluated = []
    Δe_Evaluated = []
    Qvec_plots = []

    xmaxtemp = 1000.0

    
    for i in 1:1 #length(q0vec)
        capE,capF,capH,capM,capK,capG,Lr,Qr,capE_plus_capF = inner_terms[i]
        dQ = generate_D_Q(xmaxtemp,Δl_vec[i],Δr_vec[i],capE,capF,capH,capG,capK; k1=5,k2=3, N=2000, truncate_terms=[-1;-1;-1])
        Δo,Δe = generate_inner_ratios(xmaxtemp,capE,capF,capH,capG,capK; k1=5,k2=3, N=2000, truncate_terms=[-1;-1;-1]) 

        push!(d_Q_Evaluated,dQ.(qscanvec))
        push!(Δo_Evaluated,Δo.(qscanvec))
        push!(Δe_Evaluated,Δe.(qscanvec))
        push!(Qvec_plots,qscanvec)
        print(i,"\n")
    end
    

        #pdidlle=[]
        #for i in 1:length(q0vec)
        #    if i==1
        #        pdidlle= plot(qscanvec,abs.(d_Q_Evaluated[i]))#,xaxis=:log)
        #    else
        #        pdidlle= plot!(qscanvec,abs.(d_Q_Evaluated[i]))
        #    end
        #end
        #display(pdidlle)

        Δe_Case = Δe_Evaluated[1]
        Δo_Case = Δo_Evaluated[1]
        Q_case = d_Q_Evaluated[1]
        #QvecNew_case = QvecNew

    #@load "q1p1_D_Qdata.jld2" Q_case Δo_Case Δe_Case QvecNew_case

    #plot(qscanvec,real.(Q_case))

    display(plot(qscanvec,real.(Δo_Case)))#,ylims = (-10.0,10.0)))
    plot(qscanvec,imag.(Δo_Case))

    plot(qscanvec,real.(Δe_Case),xaxis=:log,legend=false)
    plot(qscanvec,imag.(Δe_Case))


    #plot(QvecNew_case,real.(Δo_Case),ylims = (-10.0,10.0),title="Odd-solution asymptotic ratio Δo",ylabel="Δo",xlabel="Q",dpi=500)
    #png("q0_1p1_Δo")
    #plot(QvecNew_case,real.(Δo_Case),title="Odd-solution asymptotic ratio Δo",ylabel="Δo",xlabel="Q",dpi=500)
    #png("q0_1p1_Δo1")

    #plot(QvecNew_case,real.(Δe_Case),title="Even-solution asymptotic ratio Δe",ylabel="Δe",xlabel="Q",dpi=500)
    #png("q0_1p1_Δe")
    
    """
    sortperm(abs.(real.(Δo_Case));rev=true)
    vline!(QvecNew_case[sortperm(abs.(real.(Δo_Case));rev=true)[1:3]])

    Qmin1 = QvecNew_case[sortperm(abs.(real.(Δo_Case));rev=true)[1]]
    Qmin2 = QvecNew_case[sortperm(abs.(real.(Δo_Case));rev=true)[3]]

    plot(QvecNew_case,real.(Q_case),ylims = (-5*1e4,5*1e4), title = "Solving dispersion relation for growth-rate Q", xlabel="Q", ylabel="Δ(Q)",label=false)
    vline!(QvecNew_case[sortperm(abs.(Q_case))[1:1]], linewidth=2, color=:black, linestyle=:dot,dpi=500,label = "Q = 1.09")
    png("q0_1p1_D_Q_crossing")
    
    #hline!([0.0], color=:black, label = false)

    Q1 = QvecNew_case[sortperm(abs.(Q_case))[1:1]]
    Q1 = QvecNew_case[sortperm(abs.(Q_case))[1:1]]

    #@save "q1p1_D_Qdata.jld2" Q_case Δo_Case Δe_Case QvecNew_case
    """
    display(inner_terms[1])

end

if true
    df = CSV.read("/Users/sbenjamin/desktop/phd/rmatch_local/rmatch/scanq.csv", DataFrame)

    plot(df.qvals,df.odd_ratio_real,ylims = (-10.0,10.0),title="Odd-solution asymptotic ratio Δo",ylabel="Δo",xlabel="Q",dpi=500)
    png("q0_1p1_Δo_deltac")


    plot(df.qvals,df.even_ratio_real,title="Even-solution asymptotic ratio Δe",ylabel="Δe",xlabel="Q",dpi=500)
    png("q0_1p1_Δe_deltac")

    plot(df.qvals,df.even_ratio_real,title="Even-solution asymptotic ratio Δe",ylims = (-5000.0,5000.0),ylabel="Δe",xlabel="Q",dpi=500)
    png("q0_1p1_Δe_deltac_2")


end

if false #Δe Δo multi-scan
    #Checking Xmax
    qlength=1000
    QvecNew = range(-5,5,length=qlength)
    xmax = 1e7
    checkXmax_Q = Q -> checkXmax(Q,xmax,3,2,capE,capF,capH,capG,capK; convergence_tol = 10.0^(-7), truncate_terms=[2,1,1])
    checkXmax_Q_bool = Q -> checkXmax_Q(Q)[1]
    sum(checkXmax_Q_bool.(QvecNew))
    qlength

    #Cookin 
    d_Q_Evaluated_new = []
    Δo_Evaluated_new = []
    Δe_Evaluated_new = []
    Qvec_plots_new = []

    xmaxtemp = 1000.0

    for i in 1:length(q0vec)

        capE,capF,capH,capM,capK,capG,Lr,Qr,capE_plus_capF = inner_terms[i]
        dQ = generate_D_Q(xmaxtemp,Δl_vec[i],Δr_vec[i],capE,capF,capH,capG,capK; k1=5,k2=3, N=2000, truncate_terms=[-1;-1;-1])
        Δo,Δe = generate_inner_ratios(xmaxtemp,capE,capF,capH,capG,capK; k1=5,k2=3, N=2000, truncate_terms=[-1;-1;-1]) 

        push!(d_Q_Evaluated_new,dQ.(QvecNew))
        push!(Δo_Evaluated_new,Δo.(QvecNew))
        push!(Δe_Evaluated_new,Δe.(QvecNew))
        push!(Qvec_plots_new,QvecNew)
        print(i,"\n")
    end

    #@save "BigRun_q_scan_D_Qdata.jld2" d_Q_Evaluated_new Δo_Evaluated_new Δe_Evaluated_new Qvec_plots_new QvecNew



    pdidlle1=[]
    for i in 1:length(q0vec)
        if i==1
            pdidlle1= plot(QvecNew,real.(d_Q_Evaluated_new[i]),ylims = (-5*1e4,5*1e4),xlims = (0.0,2.0))#,xaxis=:log)
        else
            pdidlle1= plot!(QvecNew,real.(d_Q_Evaluated_new[i]),ylims = (-5*1e4,5*1e4),xlims = (0.0,2.0))
        end
    end
    display(pdidlle1)

    pdidlle2=[]
    for i in 1:length(q0vec)
        if i==1
            pdidlle2= plot(QvecNew,real.(Δo_Evaluated_new[i]),ylims = (-2,2))#,xaxis=:log)
        else
            pdidlle2= plot!(QvecNew,real.(Δo_Evaluated_new[i]),ylims = (-2,2))
        end
    end
    display(pdidlle2)

    pdidlle3=[]
    for i in 1:length(q0vec)
        if i==1
            pdidlle3= plot(QvecNew,real.(Δe_Evaluated_new[i]),ylims = (-1*1e4,1*1e4))#,xaxis=:log)
        else
            pdidlle3= plot!(QvecNew,real.(Δe_Evaluated_new[i]),ylims = (-1*1e4,1*1e4))
        end
    end
    display(pdidlle3)

    """
    Δe_Case = Δe_Evaluated[1]
    Δo_Case = Δo_Evaluated[1]
    Q_case = d_Q_Evaluated[1]

    plot(QvecNew,real.(Q_case))

    plot(QvecNew,real.(Δo_Case))
    plot(QvecNew,imag.(Δo_Case))

    plot(QvecNew,real.(Δe_Case))
    plot(QvecNew,imag.(Δe_Case))

    plot(QvecNew,real.(Δo_Case),ylims = (-2.0,2.0))
    
    sortperm(abs.(real.(Δo_Case));reverse=true)
    vline!(QvecNew[sortperm(abs.(real.(Δo_Case));rev=true)[1:3]])

    Qmin1 = QvecNew[sortperm(abs.(real.(Δo_Case));rev=true)[1]]
    Qmin2 = QvecNew[sortperm(abs.(real.(Δo_Case));rev=true)[3]]

    plot(QvecNew,real.(Q_case),ylims = (-5*1e4,5*1e4), title = "Solving dispersion relation for growth-rate Q", xlabel="Q", ylabel="Δ(Q)",label=false)
    vline!(QvecNew[sortperm(abs.(Q_case))[1:1]], linewidth=2, color=:black, linestyle=:dot,dpi=500,label = "Q = 1.09")
    png("q0_1p1_D_Q_crossing")
    #hline!([0.0], color=:black, label = false)

    Q1 = QvecNew[sortperm(abs.(Q_case))[1:1]]
    Q1 = QvecNew[sortperm(abs.(Q_case))[1:1]]
    """
end

if false #DEPRECATED/OLD?
    scatter(QvecNew,zero_spot,xaxis=:log)
    argmin(abs.(zero_spot))
    plot(QvecNew,abs.(d_Q_Evaluated[1]),yaxis=:log)
    plot(QvecNew,real.(Δo_Evaluated[1]))
    plot(QvecNew,real.(Δe_Evaluated[1]))


    Δl_Δr_xmin_visualiser(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, delvec; finite_pressure=false, plot_solution=true, plotwidth=10, plot_soln = true)

    Δl_zero,Δr_zero,del2 = Δl_Δr_calculator_zeroPressure(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, 1e-5; test_κ=false, integrator_reltol=10^(-20), plot_solution=false, del2 = 0.0, plotwidth=3, plot_soln = false, plot_soln_1=false)
    Δl_zero,Δr_zero,del2 = Δl_Δr_calculator_zeroPressure(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, 1e-7; test_κ=false, integrator_reltol=10^(-20), plot_solution=false, del2 = 0.0, plotwidth=3, plot_soln = false, plot_soln_1=false)

    del_min,Δl_Δr_xmin_output,abs_err = Δl_Δr_xmin(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec; big_n=20, xmin_tol = 1e-5, integrator_reltol=10^(-20), plot_solution=true, plotwidth=10, plot_soln = true, verbose = true)
    #plot(rsvec,Δ_raw_matrix[51,:])

    #Two sources of Furth fucking up:
        #Error in one of the solutions -> most likely Δr
        #DISPROVED: Need to elimimnate the pressure from the actual solver
        plot(rsvec,Δl_zero_vec)
        plot!(rsvec,Δr_zero_vec)
        plot!(rsvec,Δl_zero_vec.-Δr_zero_vec;color=:black)  


    #using PythonPlot
    #PythonPlot.contour(rsvec, delvec, Δ_raw_matrix,levels=10, color=:turbo, clabels=true, cbar=false, lw=1)


    minbound = 0.5-abs(minimum(Δ_raw_matrix))/maximum(abs.(Δ_raw_matrix))*(0.5)
    using PlotlyJS
    colorscale = [[0.0, "black"],[0.38655, "blue"], [0.499999999, "lightblue"],[0.5000001, "orangered"], [1.0, "darkred"]]

    minspot=0.4999999*(1-abs(minimum(Δ_raw_matrix))/maximum(abs.(Δ_raw_matrix)))
    maxspot=0.4999999*(1+abs(maximum(Δ_raw_matrix))/maximum(abs.(Δ_raw_matrix)))
    colorscale2 = [[0.0, "black"],[minspot, "blue"],[0.499, "turquoise"],[0.5, "white"],[0.501, "lavender"],[(maxspot+0.501)/2, "red"],[maxspot, "darkred"], [1, "black"]]

    #maximum(abs.(Δ_raw_matrix))/connum
    connum=40
    PlotlyJS.plot(PlotlyJS.contour(x = q0vec, y = delvec, z=Δ_raw_matrix, colorscale=colorscale2,
        contours_start=-maximum(abs.(Δ_raw_matrix)),contours_end=maximum(abs.(Δ_raw_matrix)),contours_size=(2*maximum(abs.(Δ_raw_matrix))/connum)),Layout(width=500, height=400,title_x=0.5,title_text= "Δ' (constant-psi approximation)",xaxis_title="q0", yaxis_title="δ (m)", yaxis_type="log"))

    #PlotlyJS.plot(PlotlyJS.contour(x = q0vec, y = delvec, z=Δ_raw_matrix, colorscale=colorscale2,
    #contours_start=-maximum(abs.(Δ_raw_matrix)),contours_end=maximum(abs.(Δ_raw_matrix)),contours_size=(2*maximum(abs.(Δ_raw_matrix))/connum)),Layout(width=500, height=400,title_x=0.5,title_text= "Δ' (constant-psi approximation)",xaxis_title="q0", yaxis_title="δ (m)", yaxis_type="log"))


    #plot(q0vec,delta_vec)
    #plot(rsvec,delta_vec)
end

##########################################################################################################################################################
#Chandra
    #Same as Furth except pressure profile is nearly* pre-set and Bt varies with r to exactly satisfy power balance (Bt on axis set).
    #Pressure is perturbed like Bt.
    #Bp not zero.
##########################################################################################################################################################
β = 0.0000000000001 
rs0 = 2.0 
R0 = 20 
ν = 1.0  
xb = 1.0 
Bp0 = 1.0 
q0 = 1.1

Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))