using MatrixEquations
using ForwardDiff
using Distances
using Plots
using Statistics
using NaNStatistics
using LinearAlgebra
using RootsAndPoles
using SpecialFunctions

##########################################################################################################################################################
#Matching given external Delta'
    #UNFINISHED
##########################################################################################################################################################

#TESTED
function contour_inner(Qmin,Qmax;args=[pi/2,-pi/2],in_res=100,out_res=100,rad_res=100)
    in_circ = Qmin.*exp.(im.*range(args[1],args[2], length=in_res))
    out_circ = Qmax.*exp.(im.*range(args[2],args[1], length=out_res))
    rad1 = range(Qmin,Qmax, length=rad_res).*exp(im*args[2])
    rad2 = range(Qmax,Qmin, length=rad_res).*exp(im*args[1])

    return vcat(in_circ,rad1,out_circ,rad2)
end

#TESTED
function generate_inner_ratios(xmax::Real,capE,capF,capH,capG,capK; k1=3,k2=2, N=2000, truncate_terms=[2;1;1]) 
    umatrix_(Q) = X -> generate_Umatrix(X,Q,k1,k2,capE,capF,capH,capG,capK;truncate_terms=truncate_terms)[1];

    EN_minus_One_odd = Q -> compute_EN_minus_One(1,N,xmax,Q,capE,capF,capH,capG,capK);
    EN_minus_One_even = Q -> compute_EN_minus_One(-1,N,xmax,Q,capE,capF,capH,capG,capK);

    Δo = Q -> alpha3_on_alpha4(umatrix_(Q),EN_minus_One_odd(Q),N,xmax,Q,capE,capF,capH,capG,capK)
    Δe = Q -> alpha3_on_alpha4(umatrix_(Q),EN_minus_One_even(Q),N,xmax,Q,capE,capF,capH,capG,capK)

    return Δo,Δe
end

#Errors sometimes -> NOT RELIABLE
function generate_inner_ratios_Xvar(capE,capF,capH,capG,capK; k1=3,k2=2, search_range=[0,6], grain=500, N=2000, truncate_terms=[2;1;1], kwargs...) 
    umatrix_(Q) = X -> generate_Umatrix(X,Q,k1,k2,capE,capF,capH,capG,capK;truncate_terms=truncate_terms)[1];

    Xmax = Q -> findXmax(Q,k1,k2,capE,capF,capH,capG,capK; search_range = search_range, verbose = false, grain=grain, truncate_terms=truncate_terms, kwargs...)[1]

    EN_minus_One_odd_(xmax) = Q -> compute_EN_minus_One(1,N,xmax,Q,capE,capF,capH,capG,capK);
    EN_minus_One_even_(xmax) = Q -> compute_EN_minus_One(-1,N,xmax,Q,capE,capF,capH,capG,capK);

    Δosub(xmax) = Q -> alpha3_on_alpha4(umatrix_(Q),EN_minus_One_odd_(xmax)(Q),N,xmax,Q,capE,capF,capH,capG,capK)
    Δesub(xmax) = Q -> alpha3_on_alpha4(umatrix_(Q),EN_minus_One_even_(xmax)(Q),N,xmax,Q,capE,capF,capH,capG,capK)

    Δo_Δe_ratios_subfunc(xmax) = Q -> (Δosub(xmax)(Q),Δesub(xmax)(Q))
    Δo_Δe_ratios = Q -> Δo_Δe_ratios_subfunc(Xmax(Q))(Q)

    return Δo_Δe_ratios
end

#TESTED
function generate_D_Q(xmax::Real,Δr,Δl,capE,capF,capH,capG,capK; k1=3,k2=2, N=2000, truncate_terms=[2;1;1]) 
    Δo,Δe = generate_inner_ratios(xmax,capE,capF,capH,capG,capK; k1=k1,k2=k2, N=N, truncate_terms=truncate_terms) 
    D_Q_indiv(Q) = det([(Δr-Δe(Q))  (Δr-Δo(Q));(Δl-Δe(Q))   -(Δl-Δo(Q))])

    return D_Q_indiv
end

function countour_Q(xmax::Real,Δr,Δl,ctour,capE,capF,capH,capG,capK; k1=3,k2=2, N=2000, truncate_terms=[2;1;1]) 
    D_Q_indiv = generate_D_Q(xmax,Δr,Δl,capE,capF,capH,capG,capK; k1=k1,k2=k2, N=N, truncate_terms=truncate_terms) 
    M_contour = D_Q_indiv.(ctour)

    return M_contour
end

#TOO SLOW...
function GRPF_Q(domain_coords,xmax::Real,Δr,Δl,capE,capF,capH,capG,capK; k1=3,k2=2, N=2000, truncate_terms=[2;1;1], plotdata=false, multithreading=false, tolerance=10^-9, maxiterations=100, maxnodes=500, tess_sizehint=5000) 
    D_Q_indiv = generate_D_Q(xmax,Δr,Δl,capE,capF,capH,capG,capK; k1=k1,k2=k2, N=N, truncate_terms=truncate_terms) 
    
    if plotdata
        Qroots, Qpoles = grpf(D_Q_indiv, domain_coords, PlotData(), GRPFParams(maxiterations, maxnodes, 3, tess_sizehint, tolerance, multithreading)) 
    else
        Qroots, Qpoles = grpf(D_Q_indiv, domain_coords, GRPFParams(maxiterations, maxnodes, 3, tess_sizehint, tolerance, multithreading)) 
    end

    return Qroots, Qpoles
end


function D_Q_Xvarying(Q,Δr,Δl,Xmax,Δo,Δe) #input Δo(Xmax),Δe(Xmax)
    xmax = Xmax(Q)
    d_q = det([(Δr-Δe(xmax)(Q)) (Δr-Δo(xmax)(Q));(Δl-Δe(xmax)(Q)) -(Δl-Δo(xmax)(Q))])

    return d_q, xmax
end

#TESTED
function countour_Q_Xvarying(Δr,Δl,ctour,capE,capF,capH,capG,capK; k1=3,k2=2, search_range=[0,6], grain=500, N=2000, truncate_terms=[2;1;1]) 
    umatrix_(Q) = X -> generate_Umatrix(X,Q,k1,k2,capE,capF,capH,capG,capK;truncate_terms=truncate_terms)[1];

    Dr = capE+capF+capH^2

    Xmax = Q -> findXmax(Q,k1,k2,capE,capF,capH,capG,capK; search_range = search_range, verbose = false, grain=grain, truncate_terms=truncate_terms)[1]

    EN_minus_One_odd_(xmax) = Q -> compute_EN_minus_One(1,N,xmax,Q,capE,capF,capH,capG,capK);
    EN_minus_One_even_(xmax) = Q -> compute_EN_minus_One(-1,N,xmax,Q,capE,capF,capH,capG,capK);

    Δo(xmax) = Q -> alpha3_on_alpha4(umatrix_(Q),EN_minus_One_odd_(xmax)(Q),N,xmax,Q,capE,capF,capH,capG,capK)
    Δe(xmax) = Q -> alpha3_on_alpha4(umatrix_(Q),EN_minus_One_even_(xmax)(Q),N,xmax,Q,capE,capF,capH,capG,capK)

    D_Q_indiv = Q -> D_Q_Xvarying(Q,Δr,Δl,Xmax,Δo,Δe)
    M_contour_raw = D_Q_indiv.(ctour)

    M_contour = [x[1] for x in M_contour_raw]
    xmaxes = [x[2] for x in M_contour_raw]

    return M_contour, xmaxes, Dr
end

#χ'(r) = Bp(r)/(2*Pi*r)
#Ψ'(r) = Bt(r)/(2*Pi*R0)
#abs(∇V) = ((2*Pi)^2*R0)*rs

#σ(r) = (Jp(r)*Bp(r) + Jt(r)*Bt(r))/(abs(B(rs))^2)
#capM = (abs(B(rs))^2/(abs(∇V)^2))*((abs(∇V)^2/abs(B(rs))^2) + (mu0*p'(rs)*(1/((2*Pi)^2*R0*rs)))^(-2)*(σ(rs)^2*abs(B(rs))^2 - (σ(rs)*abs(B(rs))^2)^2/abs(B(rs))^2))

#Λ = Ψ'(rs)^2*(χ'/Ψ')'*(1/((2*Pi)^2*R0*rs))   #(χ'/Ψ')' evaluated at rs!
#α = 2*Pi*m/Ψ'(rs))
#X0 = (ρ*capM*η^2*abs(B(rs))^4/(α^2*Λ^2*(abs(B(rs))^2/(abs(∇V)^2))^2))^(1/6)
function cylinder_root_start_Glass75(Δr,Δl,X0,capE,capF,capH,capG,capK; Vs=1, scalediff=0.1)  #NEEDS FIXING
    Dr = capE+capF+capH^2

    dr_small = abs(pi*Dr/4) < scalediff*(abs(1/(2*pi)*gamma(1/4)/gamma(3/4)*X0/Vs*(Δr+Δl))^(6/5))

    if Dr < 0  
        Δc = 1.54*(Vs/X0)*abs(Dr)^(5/6)

        if Δr+Δl > Δc
            Δsmall = (Δr+Δl)-Δc
            Qc = 0.473*abs(Dr)^(2/3)
            λ = 0.855-0.639im
            Qsmall = λ*Qc*Δsmall/Δc
            Qstart = [im*Qc+Qsmall,conj(im*Qc+Qsmall)]
        end
        scen=1
    elseif dr_small
        Qstart = [abs(1/(2*pi)*gamma(1/4)/gamma(3/4)*X0/Vs*(Δr+Δl))^(4/5)+0*im]
        scen=2
    else
        Qstart = [(pi*Dr/4)^(2/3)+0*im]
        scen=3
    end

    return Qstart, Dr, abs(pi*Dr/4)/(abs(1/(2*pi)*gamma(1/4)/gamma(3/4)*X0/Vs*(Δr+Δl))^(6/5)), scen
end