using MatrixEquations
using ForwardDiff
using Distances
using Plots
using Statistics
using NaNStatistics
using LinearAlgebra
using RootsAndPoles
using SpecialFunctions

#cd("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia")
include("InnerAsymptotics.jl")
include("InnerIntegrator.jl")
include("InnerOuterMatching.jl")

##########################################################################################################################################################
#Inner equation examples:
##########################################################################################################################################################

#Screw pinch characteristics (*Inner solution plug*)
##########################################################################################################################################################
capE = -689342.
capF = 0.912824
capH = 0.0
capM = 1.0
capK = 1.09537
capG = 27.0851
Δr=25.684/2
Δl=25.684/2
X0 = 8.39931

#Inner equation integrator, Table III scenario Glasser 1984:
##########################################################################################################################################################
capE = 0.1;
capF = 0;
capH = 0;
capG = 0;
capK = 0;
Q = 0.01;
Xmax=3.1623;

Dr = capE+capF+capH^2

U_20_20_18 = X -> generate_Umatrix(X,Q,20,20,capE,capF,capH,capG,capK)[1]; #Confirmed identical with Mathematica values
umatrix = U_20_20_18;
Amat = X -> A(X,Q,capE,capF,capH,capG,capK);
UPrime_20_20_18 = XX -> ForwardDiff.derivative(U_20_20_18, XX);
DT_20_20_18_Pos = XX -> deltaTest(XX,U_20_20_18,UPrime_20_20_18,Amat;positive=true);
DT_20_20_18_Neg = XX -> deltaTest(XX,U_20_20_18,UPrime_20_20_18,Amat;positive=false);

Xvec = 10 .^ range(0, 2, length=100);
plot([Xvec,Xvec],[DT_20_20_18_Pos.(Xvec),DT_20_20_18_Neg.(Xvec)],xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3),label=["+"    "-"])
plot!([Xmax]; seriestype = :vline,label="Xmatch = $(Xmax)")

##########################################################################################################################################################
#NEEDS FIXING, IS UNRELIABLE!
inner_ratios_capE_var2 = E -> generate_inner_ratios_Xvar(E,capF,capH,capG,capK; k1=20,k2=18, search_range=[-1,6], grain=1000, N=2000, truncate_terms=[-1;-1;-1], convergence_tol=10.0^-14) 

α_ratios_o_E0p1 = inner_ratios_capE_var2(0.1)
α_ratios_o_E0pminus1 = inner_ratios_capE_var2(-0.1)

Qvec = 10 .^ range(-5,1,length=3)
Qvec_ratiod_minus1 = []
for i in 1:length(Qvec)
    push!(Qvec_ratiod_minus1,α_ratios_o_E0pminus1(Qvec[i]))
    print(i,"\n")
end

scatter(Qvec,[real(i[1]) for i in Qvec_ratiod_minus1], xaxis=:log, ylims=(-1.0,2.0))
plot(Qvec,[real(i[2]) for i in Qvec_ratiod_minus1], xaxis=:log)

#plot(Qvec,[real(i[1]) for i in Qvec_ratiod], xaxis=:log,ylims=(-1.0,2.0))
#plot(Qvec,[real(i[2]) for i in Qvec_ratiod], xaxis=:log)


#Recreating Glass 1984 Fig2
##########################################################################################################################################################
findXmax_QE(E) = Q -> findXmax(Q,18,18,E,capF,capH,capG,capK; search_range = [-1,5], grain = 1000, convergence_tol = 10.0^(-14), verbose = false, truncate_terms=[-1,-1,-1])
#findXmax_QE(0.1)(10^-5) = 0.8768856094587429
#findXmax_QE(0.1)(10^1) = 33
#findXmax_QE(-0.1)(10^-5) = 1.17
#findXmax_QE(-0.1)(10^1) = 100.7
#findXmax_QE(-0.1)(1) = 9.33

inner_ratios_capE_sus = E -> generate_inner_ratios(100.0,E,capF,capH,capG,capK; k1=20,k2=18, N=2000, truncate_terms=[-1;-1;-1]) 

Δs_E0p1 = inner_ratios_capE_sus(0.1)
Δs_E0p06 = inner_ratios_capE_sus(0.06)
Δs_E0p02 = inner_ratios_capE_sus(0.02)
Δs_E0pm02 = inner_ratios_capE_sus(-0.02)
Δs_E0pm06 = inner_ratios_capE_sus(-0.06)
Δs_E0pm1 = inner_ratios_capE_sus(-0.1)

QvecNew = 10 .^ range(-5,1,length=3)
Qvec_ratio_pos = []
Qvec_ratio_neg = []
for i in 1:length(QvecNew)
    push!(Qvec_ratio_pos,[Δs_E0p1[1](QvecNew[i])    Δs_E0p06[1](QvecNew[i]) Δs_E0p02[1](QvecNew[i]) Δs_E0pm02[1](QvecNew[i])    Δs_E0pm06[1](QvecNew[i])    Δs_E0pm1[1](QvecNew[i])])
    push!(Qvec_ratio_neg,[Δs_E0p1[2](QvecNew[i])    Δs_E0p06[2](QvecNew[i]) Δs_E0p02[2](QvecNew[i]) Δs_E0pm02[2](QvecNew[i])    Δs_E0pm06[2](QvecNew[i])    Δs_E0pm1[2](QvecNew[i])])
    print(i,"\n")
end

scatter([QvecNew;QvecNew;QvecNew;QvecNew;QvecNew;QvecNew],[[real(i[1][1]) for i in Qvec_ratio_pos];[real(i[2][1]) for i in Qvec_ratio_pos];[real(i[3][1]) for i in Qvec_ratio_pos];[real(i[4][1]) for i in Qvec_ratio_pos];[real(i[5][1]) for i in Qvec_ratio_pos];[real(i[6][1]) for i in Qvec_ratio_pos]], xaxis=:log, ylims=(-1.0,2.0))

scatter([QvecNew;QvecNew;QvecNew;QvecNew;QvecNew;QvecNew],[[real(i[1][1]) for i in Qvec_ratio_neg];[real(i[2][1]) for i in Qvec_ratio_neg];[real(i[3][1]) for i in Qvec_ratio_neg];[real(i[4][1]) for i in Qvec_ratio_neg];[real(i[5][1]) for i in Qvec_ratio_neg];[real(i[6][1]) for i in Qvec_ratio_neg]], xaxis=:log)

#Glass 1984 Fig2 but with actual dispersion relation
##########################################################################################################################################################
QvecNew = 10 .^ range(-5,1,length=20)
generate_D_Q_E = E -> generate_D_Q(100.0,0.5,0.5,E,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1]) 

zero_spot = []
for i in 1:length(QvecNew)
    push!(zero_spot,generate_D_Q_E(0.1)(complex(QvecNew[i])))
    print(i,"\n")
end

scatter(QvecNew,zero_spot,xaxis=:log)

#
##########################################################################################################################################################
#Running Scaffidi values:
    

#Testing contour is working:
##########################################################################################################################################################
ctour = contour_inner(0.001,0.61;args=[pi/2,-pi/2],in_res=10,out_res=10,rad_res=10)
M_contour = countour_Q(100.0,0.5,0.5,ctour,0.1,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1])
ctour2 = contour_inner(0.001,10;args=[pi/2,-pi/2],in_res=10,out_res=10,rad_res=10)
M_contour2 = countour_Q(100.0,0.5,0.5,ctour2,0.1,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1])

plot(real(ctour),imag(ctour))
plot(real(M_contour),imag(M_contour))
plot(real(ctour2),imag(ctour2))
plot(real(M_contour2),imag(M_contour2))

#Testing solver is working:
##########################################################################################################################################################
D_Q = generate_D_Q(100.0,0.5,0.5,capE,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1]) 

origcoords = rectangulardomain(0.4-0.1im,0.6+0.1im,0.2/8)
origcoords2 = rectangulardomain(0.061-1im,10+1im,0.5)
output1 = GRPF_Q(origcoords,100.0,0.5,0.5,capE,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1], tolerance=1e-13)
output2 = GRPF_Q(origcoords2,100.0,0.5,0.5,capE,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1]) 


output3 = GRPF_Q(origcoords,100.0,0.5,0.5,capE,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1], tolerance=1e-13, maxnodes=1000, plotdata=true)


root1 = 0.47635638776525874 - 3.403051951345049e-10im
(ComplexF64[0.47635638843736167 - 2.1269030980874958e-11im, 0.4763563885819916 - 1.4179335483532896e-11im, 0.47635638774257155 - 1.4179335483532896e-11im, 0.4763563882955678 - 2.1269030980874958e-11im, 0.476356387597942 - 2.1269030980874958e-11im, 0.4763563878815297 - 2.1269030980874958e-11im, 0.47635638802048763 - 1.4179335483532896e-11im, 0.4763563881622816 - 1.4179335483532896e-11im], ComplexF64[0.47635638852243806 - 2.1269030980874958e-11im, 0.47635638779928924 - 1.4179335483532896e-11im, 0.4763563879410831 - 1.4179335483532896e-11im, 0.47635638838348004 - 1.4179335483532896e-11im, 0.4763563880942205 - 1.4179335483532896e-11im, 0.47635638823601434 - 1.4179335483532896e-11im, 0.47635638765182364 - 1.4179335483532896e-11im])

xvec=range(0.475,0.477,length=100)
d_Qxvec = D_Q.(xvec)
plot(xvec,real.(d_Qxvec))
argmin(abs.(real.(d_Qxvec)))
d_Qxvec[68]
abs(D_Q(0.47635638843736167 - 2.1269030980874958e-11im))

origcoords3 = rectangulardomain(0.476-0.1,0.477,0.5)

cylinder_root_start_Glass75(0.5,0.5,1,capE,capF,capH,capG,capK; Vs=1, scalediff=0.001)

##########################################################################################################################################################

#Inner equation integrator, Table IV scenario Glasser 1984:
##########################################################################################################################################################
capE = 0.1;
capF = 0.005;
capH = 0.05;
capG = 40.0;
capK = 500.0;
Q = 3.0;
Xmax=39.11;

U_20_20_18 = X -> generate_Umatrix(X,Q,20,20,capE,capF,capH,capG,capK)[1]; #Confirmed identical with Mathematica values
umatrix = U_20_20_18;
Amat = X -> A(X,Q,capE,capF,capH,capG,capK);
UPrime_20_20_18 = XX -> ForwardDiff.derivative(U_20_20_18, XX);
DT_20_20_18_Pos = XX -> deltaTest(XX,U_20_20_18,UPrime_20_20_18,Amat;positive=true);
DT_20_20_18_Neg = XX -> deltaTest(XX,U_20_20_18,UPrime_20_20_18,Amat;positive=false);

Xvec = 10 .^ range(0, 2, length=1000);
plot([Xvec,Xvec],[DT_20_20_18_Pos.(Xvec),DT_20_20_18_Neg.(Xvec)],xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3),label=["+"    "-"])
plot!([Xmax]; seriestype = :vline,label="Xmatch = $(Xmax)")

N=1600;
EnMinusOne = EN_minus_One_odd = compute_EN_minus_One(1,N,Xmax,Q,capE,capF,capH,capG,capK);
EN_minus_One_even = compute_EN_minus_One(-1,N,Xmax,Q,capE,capF,capH,capG,capK);

α_ratio_o=alpha3_on_alpha4(umatrix,EnMinusOne,N,Xmax,Q,capE,capF,capH,capG,capK;var=1)
α_ratio_e=alpha3_on_alpha4(umatrix,EN_minus_One_even,N,Xmax,Q,capE,capF,capH,capG,capK)

α_ratio_o2,α_ratio_e2 = generate_inner_ratios(Xmax,capE,capF,capH,capG,capK; k1=20,k2=20, N=N, truncate_terms=[-1;-1;-1])

#Recreating Glasser 2020 asymptotic k-scaling:
##########################################################################################################################################################
capE = -3.369*10^(-2);
capF = 2.409*10^(-3);
capH = 1.292*10^(-2);
capG = 8.950*10^1;
capK = 2.332*10^2;
Qtest = 1.234*10^(-1);
k1 = 16
k2 = 16
Q = Qtest

Up0 = X -> generate_Umatrix(X,Qtest,18,18,capE,capF,capH,capG,capK)[1]
Up1prime(X) = ForwardDiff.derivative(Up0, X)
Amat = X -> A(X,Qtest,capE,capF,capH,capG,capK)

dT = X -> deltaTest(X,Up0,Up1prime,Amat;positive=true)

X = 10 .^ range(0, 2, length=1000)
plot(X,dT.(X),xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3))

#Uprime = X -> ForwardDiff.derivative(U, [X])

Up20p20p18 = X -> generate_Umatrix(X,Qtest,20,20,capE,capF,capH,capG,capK)[1] #Confirmed identical with Mathematica values

us = []
for i=1:5

    k=4*i

    Up000 = X -> generate_Umatrix(X,Qtest,k,k,capE,capF,capH,capG,capK)[1]
    Upprime00(X) = ForwardDiff.derivative(Up000, X)
    Amat00 = X -> A(X,Qtest,capE,capF,capH,capG,capK)
    dT = X -> deltaTest(X,Up000,Upprime00,Amat00;positive=true)

    X = 10 .^ range(0, 2, length=1000)
    push!(us, dT.(X))
end
plot(X,us,xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3))



#Artificially truncating asymptotic series from Glasser 2020:
##########################################################################################################################################################
capE = -3.369*10^(-2);
capF = 2.409*10^(-3);
capH = 1.292*10^(-2);
capG = 8.950*10^1;
capK = 2.332*10^2;
Qtest = 1.234*10^(-1);
k1 = 16
k2 = 16
Q = Qtest

U_02 = X -> generate_Umatrix(X,Qtest,20,20,capE,capF,capH,capG,capK;truncate_terms=[2,1,1])[1]; #Confirmed identical with Mathematica values
U_182 = X -> generate_Umatrix(X,Qtest,20,20,capE,capF,capH,capG,capK;truncate_terms=[-1,-1,-1])[1]; #Confirmed identical with Mathematica values
Amat = X -> A(X,Q,capE,capF,capH,capG,capK);
U_02_p = XX -> ForwardDiff.derivative(U_02, XX);
U_182_p = XX -> ForwardDiff.derivative(U_182, XX);

DT_02 = XX -> deltaTest(XX,U_02,U_02_p,Amat;positive=true);
DT_182 = XX -> deltaTest(XX,U_182,U_182_p,Amat;positive=true);

Xvec = 10 .^ range(0, 2, length=1000);
plot([Xvec,Xvec],[DT_02.(Xvec),DT_182.(Xvec)],xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3),label=["0"    "full"])


#Recreating Glasser 2020 asymptotic Q-scaling:
##########################################################################################################################################################
capE = -3.369*10^(-2);
capF = 2.409*10^(-3);
capH = 1.292*10^(-2);
capG = 8.950*10^1;
capK = 2.332*10^2;
Qtest = 1.234*10^(-1);

xmaxTest_vec = []
vecoos = []
vecoos2 = []
for i=1:4
    q = 10^-2*(10^(i-1))
    spot = findXmax(q,10,10,capE,capF,capH,capG,capK;search_range = [0,2.5], verbose = false,grain=500)
    push!(xmaxTest_vec,spot[1])
    push!(vecoos,spot[2])
    push!(vecoos2,spot[3])
end
plot(vecoos2,vecoos,xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3),label=["+"    "-"])

##########################################################################################################################################################

#Recreating Glasser 2016 inner solution ratios
##########################################################################################################################################################
capF = 1.074*10^(-6);
capH = 1.801*10^(-4);
capG = 143.9;
capK = 5.191*10^5;
capM = 15.198
Qtest = 1.234*10^(-1);
k1 = 16
k2 = 16
Q = Qtest

Dr = -0.1
capE = Dr-capF-capH^2

findXmax_QE(E) = Q -> findXmax(Q,18,18,E,capF,capH,capG,capK; search_range = [0,5], grain = 2000, convergence_tol = 10.0^(-7), verbose = true, truncate_terms=[-1,-1,-1])
#findXmax_QE(0.1)(10^-4) = 1.925
#findXmax_QE(0.1)(10^3) = 22216
#findXmax_QE(-0.1)(10^-4) = 1.950
#findXmax_QE(-0.1)(10^3) = 18077.091833207232 #tol=1e-7
xmax = 100

inner_ratios_capE_sus = E -> generate_inner_ratios(xmax,E,capF,capH,capG,capK; k1=20,k2=18, N=10000, truncate_terms=[-1;-1;-1]) 

Δs_E0p1 = inner_ratios_capE_sus(0.1)
Δs_E0p06 = inner_ratios_capE_sus(0.06)
Δs_E0p02 = inner_ratios_capE_sus(0.02)
Δs_E0pm02 = inner_ratios_capE_sus(-0.02)
Δs_E0pm06 = inner_ratios_capE_sus(-0.06)
Δs_E0pm1 = inner_ratios_capE_sus(-0.1)

QvecNew = 10 .^ range(-4,3,length=100)
Qvec_ratio_pos = []
Qvec_ratio_neg = []
for i in 1:length(QvecNew)
    push!(Qvec_ratio_pos,[Δs_E0p1[1](QvecNew[i])    Δs_E0p06[1](QvecNew[i]) Δs_E0p02[1](QvecNew[i]) Δs_E0pm02[1](QvecNew[i])    Δs_E0pm06[1](QvecNew[i])    Δs_E0pm1[1](QvecNew[i])])
    push!(Qvec_ratio_neg,[Δs_E0p1[2](QvecNew[i])    Δs_E0p06[2](QvecNew[i]) Δs_E0p02[2](QvecNew[i]) Δs_E0pm02[2](QvecNew[i])    Δs_E0pm06[2](QvecNew[i])    Δs_E0pm1[2](QvecNew[i])])
    print(i,"\n")
end

scatter([QvecNew;QvecNew;QvecNew;QvecNew;QvecNew;QvecNew],[[real(i[1][1]) for i in Qvec_ratio_pos];[real(i[2][1]) for i in Qvec_ratio_pos];[real(i[3][1]) for i in Qvec_ratio_pos];[real(i[4][1]) for i in Qvec_ratio_pos];[real(i[5][1]) for i in Qvec_ratio_pos];[real(i[6][1]) for i in Qvec_ratio_pos]], xaxis=:log, ylims = (-1.8,1.7))

scatter([QvecNew;QvecNew;QvecNew;QvecNew;QvecNew;QvecNew],[[abs(i[1][1]) for i in Qvec_ratio_neg];[abs(i[2][1]) for i in Qvec_ratio_neg];[abs(i[3][1]) for i in Qvec_ratio_neg];[abs(i[4][1]) for i in Qvec_ratio_neg];[abs(i[5][1]) for i in Qvec_ratio_neg];[abs(i[6][1]) for i in Qvec_ratio_neg]], xaxis=:log, yaxis=:log)

##########################################################################################################################################################
generate_D_Q_E = E -> generate_D_Q(xmax,0.2,0.2,capE,capF,capH,capG,capK; k1=18,k2=18, N=10000, truncate_terms=[-1;-1;-1]) 

QvecNew = 10 .^ range(-4,1,length=100)
zero_spot = []
for i in 1:length(QvecNew)
    push!(zero_spot,generate_D_Q_E(-0.1)(QvecNew[i]))
    print(i,"\n")
end

scatter(QvecNew,zero_spot,xaxis=:log)
argmin(abs.(zero_spot))

##########################################################################################################################################################
ctour = contour_inner(0.02983647240283-0.003,0.02983647240283+0.003;args=[pi/2,-pi/2],in_res=100,out_res=100,rad_res=100)
M_contour = countour_Q(xmax,0.2,0.2,ctour,capE,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1])
#ctour2 = contour_inner(0.001,10;args=[pi/2,-pi/2],in_res=100,out_res=100,rad_res=100)
#M_contour2 = countour_Q(100.0,0.5,0.5,ctour2,0.1,capF,capH,capG,capK; k1=18,k2=18, N=2000, truncate_terms=[-1;-1;-1])

plot(real(ctour),imag(ctour))
plot(real(M_contour),imag(M_contour))
plot(real(ctour2),imag(ctour2))
plot(real(M_contour2),imag(M_contour2))

