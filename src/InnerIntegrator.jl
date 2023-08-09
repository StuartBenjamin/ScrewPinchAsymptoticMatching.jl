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
#Inner equation integrator (Glasser, Jardin & Tesauro 1984, finite difference method taken from section III.)
    #TESTED
##########################################################################################################################################################
U(X,Q,capE,capF,capG,capK) = 
[Q -X*Q 0;
(-X/Q) (X^2/Q) -(capE + capF)*Complex(Q)^(-2);
(-X/Q) (-(capG - capK*capE)*Q) (X^2/Q + (capG + capK*capF)*Q)]
Uprime(X,Q) =
[0 -Q 0;
-(1/Q) (2*X)/Q 0;
-(1/Q) 0 (2*X)/Q]
Udprime(Q) = 
[0 0 0;
0 2/Q 0;
0 0 2/Q]
V(Q,capH,capK) = [0 0 capH;
        -capH*Complex(Q)^(-2) 0 0;
        capH*capK*Q 0 0]

h(N,Xmax) = Xmax/(N - 1)
Xi(i,h) = (i - 1)*h

id3 = [1 0 0;0 1 0;0 0 1]

Ai(h,Xi,Q,capE,capF,capH,capG,capK) = id3 - (h/2)*V(Q,capH,capK) - (h^3/24)*(2*Uprime(Xi, Q) + U(Xi,Q,capE,capF,capG,capK)*V(Q,capH,capK) - V(Q,capH,capK)*U(Xi,Q,capE,capF,capG,capK) - V(Q,capH,capK)*V(Q,capH,capK)*V(Q,capH,capK))
Ci(h,Xi,Q,capE,capF,capH,capG,capK) = id3 + (h/2)*V(Q,capH,capK) + (h^3/24)*(2*Uprime(Xi, Q) + U(Xi,Q,capE,capF,capG,capK)*V(Q,capH,capK) - V(Q,capH,capK)*U(Xi,Q,capE,capF,capG,capK) - V(Q,capH,capK)*V(Q,capH,capK)*V(Q,capH,capK))
Bi(h,Xi,Q,capE,capF,capH,capG,capK) = -2*id3 - h^2*U(Xi,Q,capE,capF,capG,capK) - (h^4/12)*(Udprime(Q) - V(Q,capH,capK)*Uprime(Xi, Q) - V(Q,capH,capK)*V(Q,capH,capK)*U(Xi,Q,capE,capF,capG,capK) + U(Xi,Q,capE,capF,capG,capK)*U(Xi,Q,capE,capF,capG,capK))
Ai(i,N,Xmax,Q,capE,capF,capH,capG,capK) = Ai(h(N,Xmax),Xi(i,h(N,Xmax)),Q,capE,capF,capH,capG,capK)
Bi(i,N,Xmax,Q,capE,capF,capH,capG,capK) = Bi(h(N,Xmax),Xi(i,h(N,Xmax)),Q,capE,capF,capH,capG,capK)
Ci(i,N,Xmax,Q,capE,capF,capH,capG,capK) = Ci(h(N,Xmax),Xi(i,h(N,Xmax)),Q,capE,capF,capH,capG,capK)

sigmaMatrix(signVar,N,Xmax,Q,capE,capF,capH,capG,capK) = 
[(Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[1, 1] + signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[1, 1])  (Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[1, 2] - signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[1, 2])   (Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[1, 3] - signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[1, 3]);
(Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[2, 1] + signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[2, 1])   (Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[2, 2] - signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[2, 2])   (Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[2, 3] - signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[2, 3]);
(Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[3, 1] + signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[3, 1])   (Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[3, 2] - signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[3, 2])   (Ai(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[3, 3] - signVar*Ci(h(N,Xmax),0,Q,capE,capF,capH,capG,capK)[3, 3])]

function compute_EN_minus_One(signVar,N,Xmax,Q,capE,capF,capH,capG,capK)
    E1 = -inv(Bi(h(N,Xmax),0,Q,capE,capF,capH,capG,capK))*sigmaMatrix(signVar,N,Xmax,Q,capE,capF,capH,capG,capK)

    Ei=E1
    for i=2:(N-1)
    Ei = -inv(Bi(i,N,Xmax,Q,capE,capF,capH,capG,capK)+Ci(i,N,Xmax,Q,capE,capF,capH,capG,capK)*Ei)*Ai(i,N,Xmax,Q,capE,capF,capH,capG,capK)
    end

    return Ei
end

#UNTESTED
function compute_Es(signVar,N,Xmax,Q,capE,capF,capH,capG,capK) 
    E1 = -inv(Bi(h(N,Xmax),0,Q,capE,capF,capH,capG,capK))*sigmaMatrix(signVar,N,Xmax,Q,capE,capF,capH,capG,capK)

    Ei=E1

    Es = []

    push!(Es,E1)
    for i=2:(N-1)
        Ei = -inv(Bi(i,N,Xmax,Q,capE,capF,capH,capG,capK)+Ci(i,N,Xmax,Q,capE,capF,capH,capG,capK)*Ei)*Ai(i,N,Xmax,Q,capE,capF,capH,capG,capK)
        push!(Es,Ei)
    end

    return Ei, Es
end


Dmat(EnMinusOne,N,Xmax,Q,capE,capF,capH,capG,capK) = Dmat0(EnMinusOne, h(N,Xmax), V(Q,capH,capK), U(Xi(N,h(N,Xmax)),Q,capE,capF,capG,capK), Uprime(Xi(N,h(N,Xmax)), Q), Udprime(Q))
Dmat0(EnMinusOne, h, V, Un, UnPrime, UndPrime) = (1/h)*inv(id3 - (1/2)*h*V + (1/6)*h^2*(V*V + Un) - (1/24)*
h^3*(2*UnPrime + Un*V + V*Un + V*V*V))*(id3 - EnMinusOne + (1/2)*h^2*Un - (1/6)*
h^3*(UnPrime + V*Un) + (1/24)*h^4*(UndPrime + V*UnPrime + V*V*Un + Un*Un))

T3(umatrix, X) = [umatrix(X)[1, 1]/X; umatrix(X)[2, 1]; umatrix(X)[3, 1]]
T3prime(umatrix, X) = [(umatrix(X)[4, 1] - umatrix(X)[1, 1]*X^(-2)); umatrix(X)[5, 1]*X; umatrix(X)[6, 1]*X]
T4(umatrix, X) = [umatrix(X)[1, 2]/X; umatrix(X)[2, 2]; umatrix(X)[3, 2]]
T4prime(umatrix, X) = [(umatrix(X)[4, 2] - umatrix(X)[1, 2]*X^(-2)); umatrix(X)[5, 2]*X; umatrix(X)[6, 2]*X]

alpha3_on_alpha4(umatrix,EnMinusOne,N,Xmax,Q,capE,capF,capH,capG,capK;var=2) = (-(T4prime(umatrix, Xmax) - Dmat(EnMinusOne,N,Xmax,Q,capE,capF,capH,capG,capK)*T4(umatrix, Xmax))./(T3prime(umatrix, Xmax) - Dmat(EnMinusOne,N,Xmax,Q,capE,capF,capH,capG,capK)*T3(umatrix, Xmax)))[var]    


Ψ_plus(umatrix, X) = T3(umatrix, X)  #UNTESTED
Ψ_minus(umatrix, X) = T4(umatrix, X)  #UNTESTED

#ce = co+cl
#co*Δo(Q)=(ce)*Δe(Q)-cl*Δl
#co*Δo(Q)=(co+cl)*Δe(Q)-cl*Δl
#co*(Δo(Q)-Δe(Q))=cl*Δe(Q)-cl*Δl
##co = cl*(Δe(Q)-Δl)/(Δo(Q)-Δe(Q))
##ce = co+cl

co_(Δe,Δo,cl,cr) = cl*(Δe(Q)-Δl)/(Δo(Q)-Δe(Q)) #UNTESTED
ce_(Δe,Δo,cl,cr) = co(Δe,Δocl,cr)+cl  #UNTESTED

#UNTESTED
function compute_Ψs_inner(umatrix,xmatch,cl,cr,Δl,Δr,Δo,Δe,N,Xmax,Q,capE,capF,capH,capG,capK)
    Es_o = compute_Es(1,N,Xmax,Q,capE,capF,capH,capG,capK)[2] #Final one is _EN_minus_One
    Es_e = compute_Es(-1,N,Xmax,Q,capE,capF,capH,capG,capK)[2] #Final one is _EN_minus_One

    co = co_(Δe,Δo,cl,cr)
    ce = ce_(Δe,Δo,cl,cr)

    ΨL_o = (((ce*cl/co) + cr/co)/(ce+1))*Ψ_minus(umatrix, xmatch) + (((ce*cl/co)*Δl + (cr/co)*Δr)/(ce+1))*Ψ_plus(umatrix, xmatch)
    ΨL_e = ((cr*Ψ_minus(umatrix, xmatch) + cr*Δr*Ψ_plus(umatrix, xmatch)) - co*Ψ_o)/ce

    Ψ_o = []
    Ψ_e = []

    push!(Ψ_o,ΨL_o)
    push!(Ψ_e,ΨL_e)

    for i in 1:length(Es_e)
        push!(Ψ_o,Es_o[length(Es_o)-(i-1)]*Ψ_o[end])
        push!(Ψ_e,Es_e[length(Es_e)-(i-1)]*Ψ_e[end])
    end

    return reverse(Ψ_o),reverse(Ψ_e)
end

"""
function plot_Ψs_inner(N,umatrix,xmatch,cl,cr,Δl,Δr,Δo,Δe,N,Xmax,Q,capE,capF,capH,capG,capK;var=2)
    Ψ_o,Ψ_e = compute_Ψs_inner(N,umatrix,xmatch,cl,cr,Δl,Δr,Δo,Δe,N,Xmax,Q,capE,capF,capH,capG,capK)
    
    ξ_o = [i[var] for i in Ψ_o]
    ξ_e = [i[var] for i in Ψ_e]

    h0 = h(N,xmatch)
    xvec = [Xi(i,h0) for i in 1:length(Ψ_o)]

end
"""

"""
    Up1prime(xmatch)[1:6,1]

    co*[o]+ce*[e]=cr*[-]+cr*Δr*[+]
    ^have ^have   ^have  ^have ^have
                    ^have  ^have 
    co*[o]-ce*[e]=cl*[-]+cl*Δl*[+]
    ^have   ^have  ^have  ^have ^have
                    ^have  ^have 

    [e]_ξ = [(cr*[-]_ξ + cr*Δr*[+]_ξ) - co*[o]_ξ]/ce
    => [o]_ξ = (cl*[-]_ξ + cl*Δl*[+]_ξ + ce*[e]_ξ)/co
        = (cl*[-]_ξ + cl*Δl*[+]_ξ + [(cr*[-]_ξ + cr*Δr*[+]_ξ) - co*[o]_ξ]/ce)/co

    => [o]_ξ = (cl/co)*[-]_ξ + (cl/co)*Δl*[+]_ξ + (cr/(ce*co)*[-]_ξ + (cr/(ce*co))*Δr*[+]_ξ)-(1/ce)*[o]_ξ
    => [o]_ξ+(1/ce)*[o]_ξ = (cl/co)*[-]_ξ + (cl/co)*Δl*[+]_ξ + cr/(ce*co)*[-]_ξ + (cr/(ce*co))*Δr*[+]_ξ
    => [o]_ξ = ((cl/co)*[-]_ξ + (cl/co)*Δl*[+]_ξ + cr/(ce*co)*[-]_ξ + (cr/(ce*co))*Δr*[+]_ξ)/(1+1/ce)
    => [o]_ξ = ((cl/co)*[-]_ξ + cr/(ce*co)*[-]_ξ + (cl/co)*Δl*[+]_ξ + (cr/(ce*co))*Δr*[+]_ξ)/(1+1/ce)
    => [o]_ξ = (((ce*cl/co) + cr/co)/(ce+1))*[-]_ξ + (((ce*cl/co)*Δl + (cr/co)*Δr)/(ce+1))*[+]_ξ
    => [o] = (((ce*cl/co) + cr/co)/(ce+1))*[-] + (((ce*cl/co)*Δl + (cr/co)*Δr)/(ce+1))*[+]
    => [e] = ((cr*[-] + cr*Δr*[+]) - co*[o])/ce
"""