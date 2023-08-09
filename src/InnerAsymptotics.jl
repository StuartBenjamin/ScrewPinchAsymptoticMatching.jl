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
#Inner equation asymptotic solution (Glasser & Wang 2020)
    #TESTED
##########################################################################################################################################################
A0(Q,capH) = 
    [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 -Q 0 0 0 capH;
    0 1/Q 0 0 0 0;
    0 0 1/Q 0 0 0]

A1(Q,capE,capF,capH,capG,capK) = 
    [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    Q 0 0 1 0 0;
    -1/Q 0 -(capE + capF)/(Complex(Q)^2) -capH/(Complex(Q)^2) -1 0;
    -1/Q -(capG - capK*capE)*Q (capG + capK*capF)*Q capH*capK*Q 0 -1]

A2(Q,capH,capK) = 
    [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    -2 0 0 0 0 0;
    capH/(Complex(Q)^2) 0 0 0 0 0;
    -capH*capK*Q 0 0 0 0 0]

T(Q,capH) = 
    [1 0 capH*Q Complex(Q)^(5/2) capH*Q -Complex(Q)^(5/2);
    0 0 0 -Complex(Q)^(1/2) 0 Complex(Q)^(1/2);
    0 0 -Complex(Q)^(1/2) 0 Complex(Q)^(1/2) 0;
    0 1 -capH*Complex(Q)^(1/2) -Complex(Q)^2 capH*Complex(Q)^(1/2) -Complex(Q)^2;
    0 0 0 1 0 1;
    0 0 1 0 1 0]

A(X,Q,capE,capF,capH,capG,capK) = A0(Q,capH) + X^(-2)*A1(Q,capE,capF,capH,capG,capK) + X^(-4)*A2(Q,capH,capK)

Jk0(Q,capE,capF,capH,capG,capK) = 
        [0 1 0 0 0 0; 
        0 0 0 0 0 0; 
        0 0 -(1/sqrt(Complex(Q))) 0 0 0;
        0 0 0 -(1/sqrt(Complex(Q))) 0 0;
        0 0 0 0 1/sqrt(Complex(Q)) 0;
        0 0 0 0 0 1/sqrt(Complex(Q))]

Jk1(Q,capE,capF,capH,capG,capK) = 
        [(capH)  (-capH^2*capK*Complex(Q)^2)  (capH*Q*(1 + capH + capH^2*capK*Complex(Q)^(3/2) + (capG + capF*capK)*Complex(Q)^(3/2)))  (capH*Complex(Q)^(5/2)*(1 - capG + capE*capK + capH*capK*Complex(Q)^(3/2)))  (capH*Q*(1 + capH - capH^2*capK*Complex(Q)^(3/2) - (capG + capF*capK)*Complex(Q)^(3/2)))  (capH*Complex(Q)^(5/2)*(-1 + capG - capE*capK + capH*capK*Complex(Q)^(3/2)));
        (0)  (1-capH)  ((capE + capF + (-1 + capH)*capH)*sqrt(Complex(Q)))  ((-2 + capH)*Complex(Q)^2)  (-((capE + capF + (-1 + capH)*capH)*sqrt(Complex(Q))))  ((-2 + capH)*Complex(Q)^2);
        (-(1/(2*Q)))  ((capH*capK*Q)/2)  ((1/2)*(-1 - capH - capH^2*capK*Complex(Q)^(3/2) - (capG + capF*capK)*Complex(Q)^(3/2)))  (-(1/2)*Complex(Q)^(3/2)*(1 - capG + capE*capK + capH*capK*Complex(Q)^(3/2)))  ((1/2)*(-1 - capH + capG*Complex(Q)^(3/2) + capF*capK*Complex(Q)^(3/2) + capH^2*capK*Complex(Q)^(3/2)))  (-(1/2)*Complex(Q)^(3/2)*(-1 + capG - capE*capK + capH*capK*Complex(Q)^(3/2)));
        (-(1/(2*Q)))  (-(capH/(2*Complex(Q)^2)))  ((capE + capF + capH*(capH - Complex(Q)^(3/2)))/(2*Complex(Q)^(3/2)))  ((1/2)*(-1 + capH - Complex(Q)^(3/2)))  ((-capE - capF - capH*(capH + Complex(Q)^(3/2)))/(2*Complex(Q)^(3/2)))  ((1/2)*(-1 + capH + Complex(Q)^(3/2)));
        (-(1/(2*Q)))  ((capH*capK*Q)/2)  ((1/2)*(-1 - capH - capH^2*capK*Complex(Q)^(3/2) - (capG + capF*capK)*Complex(Q)^(3/2)))  (-(1/2)*Complex(Q)^(3/2)*(1 - capG + capE*capK + capH*capK*Complex(Q)^(3/2)))  ((1/2)*(-1 - capH + capG*Complex(Q)^(3/2) + capF*capK*Complex(Q)^(3/2) + capH^2*capK*Complex(Q)^(3/2)))  (-(1/2)*Complex(Q)^(3/2)*(-1 + capG - capE*capK + capH*capK*Complex(Q)^(3/2)));
        (-(1/(2*Q)))  (-(capH/(2*Complex(Q)^2)))  ((capE + capF + capH*(capH - Complex(Q)^(3/2)))/(2*Complex(Q)^(3/2)))  ((1/2)*(-1 + capH - Complex(Q)^(3/2)))  ((-capE - capF - capH*(capH + Complex(Q)^(3/2)))/(2*Complex(Q)^(3/2)))  ((1/2)*(-1 + capH + Complex(Q)^(3/2)))]

Jk2(Q,capE,capF,capH,capG,capK) = 
        [(capH^2*capK*Complex(Q)^2)  ( 0)  ( capH^3*capK*Complex(Q)^3)  ( capH^2*capK*Complex(Q)^(9/2))  ( capH^3*capK*Complex(Q)^3)  ( -capH^2*capK*Complex(Q)^(9/2));
        ((-2 + capH))  ( 0)  ( (-2 + capH)*capH*Q)  ( (-2 + capH)*Complex(Q)^(5/2))  ( (-2 + capH)*capH*Q)  ( -((-2 + capH)*Complex(Q)^(5/2)));
        (-(1/2)*capH*capK*Q)  ( 0)  ( -(1/2)*capH^2*capK*Complex(Q)^2)  ( -(1/2)*capH*capK*Complex(Q)^(7/2))  ( -(1/2)*capH^2*capK*Complex(Q)^2)  ( (1/2)*capH*capK*Complex(Q)^(7/2));
        ( capH/(2*Complex(Q)^2))  ( 0)  ( capH^2/(2*Q))  ( (capH*sqrt(Complex(Q)))/2)  ( capH^2/(2*Q))  ( -((capH*sqrt(Complex(Q)))/2));
        ( -(1/2)*(capH*capK*Q))  ( 0)  ( -(1/2)*capH^2*capK*Complex(Q)^2)  ( -(1/2)*capH*capK*Complex(Q)^(7/2))  ( -(1/2)*capH^2*capK*Complex(Q)^2)  ( (1/2)*capH*capK*Complex(Q)^(7/2));
        ( capH/(2*Complex(Q)^2))  ( 0)  ( capH^2/(2*Q))  ( (capH*sqrt(Complex(Q)))/2)  ( capH^2/(2*Q))  ( -((capH*sqrt(Complex(Q)))/2))]

function Jk(k::Int,Q,capE,capF,capH,capG,capK)
    if k==0
        return Jk0(Q,capE,capF,capH,capG,capK)
    elseif k==1
        return Jk1(Q,capE,capF,capH,capG,capK) 
    elseif k==2 
        return Jk2(Q,capE,capF,capH,capG,capK)
    else return [0 0 0 0 0 0;
                0 0 0 0 0 0; 
                0 0 0 0 0 0; 
                0 0 0 0 0 0; 
                0 0 0 0 0 0; 
                0 0 0 0 0 0]
    end
end

Jk11(k::Int,Q,capE,capF,capH,capG,capK) = Jk(k,Q,capE,capF,capH,capG,capK)[1:2,1:2]
Jk22(k::Int,Q,capE,capF,capH,capG,capK) = Jk(k,Q,capE,capF,capH,capG,capK)[3:6,3:6]
Jk12(k::Int,Q,capE,capF,capH,capG,capK) = Jk(k,Q,capE,capF,capH,capG,capK)[1:2,3:6]
Jk21(k::Int,Q,capE,capF,capH,capG,capK) = Jk(k,Q,capE,capF,capH,capG,capK)[3:6,1:2]

#Pk21[Q_, k_] := LyapunovSolve[Jk22[Q, 0], -Jk11[Q, 0], -Kk21[Q, k]]
#Pk21[Q_, 0] := IdentityMatrix[6][[3 ;; 6, 1 ;; 2]]
#Bk11[Q_, k_] := Jk11[Q, k] + Sum[Dot[Jk12[Q, k - l], Pk21[Q, l]], {l, 1, k - 1}]
#Kk21[Q_, k_] :=  Jk21[Q, k] + 2*(k - 1)*Pk21[Q, k - 1] +  Sum[Dot[Jk22[Q, k - l], Pk21[Q, l]] - Dot[Pk21[Q, l], Bk11[Q, k - l]], {l, 1, k - 1}]

function recursive_1(ktot::Int,Q,capE,capF,capH,capG,capK)
    Bk11s = []
    Pk21s = []
    Kk21s = []

    push!(Bk11s, Jk11(0,Q,capE,capF,capH,capG,capK))
    push!(Pk21s, [0 0;0 0;0 0;0 0])

    for k=1:ktot
        if k==1
            push!(Bk11s, Jk11(1,Q,capE,capF,capH,capG,capK))
            push!(Kk21s, Jk21(1,Q,capE,capF,capH,capG,capK))
        else
            internal_sum_1 = [0 0;0 0]
            internal_sum_2 = [0 0;0 0;0 0;0 0]
            for l=1:k-1
                internal_sum_1 += Jk12(k-l,Q,capE,capF,capH,capG,capK)*Pk21s[l+1]
                internal_sum_2 += (Jk22(k-l,Q,capE,capF,capH,capG,capK)*Pk21s[l+1]-Pk21s[l+1]*Bk11s[k-l+1])
            end
            push!(Bk11s, Jk11(k,Q,capE,capF,capH,capG,capK)+internal_sum_1)
            push!(Kk21s, Jk21(k,Q,capE,capF,capH,capG,capK)+2*(k-1)*Pk21s[k-1+1]+internal_sum_2)
        end
        push!(Pk21s,sylvc(Jk22(0,Q,capE,capF,capH,capG,capK),  -Jk11(0,Q,capE,capF,capH,capG,capK),  -Kk21s[end]))
    end

    return Bk11s, Pk21s
end

function p62(X,Pk21s;truncate_terms=-1)
    if (truncate_terms==-1)
        imax = length(Pk21s)
    else
        imax = min(truncate_terms,length(Pk21s))
    end
    
    p62mat = [1 0;0 1;0 0;0 0;0 0;0 0]
    for i in 1:imax
        p62mat += X^(-2*(i-1))*vcat([0 0;0 0],Pk21s[i])
    end
    return p62mat
end

function Bmatrix(X,Bk11s)
    Bmat = [0 0;0 0]
    for i in 1:length(Bk11s)
        Bmat += X^(-2*(i-1))*Bk11s[i]
    end
    return Bmat
end

#Kk[Q_, k_] :=  Bk11[Q, k] + 2*(k - 1)*Qk[Q, k - 1] + Sum[Dot[Bk11[Q, k - l], Qk[Q, l]] - Dot[Qk[Q, l], Ck[Q, k - l]], {l, 1, k - 1}]
#Qk[Q_, k_] := -{{0, 0}, {Kk[Q, k][[1, 1]], Kk[Q, k][[1, 2]]}} 
#Ck[Q_, k_] := {{0, 0}, {Kk[Q, k][[2, 1]], Kk[Q, k][[1, 1]] + Kk[Q, k][[2, 2]]}}

function recursive_2(ktot::Int,Bk11s,Q,capE,capF,capH,capG,capK)
    if (ktot+1) > length(Bk11s)
        ktot = length(Bk11s)-1
    end
    
    Kks = []
    Qks = []
    Cks = []

    push!(Qks, [1 0;0 1])
    push!(Cks, Bk11s[1])

    for k=1:ktot
        if k==1
            push!(Kks, Bk11s[k+1])
        else
            internal_sum = [0 0;0 0]
            for l=1:k-1
                internal_sum += (Bk11s[k-l+1]*Qks[l+1]-Qks[l+1]*Cks[k-l+1])
            end
            push!(Kks, Bk11s[k+1] + 2*(k - 1)*Qks[k-1+1] + internal_sum)
        end

        push!(Qks,-1*[0 0;Kks[end][1,1] Kks[end][1,2]])
        push!(Cks, [0 0;Kks[end][2,1] (Kks[end][1,1]+Kks[end][2,2])])
    end

    return Qks, Cks
end

function Qmatrix(X,Qks;truncate_terms=-1)
    if (truncate_terms==-1)
        imax = length(Qks)
    else
        imax = min(truncate_terms,length(Qks))
    end

    Qmat = [0 0;0 0]
    for i in 1:imax
        Qmat += X^(-2*(i-1))*Qks[i]
    end
    return Qmat
end

S(X) = [1 0; 0 1/(X^2)]
Sinv(X) = [1 0; 0 X^2]
Sprime(X) = [0 0; 0 -2/(X^3)]

#Dk[Q_, 0] := {{0, 1}, {Ck[Q, 2][[2, 1]], 3}}
#Dk[Q_, k_] := {{0, 0}, {Ck[Q, k + 2][[2, 1]], Ck[Q, k + 1][[2, 2]]}}

function Generate_Dks(Cks)
    Dks = []
    push!(Dks, [0 1;Cks[2+1][2,1] 3])

    for k in 1:(length(Cks)-3)
        push!(Dks, [0 0;Cks[k+2+1][2,1] Cks[k+1+1][2,2]])
    end

    return Dks
end

Di(capE,capF,capH) = capE + capF + capH - 1/4

rp(capE,capF,capH) = (3/2) + sqrt(-Di(capE,capF,capH))
rm(capE,capF,capH) = (3/2) - sqrt(-Di(capE,capF,capH))

Rmat(capE,capF,capH) = [rp(capE,capF,capH) 0;0 rm(capE,capF,capH)]
Y0(capE,capF,capH) = [1 1;rp(capE,capF,capH) rm(capE,capF,capH)]
Y0Inv(capE,capF,capH) = [rm(capE,capF,capH)/(rm(capE,capF,capH) - rp(capE,capF,capH)) -(1/(rm(capE,capF,capH) - rp(capE,capF,capH))); -(rp(capE,capF,capH)/(rm(capE,capF,capH) - rp(capE,capF,capH))) 1/(rm(capE,capF,capH) - rp(capE,capF,capH))]

function Generate_Eks(Dks,Y0,Y0Inv)
    Eks = []

    for i in Dks
        push!(Eks, Y0Inv*i*Y0)
    end

    return Eks
end

#Zk[Q_, 0] := IdentityMatrix[2]
#Zk[Q_, k_] :=  LyapunovSolve[-R, (R - 2*k*IdentityMatrix[2]), Sum[Dot[Ek[Q, l], Zk[Q, k - l]], {l, 1, k}]]

function Generate_Zks(Eks,Rmat0)
    Zks = []

    push!(Zks, [1 0;0 1])

    for k=1:(length(Eks)-1)
        internal_sum=[0 0;0 0]

        for l = 1:k
            internal_sum+= Eks[l+1]*Zks[k-l+1]
        end

        push!(Zks, sylvc(-Rmat0, (Rmat0-2*k*[1 0;0 1]), internal_sum))
    end

    return Zks
end

function Zmatrix(X,Zks,capE,capF,capH;truncate_terms=-1) #CHECK 
    if (truncate_terms==-1)
        imax = length(Zks)
    else
        imax = min(truncate_terms,length(Zks))
    end
    
    Zmat0 = [X^rp(capE,capF,capH) 0;0 X^rm(capE,capF,capH)]
    Zmat = [0 0;0 0]

    for i in 1:imax
        Zmat += X^(-2*(i-1))*Zks[i]*Zmat0
    end
    return Zmat
end

function ZmatrixN2(X,R,Zks,capE,capF,capH) 
    Zmat = [0 0;0 0]
    spot = []

    for i in 1:length(Zks)
        push!(spot,Zks[i]*(exp(log(X)*(R - 2*(i-1)*[1 0;0 1]))))
        Zmat += Zks[i]*(exp(log(X)*(R - 2*(i-1)*[1 0;0 1])))
    end
    return Zmat, spot
end

Ymatrix(X,Zks,capE,capF,capH;truncate_terms=-1) = Y0(capE,capF,capH)*Zmatrix(X,Zks,capE,capF,capH;truncate_terms=truncate_terms)

function generate_Umatrix(X,Q,k1::Int,k2::Int,capE,capF,capH,capG,capK; truncate_terms=[-1;-1;-1])

    Bk11s, Pk21s = recursive_1(k1,Q,capE,capF,capH,capG,capK)

    Qks, Cks = recursive_2(k2,Bk11s,Q,capE,capF,capH,capG,capK)

    Dks = Generate_Dks(Cks)
    Eks = Generate_Eks(Dks,Y0(capE,capF,capH),Y0Inv(capE,capF,capH))
    Zks = Generate_Zks(Eks,Rmat(capE,capF,capH))

    Tmat = T(Q,capH)
    p62mat = p62(X,Pk21s;truncate_terms=truncate_terms[1])
    Qmat = Qmatrix(X,Qks;truncate_terms=truncate_terms[2])
    Smat = S(X)
    Ymat = Ymatrix(X,Zks,capE,capF,capH;truncate_terms=truncate_terms[3])

    return Tmat*p62mat*Qmat*Smat*Ymat, [length(Pk21s),length(Qks),length(Zks)], [Pk21s,Qks,Zks], [Tmat,p62mat,Qmat,Smat,Ymat]
end
#P62[Q_, k_, X_] :=  IdentityMatrix[6][[1 ;; 6, 1 ;; 2]] +   Sum[X^(-2*kprime)*ArrayFlatten[{{{0, 0}, {0, 0}}, Pk21[Q, kprime]}, 1], {kprime, 1, k}]
#B[Q_, k_] := Sum[X^(-2*kprime)*Bk11[Q, kprime], {kprime, 0, k}]

function deltaTest(X,Up1,Up1prime,A1;positive=true)
    if positive
        d1 = Up1prime(X)[1:6,1]
        d2 = (X*A1(X)*Up1(X))[1:6,1]
    else
        d1 = Up1prime(X)[1:6,2]
        d2 = (X*A1(X)*Up1(X))[1:6,2]
    end
    d0 = [0 0 0 0 0 0]

    return euclidean(d1,d2)/max(euclidean(d1,d0),euclidean(d2,d0))
end
function dTfull(X, Upfn,Amat)
    Uprime(XX) = ForwardDiff.derivative(Upfn, XX)

    return [deltaTest(X,Upfn,Uprime,Amat;positive=true) deltaTest(X,Upfn,Uprime,Amat;positive=false)]
end

#Errors sometimes, not very smart..., should use secant method too
function findXmax(Q,k1::Int,k2::Int,capE,capF,capH,capG,capK; search_range = [0,2], grain = 1000, convergence_tol = 10.0^(-7), verbose = false, truncate_terms=[2,1,1])
    umatrix = X -> generate_Umatrix(X,Q,k1,k2,capE,capF,capH,capG,capK;truncate_terms=truncate_terms)[1];
    Amat = X -> A(X,Q,capE,capF,capH,capG,capK);
    uPrime = X -> ForwardDiff.derivative(umatrix, X);
    DT_Pos = XX -> deltaTest(XX,umatrix,uPrime,Amat;positive=true);
    DT_Neg = XX -> deltaTest(XX,umatrix,uPrime,Amat;positive=false);

    Xvec = 10 .^ range(search_range[1],search_range[2], length=grain);
    DT_vec_pos = DT_Pos.(Xvec)
    DT_vec_neg = DT_Neg.(Xvec)

    Xmax_arg_pos = nanargmin(abs.(DT_vec_pos .- convergence_tol))
    Xmax_arg_neg = nanargmin(abs.(DT_vec_neg .- convergence_tol))
    
    Xmax_arg = max(Xmax_arg_pos,Xmax_arg_neg)

    if Xmax_arg==length(Xvec)
        @warn "Series hasn't converged by X = $(Xvec[end]). Increase magnitude of search_range."
    end

    if Xmax_arg==1
        @warn "Series has already converged by X = $(Xvec[1]). Decrease magnitude of search_range."
    end

    if verbose
        print(Xmax_arg_pos," ",Xmax_arg_neg,"\n")
        p1 = plot([Xvec,Xvec],[DT_Pos.(Xvec),DT_Neg.(Xvec)],xaxis=:log, yaxis=:log,ylims=(10^-15,10^0.3),label=["+"    "-"])
        p1 = plot!([Xvec[Xmax_arg]]; seriestype = :vline,label="Xmatch = $(Xmax)")
        display(p1)
    end

    return Xvec[Xmax_arg],[DT_vec_pos[Xmax_arg],DT_vec_neg[Xmax_arg]],[DT_vec_pos,Xvec]
end 

#returns 1 if converged given Q, xmax, and 0 if not
function checkXmax(Q,xmax,k1::Int,k2::Int,capE,capF,capH,capG,capK; convergence_tol = 10.0^(-7), truncate_terms=[2,1,1])
    umatrix = X -> generate_Umatrix(X,Q,k1,k2,capE,capF,capH,capG,capK;truncate_terms=truncate_terms)[1];
    Amat = X -> A(X,Q,capE,capF,capH,capG,capK);
    uPrime = X -> ForwardDiff.derivative(umatrix, X);
    DT_Pos = XX -> deltaTest(XX,umatrix,uPrime,Amat;positive=true);
    DT_Neg = XX -> deltaTest(XX,umatrix,uPrime,Amat;positive=false);

    #Xvec = 10 .^ range(search_range[1],search_range[2], length=grain);
    DT_vec_pos = DT_Pos(xmax)
    DT_vec_neg = DT_Neg(xmax)

    if DT_vec_pos>convergence_tol || DT_vec_neg>convergence_tol
        return false, maximum([DT_vec_pos,DT_vec_neg]), convergence_tol
    else
        return true, maximum([DT_vec_pos,DT_vec_neg]), convergence_tol
    end
end 