#cd("/Users/sbenjamin/Desktop/PHD/Cylindrical Delta Prime Widget/Screw pinch asymptotic matching in Julia")
include("ChandraScrewPinchEquilibrium.jl")
include("OuterIntegrator.jl")

############################################################################################################################################################
#Glasser -> Don't use this since I don't understand, see line 138 below.
############################################################################################################################################################
    #Turning cylinder into Hamada coordinates:
        #V=(2*pi)^2*R0*r^2/2, Grad V = ((2*pi)^2*R0)*r rhat
        #thetaG=theta/(2*pi), Grad thetaG = 1/(2*pi*r) thetahat
        #zeta=z/(2*pi*R0), Grad zeta = 1/(2*pi*R0) zhat 
    ∇V_(R0) = r -> ((2*pi)^2*R0)*r #Evaluated at rs in Glasser '75
    drdV(R0) = r -> 1/∇V_(R0)(r) #(1/((2*pi)^2*R0*r)), needed for derivatives w.r.t V via the chain rule
    B_(Bt,Bp) = r -> sqrt(Bt(r)^2 + Bp(r)^2)
    sig_(Bt,Bp,Jt,Jp) = r -> (Jp(r)*Bp(r) + Jt(r)*Bt(r))/(B_(Bt,Bp)(rs)^2)

    ############################################################################################################################################################
    Ψprime_(Bt,R0) = r -> Bt(r)/(2*pi*R0)
    Ψdprime_(Bt,R0) = r -> ForwardDiff.derivative(Ψprime_(Bt,R0),r)*drdV(R0)(r)
    χprime_(Bp) = r -> Bp(r)/(2*pi*r)
    χdprime_(Bp,R0) = r -> ForwardDiff.derivative(χprime_(Bp),r)*drdV(R0)(r)

    ############################################################################################################################################################
    Jprime_(Jp) = r -> Jp(r)/(2*pi*r)
    Iprime_(Jt,R0) = r -> Jt(r)/(2*pi*R0)

    ############################################################################################################################################################
    ι_GlassOn2Pi(Bt,Bp,R0) = r -> χprime_(Bp)(r)/Ψprime_(Bt,R0)(r) #Should be equal to n/m when r = rs
    ι_GlassOn2Pi_prime(Bt,Bp,R0) = r -> ForwardDiff.derivative(ι_GlassOn2Pi(Bt,Bp,R0),r)*drdV(R0)(r) 
    ι_On2Pi(Bt,Bp,k) = r -> Bp(r)/(k*r*Bt(r)) #Should be the same... #TESTED

    ############################################################################################################################################################
    Λ_(Bt,Bp,R0) = r -> Ψprime_(Bt,R0)(r)^2*ι_GlassOn2Pi_prime(Bt,Bp,R0)(r) #Evaluated at rs in Glasser '75
    Λ2_(Bt,Bp,R0) = r -> Ψprime_(Bt,R0)(r)*χdprime_(Bp,R0)(r) - χprime_(Bp)(r)*Ψdprime_(Bt,R0)(r)  #Should be the same... #TESTED

    ############################################################################################################################################################
    α1(Bt,R0,m) = r -> 2*pi*m/Ψprime_(Bt,R0)(r)
    α2(Bp,n) = r -> 2*pi*n/χprime_(Bp)(r) #Should be the same (at least at rs) #TESTED
    α_(Bp,n,rs) = 2*pi*n/χprime_(Bp)(rs)

    #To clarify:
        #γ -> the ratio of specific heats (can't be 1,makes a bunch of singularities in the equations)(defaults to 5/3 in M3D-C1 https://m3dc1.pppl.gov/NEWDOC-latest.pdf)
            #See equation 4, page 875, 876 in Glasser '75 -> comes straight from equation 16, Johnson and Greene, 1967.
            #γ is defined in Johnson and Greene '67 in equation 9, which is a simplification of equation 4 from Coppi 1966.
                #Coppi 1966's equations 1-7 are identical equilibrium to equations 6-12 in Johnson and Greene, save for two equations...
                    #equation 1, where 2nd order spatial variations in flow are killed, and equation 4, where η|J^2| and second order variations in temperature and flow are killed.
                    #Why on earth is η|J|^2 killed (η being resistivity)?? Idk. Regardless, equations work (scaffidi)
            #Resistive inner MHD equations: "The major limitations in this model are the neglect of finite gyration radius, viscosity and thermal conductivity terms and the use of scalar pressure and resistivity. - Greene and Johnson '67"
        #Vs, NEEDS SORTING OUT, SEE GLASSER '75 equ 81
        #ρ, η needed
    function generateInnerTerms_ScrewPinch_Glass75(Bt, Bp, Jt, Jp, p, dpdr, R0, rs,  ρ, η; Vs=1.0, γ=5/3)
        ∇V = ∇V_(R0)(rs)
        B = B_(Bt,Bp)(rs)
        sig = sig_(Bt,Bp,Jt,Jp)(rs)

        Ψprime = Ψprime_(Bt,R0)(rs)
        Ψdprime = Ψdprime_(Bt,R0)(rs)
        χprime = χprime_(Bp)(rs)
        χdprime= χdprime_(Bp,R0)(rs)

        Jprime = Jprime_(Jp)(rs)
        Iprime = Iprime_(Jt,R0)(rs)

        ι_on2Pi = ι_GlassOn2Pi(Bt,Bp,R0)(rs)
        Λ = Λ_(Bt,Bp,R0)(rs)
        α=α_(Bp,n,rs)

        #mu0 not needed from dpdr, needed for p
        capE = ((B^2/(∇V^2))/(Λ^2))*(Jprime*Ψdprime - Iprime*χdprime + Λ*(sig));
        capF = ((B^2/(∇V^2))/(Λ^2))*(dpdr(rs)*drdV(R0)(rs))^2*(1/B^2); #capF = ((B^2/(∇V^2))/(Λ^2))*(sig^2*B^2/(∇V^2) - (sig*B^2/(∇V^2))^2/(B^2/(∇V^2)) + (dpdr(rs)*drdV(R0)(rs))^2*(1/B^2));
        capH = 0.0; #capH = ((B^2/(∇V^2))/Λ)*(((sig*B^2/(∇V^2)))/(B^2/(∇V^2)) - sig);
        capM = 1.0; #capM = (B^2/(∇V^2))*((∇V^2/B^2) + (dpdr(rs)*drdV(R0)(rs))^(-2)*(sig^2*B^2 - (sig*B^2)^2/B^2)); 
        capK = 1/capF; #capK = (Λ^2/(capM*(dpdr(rs)*drdV(R0)(rs))^2))*B^2/(B^2/(∇V^2))
        capG = B^2/(capM*γ*mu0*p(rs))

        X0 = (ρ*capM*η^2*B^4/(α^2*Λ^2*(B^2/(∇V^2))^2))^(1/6)
        Q0 = (η*α^2*Λ^2*B^2/(ρ*capM*(B^2/(∇V^2))))^(1/3);

        return capE,capF,capH,capM,capK,capG,X0,Q0,(capE+capF)
    end

############################################################################################################################################################
#Scaffidi
############################################################################################################################################################
#η_diff is a diffusion coefficient (Scaffidi pg18))
    function generateInnerTerms_ScrewPinch_Scaffidi(Bp, q, p, dpdr, rs, k, η_diff; γ=5/3, ρ=1.0)
        qprime = ForwardDiff.derivative(q,rs)
        B = B_(Bt,Bp)(rs)

        Lr = (ρ*η_diff^2*rs^2/(n^2*Bp(rs)^2*qprime^2))^(1/6)
        
        capF = dpdr(rs)^2*rs^2*k^2/(Bp(rs)^4*qprime^2); #This the same
        capE = -2*dpdr(rs)*rs*k^2/(Bp(rs)^2*qprime^2) - capF; #This different
        capE_plus_capF = -2*dpdr(rs)*rs*k^2/(Bp(rs)^2*qprime^2)
        capH = 0.0; #capH = ((B^2/(∇V^2))/Λ)*(((sig*B^2/(∇V^2)))/(B^2/(∇V^2)) - sig);
        capM = 1.0; #capM = (B^2/(∇V^2))*((∇V^2/B^2) + (dpdr(rs)*drdV(R0)(rs))^(-2)*(sig^2*B^2 - (sig*B^2)^2/B^2)); 
        capK = 1/capF; #capK = (Λ^2/(capM*(dpdr(rs)*drdV(R0)(rs))^2))*B^2/(B^2/(∇V^2))
        capG = B^2/(capM*γ*mu0*p(rs))

        Qr = η_diff/(Lr^2) #See Scaffidi page 18 for more detail

        return capE,capF,capH,capM,capK,capG,Lr,Qr,capE_plus_capF
    end

############################################################################################################################################################
#Comparison/Discussion
############################################################################################################################################################

    #Different:
        #capE = ((B^2/(∇V^2))/(Λ^2))*(Jprime*Ψdprime - Iprime*χdprime + Λ*(sig));
        #capE = -2*dpdr(rs)*rs*k^2/(Bp(rs)^2*qprime^2) - capF;

        #X0 = (ρ*capM*η^2*B^4/(α^2*Λ^2*(B^2/(∇V^2))^2))^(1/6)
        #Lr = (ρ*η_diff^2*rs^2/(n^2*Bp(rs)^2*qprime^2))^(1/6)
            #^Makes sense these two different since their units are different

    #Same:
        #capF = dpdr(rs)^2*rs^2*k^2/(Bp(rs)^4*qprime^2) 
        #capF = ((B^2/(∇V^2))/(Λ^2))*(dpdr(rs)*drdV(R0)(rs))^2*(1/B^2)

        #Qr = η_diff/((ρ*η_diff^2*rs^2/(n^2*Bp(rs)^2*qprime^2))^(2/6))
        #Q0 = (η_diff*α^2*Λ^2*B^2/(ρ*capM*(B^2/(∇V^2))))^(1/3);

    #Discussion:
        #capF and capE should have the same units
            #Double check against Scaffidi/Equations directly
                #If this is the case, don't see why one pair should be different and the other the same...

#Two options
    #1: Just throw in the predictions into the growth rate... [I'm doing this for now]
        #Ask what density he used...
        #Furth w finite pressure...
    #2: FIND OUT UNITS OF the equations [Will do soon], convert expressions between one another
            #Convert the units of one expression to the other, see if they converge...
                #x is dimensionless
            #eta is a perturbed displacement. If its in units of x, it's dimensionless
            #First check Glasser repeat terms are behaving correctly
    #3: Figure out how to match the asymptotic forms when one's in Psi (prop to r^2) and the other's in r... I don't actually know...
        #Revert to Scaffidi I reckon. He knows whats up

############################################################################################################################################################
############################################################################################################################################################
#Testing:

if false
    q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
    rs0 = 2 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt) 
    R0 = 20 #Controls ratio of Bt to Bp
    ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
    xb = 1 #Widen device to increase q at edge (peaks current more stongly)
    Bp0 = 1.0 #Size of Bp (ratio of Bp to Bt unchanged)
    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=true, print_mathematica_inputs=true, plotrvec = range(0.000001,xb*rs0,200))

    m = 2;
    n = 1;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.001;

    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]

    η_diff = 1e-5
    η = 1e-5
    ρ = 1.0 #10^20*2.5*1.66054e-27

    generateInnerTerms_ScrewPinch_Scaffidi(Bp, q, p, dpdr, rs, k, η_diff; γ=5/3, ρ=1.0)
    generateInnerTerms_ScrewPinch_Glass75(Bt, Bp, Jt, Jp, p, dpdr, R0, rs,  ρ, η; Vs=1.0, γ=5/3)

    #capE's different (why this and none of the others?)
    #X0's/Lr different (makes sense)

    ############################################################################################################################################################


    ############################################################################################################################################################
end