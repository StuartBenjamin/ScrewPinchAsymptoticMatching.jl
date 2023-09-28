
############################################################################################################################################################
#Test, plotting solutions and increasing convergence of higher order expansions
############################################################################################################################################################

if true  #THIS IS WORKING
    β = 0.0005
    rs0 = 0.3
    R0 = 20 
    ν = 1.0  
    xb = 1.0 
    Bp0 = 1.0 
    q0 = 1.1

    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

    m = 2;
    n = 1;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.001;

    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]

    nmax = 5
    del=0.01
    Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, 20, del; integrator_reltol=10^(-20), plot_solution=true, plotwidth=10, plot_soln = true)


    plots = []
    del = 0.001
    del2= 0.001
    res = 1000
    plotwidth = 50
    closevec_in = range(rs-plotwidth*del2,rs-del2,length=res)
    closevec_out = range(rs+del2,rs+plotwidth*del2,length=res)

    imax = 10
    for i in 0:imax
        nmax = i*4

        ξs_plus, ξs_minus, σ_plus, σ_minus, f2, g0 = ξseries_plus_minus(f(k, m, Bt, Bp, k0_), g(k, m, dpdr, Bt, Bp, k0_), rs, nmax; ξ0 = 1.0)
        ξ_plus, ξ_minus, σ_plus, σ_minus = ξ_plus_minus(f(k, m, Bt, Bp, k0_), g(k, m, dpdr, Bt, Bp, k0_), rs, rs0, nmax; ξ0 = 1.0)
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, del; g=g_, integrator_reltol=1e-15, plot_solution=false)

        ξ_plus_prime(r) = ForwardDiff.derivative(ξ_plus,r)
        ξ_minus_prime(r) = ForwardDiff.derivative(ξ_minus,r)

        zl, zr, del2, outScale2 = zl_zr_calculator(solin, solout, inScale, outScale, del, k, m, Bt, Bp, rs; del2 = del2)

        Δl = (zl*ξ_minus_prime(rs-del2)-ξ_minus(rs-del2))/(ξ_plus(rs-del2)-zl*ξ_plus_prime(rs-del2))
        Δr = (zr*ξ_minus_prime(rs+del2)-ξ_minus(rs+del2))/(ξ_plus(rs+del2)-zr*ξ_plus_prime(rs+del2))

        ############################################################################################################################################################
        ξl(r) = ξ_minus(r)+Δl*ξ_plus(r)
        ξr(r) = ξ_minus(r)+Δr*ξ_plus(r)

        psi_l_raw(r) = ξl(r)*F(k, m, Bt, Bp)(r)
        psi_r_raw(r) = ξr(r)*F(k, m, Bt, Bp)(r)

        cl = solin(rs-del2)[2]*inScale/psi_l_raw(rs-del2)
        cr = solout(rs+del2)[2]*outScale2/psi_r_raw(rs+del2)

        psi_l(r) = cl*psi_l_raw(r)
        psi_r(r) = cr*psi_r_raw(r)

        closevec_in = range(rs-plotwidth*del2,rs-del2,length=res)
        closevec_out = range(rs+del2,rs+plotwidth*del2,length=res)

        if i == 0
            #p1 = plot([closevec_in  closevec_out],[[solin(i)[2]*inScale for i in closevec_in]  [solout(i)[2]*outScale2 for i in closevec_out]], label = "Psi")
            p1 = plot([closevec_in  closevec_out],[[psi_l(i) for i in closevec_in]  [psi_r(i) for i in closevec_out]], xlims = (rs-plotwidth*del2,rs+plotwidth*del2), linewidth=2, label=false)#["Psi expansion O($(nmax))"   "two"], linewidth=10) #label=false
        end
        p1 = plot!([closevec_in  closevec_out],[[psi_l(i) for i in closevec_in]  [psi_r(i) for i in closevec_out]], xlims = (rs-plotwidth*del2,rs+plotwidth*del2), linewidth=2, label=false)#"Psi expansion O($(nmax))", linewidth=10) #label=false

        if i == imax
            plot!([closevec_in],[solin(i)[2]*inScale for i in closevec_in], label = false, color=:black)
            plot!([closevec_out],[solout(i)[2]*outScale2 for i in closevec_out], label = "Psi numerical", color=:black)
            vline!([rb], label = "Edge of device")
            display(p1)
        end
        #plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=5)
        #p0 = plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=5)
        #push!(plots,plot_full_Psis_w_frobenius(solin, solout, inScale, outScale2, psi_l, psi_r, del2, nmax, rs; plotwidth=5))
    end
############################################################################################################################################################
#Checking del-dependence: THIS IS HOW TO CHECK YOU'RE IN XMIN
############################################################################################################################################################
    β = 0.1 #Core beta
    ν = 1   #Current shape modulation parameter
    rs0 = 1.3   #radial normalisation unit, setting of singular layer (Chandra)
    xb = 2.0    #wall location in normalised units
    rb = xb*rs0 #wall physical location
    R0 = 20

    Bp = Bp_Chand(rs0,R0,ν)
    Bt = Bt_Chand(β,rs0,R0,ν,xb)
    dpdr = dpdr_Chand(β,rs0,R0,ν,xb)
    #q_(Bt,Bp,R0)

    m = 2;
    n = 1;
    k_(n,R0) = n/R0;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.000001;

    rs, rs_plot = find_rs(q_(Bt,Bp,R0),m,n,rb,q0,ν,rs0)
    rs = rs[1]
    nmax = 5

    delvec = 10 .^ range(-7,-0.5, length = 10)

    Δl_Δr_calc1(nmax) = del -> Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false)

    Δvecs = []
    delvecs = []
    pdel = nothing
    for j = 0:10
        Δvec = []
        for i in 1:length(delvec)
            Δl,Δr,del2 = Δl_Δr_calc1(j)(delvec[i])
            if del2 != delvec[i]
                @error("cat in the hat")
            end
            push!(Δvec,(Δl+Δr))
        end
        j==0 ? (pdel = plot(delvec, Δvec, xaxis=:log, legend=:bottomleft, label = "O($(j))", title = "Δl + Δr as a function of matching width", xlabel="del", ylabel = "Δl + Δr")) : (pdel = plot!(delvec, Δvec, xaxis=:log, legend=:bottomleft, label = "O($(j))", title = "Δl + Δr as a function of matching width", xlabel="del", ylabel = "Δl + Δr"))

        push!(Δvecs,Δvec)
        push!(delvecs,delvec)
    end
    display(pdel)
    plot(delvecs,Δvecs, xaxis=:log, legend=:bottomleft)

############################################################################################################################################################
#Finding xmin
############################################################################################################################################################
    
    q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
    rs0 = 2 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt) 
    R0 = 20 #Controls ratio of Bt to Bp
    ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
    xb = 1 #Widen device to increase q at edge (peaks current more stongly)
    Bp0 = 1.0 #Size of Bp (ratio of Bp to Bt unchanged)

    m = 2;
    n = 1;
    k = k_(n,R0)
    c0 = 1;
    r0 = 0.000001;
    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=true, print_mathematica_inputs=true, plotrvec = range(0.000001,xb*rs0,200))
    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]

    delvec = 10 .^ range(-7.2,-0.5, length = 10)
    nmax=0
    #xmin, Δl_Δr_del2_tuple, abs_err = Δl_Δr_xmin_deprecated(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec;  xmin_tol = 1e-3, integrator_reltol=10^(-20), plot_solution=true, plotwidth=10, plot_soln = false)
    xmin_new, Δl_Δr_del2_tuple_new, abs_err = Δl_Δr_xmin(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, delvec;  xmin_tol = 1e-3, integrator_reltol=10^(-20), plot_solution=true, plotwidth=10, plot_soln = true)

    #xmin_new, Δl_Δr_del2_tuple_new = Δl_Δr_xmin_visualiser(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, delvec; plot_solution=true, finite_pressure=true)

    big_n=20
    xmin_tol = 1e-5
    integrator_reltol=10^(-20)
    plot_solution=true
    del2 = 0.0
    plotwidth=10
    plot_soln = false
    verbose = false

############################################################################################################################################################
#Chandra Beta pressure-dependence:
#Do the full match (checking correct defs of Δl, Δr) then come back to this
############################################################################################################################################################
    ν = 1   #Current shape modulation parameter
    rs0 = 1.0   #radial normalisation unit, setting of singular layer (Chandra)
    xb = 2.0    #wall location in normalised units
    rb = xb*rs0 #wall physical location
    R0 = 20

    βvec = []
    delta_vec = []
    for i in 1:76
        print(i,"\n")
        β = i*0.01 #Core beta
        push!(βvec, β)

        Bp = Bp_Chand(rs0,R0,ν)
        Bt = Bt_Chand(β,rs0,R0,ν,xb)
        dpdr = dpdr_Chand(β,rs0,R0,ν,xb)
        #q_(Bt,Bp,R0)

        m = 2;
        n = 1;
        k_(n,R0) = n/R0;
        k = k_(n,R0)
        qtest = m/n;
        c0 = 1;
        r0 = 0.000001;

        rs, rs_plot = find_rs(q_(Bt,Bp,R0),m,n,rb,q0,ν,rs0)
        rs = rs[1]

        nmax = 20
        del = 10^-7

        Δl,Δr,del2 = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false)

        push!(delta_vec, Δl+Δr)
    end

    plot(βvec,delta_vec)

############################################################################################################################################################
#Nishimura Scaffidi Beta pressure-dependence -> works!!
#Do the full match (checking correct defs of Δl, Δr) then come back to this
############################################################################################################################################################
    ν = 1   #Current shape modulation parameter
    rs0 = 1.0   #radial normalisation unit, setting of singular layer (Chandra)
    xb = 2.0    #wall location in normalised units
    rb = xb*rs0 #wall physical location
    R0 = 20

    βvec = []
    delta_vec = []
    for i in 1:49
        print(i,"\n")
        β = i*0.01 #Core beta
        push!(βvec, β)

        Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Scaffidi_Equil(β,rs0,R0,ν,xb; Bp0=0.0, q0=1.0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

        m = 2;
        n = 1;
        k_(n,R0) = n/R0;
        k = k_(n,R0)
        qtest = m/n;
        c0 = 1;
        r0 = 0.000001;

        rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
        rs = rs[1]

        nmax = 20
        del = 10^-7

        Δl,Δr,del2 = Δl_Δr_calculator(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, nmax, del; integrator_reltol=10^(-20), plot_solution=false, plotwidth=10, plot_soln = false)

        push!(delta_vec, Δl+Δr)
    end

    plot(βvec,delta_vec)
end

############################################################################################################################################################
#Zero-pressure expansion testing
############################################################################################################################################################
if true
    ############################################################################################################################################################
    #Basic Test -> WORKING
    ############################################################################################################################################################
    q0 = 1.3 #Moves qstart, q still increases by the same proportion regardless of where you place it
    rs0 = 2.0 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt) 
    R0 = 20 #Controls ratio of Bt to Bp
    ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
    xb = 1.0 #Widen device to increase q at edge (peaks current more stongly)
    Bp0 = 1.0 #Size of Bp (ratio of Bp to Bt unchanged)
    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=true, print_mathematica_inputs=true, plotrvec = range(0.000001,xb*rs0,200))
    m = 2;
    n = 1;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.00001;

    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]

    κ= κ_(Bp, Bt, dpdr, k, m, rs; g1=gzero,test=true)

    nmax=2
    del = 0.005

        Δl_Δr_calculator_zeroPressure(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, 2*del; test_κ=true, integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=3, plot_soln = true, plot_soln_1=false)

        ψ_plus, ψ_minus, ψ_plus_prime, ψ_minus_prime = ψ_plus_minus_zeroPressure(Bp, Bt, dpdr, k, m, rs; test_κ=true)
        solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, 2*del; g=gzero, integrator_reltol=10^(-20), plot_solution=true)

        #el2 = 2*del
        Δl_Δr_calculator_zeroPressure(rs, ψ_plus, ψ_minus, ψ_plus_prime, ψ_minus_prime, solin, solout, inScale, outScale, 2*del, 2*del; test_κ=true, plot_solution=true, plotwidth=3, plot_soln = true, plot_soln_1=false)


        delvec = 10 .^ range(-9.0,-0.5, length = 500)
        xmin_new, Δl_Δr_del2_tuple_new = Δl_Δr_xmin_visualiser(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, delvec; plot_solution=true, finite_pressure=false)

        Δl_Δr_calculator_zeroPressure(Bp, Bt, dpdr, k, m, r0, rs, rb, rs0, 0.00048907870941395; test_κ=true, integrator_reltol=10^(-20), plot_solution=true, del2 = 0.0, plotwidth=10, plot_soln = true, plot_soln_1=false)

    ############################################################################################################################################################
    #Beta_match
    ############################################################################################################################################################





end
######

