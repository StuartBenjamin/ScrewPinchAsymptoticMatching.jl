##########################################################################################################################################################
#Tests:
if true #Checking functions are working properly
    β = 0.0000001 
    rs0 = 1.0 
    R0 = 20 
    ν = 1.0  
    xb = 2.0 
    Bp0 = 0.0 
    q0 = 1.0

    Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Chandra_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

    m = 2;
    n = 1;
    k = k_(n,R0)
    qtest = m/n;
    c0 = 1;
    r0 = 0.001;

    rs, rs_plot = find_rs(q,m,n,rb,q0,ν,rs0)
    rs = rs[1]

    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, 0.0001; integrator_reltol=10^(-20), plot_solution=true)
    solin, solout, inScale, outScale, del = Psi_w_scales(Bp, Bt, dpdr, k, m, r0, rs, rb, 0.0001; g=gzero, integrator_reltol=10^(-20), plot_solution=true)
    delta_prime_raw, solin, solout, inScale, outScale, del =raw_delta_prime(Bp, Bt, dpdr, n/R0, m, r0, rs, rb, 0.01; del2=0.0, integrator_reltol=10^(-20));
    delta_prime_raw, solin, solout, inScale, outScale, del =raw_delta_prime(Bp, Bt, dpdr, n/R0, m, r0, rs, rb, 0.01; g=gzero, del2=0.0, integrator_reltol=10^(-20));
    raw_delta_prime(solin, solout, inScale, outScale, del; del2=0.02)
end
if true #Sorting out different ways to write g
    β = 0.1 
    rs0 = 2.0 
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
    rvec = range(0.000001,rb,200)

    gCop(k, n, m, dpdr, Bt, Bp) = r -> (2*n^2*k^2*r^2/(n^2*k^2*r^2+m^2))*dpdr(r) + (1/r)*(n*k*r*Bt(r)-m*Bp(r))^2*(n^2*k^2*r^2+m^2-1)/(n^2*k^2*r^2+m^2)+(2*n^2*k^2*r)/((n^2*k^2*r^2+m^2)^2)*(n^2*k^2*r^2*Bt(r)^2-m^2*Bp(r)^2)  #From Coppi '66
    g_old(k, m, dpdr, Bt, Bp) = r -> ((m^2 - 1)*r*F(k, m, Bt, Bp)(r)^2 + k^2*r^2*(2*dpdr(r) + r*F(k, m, Bt, Bp)(r)^2 + F(k, m, Bt, Bp)(r)*2*(k*r*Bt(r) - m*Bp(r))/(k^2*r^2 + m^2)))/(k^2*r^2 + m^2)  #From Furth '73

    #gCop is identical to g_ from Ideal MHD. g_old (from Furth) is very similar (small difference on order of <1%).

    p2b = plot([rvec,rvec,rvec],[g_(k, m, dpdr, Bt, Bp).(rvec),g_old(k, m, dpdr, Bt, Bp).(rvec),gCop(k, n, m, dpdr, Bt, Bp).(rvec)];label=["g_"   "g_old" "gCop"])
    plot(p1b,p2b)
    p3b = plot([rvec,rvec],[(g_(k, m, dpdr, Bt, Bp).(rvec)-gCop(k, n, m, dpdr, Bt, Bp).(rvec))/g_(k, m, dpdr, Bt, Bp)(2.0),(g_(k, m, dpdr, Bt, Bp).(rvec)-g_old(k, m, dpdr, Bt, Bp).(rvec))/g_(k, m, dpdr, Bt, Bp)(2.0)];label=["g_ - gCop (normalised)"   "g_ - g_old (normalised)"])

    ##########################################################################################################################################################
    #This below (k0_big,g_big) is wrong, from an old error I made
    k0_big(k,m) = r -> k^2+m^2/r^2 #wrong, an old error I made
    g_big(k, m, dpdr, Bt, Bp) = r -> g(k, m, dpdr, Bt, Bp, k0_big)(r) #wrong, see above
    p1b = plot([rvec,rvec,rvec,rvec],[g_(k, m, dpdr, Bt, Bp).(rvec),g_old(k, m, dpdr, Bt, Bp).(rvec),gCop(k, n, m, dpdr, Bt, Bp).(rvec),g_big(k, m, dpdr, Bt, Bp).(rvec)];label=["g_"   "g_old" "gCop"  "g_big"])
    ##########################################################################################################################################################

    plot(p1b,p2b,p3b)
end