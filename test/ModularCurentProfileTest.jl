#Plotting & Testing
############################################################################################################################################################
############################################################################################################################################################
    
#Jt = CubicSpline([0.0,0.5,1.0,1.5,2.0],[8*10^4,6.4*10^4,4*10^4,1.0*10^4,3.0*10^4])
#ppp = CubicSpline([0.0,0.5,1.0,1.5,2.0],[4.0*10^4,3.8*10^4,3.0*10^4,1.0*10^4,0.0*10^4])
#plot(plotrvec,1 ./ Bt.(plotrvec))

#plot(plotrvec, gradspln(p).(plotrvec))
#plot(plotrvec, gradspln(Jt).(plotrvec))
#plot!(plotrvec, forwarddiff(Jt).(plotrvec))
#plot(plotrvec, forwarddiff(Jt,2).(plotrvec))
#plot(plotrvec, gradspln(Jt).(plotrvec))
#plot(plotrvec, Bp.(plotrvec).*Jt.(plotrvec))

#WORKING ON CUBIC SPLINE TESTING:
r_raw = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
r_raw_double = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
r_inf = [0.0,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 
0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 
0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 
0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 
0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 
0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 
0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 
0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1., 
1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 
1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 
1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 
1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43, 1.44, 
1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 
1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 
1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 
1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 
1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0]



Jt_raw=[79577.5, 73573.8, 59139., 43024.2, 29587.1, 19894.4, 13366.3, 9082.53, 6279., 4426.48, 3183.1]
Jt_raw_double=[79577.5, 78009.5, 73573.8, 66978.8, 59139., 50929.6, 43024.2, 
35844.1, 29587.1, 24290.3, 19894.4, 16293.2, 13366.3, 10997.3, 
9082.53, 7533.96, 6279., 5258.85, 4426.48, 3744.45, 3183.1]
Jt_inf = [79577.5, 79561.6, 79513.8, 79434.4, 79323.4, 79181.1, 79007.6, 
78803.3, 78568.6, 78303.8, 78009.5, 77686.1, 77334.2, 76954.4, 
76547.4, 76113.8, 75654.4, 75169.9, 74661.1, 74128.8, 73573.8, 
72997.2, 72399.6, 71782.1, 71145.5, 70490.8, 69818.9, 69130.8, 
68427.5, 67709.8, 66978.8, 66235.3, 65480.5, 64715.1, 63940.1, 
63156.4, 62365., 61566.7, 60762.3, 59952.8, 59139., 58321.7, 57501.6, 
56679.6, 55856.3, 55032.6, 54209., 53386.3, 52565.1, 51746., 50929.6, 
50116.4, 49307.1, 48502., 47701.7, 46906.7, 46117.3, 45334., 44557.1, 
43787.1, 43024.2, 42268.7, 41520.9, 40781.2, 40049.6, 39326.5, 
38612.1, 37906.5, 37209.8, 36522.3, 35844.1, 35175.2, 34515.8, 
33865.9, 33225.6, 32594.9, 31974., 31362.7, 30761.1, 30169.3, 
29587.1, 29014.6, 28451.8, 27898.6, 27354.9, 26820.8, 26296.1, 
25780.8, 25274.8, 24778., 24290.3, 23811.7, 23342., 22881.2, 22429.1, 
21985.7, 21550.8, 21124.4, 20706.2, 20296.3, 19894.4, 19500.4, 
19114.3, 18736., 18365.2, 18001.9, 17645.9, 17297.2, 16955.6, 
16620.9, 16293.2, 15972.1, 15657.7, 15349.8, 15048.2, 14752.9, 
14463.8, 14180.7, 13903.4, 13632., 13366.3, 13106.1, 12851.4, 12602., 
12357.9, 12118.9, 11884.9, 11655.9, 11431.6, 11212.2, 10997.3, 
10786.9, 10581., 10379.5, 10182.2, 9989.01, 9799.92, 9614.81, 9433.6, 
9256.2, 9082.53, 8912.51, 8746.06, 8583.1, 8423.56, 8267.35, 8114.42, 
7964.69, 7818.08, 7674.52, 7533.96, 7396.32, 7261.55, 7129.57, 
7000.32, 6873.75, 6749.8, 6628.4, 6509.5, 6393.05, 6279., 6167.27, 
6057.84, 5950.64, 5845.63, 5742.76, 5641.98, 5543.24, 5446.5, 
5351.72, 5258.85, 5167.85, 5078.68, 4991.3, 4905.66, 4821.74, 
4739.48, 4658.87, 4579.85, 4502.4, 4426.48, 4352.06, 4279.1, 4207.57, 
4137.44, 4068.69, 4001.27, 3935.17, 3870.35, 3806.79, 3744.45, 
3683.32, 3623.36, 3564.55, 3506.87, 3450.29, 3394.79, 3340.34, 
3286.92, 3234.52, 3183.1]

Jt=CubicSpline(r_raw,Jt_raw)
Jt_double=CubicSpline(r_raw_double,Jt_raw_double)
Ji_inf=CubicSpline(r_inf,Jt_inf)

Jtbig = CubicSpline(rbig,Jtbig_vals)

Ji_inf_negative=CubicSpline(r_inf,Jt_inf.-11212.2)


_checkbounds(Jt, 0.0,2.0;debug=true)
_checkbounds(Jt_double, 0.0,2.0;debug=true)
_checkbounds(Ji_inf, 0.0,2.0;debug=true)

integral(Jt,0.0,2.0;debug=false)
integral(Jt_double,0.0,2.0;debug=false)
integral(Ji_inf,0.0,2.0;debug=false)
integral(Ji_inf_negative,0.0,2.0;debug=false)

integral(Jtbig,0.0,2.0;debug=false)

rJt_integral_components(Ji_inf)

rJt_components = rJt_integral_components(Jt)
rJt_components_inf = rJt_integral_components(Ji_inf)
rJt_compbig= rJt_integral_components(Jtbig)

rJt_integral(rJt_components,length(Jt.xs)-1)
rJt_integral(rJt_components_inf,length(Ji_inf.xs)-1)
rJt_integral(rJt_compbig,length(Jtbig.xs)-1)

Jt_on_r_integral = sum(Jt_on_r_integral_components(Jt)[1:end]) #from 0.2 to 1.0
Jt_on_r_integral_inf = sum(Jt_on_r_integral_components(Ji_inf)[20:end])
Jt_on_r_integral_big = sum(Jt_on_r_integral_components(Jtbig)[2000:end])

plot(rbig[2001:end],Jt.(rbig[2001:end])./rbig[2001:end],xlims=(0.0,2.0))
plot!(rbig[2001:end],Ji_inf.(rbig[2001:end])./rbig[2001:end],xlims=(0.0,2.0))
plot!(rbig[2001:end],Jtbig.(rbig[2001:end])./rbig[2001:end],xlims=(0.0,2.0))

y0test = [0.0,pi/8,0.42*0.5*pi,0.7*pi,pi]
xvec = range(0.0,y0test[end];step=y0test[end]/1000)



tspline=CubicSpline(y0test,sin.(y0test))

evpol(xref,ps) = r -> evalpoly(r-xref,ps)


plot(xvec,evpol(tspline.xs[1],tspline.ps[1,:]).(xvec), ylims = (0.0,2.0))
plot!(xvec,evpol(tspline.xs[2],tspline.ps[2,:]).(xvec),ylims = (0.0,2.0))
plot!(xvec,evpol(tspline.xs[3],tspline.ps[3,:]).(xvec),ylims = (0.0,2.0))
plot!(xvec,evpol(tspline.xs[4],tspline.ps[4,:]).(xvec),ylims = (0.0,2.0))
plot!(xvec,tspline.(xvec),color=:black)
scatter!(y0test,sin.(y0test))

#Jt_on_r_integral_components IS WORKING. Just remember, component 1 is from x2 to x3, etc

f1 = sum(full_integral_components(Jt))
f2 = sum(full_integral_components(Ji_inf))
f3 = sum(full_integral_components(Jtbig))

xve = range(0.0,0.2;step=0.2/1000)
#plot(xve,evpol(Jt.xs[1],Jt.ps[1,:]).(xve))

rjt_coeff = _int_rJt(Jt.ps[1,:],Jt.xs[1])
fullint_coeff = _int_full(Jt.ps[1,:],Jt.xs[1])

plot(xve,Jt.(xve))
plot(xve,evpol(0.0,rjt_coeff).(xve))
plot(xve,evpol(0.0,full_coeff).(xve))
plot!(xve,Jt.(xve[2:end]),color=:black)

plot(xve,evpol(0.0,fullint_coeff).(xve))

plot(xve,(evpol(0.0,fullint_coeff).(xve)).*(75.0453/evpol(0.0,fullint_coeff)(xve[end])))

Bp_spl(Jt) = r -> Bp_Spline(Jt, r)

plot(xve,Bp_spl(Jt).(xve))

sum(initialise_internalInt_Spline(Jt).*mu0)
sum(initialise_internalInt_Spline(Ji_inf).*mu0)
sum(initialise_internalInt_Spline(Jtbig).*mu0)

intintcomps,rJt_components=initialise_internalInt_Spline(Jt)
intintcomps_inf,rJt_components_inf = initialise_internalInt_Spline(Ji_inf)
intintcomps_big,rJt_components_big = initialise_internalInt_Spline(Jtbig)

internalInt_Spline(Jt, intintcomps, rJt_components, 2.0)
internalInt_Spline(Ji_inf, intintcomps_inf, rJt_components_inf, 2.0)

p(r) = 39788.7*(1.0 - ((r/1.0)/2.0)^2)

Bp1 = Bp_Spline(Jt; Bp_ref=0.0, r_ref=0.0)
Bp2 = Bp_Spline(Ji_inf; Bp_ref=0.0, r_ref=0.0)

plot(xvec,Bp2.(xvec))

Bt1,internalint1 = Bt_Spline(Jt,p,1.0)

plot(xvec,mu0.*internalint1.(xvec))


plot(xvec,Bt1.(xvec))


#Testing
Jt_validation = r -> (1/mu0)*(1/r)*ForwardDiff.derivative(rr -> rr*Bp1(rr), r)
Jp_validation = r -> -(1/mu0)*ForwardDiff.derivative(Bt1, r)
dpdr_validation = r -> ForwardDiff.derivative(p, r)


force_balance = r -> (Jp_validation(r)*Bt1(r)-Jt_validation(r)*Bp1(r)-dpdr_validation(r))/(p(0.0)^2)
force_balance_JB = r -> (Jp_validation(r)*Bt1(r)-Jt_validation(r)*Bp1(r))
force_balance_p = r -> (-dpdr_validation(r))

plot(xvec[1:end],force_balance_JB.(xvec[1:end]))
plot!(xvec[1:end],force_balance_p.(xvec[1:end]))
plot!(xvec[1:end],p.(xvec[1:end]))

plot(xvec[1:end],force_balance.(xvec[1:end]))

#We provide Jt, we provide p. Where does q come into it??


#using Conda
#using PyCall
#Conda.add("scipy")
#@pyimport scipy
####