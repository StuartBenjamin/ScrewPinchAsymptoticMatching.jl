#Testing!!! 
##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################
#FURTH
    #Bt doesn't vary w r, Bp is zero.
    #There's finite pressure to balance Bp. 
    #If the ratio of Bp to Bt is small, beta is small (controlled by rs0).
##########################################################################################################################################################
#Furth is working perfectly. This replicated Cesar's Modular Pinch Furth mathematica script
q0 = 1.1 #Moves qstart, q still increases by the same proportion regardless of where you place it
rs0 = 2 #Make smaller to increase current in core while leaving q-profile unchanged (stronger Bt)
R0 = 20 #Controls ratio of Bt to Bp
ν = 1  #Increase ν to increase total current by widening peak. Leaves peak magnitude unchanged
xb = 1 #Widen device to increase q at edge (peaks current more stongly)
Bp0 = 1.0 #Size of Bp (ratio of Bp to Bt unchanged)
    #Fix total current, change shape (figure out how to do that)
    #Keep the same, change current magnitude
Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Furth_Equil(q0,rs0,R0,ν,xb; Bp0=Bp0, plot_equil=true, print_mathematica_inputs=true, plotrvec = range(0.000001,xb*rs0,200))

r2 = 1.23456
Bp(r2)
Bt(r2)
q(r2)
dpdr(r2)
p(r2)
Jt(r2)
Jp(r2)
rb

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

r2 = 1.2345
Bp(r2)
Bt(r2)
q(r2)
dpdr(r2)
p(r2)
Jt(r2)
Jp(r2)
rb

##########################################################################################################################################################
#Scaffidi
    #Same as Chandra except Bt at the boundary is set, and so the variation occurs towards the axis.
    #Pressure is perturbed like Bt, deviates towards the axis.
    #Bp not zero.
##########################################################################################################################################################
β = 0.0000000000001 
rs0 = 2.0 
R0 = 20 
ν = 1.0  
xb = 1.0 
Bp0 = 1.0# #Size of Bp (ratio of Bp to Bt unchanged)
q0 = 1.1

Bp,Bt,q,dpdr,p,Jt,Jp,rb,outerp6 = Scaffidi_Equil(β,rs0,R0,ν,xb; Bp0=Bp0, q0=q0, plot_equil=true, print_mathematica_inputs=false, plotrvec = range(0.000001,xb*rs0,200))

r2 = 0.12345
Bp(r2)
Bt(r2)
q(r2)
dpdr(r2)
p(r2)
Jt(r2)
Jp(r2)
rb

##########################################################################################################################################################
##########################################################################################################################################################

passed = true