#########################################################################################
#General Spline Control:
#########################################################################################

"
Fails if sign of gradient changes more than 'max' times.
"
function num_inflection_points(f,rb;max=2)
    fderiv(r) = ForwardDiff.derivative(f,r)
    fdoublederiv(r) = ForwardDiff.derivative(fderiv,r)

    if length(find_zeros(fdoublederiv,0.0,rb))>max
        return false
    end    
    return true
end

"
Fails if sign of gradient changes more than 'max' times. Coarse grained to spline inputs.
"
function num_inflection_points_coarse(spline::CubicSpline;max=2)
    x_intervals = [spline.xs[i+1]-spline.xs[i] for i in 1:(length(spline.xs)-1)]
    y_intervals = [spline.ys[i+1]-spline.ys[i] for i in 1:(length(spline.xs)-1)]
    grads = y_intervals./x_intervals

    count=0
    for i in 1:(length(grads)-1)
        if sign(grads[i])!=sign(grads[i+1])
            count+=1
        end
    end

    if count>max
        return false
    end  
    return true
end

"
Fails if gradmax is equalled or exceeded, save at 'n_exceeds' locations.
"
function gradmax_exceeded(f,rb,gradmax; n_exceeds=0, abs_val=false)
    fderiv(r) = ForwardDiff.derivative(f,r)

    if !abs_val
        zeros = find_zeros(r -> fderiv(r)-gradmax,0.0,rb)
    else
        zeros = find_zeros(r -> abs(fderiv(r))-gradmax,0.0,rb)
    end

    if length(zeros)>n_exceeds || (n_exceeds==0 && (fderiv(rb/2)>gradmax))
        return false
    else
        return true
    end
end

"
Fails if gradmin is equalled or exceeded (lower), save at 'n_exceeds' locations.
"
function gradmin_exceeded(f,rb,gradmin; n_exceeds=0, abs_val=false)
    fderiv(r) = ForwardDiff.derivative(f,r)

    if !abs_val
        zeros = find_zeros(r -> fderiv(r)-gradmin,0.0,rb)
    else
        zeros = find_zeros(r -> abs(fderiv(r))-gradmin,0.0,rb)
    end

    if length(zeros)>n_exceeds || (n_exceeds==0 && (fderiv(rb/2)<gradmin))
        return false
    else
        return true
    end
end

"
Fails if gradmax is equalled or exceeded. Coarse grained to spline inputs.
"
function gradmax_exceeded_coarse(spline::CubicSpline,gradmax; abs_val=false)
    x_intervals = [spline.xs[i+1]-spline.xs[i] for i in 1:(length(spline.xs)-1)]
    y_intervals = [spline.ys[i+1]-spline.ys[i] for i in 1:(length(spline.xs)-1)]
    grads = y_intervals./x_intervals

    if !abs_val
        maximum(grads) > gradmax && (return false)
    else
        maximum(abs.(grads)) > gradmax && (return false)
    end

    return true
end

"
Fails if max is equalled or exceeded, save at 'n_exceeds' locations.
"
function max_exceeded(f,rb,fmax; n_exceeds=0)
    zeros = find_zeros(r -> f(r)-fmax,0.0,rb)

    if length(zeros)>n_exceeds || (n_exceeds==0 && (f(rb/2)>fmax))
        return false
    else
        return true
    end
end

"
Fails if max is equalled or exceeded. Coarse grained to spline inputs.
"
function max_exceeded_coarse(spline::CubicSpline,fmax)
    if maximum(spline.ys)>fmax
        return false
    else
        return true
    end
end

"
Fails if min is equalled or gone beneath, save at 'n_exceeds' locations.
"
function min_exceeded(f,rb,fmin; n_exceeds=0)
    zeros = find_zeros(r -> f(r)-fmin,0.0,rb)

    if length(zeros)>n_exceeds || (n_exceeds==0 && (f(rb/2)<fmin))
        return false
    else
        return true
    end
end

"
Fails if min is equalled or gone beneath, save at 'n_exceeds' locations. Coarse grained to spline inputs.
"
function min_exceeded_coarse(spline::CubicSpline,fmin)
    if minimum(spline.ys)<fmin
        return false
    else
        return true
    end
end

"
Fails if sign of gradient is positive. Coarse grained to spline inputs.
"
function monotonic_coarse(spline::CubicSpline)
    x_intervals = [spline.xs[i+1]-spline.xs[i] for i in 1:(length(spline.xs)-1)]
    y_intervals = [spline.ys[i+1]-spline.ys[i] for i in 1:(length(spline.xs)-1)]
    grads = y_intervals./x_intervals

    for i in grads
        sign(i)==1.0 && (return false)
    end

    return true
end

function spline_narrow_max_(max_ref, rb, maxgrad_width; area_integrated_max=nothing)
    if !(area_integrated_max isa Nothing)
        max_ref=area_integrated_max/(pi*rb^2)
    elseif max_ref isa Nothing
        throw("Need either max_ref or area_integrated_max defined.")
    end

    maxgrad = max_ref/(rb*maxgrad_width)

    return maxgrad
end

#########################################################################################
#Specific Spline Control:
#########################################################################################

function total_plasma_current(Jt::CubicSpline)
    return total_plasma_current(Jt,Jt.xs[end])
end

function total_plasma_current(Jt::Union{Function,CubicSpline},rb; use_num_integrate=false)
    if (Jt isa CubicSpline) && !use_num_integrate
        Jtot = 2*pi*rJt_integral(Jt, rb, 0.0)
    else
        Jtot = 2*pi*quadgk(r->r*Jt(r), 0.0, rb, rtol=Jt(0.000000001)*1e-15)[1]
    end

    return Jtot
end

#########################################################################################
#Generating Cleaners:
#########################################################################################

function generate_cleaners(rb; max_ref = nothing, area_integrated_max=nothing, 
        Jtot_range=nothing, use_coarse=true, use_fine=true, maxval = nothing, minval=nothing,
        wobbles = 3, monotonic=false, maxgrad_width = rb/6)

    cleaner_functions = []
    
    monotonic_(f) = gradmax_exceeded(f,rb,0.0)
    wobbles_(f) = num_inflection_points(f,rb;max=wobbles)
    wobbles_coarse(f) = num_inflection_points_coarse(f;max=wobbles)

    if !(Jtot_range isa Nothing)
        function Jt_tot_(f)
            if (Jtot_range[1]<total_plasma_current(f,rb)<Jtot_range[2]) || (Jtot_range[2]<total_plasma_current(f,rb)<Jtot_range[1])
                return true
            end 
            return false
        end
    end
    
    if use_coarse
        #1
        if !(maxval isa Nothing)
            push!(cleaner_functions,f -> max_exceeded_coarse(f,maxval))
        end

        if !(minval isa Nothing)
            push!(cleaner_functions,f -> min_exceeded_coarse(f,minval))
        end

        #2
        monotonic && (push!(cleaner_functions,monotonic_coarse)) 
        push!(cleaner_functions,f -> gradmax_exceeded_coarse(f,spline_narrow_max_(max_ref, rb, maxgrad_width;area_integrated_max=area_integrated_max); abs_val=true))
                                    
        #3
        wobbles>=0 && (push!(cleaner_functions,wobbles_coarse))
    end

    #0 Integration
    !(Jtot_range isa Nothing) && push!(cleaner_functions,Jt_tot_)

    if use_fine
        #1 Function eval
        if !(maxval isa Nothing)
            push!(cleaner_functions,f -> max_exceeded(f,rb,maxval))
        end

        if !(minval isa Nothing)
            push!(cleaner_functions,f -> min_exceeded(f,rb,minval))
        end

        #2 Derivative eval
        monotonic && (push!(cleaner_functions,monotonic_)) 
        push!(cleaner_functions,f -> gradmax_exceeded(f,rb,spline_narrow_max_(max_ref, rb, maxgrad_width;area_integrated_max=area_integrated_max);abs_val=true))

        #3 Double-derivative eval
        wobbles>=0 && (push!(cleaner_functions,wobbles_))
    end

    return cleaner_functions
end

#########################################################################################
#Random Spline Generation:
#########################################################################################

function random_in_range(valmax,valmin,numvals)
    return valmin.+(valmax-valmin).*rand(Float64, numvals)
end

function random_Int_in_range(intmax,intmin,numvals)
    return rand(intmin:intmax,numvals)
end

function randomJt(Jtotmax, numJts, rb; knotmax=7, knotmin=5, minJt=0.0, 
    plot_profs=false, J0bounds=nothing, Jedgebounds=nothing)

    knotnums = random_Int_in_range(knotmax,knotmin,numJts)

    splines_vals = random_in_range(Jtotmax/(2*pi),minJt,(numJts,knotmax))

    if J0bounds isa Nothing
        start_vals = splines_vals[:,1]
    elseif length(J0bounds)==1
        start_vals=J0bounds[1].*ones(numJts)
    elseif J0bounds isa Real
        start_vals=J0bounds.*ones(numJts)
    else
        start_vals = random_in_range(maximum(J0bounds),minimum(J0bounds),numJts)
    end

    if Jedgebounds isa Nothing
        end_vals = splines_vals[:,end]
    elseif length(Jedgebounds)==1
        end_vals=Jedgebounds[1].*ones(numJts)
    elseif Jedgebounds isa Real
        end_vals=Jedgebounds.*ones(numJts)
    else
        end_vals = random_in_range(maximum(Jedgebounds),minimum(Jedgebounds),numJts)
    end

    Jt_vec = CubicSpline[]

    for i in 1:numJts
        Jt = CubicSpline(range(0.0,rb;length=knotnums[i]),vcat(start_vals[i],splines_vals[i,2:(knotnums[i]-1)],end_vals[i]))
        push!(Jt_vec, Jt)
    end

    plot_profs && plot_profiles(Jt_vec, rb)
    
    return Jt_vec
end

function random_pressure(p0bounds, num_splines, rb; knotmax=12, knotmin=5, maxp=nothing, maxp_mult=1.0, minp=0.0, plot_profs=false,  pedgebounds=0.0)
    knotnums = random_Int_in_range(knotmax,knotmin,num_splines)

    if length(p0bounds)==1
        start_vals=p0bounds[1].*ones(num_splines)
        maxp0=p0bounds[1]
    elseif p0bounds isa Real
        start_vals=p0bounds.*ones(num_splines)
        maxp0=p0bounds
    else
        start_vals = random_in_range(maximum(p0bounds),minimum(p0bounds),num_splines)
        maxp0=maximum(p0bounds)
    end

    if maxp isa Nothing
        maxp = maxp_mult*maxp0
    end

    splines_vals = random_in_range(maxp,minp,(num_splines,knotmax))

    if pedgebounds isa Nothing
        end_vals = splines_vals[:,end]
    elseif length(pedgebounds)==1
        end_vals=pedgebounds[1].*ones(num_splines)
    elseif pedgebounds isa Real
        end_vals=pedgebounds.*ones(num_splines)
    else
        end_vals = random_in_range(maximum(pedgebounds),minimum(pedgebounds),num_splines)
    end

    p_vec = CubicSpline[]

    for i in 1:num_splines
        pspl = CubicSpline(range(0.0,rb;length=knotnums[i]),vcat(start_vals[i],splines_vals[i,2:(knotnums[i]-1)],end_vals[i]))
        push!(p_vec, pspl)
    end

    plot_profs && plot_profiles(p_vec, rb)
    
    return p_vec
end

#########################################################################################
#Clean Spline Generation and plotting:
#########################################################################################

function clean_vec(in_vec,cleaner_functions)
    valid = trues(length(in_vec))
    valid_inds = 1:length(in_vec)

    for cleaner_func in cleaner_functions
        valid_inds = filter(x->cleaner_func(in_vec[x]),valid_inds) 
    end

    return in_vec[valid_inds]
end

function gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean,cleaner_functions; maxbatches=100, verbose=true, kwargs...)
    clean_counter = 0
    Jt_vec = CubicSpline[]
    i=1

    while clean_counter < num_clean && (i<(maxbatches+1))
        verbose && print("Batch $(i) of $(batch_size) random splines being generated. Valid count = $(clean_counter).\n")
        Jt_temps = clean_vec(randomJt(Jtotmax, batch_size, rb; kwargs...),cleaner_functions)

        clean_counter += length(Jt_temps)

        append!(Jt_vec,  Jt_temps)
        i+=1
    end
    print("$(clean_counter) of $(i*batch_size) random splines meet requirements.\n")

    return Jt_vec
end

function gen_n_clean_pressure_profiles(p0bounds,rb,batch_size,num_clean,cleaner_functions; maxbatches=100, verbose=true,  kwargs...)
    clean_counter = 0
    pressure_vec = CubicSpline[]
    i=1

    while clean_counter < num_clean && (i<(maxbatches+1))
        verbose && print("Batch $(i) of $(batch_size) random splines being generated. Valid count = $(clean_counter).\n")
        p_temps = clean_vec(random_pressure(p0bounds, batch_size, rb; kwargs...),cleaner_functions)

        clean_counter += length(p_temps)

        append!(pressure_vec,  p_temps)
        i+=1
    end
    print("$(clean_counter) of $(i*batch_size) random splines meet requirements.\n")

    return pressure_vec
end

function plot_profiles(Jt_vec, rb; p_vec=nothing, plotrvec = range(0.000001,rb,200), ylims=nothing)
    plotvec=[]
          
    for i in 1:length(Jt_vec)
        Jt=Jt_vec[i]
        p5=plot(plotrvec,Jt.(plotrvec),title = "Toroidal current density",xlabel="r (m)",ylabel="Amps/m^2",label=false,ylims=ylims)
        p5=scatter!(Jt.xs,Jt.(Jt.xs),title = "Toroidal current",xlabel="r (m)",ylabel="Amps/m^2",label=false,ylims=ylims)
        push!(plotvec,p5)
    end

    if length(Jt_vec)==1
        display(plot(plotvec[1]))
    elseif length(Jt_vec)==2
        display(plot(plotvec[1],plotvec[2]))
    elseif length(Jt_vec)==3
        display(plot(plotvec[1],plotvec[2],plotvec[3]))
    elseif length(Jt_vec)==4
        display(plot(plotvec[1],plotvec[2],plotvec[3],plotvec[4]))
    elseif length(Jt_vec)>=5
        display(plot(plotvec[1],plotvec[2],plotvec[3],plotvec[4],plotvec[5]))
    end
end

function gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean; maxbatches=100, verbose=true, 
    max_ref = nothing, area_integrated_max=nothing, Jtot_range=nothing, use_coarse=true, 
    use_fine=true, maxval = nothing, minval=nothing, wobbles = 3, monotonic=false, maxgrad_width = rb/6, kwargs...)

    return gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean,generate_cleaners(rb; max_ref = max_ref, 
        area_integrated_max=area_integrated_max, Jtot_range=Jtot_range, use_coarse=use_coarse, use_fine=use_fine, maxval = maxval, minval=minval,
        wobbles = wobbles, monotonic=monotonic, maxgrad_width = maxgrad_width); maxbatches=maxbatches,verbose=verbose, kwargs...)
end

function axis_beta(p,Bt0)
    pm0=(1/(2*mu0))*(Bt0^2)
    return p(0.0)/pm0
end

#########################################################################################
#Equilibrium generation with input Jts, pressure_profs
#########################################################################################

struct Equilibrium
    Bp::Function
    Bt::Function
    q::Function
    dpdr::Function
    p::Union{Function,CubicSpline}
    Jt::Union{Function,CubicSpline}
    Jp::Union{Function}
    rb::Number
    R0::Number
    rs::Union{Number,Nothing}
    rs0::Number #Deprecated re-scaling of r, set to 1.0 (defined for Chandra equilibria but not splines).
end

struct ΔprimeScrew
    rs::Number
    m::Union{Number,Nothing}
    n::Union{Number,Nothing}
    Δl::Union{Number,Nothing}
    Δr::Union{Number,Nothing}
    Δprime::Union{Number,Nothing}
    del::Union{Number,Nothing}
    nmax::Union{Number,Nothing}
    Δlzero::Union{Number,Nothing}
    Δrzero::Union{Number,Nothing}
    Δprimezero::Union{Number,Nothing}
    delzero::Union{Number,Nothing}
end

struct ResistiveEquilibrium
    equilibrium::Equilibrium
    Δprimes::Union{ΔprimeScrew,AbstractArray{ΔprimeScrew}}
end

function gen_equilibria_Jts(Jts,pressure_profs; Bt0=10, R0=3, dpdr_vec=nothing, dpdr=nothing, rs0=1.0, qtest=nothing)
    function dpdr0(i)
        return nothing
    end

    # Function currently set up to either take one pressure profile (w. optional dpdr), or a vector of pressure profiles (w. optional equal length vector of dpdrs)
    # Not yet able to take non-zero, non-equal length vectors of pressure and Jt, for example to cycle through 3 characteristic pressure proiles
    if (dpdr isa Nothing) && (dpdr_vec isa Nothing)
        function dpdr0(i)
            return nothing
        end
    elseif !(dpdr_vec isa Nothing) && (dpdr isa Nothing)
        @assert length(dpdr_vec)==length(pressure_profs)
        function dpdr0(i)
            return dpdr_vec[i]
        end
    elseif (dpdr_vec isa Nothing) && !(dpdr isa Nothing)
        function dpdr0(i)
            return dpdr
        end
    else
        error("Both dpdr and dpdr_vec defined in gen_equilibria - which one should we use?")
    end

    equilibria = Equilibrium[]

    if !(pressure_profs isa Union{Function,CubicSpline}) && length(pressure_profs)!=1
        @assert length(Jts)==length(pressure_profs)

        for i in 1:length(Jts)
            Bp,Bt,q,dpdr,p,Jts[i],Jp,Jts[i].xs[end],outerp6,Jp2 = Spline_Equil(Jts[i],pressure_profs[i],Bt0,R0; dpdr = dpdr0(i), plot_equil=false, print_mathematica_inputs=false)
            push!(equilibria,Equilibrium(Bp,Bt,q,dpdr,p,Jts[i],Jp,Jts[i].xs[end],R0,find_rs(q,0,1,Jts[i].xs[end];qtest=qtest,verbose=false),rs0))
        end
    else
        for i in 1:length(Jts)
            Bp,Bt,q,dpdr,p,Jts[i],Jp,Jts[i].xs[end],outerp6,Jp2 = Spline_Equil(Jts[i],pressure_profs,Bt0,R0; dpdr = dpdr0(i), plot_equil=false, print_mathematica_inputs=false)
            push!(equilibria,Equilibrium(Bp,Bt,q,dpdr,p,Jts[i],Jp,Jts[i].xs[end],R0,find_rs(q,0,1,Jts[i].xs[end];qtest=qtest,verbose=false),rs0))
        end
    end

    return equilibria
end

function gen_clean_equilibria(Jts,pressure_profs; Bt0=10, R0=3, dpdr_vec=nothing, dpdr=nothing, rs0=1.0, r0=1e-2,
    max_beta=0.1, qtest=nothing, maxq=nothing, minq=nothing, qedgerange=nothing, shear_max=nothing, shear_min=nothing, ideal_mhd=true, m1ncap = 7, m0ncap=3, nmax=8, del=1e-5, integrator_reltol=10^(-20), integrator_reltol_no_rs=1e-5, verbose=false, ideal_verbose=false, verify_sols=false, case=0, ignore_Suydam=false, return_psi_small=false)

    equilibria = gen_equilibria_Jts(Jts,pressure_profs; Bt0=Bt0, R0=R0, dpdr_vec=dpdr_vec, dpdr=dpdr, rs0=rs0, qtest=qtest)
    cleaners = equilibrium_cleaners(; r0=r0, max_beta=max_beta, qtest=qtest, maxq=maxq, minq=minq, qedgerange=qedgerange, shear_max=shear_max, shear_min=shear_min, ideal_mhd=ideal_mhd, m1ncap = m1ncap, m0ncap=m0ncap, nmax=nmax, del=del, integrator_reltol=integrator_reltol, integrator_reltol_no_rs=integrator_reltol_no_rs, verbose=verbose, ideal_verbose=ideal_verbose, verify_sols=verify_sols, case=case, ignore_Suydam=ignore_Suydam, return_psi_small=return_psi_small)

    return clean_vec(equilibria,cleaners)
end

function gen_n_clean_equilibria(Jtotmax,p,rb,batch_size,num_clean; 
                            maxbatches=100, maxJtbatches=1, verbose=true, Bt0=10, R0=3, dpdr_vec=nothing, 
                            dpdr=nothing, max_beta=0.1, qtest=nothing, maxq=nothing, minq=nothing, 
                            qedgerange=nothing, shear_max=nothing, shear_min=nothing,
                            ideal_mhd=true, m1ncap = 7, m0ncap=3, rs0=1.0, r0=1e-2, nmax=8, del=1e-5, integrator_reltol=10^(-20), integrator_reltol_no_rs=1e-5, ideal_verbose=false, verify_sols=false, case=0, ignore_Suydam=false, return_psi_small=false,
                            kwargs...)

    clean_counter = 0
    equilibria_vec = Equilibrium[]
    i=1

    while clean_counter < num_clean && (i<(maxbatches+1))
        verbose && print("Batch $(i) of $(batch_size) random equilibria being generated. Valid count = $(clean_counter).\n")
        Jt_temps = gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean; maxbatches=maxJtbatches, verbose=false, kwargs...)

        equil_temps = gen_clean_equilibria(Jt_temps,p; Bt0=Bt0, R0=R0, dpdr_vec=dpdr_vec, 
                                            dpdr=dpdr, rs0=rs0, r0=r0, max_beta=max_beta, qtest=qtest, maxq=maxq, minq=minq, 
                                            qedgerange=qedgerange, shear_max=shear_max, shear_min=shear_min,ideal_mhd=ideal_mhd, m1ncap = m1ncap, m0ncap=m0ncap, nmax=nmax, del=del, integrator_reltol=integrator_reltol, integrator_reltol_no_rs=integrator_reltol_no_rs, verbose=verbose, ideal_verbose=ideal_verbose, verify_sols=verify_sols, case=case, ignore_Suydam=ignore_Suydam, return_psi_small=return_psi_small)

        clean_counter += length(equil_temps)

        append!(equilibria_vec,  equil_temps)
        i+=1
    end

    print("$(clean_counter) of $(i*batch_size) randomly generated equilibria meet requirements.\n")

    return equilibria_vec
end

#########################################################################################
#Plotting Equilibria
#########################################################################################

function plot_equil(equilibrium::Equilibrium; plotrvec = range(0.000001,equilibrium.rb,200))
    p1=plot(plotrvec,equilibrium.Bp.(plotrvec),title = "Bp in Teslas",xlabel="r (m)",ylabel="T",label=false)
    p2=plot(plotrvec,equilibrium.Bt.(plotrvec),title = "Bt in Teslas",label=false,ylims=(0.0,2*equilibrium.Bt(equilibrium.rb)),xlabel="r (m)",ylabel="T")
    p3=plot(plotrvec,equilibrium.q.(plotrvec),title = "q",label=false,xlabel="r (m)")
    p4=plot(plotrvec,local_beta(equilibrium.p,equilibrium.Bt,equilibrium.Bp).(plotrvec),title = "Local Plasma β",xlabel="r (m)",label=false)
    p5=plot(plotrvec,equilibrium.Jt.(plotrvec),title = "Toroidal current density",xlabel="r (m)",ylabel="Amps/m^2",label=false)
    p6=plot(plotrvec,equilibrium.Jp.(plotrvec),title = "Poloidal current density",xlabel="r (m)",ylabel="Amps/m^2",label=false)

    outerp6 = plot(p1,p2,p3,p4,p5,p6)
    display(outerp6)
end

function plot_equil_short(equilibrium::Equilibrium; plotrvec = range(0.000001,equilibrium.rb,200),guidefontsize=3,titlefontsize=3,tickfontsize=2,kwargs...)
    p1=plot(plotrvec,equilibrium.Bp.(plotrvec),title = "Bp in Teslas",xlabel="r (m)",ylabel="T",label=false, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
    p2=plot(plotrvec,equilibrium.Bt.(plotrvec),title = "Bt in Teslas",label=false,ylims=(0.0,2*equilibrium.Bt(equilibrium.rb)),xlabel="r (m)",ylabel="T", titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
    p3=plot(plotrvec,equilibrium.q.(plotrvec),title = "q",label=false,xlabel="r (m)", titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
    p4=plot(plotrvec,local_beta(equilibrium.p,equilibrium.Bt,equilibrium.Bp).(plotrvec),title = "Local Plasma β",xlabel="r (m)",label=false, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

    display(plot(p1,p2,p3,p4))

    return p1,p2,p3,p4
end

function plot_equil(equilibria::AbstractArray{Equilibrium}; case=0, plotrvec = range(0.000001,equilibria[1].rb,200), titlefontsize=5, guidefontsize=5, tickfontsize=2, kwargs...)
    plots=[]

    if length(equilibria)==1 || case==1
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            layout = (1, 4)))
    elseif length(equilibria)==2 || case==2
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            layout = (2, 4)))
    elseif length(equilibria)==3 || case==3
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        c1,c2,c3,c4=plot_equil_short(equilibria[3];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            c1,c2,c3,c4,
                            layout = (3, 4)))
    elseif length(equilibria)==4 || case==4
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        c1,c2,c3,c4=plot_equil_short(equilibria[3];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        d1,d2,d3,d4=plot_equil_short(equilibria[4];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            c1,c2,c3,c4,
                            d1,d2,d3,d4,
                            layout = (4, 4)))
    elseif length(equilibria)>=5 || case==5
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        c1,c2,c3,c4=plot_equil_short(equilibria[3];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        d1,d2,d3,d4=plot_equil_short(equilibria[4];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
        e1,e2,e3,e4=plot_equil_short(equilibria[5];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            c1,c2,c3,c4,
                            d1,d2,d3,d4,
                            e1,e2,e3,e4,
                            layout = (5, 4)))
    end

    return pdog
end

function local_beta_(equilibrium::Equilibrium)
    return local_beta(equilibrium.p,equilibrium.Bt,equilibrium.Bp)
end

function pm_(equilibrium::Equilibrium)
    return r -> (1/(2*mu0))*(equilibrium.Bp(r)^2 + equilibrium.Bt(r)^2)
end

function plot_pm(equilibrium)
    plotrvec = range(0.000001,equilibrium.rb,200)
    display(plot(plotrvec,pm_(equilibrium).(plotrvec)))
end

function plot_local_beta(equilibrium)
    plotrvec = range(0.000001,equilibrium.rb,200)
    display(plot(plotrvec,local_beta_(equilibrium).(plotrvec)))
end

function plot_Suydam(equilibrium)
    test_Suydam(equilibrium.Bt, equilibrium.q, equilibrium.dpdr, equilibrium.rb; plotresults=true)
end

#########################################################################################
#Generate equilibrium cleaners
#########################################################################################

function qedge_test(equilibrium::Equilibrium,qedgerange)
    if minimum(qedgerange) < equilibrium.q(equilibrium.rb) <maximum(qedgerange)
        return true
    end
    
    return false
end

function beta_limit(equilibrium::Equilibrium, max_beta)
    pm(r) = (1/(2*mu0))*(equilibrium.Bp(r)^2 + equilibrium.Bt(r)^2)
    local_beta(r) = equilibrium.p(r)/pm(r)

    return max_exceeded(local_beta,equilibrium.rb,max_beta)
end

function q_one_rs(equilibrium::Equilibrium, qtest)
    num_rss = length(find_zeros(r->equilibrium.q(r)-qtest,0.0,equilibrium.rb))

    if num_rss!=1
        return false
    else
        return true
    end
end

function equilibrium_cleaners(; r0=1e-2, max_beta=0.1, qtest=nothing, maxq=nothing, minq=nothing, qedgerange=nothing, shear_max=nothing, shear_min=nothing, ideal_mhd=true, m1ncap = 7, m0ncap=3, nmax=8, del=1e-5, integrator_reltol=10^(-20), integrator_reltol_no_rs=1e-5, verbose=false, rs_verbose=false, ideal_verbose=false, verify_sols=false, case=0, ignore_Suydam=false, return_psi_small=false)
    cleaner_functions = Function[]

    if !(max_beta isa Nothing)
        push!(cleaner_functions, equilibrium -> beta_limit(equilibrium,max_beta))
    end

    if !(maxq isa Nothing)
        push!(cleaner_functions, equilibrium -> max_exceeded(equilibrium.q,equilibrium.rb,maxq))
    end
    if !(minq isa Nothing)
        push!(cleaner_functions, equilibrium -> min_exceeded(equilibrium.q,equilibrium.rb,minq))
    end

    if !(qedgerange isa Nothing)
        push!(cleaner_functions, equilibrium -> qedge_test(equilibrium,qedgerange))
    end

    if !(qtest isa Nothing)
        push!(cleaner_functions, equilibrium -> q_one_rs(equilibrium,qtest))
    end

    if !(shear_max isa Nothing) && (shear_min isa Nothing)
        push!(cleaner_functions, equilibrium -> gradmax_exceeded(equilibrium.q,equilibrium.rb,shear_max;abs_val=true))
    elseif !(shear_max isa Nothing) && !(shear_min isa Nothing)
        push!(cleaner_functions, equilibrium -> gradmax_exceeded(equilibrium.q,equilibrium.rb,shear_max;abs_val=false))
        push!(cleaner_functions, equilibrium -> gradmin_exceeded(equilibrium.q,equilibrium.rb,shear_min;abs_val=false))
    elseif (shear_max isa Nothing) && !(shear_min isa Nothing)
        push!(cleaner_functions, equilibrium -> gradmin_exceeded(equilibrium.q,equilibrium.rb,shear_min;abs_val=false))
    end

    if ideal_mhd
        push!(cleaner_functions,equilibrium -> test_ideal_stability(equilibrium.Bp, equilibrium.Bt, equilibrium.dpdr, equilibrium.R0, r0*equilibrium.rb, equilibrium.rb; m1ncap = m1ncap, m0ncap=m0ncap, rs0=equilibrium.rs0, nmax=nmax, del=del, integrator_reltol=integrator_reltol, integrator_reltol_no_rs=integrator_reltol_no_rs, verbose=ideal_verbose, rs_verbose=rs_verbose, verify_sols=verify_sols, case=case, ignore_Suydam=ignore_Suydam, return_psi_small=return_psi_small))
    end

    return cleaner_functions
end

#########################################################################################
#Run stability codes 
#########################################################################################

function run_Δl_Δr_calculator(equilibria, m, n, r0, nmax, del; integrator_reltol=1e-20, plot_soln_equil=false, report_err=false, run_zero_pressure=false)
    outequils = ResistiveEquilibrium[]
    inds = []

    k = k_(n, equilibria[1].R0)

    for (io,i) in enumerate(equilibria)    
        
        try                                    
            Δl,Δr,del2 = Δl_Δr_calculator(i.Bp, i.Bt, i.dpdr, k, m, r0, i.rs, i.rb, i.rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)

            if run_zero_pressure
                Δlzero, Δrzero, delzero = Δl_Δr_calculator_zeroPressure(i.Bp, i.Bt, i.dpdr, k, m, r0, i.rs, i.rb, i.rs0, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)

                push!(outequils,ResistiveEquilibrium(i,ΔprimeScrew(i.rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,Δlzero,Δrzero,Δlzero-Δrzero,delzero)))
            else
                push!(outequils,ResistiveEquilibrium(i,ΔprimeScrew(i.rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,nothing,nothing,nothing,nothing)))
            end
            push!(inds,io)
        catch
            if report_err
                Δl_Δr_calculator(i.Bp, i.Bt, i.dpdr, k, m, r0, i.rs, i.rb, i.rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)
                run_zero_pressure && Δl_Δr_calculator_zeroPressure(i.Bp, i.Bt, i.dpdr, k, m, r0, i.rs, i.rb, i.rs0, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)
                run_zero_pressure && ResistiveEquilibrium(i,ΔprimeScrew(i.rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,Δlzero,Δrzero,Δlzero-Δrzero,delzero))
                ResistiveEquilibrium(i,ΔprimeScrew(i.rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,nothing,nothing,nothing,nothing))
            end
            continue
        end
    end

    return outequils, inds
end

#########################################################################################
#Old Testing
#########################################################################################

if false
    #Jt_vec=randomJt(Jtotmax, Jtotmin, numJts, rb);

    Jt_vec = clean_vec(randomJt(Jtotmax, 5000, rb;J0bounds=5e5, Jedgebounds=1e5),generate_cleaners(rb ; Jtot_range = [Jtotmax, Jtotmin], maxgrad_width = rb/5,monotonic=false))
    plot_profiles(Jt_vec, rb;ylims = (0.0,1.5*Jtotmax/(pi*rb^2)))
    display(length(Jt_vec))

    Jtotmax = 1.5*total_plasma_current(Jt,rb)
    Jtotmin = 0.5*total_plasma_current(Jt,rb)
    Jtotrange=[Jtotmax, Jtotmin]

    #Jts = gen_n_clean_Jts(Jtotmax,rb,1000,10,generate_cleaners(rb, [Jtotmax, Jtotmin];Jt_maxgrad_width = rb/5,Jtmonotonic=true); J0bounds=5e5, Jedgebounds=1e5)
    Jts = gen_n_clean_Jts(Jtotmax,rb,5000,80,generate_cleaners(rb ; Jtot_range = Jtotrange,maxgrad_width = rb/100,use_fine=false,use_coarse=true, monotonic=false, wobbles=1, area_integrated_max=maximum(Jtotrange)); 
                            maxbatches=10, J0bounds=[0.8*Jt(0.0),1.2*Jt(0.0)], Jedgebounds=[1.2*Jt(2.0),0.8*Jt(2.0)])
    plot_profiles(Jts, rb;ylims = (0.0,1.4*Jt(0.0)))


    ps = gen_n_clean_pressure_profiles(1e20,rb,500,10,generate_cleaners(rb, [Jtotmax, Jtotmin]; maxgrad_width = rb/100,use_fine=false,use_coarse=true, monotonic=false, wobbles=1);
                            maxbatches=10, J0bounds=[0.8*Jt(0.0),1.2*Jt(0.0)], Jedgebounds=[1.2*Jt(2.0),0.8*Jt(2.0)])

    plot_profiles(randomJt(Jtotmax, 1000, rb;J0bounds=5e5, Jedgebounds=1e5),rb)

    Jtotmax = 1.5*total_plasma_current(Jt,rb)
    Jtotmin = 0.5*total_plasma_current(Jt,rb)
    Jtotrange=[Jtotmax, Jtotmin]


    batch_size=1000
    num_clean=20

    Jts = gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean; maxbatchesm=100, verbose=true, 
        max_ref = nothing, area_integrated_max=Jtotmax, Jtot_range=Jtotrange, use_coarse=true, 
        use_fine=true, maxval = nothing, wobbles = 3, monotonic=true, maxgrad_width = rb/6)
    plot_profiles(Jts, rb;ylims = (0.0,1.4*Jt(0.0)))

    gen_n_clean_equilibria(Jtotmax,p,rb,500,5; 
        maxbatches=10, verbose=true, Bt0=10, R0=3, dpdr_vec=nothing, 
        dpdr=nothing, max_beta=0.1, maxq=nothing, minq=nothing, 
        qedgerange=nothing, shear_max=nothing, shear_min=nothing,
        max_ref = nothing, area_integrated_max=Jtotmax, Jtot_range=Jtotrange, use_coarse=true, 
        use_fine=true, maxval = nothing, wobbles = 3, monotonic=true, maxgrad_width = rb/6
        )

    #Ideas:
        #Make a rough version of gradmax_exceeded, max_exceeded, num_wobbles and monotonic (should be way faster)
end