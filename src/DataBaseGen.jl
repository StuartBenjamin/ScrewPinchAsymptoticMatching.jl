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
    zeros = find_zeros(x -> f(r)-fmax,0.0,rb)

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
    zeros = find_zeros(x -> f(r)+fmin,0.0,rb)

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
        Jtot_range=nothing, use_coarse=true, use_fine=true, maxval = nothing, 
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

        #2
        monotonic && (push!(cleaner_functions,monotonic_coarse)) 
        push!(cleaner_functions,f -> gradmax_exceeded_coarse(f,spline_narrow_max_(nothing, rb, maxgrad_width;area_integrated_max=area_integrated_max); abs_val=true))
                                    
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

        #2 Derivative eval
        monotonic && (push!(cleaner_functions,monotonic_)) 
        push!(cleaner_functions,f -> gradmax_exceeded(f,rb,spline_narrow_max_(nothing, rb, maxgrad_width;area_integrated_max=area_integrated_max);abs_val=true))

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

function gen_n_clean_pressure_profiles(P0bounds,rb,batch_size,num_clean,cleaner_functions; maxbatches=100, verbose=true, kwargs...)
    clean_counter = 0
    pressure_vec = CubicSpline[]
    i=1

    while clean_counter < num_clean && (i<(maxbatches+1))
        verbose && print("Batch $(i) of $(batch_size) random splines being generated. Valid count = $(clean_counter).\n")
        p_temps = clean_vec(random_pressure(p0bounds, num_splines, rb; kwargs...),cleaner_functions)

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
    use_fine=true, maxval = nothing, wobbles = 3, monotonic=false, maxgrad_width = rb/6, kwargs...)

    return gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean,generate_cleaners(rb; max_ref = max_ref, 
        area_integrated_max=area_integrated_max, Jtot_range=Jtot_range, use_coarse=use_coarse, use_fine=use_fine, maxval = maxval, 
        wobbles = wobbles, monotonic=monotonic, maxgrad_width = maxgrad_width); maxbatches=maxbatches,verbose=verbose, kwargs...)
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
    Jt::Union{Function,CubicSpline,Nothing}
    Jp::Union{Function,Nothing}
    rb::Number
    R0::Number
    rs::Union{Number,Nothing}
end

function gen_equilibria_Jts(Jts,pressure_prof; Bt0=10, R0=3, dpdr_vec=nothing, dpdr=nothing)

    # Function currently set up to either take one pressure profile (w. optional dpdr), or a vector of pressure profiles (w. optional equal length vector of dpdrs)
    # Not yet able to take non-zero, non-equal length vectors of pressure and Jt, for example to cycle through 3 characteristic pressure proiles
    if (dpdr isa Nothing) &&  (dpdr_vec isa Nothing)
        dpdr0(i)=nothing
    elseif !(dpdr_vec isa Nothing) && (dpdr isa Nothing)
        @assert length(dpdr_vec)=length(pressure_profs)
        dpdr0(i)=dpdr_vec[i]
    elseif (dpdr_vec isa Nothing) && !(dpdr isa Nothing)
        dpdr0(i)=dpdr
    else
        error("Both dpdr and dpdr_vec defined in gen_equilibria - which one should we use?")
    end

    equilibria = Equilibrium[]

    if !(pressure_profs isa Union{Function,CubicSpline}) && length(pressure_profs)!=1
        @assert length(Jts)=length(pressure_profs)

        for i in 1:length(Jts)
            Bp,Bt,q,dpdr,p,Jt,Jp,Jt.xs[end],outerp6,Jp2 = Spline_Equil(Jts[1],pressure_profs[i],Bt0,R0; dpdr = dpdr0(i), plot_equil=false, print_mathematica_inputs=false)
            push!(equilibria,Equilibrium(Bp,Bt,q,dpdr,p,Jt,Jp,Jt.xs[end],R0,nothing))
        end
    else
        for i in 1:length(Jts)
            Bp,Bt,q,dpdr,p,Jt,Jp,Jt.xs[end],outerp6,Jp2 = Spline_Equil(Jts[1],pressure_profs,Bt0,R0; dpdr = dpdr0(i), plot_equil=false, print_mathematica_inputs=false)
            push!(equilibria,Equilibrium(Bp,Bt,q,dpdr,p,Jt,Jp,Jt.xs[end],R0,nothing))
        end
    end

    return equilibria
end

function gen_clean_equilibria(Jts,pressure_profs; Bt0=10, R0=3, dpdr_vec=nothing, dpdr=nothing, 
    max_beta=0.1, maxq=nothing, minq=nothing, qedgerange=nothing, shear_max=nothing, shear_min=nothing)

    equilibria = gen_equilibria_(Jts,pressure_profs; Bt0=Bt0, R0=R0, dpdr_vec=dpdr_vec, dpdr=dpdr)
    cleaners = equilibrium_cleaners(; max_beta=max_beta, maxq=maxq, minq=minq, qedgerange=qedgerange, shear_max=shear_max, shear_min=shear_min)

    return clean_vec(equilibria,cleaners)
end

function gen_n_clean_equilibria(Jtotmax,p,rb,batch_size,num_clean; 
                            maxbatches=100, maxJtbatches=1, verbose=true, Bt0=10, R0=3, dpdr_vec=nothing, 
                            dpdr=nothing, max_beta=0.1, maxq=nothing, minq=nothing, 
                            qedgerange=nothing, shear_max=nothing, shear_min=nothing,
                            kwargs...)

    clean_counter = 0
    equilibria_vec = Equilibrium[]
    i=1

    while clean_counter < num_clean && (i<(maxbatches+1))
        verbose && print("Batch $(i) of $(batch_size) random equilibria being generated. Valid count = $(clean_counter).\n")
        Jt_temps = gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean; maxbatches=maxJtbatches, verbose=false, kwargs...)

        equil_temps = gen_clean_equilibria(Jt_temps,p; Bt0=Bt0, R0=R0, dpdr_vec=dpdr_vec, 
                                            dpdr=dpdr, max_beta=max_beta, maxq=maxq, minq=minq, 
                                            qedgerange=qedgerange, shear_max=shear_max, shear_min=shear_min)

        clean_counter += length(equil_temps)

        append!(equilibria_vec,  equil_temps)
        i+=1
    end

    print("$(clean_counter) of $(i*batch_size) randomly generated equilibria meet requirements.\n")

    return equilibria_vec
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
    pm = r -> (1/(2*mu0))*(equilibrium.Bp(r)^2 + equilibrium.Bt(r)^2)
    local_beta(r) = r -> equilibrium.p(r)/pm(r)

    return max_exceeded(local_beta,equilibrium.rb,max_beta)
end

function equilibrium_cleaners(; max_beta=0.1, maxq=nothing, minq=nothing, qedgerange=nothing, shear_max=nothing, shear_min=nothing)
    cleaner_functions = Function[]

    if !(max_beta isa Nothing)
        push!(cleaner_functions, equilibrium -> qedge_test(equilibrium,qedgerange))
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

    if !(shear_max isa Nothing) && (shear_min isa Nothing)
        push!(cleaner_functions, equilibrium -> gradmax_exceeded(equilibrium.q,equilibrium.rb,shear_max;abs_val=true))
    elseif !(shear_max isa Nothing) && !(shear_min isa Nothing)
        push!(cleaner_functions, equilibrium -> gradmax_exceeded(equilibrium.q,equilibrium.rb,shear_max;abs_val=false))
        push!(cleaner_functions, equilibrium -> gradmin_exceeded(equilibrium.q,equilibrium.rb,shear_min;abs_val=false))
    elseif (shear_max isa Nothing) && !(shear_min isa Nothing)
        push!(cleaner_functions, equilibrium -> gradmin_exceeded(equilibrium.q,equilibrium.rb,shear_min;abs_val=false))
    end

    return cleaner_functions
end






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

Jts = gen_n_clean_Jts(Jtotmax,rb,batch_size,num_clean; maxbatches=100, verbose=true, 
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