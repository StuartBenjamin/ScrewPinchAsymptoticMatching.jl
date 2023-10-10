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
        wobbles = 3, monotonic=false, maxgrad_width = rb/6, diagnose=false)

    cleaner_functions = []
    labels = []
    
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
            push!(labels,"maxval coarse")
        end

        if !(minval isa Nothing)
            push!(cleaner_functions,f -> min_exceeded_coarse(f,minval))
            push!(labels,"minval coarse")
        end

        #2
        monotonic && (push!(cleaner_functions,monotonic_coarse); push!(labels,"monotonic coarse")) 
        push!(cleaner_functions,f -> gradmax_exceeded_coarse(f,spline_narrow_max_(max_ref, rb, maxgrad_width;area_integrated_max=area_integrated_max); abs_val=true))
        push!(labels,"gradmax exceeded coarse")
        
        #3
        wobbles>=0 && (push!(cleaner_functions,wobbles_coarse); push!(labels,"wobbles coarse"))  
    end

    #0 Integration
    !(Jtot_range isa Nothing) && (push!(cleaner_functions,Jt_tot_); push!(labels,"Jtot range check"))

    if use_fine
        #1 Function eval
        if !(maxval isa Nothing)
            push!(cleaner_functions,f -> max_exceeded(f,rb,maxval))
            push!(labels,"Maxval fine")
        end

        if !(minval isa Nothing)
            push!(cleaner_functions,f -> min_exceeded(f,rb,minval))
            push!(labels,"Minval fine")
        end

        #2 Derivative eval
        monotonic && (push!(cleaner_functions,monotonic_); push!(labels,"Monotonic fine")) 
        push!(cleaner_functions,f -> gradmax_exceeded(f,rb,spline_narrow_max_(max_ref, rb, maxgrad_width;area_integrated_max=area_integrated_max);abs_val=true))
        push!(labels,"gradmax exceeded fine")

        #3 Double-derivative eval
        wobbles>=0 && (push!(cleaner_functions,wobbles_); push!(labels,"wobbles fine"))
    end

    diagnose && display(labels)

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

function randomJt(Jt_upper_bound, numJts, rb; knotmax=10, knotmin=5, minJt=0.0, 
    plot_profs=false, J0bounds=nothing, Jedgebounds=nothing, useDirichlet=true, dirichlet_alpha=10)

    knotnums = random_Int_in_range(knotmax,knotmin,numJts)

    splines_vals = random_in_range(Jt_upper_bound,minJt,(numJts,knotmax))

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
        if !useDirichlet
            Jt = CubicSpline(range(0.0,rb;length=knotnums[i]),vcat(start_vals[i],splines_vals[i,2:(knotnums[i]-1)],end_vals[i]))
        else
            x0 = 0.0
            xend = rb
            xintermediate = rb.*cumsum(rand(Dirichlet(knotnums[i]-1,dirichlet_alpha)))
            xs = vcat(x0,xintermediate[1:end-1],xend)

            Jt = CubicSpline(xs,vcat(start_vals[i],splines_vals[i,2:(knotnums[i]-1)],end_vals[i]))
        end
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

function clean_vec(in_vec,cleaner_functions;diagnose=false)
    valid = trues(length(in_vec))
    valid_inds = 1:length(in_vec)

    if diagnose
        num_fails = zeros(Int,length(cleaner_functions))
        for (io,cleaner_func) in enumerate(cleaner_functions)
            num_fails[io] += length(filter(x->!cleaner_func(x),in_vec))
        end
        display(num_fails)
    end
    for cleaner_func in cleaner_functions
        valid_inds = filter(x->cleaner_func(in_vec[x]),valid_inds) 
    end

    return in_vec[valid_inds]
end

function gen_n_clean_Jts(Jt_upper_bound,rb,batch_size,num_clean,cleaner_functions; maxbatches=100, verbose=true, kwargs...)
    clean_counter = 0
    Jt_vec = CubicSpline[]
    i=1

    while clean_counter < num_clean && (i<(maxbatches+1))
        verbose && print("Batch $(i) of $(batch_size) random splines being generated. Valid count = $(clean_counter).\n")
        Jt_temps = clean_vec(randomJt(Jt_upper_bound, batch_size, rb; kwargs...),cleaner_functions)

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
        p5=plot(plotrvec,Jt.(plotrvec),title = "Toroidal current density",xlabel="r (m)",ylabel="\$\\textrm{A/m}^{2}\$",label=false,ylims=ylims)
        p5=scatter!(Jt.xs,Jt.(Jt.xs),title = "Toroidal current",xlabel="r (m)",ylabel="\$\\textrm{A/m}^{2}\$",label=false,ylims=ylims)
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

function gen_n_clean_Jts(Jt_upper_bound,rb,batch_size,num_clean; maxbatches=100, verbose=true, 
    max_ref = nothing, area_integrated_max=nothing, Jtot_range=nothing, use_coarse=true, 
    use_fine=true, maxval = nothing, minval=nothing, wobbles = 3, monotonic=false, maxgrad_width = rb/6, kwargs...)

    return gen_n_clean_Jts(Jt_upper_bound,rb,batch_size,num_clean,generate_cleaners(rb; max_ref = max_ref, 
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
    dpdr::Function #Scaled by mu0
    p::Union{Function,CubicSpline}
    Jt::Union{Function,CubicSpline}
    Jp::Union{Function}
    rb::Number
    R0::Number
    rs::Union{Number,Nothing}
    rs0::Number #Deprecated re-scaling of r, set to 1.0 (defined for Chandra equilibria but not splines).
end

function Base.show(io::IO, equilibrium::Equilibrium) 
    if equilibrium.Jt isa CubicSpline
        print(io, "Resistive equilibrium, Jtot = ",round(total_plasma_current(equilibrium.Jt);sigdigits=3),".")
    else
        print(io, "Resistive equilibrium, Jtot = ",round(total_plasma_current(equilibrium.Jt,equilibrium.rb);sigdigits=3),".")
    end
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

function Base.show(io::IO, Δprimes::ΔprimeScrew) 
    print(io, "Resistive $(Δprimes.m)/$(Δprimes.n) mode analysis")
end

struct ResistiveEquilibrium
    equilibrium::Equilibrium
    Δprimes::Union{ΔprimeScrew,AbstractArray{ΔprimeScrew}}
end

function Base.show(io::IO, req::ResistiveEquilibrium) 
    if req.Δprimes isa AbstractArray
        mnlist = ""
        for i in req.Δprimes
            mnlist = mnlist*"$(i.m)/$(i.n) "
        end

        print(io, "Resistive equilibrium, Jtot = ",round(total_plasma_current(req.equilibrium.Jt);sigdigits=3),". Modes analysed = ",mnlist)
    elseif req.Δprimes isa ΔprimeScrew
        print(io, "Resistive equilibrium, Jtot = ",round(total_plasma_current(req.equilibrium.Jt);sigdigits=3),". Modes analysed = ","$(req.Δprimes.m)/$(req.Δprimes.n) ")
    else 
        print(io, "Resistive equilibrium, Jtot = ",round(total_plasma_current(req.equilibrium.Jt);sigdigits=3),".")
    end
end

function Dr(requilibrium::ResistiveEquilibrium)
    if requilibrium.Δprimes isa ΔprimeScrew
        qprime = ForwardDiff.derivative(requilibrium.equilibrium.q,requilibrium.Δprimes.rs)
        return -2*(1/mu0)*requilibrium.equilibrium.dpdr(requilibrium.Δprimes.rs)*requilibrium.Δprimes.rs*requilibrium.Δprimes.n^2/(requilibrium.equilibrium.Bp(requilibrium.Δprimes.rs)^2*requilibrium.equilibrium.R0^2*qprime^2)
    else 
        Drs = Real[]
        for i in requilibrium.Δprimes
            qprime = ForwardDiff.derivative(requilibrium.equilibrium.q,i.rs)
            push!(Drs, -2*(1/mu0)*requilibrium.equilibrium.dpdr(i.rs)*i.rs*i.n^2/(requilibrium.equilibrium.Bp(i.rs)^2*requilibrium.equilibrium.R0^2*qprime^2))
        end

        return Drs
    end
end

function print_equil_data_(equil_series::Union{AbstractArray{ResistiveEquilibrium},AbstractArray{Equilibrium}}; directoryprefactor="", kwargs...)
    path=pwd()
    for i in 1:length(equil_series)
        mkdir(directoryprefactor*"_$(i)")
        cd(directoryprefactor*"_$(i)")

        print_equil_data_(equil_series[i]; kwargs...)
        cd(path)
    end
end


function print_equil_data_(equilibrium::Union{ResistiveEquilibrium,Equilibrium}; kwargs...)
    if equilibrium isa ResistiveEquilibrium
        return print_equil_data_(equilibrium.equilibrium.Jt,equilibrium.equilibrium.Bt,equilibrium.equilibrium.p; kwargs...)
    else
        return print_equil_data_(equilibrium.Jt,equilibrium.Bt,equilibrium.p; kwargs...)
    end
end

function print_equil_data_(Jt, Bt, p; rvec=nothing, filename_prefactor="", destination="")
    if rvec isa Nothing
        throw("Define vector of r-values 'rvec' to print equilibrium data on.")
    end

    Jtvec = Jt.(rvec)
    Btvec = Bt.(rvec)
    pvec = p.(rvec)

    open(string(destination,filename_prefactor,"profile_j.txt"), "w") do io
        #writedlm(io, [rvec Jtvec], ' ')
        for (ioo,o) in enumerate(rvec)
            @printf(io, "%f %f\n", o,Jtvec[ioo])
        end
    end
    open(string(destination,filename_prefactor,"profile_f.txt"), "w") do io
        #writedlm(io, [rvec Btvec], ' ')
        for (ioo,o) in enumerate(rvec)
            @printf(io, "%f %f\n", o,Btvec[ioo])
        end
    end
    open(string(destination,filename_prefactor,"profile_p.txt"), "w") do io
        #writedlm(io, [rvec pvec], ' ')
        for (ioo,o) in enumerate(rvec)
            @printf(io, "%f %f\n", o,pvec[ioo])
        end
    end

    return Jtvec,Btvec,pvec
end


##

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

function gen_n_clean_equilibria(Jt_upper_bound,p,rb,batch_size,num_clean; 
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
        Jt_temps = nothing

        try
            Jt_temps = gen_n_clean_Jts(Jt_upper_bound,rb,batch_size,num_clean; maxbatches=maxJtbatches, verbose=false, kwargs...)
        catch
            print("Batch $(i) Jt generation failed\n.")
            i+=1
            continue
        end
        equil_temps = nothing

        try
            equil_temps = gen_clean_equilibria(Jt_temps,p; Bt0=Bt0, R0=R0, dpdr_vec=dpdr_vec, 
                                                dpdr=dpdr, rs0=rs0, r0=r0, max_beta=max_beta, qtest=qtest, maxq=maxq, minq=minq, 
                                                qedgerange=qedgerange, shear_max=shear_max, shear_min=shear_min,ideal_mhd=ideal_mhd, m1ncap = m1ncap, m0ncap=m0ncap, nmax=nmax, del=del, integrator_reltol=integrator_reltol, integrator_reltol_no_rs=integrator_reltol_no_rs, verbose=verbose, ideal_verbose=ideal_verbose, verify_sols=verify_sols, case=case, ignore_Suydam=ignore_Suydam, return_psi_small=return_psi_small)

            clean_counter += length(equil_temps)

            append!(equilibria_vec,  equil_temps)
            i+=1
        catch
            print("Batch $(i) equilibria generation and cleaning failed\n.")
            i+=1
        end
    end

    print("$(clean_counter) of $(i*batch_size) randomly generated equilibria meet requirements.\n")

    return equilibria_vec
end

#########################################################################################
#Plotting Equilibria
#########################################################################################

function plot_equil(equilibrium::Equilibrium; plotrvec = range(0.000001,equilibrium.rb,200), Jt_ref = nothing, titlefontsize=12, guidefontsize=8, tickfontsize=6)
    p1=plot(plotrvec,equilibrium.Bp.(plotrvec),title = "Bp",xlabel="r (m)",ylabel="T",label=false,titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)
    p2=plot(plotrvec,equilibrium.Bt.(plotrvec),title = "Bt",label=false,ylims=(0.0,2*equilibrium.Bt(equilibrium.rb)),xlabel="r (m)",ylabel="T",titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)
    p3=plot(plotrvec,equilibrium.q.(plotrvec),title = "q",label=false,xlabel="r (m)",titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)
    p4=plot(plotrvec,local_beta(equilibrium.p,equilibrium.Bt,equilibrium.Bp).(plotrvec),title = "Plasma β",xlabel="r (m)",label=false,titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)

    if Jt_ref isa Nothing
        p5=plot(plotrvec,equilibrium.Jt.(plotrvec),title = "Jt density",xlabel="r (m)",ylabel="\$\\textrm{A/m}^{2}\$",label=false, ylims=nothing,titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)
    else
        p5=plot(plotrvec,equilibrium.Jt.(plotrvec),title = "Jt density",xlabel="r (m)",ylabel="\$\\textrm{A/m}^{2}\$",label=false, ylims = (0.0,Jt_ref),titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)
    end
    p6=plot(plotrvec,equilibrium.Jp.(plotrvec),title = "Jp density",xlabel="r (m)",ylabel="\$\\textrm{A/m}^{2}\$",label=false,titlefontsize=titlefontsize, guidefontsize=guidefontsize,tickfontsize=tickfontsize)

    outerp6 = plot(p1,p2,p3,p4,p5,p6)
    display(outerp6)
end

function plot_equil_short(equilibrium::Equilibrium; plotrvec = range(0.000001,equilibrium.rb,200),guidefontsize=3,titlefontsize=3,tickfontsize=2, ylims=nothing, kwargs...)
    p1=plot(plotrvec,equilibrium.Bp.(plotrvec),title = "Bp",xlabel="r (m)",ylabel="T",label=false, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
    p2=plot(plotrvec,equilibrium.Bt.(plotrvec),title = "Bt",label=false,ylims=(0.0,2*equilibrium.Bt(equilibrium.rb)),xlabel="r (m)",ylabel="T", titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
    p3=plot(plotrvec,equilibrium.q.(plotrvec),title = "q",label=false,xlabel="r (m)", titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)
    p4=plot(plotrvec,local_beta(equilibrium.p,equilibrium.Bt,equilibrium.Bp).(plotrvec),title = "Local Plasma β",xlabel="r (m)",label=false, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize, kwargs...)

    display(plot(p1,p2,p3,p4))

    return p1,p2,p3,p4
end

function plot_equil(equilibrium::ResistiveEquilibrium;kwargs...)
    return plot_equil(equilibrium.equilibrium;kwargs...)
end

function plot_equil(equilibria::AbstractArray{Equilibrium}; case=0, Jt_ref=nothing, plotrvec = range(0.000001,equilibria[1].rb,200), titlefontsize=6, guidefontsize=6, tickfontsize=3, kwargs...)
    plots=[]

    if Jt_ref isa Nothing
        ylims = nothing
    else
        ylims = (0.0,Jt_ref)
    end

    if length(equilibria)==1 || case==1
        return plot_equil(equilibria[1];plotrvec=plotrvec, Jt_ref=Jt_ref)
    elseif length(equilibria)==2 || case==2
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            layout = (2, 4)))
    elseif length(equilibria)==3 || case==3
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        c1,c2,c3,c4=plot_equil_short(equilibria[3];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            c1,c2,c3,c4,
                            layout = (3, 4)))
    elseif length(equilibria)==4 || case==4
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        c1,c2,c3,c4=plot_equil_short(equilibria[3];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        d1,d2,d3,d4=plot_equil_short(equilibria[4];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            c1,c2,c3,c4,
                            d1,d2,d3,d4,
                            layout = (4, 4)))
    elseif length(equilibria)>=5 || case==5
        a1,a2,a3,a4=plot_equil_short(equilibria[1];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        b1,b2,b3,b4=plot_equil_short(equilibria[2];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        c1,c2,c3,c4=plot_equil_short(equilibria[3];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        d1,d2,d3,d4=plot_equil_short(equilibria[4];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)
        e1,e2,e3,e4=plot_equil_short(equilibria[5];plotrvec=plotrvec, titlefontsize=titlefontsize, guidefontsize=guidefontsize, tickfontsize=tickfontsize,  ylims=ylims, kwargs...)

        pdog=display(plot(a1,a2,a3,a4,
                            b1,b2,b3,b4,
                            c1,c2,c3,c4,
                            d1,d2,d3,d4,
                            e1,e2,e3,e4,
                            layout = (5, 4)))
    end

    return pdog
end

function plot_equil(equilibria::AbstractArray{ResistiveEquilibrium}; kwargs...)
    return plot_equil([i.equilibrium for i in equilibria];kwargs...)
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

function plot_current_profiles(equilibria::AbstractArray{Equilibrium};plotrvec = range(0.000001,equilibria[1].rb,200),rss=nothing)
    p1 = nothing
    for (io,i) in enumerate(equilibria)
        if io==1 
            p1 = plot(plotrvec,i.Jt.(plotrvec);label=false)
            !(rss isa Nothing) && display(vline!([rss[io]];label=false))
        else
            p1 = plot!(plotrvec,i.Jt.(plotrvec);title = "Toroidal current profiles",label=false, xlabel="r (m)", ylabel="Jt density (\$\\textrm{A/m}^{2}\$)")
            !(rss isa Nothing) && display(vline!([rss[io]];label=false))
        end
    end
    display(p1)

    return p1
end

function plot_current_profiles(i::ResistiveEquilibrium;plotrvec = range(0.000001,i.equilibrium.rb,200),m_ind=nothing)
    p1 = nothing

    p1 = plot(plotrvec,i.equilibrium.Jt.(plotrvec);label=false)
    !(m_ind isa Nothing) && (p1=vline!([i.Δprimes[m_ind].rs];label=false))
    display(p1)

    return p1
end

function plot_current_profiles(equilibria::AbstractArray{ResistiveEquilibrium}; m_ind=nothing, n=1, plotrvec = range(0.000001,equilibria[1].equilibrium.rb,200))
    rss=nothing
    if !(m_ind isa Nothing)
        if length(equilibria[1].Δprimes) > 1
            rss = [i.Δprimes[m_ind].rs for i in equilibria]
        else
            rss = [i.Δprimes.rs for i in equilibria]
        end
    end

    return plot_current_profiles([i.equilibrium for i in equilibria];plotrvec = plotrvec, rss=rss)
end  

function plot_q_profiles(equilibria::AbstractArray{ResistiveEquilibrium}; m_ind=nothing, n=1, plotrvec = range(0.000001,equilibria[1].equilibrium.rb,200))
    rss=nothing
    if !(m_ind isa Nothing)
        if length(equilibria[1].Δprimes) > 1
            rss = [i.Δprimes[m_ind].rs for i in equilibria]
        else
            rss = [i.Δprimes.rs for i in equilibria]
        end
    end

    p1 = nothing
    for (io,i) in enumerate(equilibria)
        if io==1 
            p1 = plot(plotrvec,i.equilibrium.q.(plotrvec);label=false)
            !(rss isa Nothing) && display(vline!([rss[io]];label=false))
        else
            p1 = plot!(plotrvec,i.equilibrium.q.(plotrvec);title = "q profiles",label=false, xlabel="r (m)")
            !(rss isa Nothing) && display(vline!([rss[io]];label=false))
        end
    end
    display(p1)

    return p1
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

function run_Δl_Δr_calculator(equilibria, m::Int, n::Int, r0, nmax, del; integrator_reltol=1e-20, plot_soln_equil=false, report_err=false, run_zero_pressure=false)
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
            print("Δl_Δr_calculation failed\n.")
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


function run_Δl_Δr_calculator(equilibria, ms::AbstractArray{Int}, n::Int, r0, nmax, del; integrator_reltol=1e-20, plot_soln_equil=false, report_err=false, run_zero_pressure=false)
    outequils = ResistiveEquilibrium[]
    inds = []

    k = k_(n, equilibria[1].R0)

    for (io,i) in enumerate(equilibria)    
        try     
            Δprimes = ΔprimeScrew[]
            
            for m in ms
                rs = find_rs(i.q,m,n,i.Jt.xs[end];verbose=false)
                if rs > 0
                    Δl,Δr,del2 = Δl_Δr_calculator(i.Bp, i.Bt, i.dpdr, k, m, r0, rs, i.rb, i.rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)

                    if run_zero_pressure
                        Δlzero, Δrzero, delzero = Δl_Δr_calculator_zeroPressure(i.Bp, i.Bt, i.dpdr, k, m, r0, rs, i.rb, i.rs0, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)
                        push!(Δprimes,ΔprimeScrew(rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,Δlzero,Δrzero,Δlzero-Δrzero,delzero))
                    else
                        push!(Δprimes,ΔprimeScrew(rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,nothing,nothing,nothing,nothing))
                    end
                end
            end

            push!(outequils,ResistiveEquilibrium(i,Δprimes))
            push!(inds,io)
        catch
            print("Δl_Δr_calculation failed\n.")
            if report_err  #same code to identify the error
                Δprimes = ΔprimeScrew[]
            
                for m in ms
                    rs = find_rs(i.q,m,n,i.Jt.xs[end];verbose=false)
                    if rs > 0
                        Δl,Δr,del2 = Δl_Δr_calculator(i.Bp, i.Bt, i.dpdr, k, m, r0, rs, i.rb, i.rs0, nmax, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)
    
                        if run_zero_pressure
                            Δlzero, Δrzero, delzero = Δl_Δr_calculator_zeroPressure(i.Bp, i.Bt, i.dpdr, k, m, r0, rs, i.rb, i.rs0, del; integrator_reltol=integrator_reltol, plot_solution=plot_soln_equil, plot_soln = plot_soln_equil)
                            push!(Δprimes,ΔprimeScrew(rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,Δlzero,Δrzero,Δlzero-Δrzero,delzero))
                        else
                            push!(Δprimes,ΔprimeScrew(rs,m,n,Δl,Δr,Δl+Δr,del2,nmax,nothing,nothing,nothing,nothing))
                        end
                    end
                end
    
                #push!(outequils,ResistiveEquilibrium(i,Δprimes))
                #push!(inds,io)
            end
            continue
        end
    end

    return outequils, inds
end

#########################################################################################
#Gen euilibria and run stability codes
#########################################################################################

function gen_equil_run_Δl_Δr(f_gen_equilibria::Function, batch_size, num_clean, ms::AbstractArray{Int}, n, r0, nmax, del; filename = nothing, loop_big_batch=20, path=nothing, kwargs...)
    cd(path)
    resistive_equils_store = []
    equils_and_inds = []

    for i in 1:loop_big_batch
        equils_loop=f_gen_equilibria(batch_size,num_clean)

        resistive_equils, inds = run_Δl_Δr_calculator(equils_loop, ms, n, r0, nmax, del; kwargs...)

        push!(resistive_equils_store, resistive_equils)
        push!(equils_and_inds, (equils_loop, inds))

        if !(filename isa Nothing)
            @save filename resistive_equils_store equils_and_inds
        end
    end

    return resistive_equils_store, equils_and_inds
end


#########################################################################################
#Well Analysis
#########################################################################################

function rs_near_local_minmax(res_euil::AbstractArray,min_closeness,sep_from_edge; kwargs...)
    sutblefunc(i) = rs_near_local_minmax(i,min_closeness,sep_from_edge; kwargs...)
    inds = Int[]
    well_equils = ResistiveEquilibrium[]
    norm_distances_to_well = Number[]
    well_widths = Number[]
    well_gradients = Number[]
    Jt_double_derivatives = Number[]
    wellhill_closenesses = Number[]


    for (io,o) in enumerate(res_euil)
        tempresult = sutblefunc(o)
        if !(tempresult==false) 
            norm_distance_to_wellhill, wellhill_width, well_gradient, Jtdd, wellhill_closeness = tempresult
            push!(inds,io)
            push!(well_equils,o)
            push!(norm_distances_to_well,norm_distance_to_wellhill)
            push!(well_widths,wellhill_width)
            push!(well_gradients,well_gradient)
            push!(Jt_double_derivatives,Jtdd)
            push!(wellhill_closenesses,wellhill_closeness)
        end
    end

    if length(inds) == 0
        @warn "No wells/hills found"
        return nothing,nothing,nothing,nothing,nothing,nothing,nothing 
    end

    return well_equils, norm_distances_to_well, well_widths, well_gradients, Jt_double_derivatives, wellhill_closenesses, inds
end

function rs_near_local_minmax(res_euil,min_closeness,sep_from_edge; well=true, m_ind=1,return_wellhill_info=false)
    return rs_near_local_minmax(res_euil.equilibrium.Jt,res_euil.equilibrium.rb,res_euil.Δprimes[m_ind].rs,min_closeness,sep_from_edge; well=well, return_wellhill_info=return_wellhill_info)
end

function rs_near_local_minmax(Jt,rb,rs,min_closeness,sep_from_edge; well=true, return_wellhill_info=false)
    Jtd(r) = ForwardDiff.derivative(Jt,r)
    Jtdd(r) = ForwardDiff.derivative(Jtd,r)

    local_minmax = []
    local_minmax_ids = Int[]
    zero_gradients = find_zeros(r -> Jtd(r),0.0,rb)

    if length(zero_gradients) == 0
        return false
    end

    for (io,i) in enumerate(zero_gradients)
        if well
            if Jtdd(i) > 0
                push!(local_minmax,i)
                push!(local_minmax_ids,io)
            end
        else
            if Jtdd(i) < 0
                push!(local_minmax,i)
                push!(local_minmax_ids,io)
            end
        end
    end

    if length(local_minmax) == 0
        return false
    end

    norm_distance_to_wellhill = minimum(abs.(local_minmax.-rs))/rb
    closest_wellhill = local_minmax[argmin(abs.(local_minmax.-rs))]
    closest_wellhill_id = local_minmax_ids[argmin(abs.(local_minmax.-rs))]

    if closest_wellhill_id==1
        LHS_dist = closest_wellhill
    else
        LHS_dist=abs(zero_gradients[closest_wellhill_id-1]-zero_gradients[closest_wellhill_id])
    end
    if closest_wellhill_id==length(zero_gradients)
        RHS_dist = rb-closest_wellhill
    else
        RHS_dist=abs(zero_gradients[closest_wellhill_id+1]-zero_gradients[closest_wellhill_id])
    end

    wellhill_width = minimum([LHS_dist,RHS_dist])
    wellhill_closeness = minimum(abs.(local_minmax.-rs))/wellhill_width

    if well
        min_edge_val = minimum([Jt(maximum([closest_wellhill-wellhill_width,1e-15])),Jt(minimum([closest_wellhill+wellhill_width,rb]))])
    else
        max_edge_val = maximum([Jt(maximum([closest_wellhill-wellhill_width,1e-15])),Jt(minimum([closest_wellhill+wellhill_width,rb]))])
    end

    if well
        wellhill_gradient = (abs(min_edge_val-Jt(closest_wellhill))/wellhill_width)/Jt(closest_wellhill)
    else
        wellhill_gradient = (abs(max_edge_val-Jt(closest_wellhill))/wellhill_width)/Jt(closest_wellhill)
    end

    #Are you in a well/on a hill?
    #Outside bounds
    (norm_distance_to_wellhill > wellhill_width) && (return false)
    #Too far up side/down side
    if well 
        if (Jt(rs) > minimum([Jt(maximum([rs-LHS_dist,1e-15])),Jt(minimum([rs+RHS_dist,rb]))]))
            return false
        end
    else
        if (Jt(rs) < maximum([Jt(maximum([rs-LHS_dist,1e-15])),Jt(minimum([rs+RHS_dist,rb]))]))
            return false
        end
    end

    if !return_wellhill_info
        if norm_distance_to_wellhill < min_closeness && (abs(closest_wellhill-rb) > sep_from_edge*rb)
            return true
        end
        return false
    end

    if return_wellhill_info
        return norm_distance_to_wellhill, wellhill_width, wellhill_gradient, Jtdd(closest_wellhill), wellhill_closeness
    end

    return false
end