"""
Integration of spline from point x1 to point x2. Assume both x1 and x2 are in the range of interpolated data.
"""
function integral(spline::CubicSpline, x1::Real, x2::Number)
    idx1,idx2 = _checkbounds(spline, x1, x2)

    # Both points in same interval
    if idx1==idx2
        xref = spline.xs[idx1]
        old_coeff = spline.ps[idx1,:]

        int_coeff = _intpoly(old_coeff,xref)

        x1_val = @evalpoly(x1-xref, int_coeff...) #correct to a constant
        x2_val = @evalpoly(x2-xref, int_coeff...) #correct to a constant

        return x2_val-x1_val #constant has no effect
    end

    # Assume idx2 > idx2
    signvar = 1
    # If not, relabel x1, idx1 and x2, idx2, and swap sign of final answer:
    if idx2<idx1 
        hold = idx2
        idx2 = idx1
        idx1 = hold

        hold2 = x2
        x2 = x1
        x1 = hold2

        signvar = -1
    end

    # Evaluating contribution from idx1 
    idx1_int_coeff = _intpoly(spline.ps[idx1,:],spline.xs[idx1])

    idx1_rboundary = spline.xs[idx1+1]

    idx1_val1 = @evalpoly(x1-spline.xs[idx1], idx1_int_coeff...)
    idx1_val2 = @evalpoly(idx1_rboundary-spline.xs[idx1], idx1_int_coeff...)

    idx1_int = idx1_val2 - idx1_val1

    # Evaluating contribution from idx2
    idx2_int_coeff = _intpoly(spline.ps[idx2,:],spline.xs[idx2])

    idx2_lboundary = spline.xs[idx2]

    idx2_val1 = @evalpoly(idx2_lboundary-spline.xs[idx2], idx2_int_coeff...) 
    idx2_val2 = @evalpoly(x2-spline.xs[idx2], idx2_int_coeff...)

    idx2_int = idx2_val2 - idx2_val1

    # Evaluating contributions from full-width integrated intervals
    full_widths = idx2-idx1-1 #number of full-width integration intervals

    runsum = 0.0
    for i in 1:full_widths
        idx_int_coeff = _intpoly(spline.ps[idx1+i,:],spline.xs[idx1+i])
    
        idx_val1 = @evalpoly(spline.xs[idx1+i]-spline.xs[idx1+i], idx_int_coeff...)
        idx_val2 = @evalpoly(spline.xs[idx1+i+1]-spline.xs[idx1+i], idx_int_coeff...)

        runsum += (idx_val2-idx_val1)
    end

    return signvar*(idx1_int+idx2_int+runsum)
end

"""
Calculates the coefficients for the integrated polynomial (assumes polynomial is of degree 3 or less).
"""
function _intpoly(coeff::AbstractArray{<:Real,1},xref::Real)
    # Integrated polynomial uses same shifted form (x-xref)^n so @evalpoly can be used
    coeff = coeff[1:end] ./ collect(1:(length(coeff)))

    # New constant is required to write the integrated polynomial in the above form
    newconst = -(-coeff[1]*xref + (coeff[2]*xref^2)/2 - (coeff[3]*xref^3)/3 + (coeff[4]*xref^4)/4)

    return vcat(newconst,coeff)
end

"""
Checks the integral bounds x1 and x2 lie within spline interpolation range. No integration into extrapolated region at this point.
Also returns x1 and x2's corresponding spline intervals.
"""
function _checkbounds(spline::CubicSpline, x1::Real, x2::Number)
    @assert x1 != x2

    # Checking intervals of both points
    # Does x1 sit on upper interval border?
    if x1 == spline.xs[end]
        idx1 = length(spline.xs)-1

    # Find the corresponding interval for x1
    else

        # Find the interval index
        idx1 = _binary_search_interval(spline.xs, x1)
    end
    # Does x2 sit on upper interval border?
    if x2 == spline.xs[end]
        idx2 = length(spline.xs)-1

    # Find the corresponding interval for x2
    else

        # Find the interval index
        idx2 = _binary_search_interval(spline.xs, x2)
    end

    # Making sure x1, x2 in boundaries. No integration in the extraplolated region included at this point.
    if idx1 == 0
        throw("x1 too small ($x1 < $(spline.xs[1]))")
    elseif idx1 == length(spline.xs)
        throw("x1 too big ($x1 > $(spline.xs[end]))")
    end
    if idx2 == 0
        throw("x2 too small ($x2 < $(spline.xs[1]))")
    elseif idx2 == length(spline.xs)
        throw("x2 too big ($x2 > $(spline.xs[end]))")
    end

    return idx1,idx2
end

"""
Calculates the coefficients for the integral of j*Jt, where Jt is a cubic polynomial of degree 3 or less.
Constant of integration is set to 0. 
"""
function _int_rJt(coeff::AbstractArray{<:Real,1},xref::Real)
    p1=coeff[1]
    p2=coeff[2]
    p3=coeff[3]
    p4=coeff[4]

    np = zeros(Real, 6)

    np[6] = p4/5
    np[5] = p3/4 - (3*p4*xref)/4
    np[4] = p2/3 - (2*p3*xref)/3 + p4*xref^2
    np[3] = p1/2 - (p2*xref)/2 + (p3*xref^2)/2 - (p4*xref^3)/2
    np[2] = 0.0
    np[1] = 0.0

    return np
end

"""
Calculates the coefficients for the integral of Jt/r, where Jt is a cubic polynomial of degree 3 or less.
Constant of integration is set to 0. Don't evaluate this from r=0. 
"""
function _int_Jt_on_r(coeff::AbstractArray{<:Real,1},xref::Real)
    p1=coeff[1]
    p2=coeff[2]
    p3=coeff[3]
    p4=coeff[4]

    np = zeros(Real, 4)

    np[4] = p4/3
    np[3] = (1/2)*(p3 - 3*p4*xref)
    np[2] = p2 - 2*p3*xref + 3*p4*xref^2
    np[1] = 0.0

    log_coeff = p1 - p2*xref + p3*xref^2 - p4*xref^3

    return np,log_coeff
end

"""
Calculates the coefficients for the integral of I[Jtr]Jt/r, where Jt is a cubic polynomial of degree 3 or less.
Constant of integration is set to 0. 
"""
function _int_full(coeff::AbstractArray{<:Real,1},xref::Real)
    p1=coeff[1]
    p2=coeff[2]
    p3=coeff[3]
    p4=coeff[4]

    np = zeros(Real, 9)

    np[9] = p4^2/40
    np[8] = (9*p3*p4)/140 - (27*p4^2*xref)/140
    np[7] = p3^2/24 + (4*p2*p4)/45 - (77*p3*p4*xref)/180 + (77*p4^2*xref^2)/120
    np[6] = (7*p2*p3)/60 + (7*p1*p4)/50 - (7*p3^2*xref)/30 - (49*p2*p4*xref)/100 + (119*p3*p4*xref^2)/100 - (119*p4^2*xref^3)/100
    np[5] = p2^2/12 + (3*p1*p3)/16 - (25*p2*p3*xref)/48 - (9*p1*p4*xref)/16 + (25*p3^2*xref^2)/48 + (17/16)*p2*p4*xref^2 - (7/4)*p3*p4*xref^3 + (21*p4^2*xref^4)/16
    np[4] = (5*p1*p2)/18 - (5*p2^2*xref)/18 - (5*p1*p3*xref)/9 + 
    (5/6)*p2*p3*xref^2 + (5/6)*p1*p4*xref^2 - (5*p3^2*xref^3)/9 - 
    (10/9)*p2*p4*xref^3 + (25/18)*p3*p4*xref^4 - (5*p4^2*xref^5)/6
    np[3] = (p1^2)/4 - (p1*p2*xref)/2 + (p2^2*xref^2)/4 + (1/2)*p1*p3*xref^2 - 
    (1/2)*p2*p3*xref^3 - (1/2)*p1*p4*xref^3 + (p3^2*xref^4)/4 + 
    (1/2)*p2*p4*xref^4 - (1/2)*p3*p4*xref^5 + (p4^2*xref^6)/4
    np[2] = 0.0
    np[1] = 0.0

    return np
end

"""
Calculates the integration of rJt from point xs[i] to xs[i+1] for all intervals in the Jt spline. Sign is for integration from xs[1] to xs[end].
"""
function rJt_integral_components(Jt::CubicSpline)
    rJt_components = zeros(Real,length(Jt.xs)-1)
    for i in 1:(length(Jt.xs)-1)
        idx_int_coeff = _int_rJt(Jt.ps[i,:],Jt.xs[i])
    
        idx_val1 = @evalpoly(Jt.xs[i], idx_int_coeff...)
        idx_val2 = @evalpoly(Jt.xs[i+1], idx_int_coeff...)

        rJt_components[i] = idx_val2-idx_val1
    end

    return rJt_components
end

"""
Calculates the integration of Jt/r from point xs[i] to xs[i+1] for nearly all intervals in the Jt spline (integration limited to i >= 2 due to divergence at x = 0.0). 
Sign is for integration from xs[2] to xs[end].
"""
function Jt_on_r_integral_components(Jt::CubicSpline)
    Jt_on_r_components = zeros(Real,length(Jt.xs)-2)
    for i in 2:(length(Jt.xs)-1)
        idx_int_coeff,log_coeff = _int_Jt_on_r(Jt.ps[i,:],Jt.xs[i])
    
        idx_val1 = @evalpoly(Jt.xs[i], idx_int_coeff...) + log_coeff*log(Jt.xs[i])
        idx_val2 = @evalpoly(Jt.xs[i+1], idx_int_coeff...) + log_coeff*log(Jt.xs[i+1])

        Jt_on_r_components[i-1] = idx_val2-idx_val1
    end

    return Jt_on_r_components
end

"""
Calculates the integration of I[Jtr]Jt/r from point xs[i] to xs[i+1] for all intervals in the Jt spline. Sign is for integration from xs[i] to xs[end].
"""     
function full_integral_components(Jt::CubicSpline)   
    full_components = zeros(Real,length(Jt.xs)-1)
    for i in 1:(length(Jt.xs)-1)
        idx_int_coeff = _int_full(Jt.ps[i,:],Jt.xs[i])
    
        idx_val1 = @evalpoly(Jt.xs[i], idx_int_coeff...)
        idx_val2 = @evalpoly(Jt.xs[i+1], idx_int_coeff...) 

        full_components[i] = idx_val2-idx_val1
    end

    return full_components
end

"""
Integration of rJt from knot xs[i] to knot xs[j]. Requires running of rJt_integral_components first.
Defaults to i = 1 (integration from inner part of spline). 
"""
function rJt_integral(rJt_components::AbstractArray{<:Real,1}, j::Int; i::Int=1)
    i==j && (return 0.0)
    # Making sure i, j in boundaries.
    if i <= 0
        throw("i too small ($i <= 0)")
    elseif i > length(rJt_components)
        throw("i too big ($i > $(length(rJt_components)))")
    end
    if j <= 0
        throw("j too small ($j <= 0)")
    elseif j > length(rJt_components)
        throw("j too big ($j > $(length(rJt_components)))")
    end

    if j>i
        return sum(rJt_components[i:(j-1)])
    else
        return -sum(rJt_components[j:(i-1)])
    end
end

"""
Integration of rJt from point x1 to point x2. Assume both x1 and x2 are in the range of interpolated data.
"""
function rJt_integral(Jt::CubicSpline, x2::Number, x1::Real)
    idx1,idx2 = _checkbounds(Jt, x1, x2)

    # Both points in same interval
    if idx1==idx2
        xref = Jt.xs[idx1]
        old_coeff = Jt.ps[idx1,:]

        int_coeff = _int_rJt(old_coeff,xref)

        x1_val = @evalpoly(x1, int_coeff...) #correct to a constant
        x2_val = @evalpoly(x2, int_coeff...) #correct to a constant

        return x2_val-x1_val #constant has no effect
    end

    # Assume idx2 > idx2
    signvar = 1
    # If not, relabel x1, idx1 and x2, idx2, and swap sign of final answer:
    if idx2<idx1 
        hold = idx2
        idx2 = idx1
        idx1 = hold

        hold2 = x2
        x2 = x1
        x1 = hold2

        signvar = -1
    end

    # Evaluating contribution from idx1 
    idx1_int_coeff = _int_rJt(Jt.ps[idx1,:],Jt.xs[idx1])

    idx1_rboundary = Jt.xs[idx1+1]

    idx1_val1 = @evalpoly(x1, idx1_int_coeff...)
    idx1_val2 = @evalpoly(idx1_rboundary, idx1_int_coeff...)

    idx1_int = idx1_val2 - idx1_val1

    # Evaluating contribution from idx2
    idx2_int_coeff = _int_rJt(Jt.ps[idx2,:],Jt.xs[idx2])

    idx2_lboundary = Jt.xs[idx2]

    idx2_val1 = @evalpoly(idx2_lboundary, idx2_int_coeff...) 
    idx2_val2 = @evalpoly(x2, idx2_int_coeff...)

    idx2_int = idx2_val2 - idx2_val1

    # Evaluating contributions from full-width integrated intervals
    full_widths = idx2-idx1-1 #number of full-width integration intervals

    runsum = 0.0
    for i in 1:full_widths
        idx_int_coeff = _int_rJt(Jt.ps[idx1+i,:],Jt.xs[idx1+i])
    
        idx_val1 = @evalpoly(Jt.xs[idx1+i], idx_int_coeff...)
        idx_val2 = @evalpoly(Jt.xs[idx1+i+1], idx_int_coeff...)

        runsum += (idx_val2-idx_val1)
    end

    return signvar*(idx1_int+idx2_int+runsum)
end

function _rJt_int_comp(Jt::CubicSpline,i::Real)
    idx_int_coeff = _int_rJt(Jt.ps[i,:],Jt.xs[i])

    return @evalpoly(Jt.xs[i], idx_int_coeff...)
end

"""
Calculates all standard terms used in the Bt spline internal integral - ie the integral of Bp*Jt.
"""     
function initialise_internalInt_Spline(Jt::CubicSpline)   
    num_intervals = length(Jt.xs)-1

    rJt_components = rJt_integral_components(Jt)
    Jt_on_r_components = Jt_on_r_integral_components(Jt)
    full_integral_comps = full_integral_components(Jt)   

    internalInt_components = zeros(Real,length(Jt.xs)-1)

    internalInt_components[1] = (-_rJt_int_comp(Jt,1))+full_integral_comps[1]
    for i=2:num_intervals
        internalInt_components[i] = Jt_on_r_components[i-1]*(rJt_integral(rJt_components, i) - _rJt_int_comp(Jt,i))+full_integral_comps[i]
    end

    return internalInt_components, rJt_components
end

"""
Calculates Bp from toroidal current Jt.
"""   
function Bp_Spline(Jt::CubicSpline; Bp_ref::Real=0.0, r_ref::Real=0.0)
    function Bp_at_r(r) 
        r==0.0 && (return 0.0)
        (1/r)*(mu0*rJt_integral(Jt, r, r_ref)+Bp_ref*r_ref)
    end
    return Bp_at_r
end

"""
Calculates the integral of Bp*Jt from r = 0.0 to r = x2. 
"""     
function internalInt_Spline(spline::CubicSpline, internalInt_components::AbstractArray{<:Real,1}, rJt_components::AbstractArray{<:Real,1}, x2::Number; x1::Real=0.0) 
    @assert x1 == 0.0

    if x1==x2   
        return 0.0
    end

    idx1,idx2 = _checkbounds(spline, x1, x2)

    # Both points in same interval ()
    if idx1==idx2 && idx1==1

        idx_int_coeff = _int_full(spline.ps[idx1,:],spline.xs[idx1])
    
        idx_val1 = @evalpoly(spline.xs[idx1], idx_int_coeff...)
        idx_val2 = @evalpoly(x2, idx_int_coeff...) 

        return mu0*(idx_val2-idx_val1)
    end


    # Assume idx2 > idx2
    signvar = 1

    # Evaluating contribution from intervals up to idx2
    int_val = 0.0
    int_val += sum(internalInt_components[1:(idx2-1)])

    # Evaluating contribution from idx2: log component
    idx_int_coeff,log_coeff = _int_Jt_on_r(spline.ps[idx2,:],spline.xs[idx2])
    
    idxlog_val1 = @evalpoly(spline.xs[idx2], idx_int_coeff...) + log_coeff*log(spline.xs[idx2])
    idxlog_val2 = @evalpoly(x2, idx_int_coeff...) + log_coeff*log(x2)

    int_val += (idxlog_val2-idxlog_val1)*(rJt_integral(rJt_components, idx2) - _rJt_int_comp(spline,idx2))

    # Evaluating contribution from idx2: full component
    idx_int_coeff = _int_full(spline.ps[idx2,:],spline.xs[idx2])
    
    idxfull_val1 = @evalpoly(spline.xs[idx2], idx_int_coeff...)
    idxfull_val2 = @evalpoly(x2, idx_int_coeff...) 

    int_val += (idxfull_val2-idxfull_val1)

    return int_val*mu0
end

"""
Calculates Bt from toroidal current Jt, pressure p and toroidal magnetic field on axis Bt0.
"""     
function Bt_Spline(Jt::CubicSpline,p::Union{Function,CubicSpline},Bt0::Real)
    internalInt_components,rJt_components = initialise_internalInt_Spline(Jt)   
    internalInt_Spline_Func = r -> internalInt_Spline(Jt, internalInt_components, rJt_components, r; x1=0.0)

    function Bt(r) 
        sqrt(Bt0^2 - 2*mu0*(p(r) - p(0.0)) - 2*mu0*internalInt_Spline_Func(r))
    end

    return Bt
end

"""
Calculates all equilibrium functions from an input Jt spline and input pressure function/spline. 
If an analytic pressure is used, input the analytical derivative (scaled by mu0) in the variable dpdr.
"""
function Spline_Equil(Jt::CubicSpline,p::Union{Function,CubicSpline},Bt0::Real,R0::Real; dpdr = nothing, plot_equil=true, print_mathematica_inputs=true, plotrvec = range(0.000001,Jt.xs[end],200), rvec=range(0.02,Jt.xs[end],step=0.02), kwargs...)
    Bp = Bp_Spline(Jt)
    Bt = Bt_Spline(Jt,p,Bt0)

    q = r -> r*Bt(r)/(R0*Bp(r))
    pm = r -> (1/(2*mu0))*(Bp(r)^2 + Bt(r)^2)
    local_beta(r) = r -> p(r)/pm(r)

    Jp = r -> -(1/(mu0))*ForwardDiff.derivative(Bt,r)
    if p isa CubicSpline
        function Jp2(r)
            return (-(1/(mu0))*(1/2)*(1/Bt(r))*(- 2*mu0*gradient(p,r) - 2*mu0*(Bp(r)*Jt(r))))
        end
    else
        Jp2 = nothing
    end

    if print_mathematica_inputs
        print_equil_data(Jt,Bt,p; rvec=rvec,kwargs...)
    end
    if plot_equil
        p1=plot(plotrvec,Bp.(plotrvec),title = "Bp in Teslas",xlabel="r (m)",ylabel="T",label=false)
        p2=plot(plotrvec,Bt.(plotrvec),title = "Bt in Teslas",label=false,ylims=(0.0,2*Bt(Jt.xs[end])),xlabel="r (m)",ylabel="T")
        p3=plot(plotrvec,q.(plotrvec),title = "q",label=false,xlabel="r (m)")
        p4=plot(plotrvec,local_beta.(plotrvec),title = "Local Plasma β",xlabel="r (m)",label=false)
        p isa CubicSpline && (p4=scatter!(p.xs,local_beta.(p.xs),title = "Local Plasma β",xlabel="r (m)",label=false))
        p5=plot(plotrvec,Jt.(plotrvec),title = "Toroidal current density",xlabel="r (m)",ylabel="Amps/m^2",label=false)
        p5=scatter!(Jt.xs,Jt.(Jt.xs),title = "Toroidal current",xlabel="r (m)",ylabel="Amps/m^2",label=false)
        p6=plot(plotrvec,Jp.(plotrvec),title = "Poloidal current",xlabel="r (m)",ylabel="Amps (?)",label=false)
        p isa CubicSpline && (p6=plot!(plotrvec,Jp2.(plotrvec),title = "Poloidal current",xlabel="r (m)",ylabel="Amps (?)",label=false))

        outerp6 = plot(p1,p2,p3,p4,p5,p6)
        display(outerp6)
    else
        outerp6=nothing
    end

    if p isa CubicSpline
        dpdr = r -> mu0*gradient(p,r)
    elseif dpdr isa Nothing
        throw("Need to manually define dpdr")
    end
    

    return Bp,Bt,q,dpdr,p,Jt,Jp,Jt.xs[end],outerp6,Jp2
end

diff_spln(spln) = r -> gradient(spln,r)
forwarddiff_spln(spln) = r -> ForwardDiff.derivative(r0 -> spln(r0), r)

function forwarddiff_spln(spln::CubicSpline,n)
    return forwarddiff_spln(r -> spln(r),n)
end
function forwarddiff_spln(spln::Function,n)
    if n>0
        return forwarddiff_spln(r -> ForwardDiff.derivative(r0 -> spln(r0),r), n-1)
    end
    
    return spln
end

function find_rs_Optim(q,m,n,rb;verbose=true)
    qtest = m/n
    f1 = x -> abs(q(first(x))-qtest)

    res1 = optimize(f1, 0.0, rb,  [rb/10])
    res2 = optimize(f1, 0.0, rb,  [rb/2])
    res3 = optimize(f1, 0.0, rb,  [9*rb/10])
    #optimize(f1,[q_Furth_find_rs(q0,ν,rs0)(qtest)],LBFGS())
    #res = optimize(ftest,xb/2,rb)

    if (abs(res1.minimizer[1]-res2.minimizer[1]) > 1e-10*rb) || (abs(res3.minimizer[1]-res2.minimizer[1]) > 1e-5*rb)
        print("Multiple rs values detected: \n") 
        print("rs = {$(res1.minimizer[1]), $(res2.minimizer[1]), $(res3.minimizer[1])}")
    end

    res = res2

    if Optim.converged(res) #&& !(abs(res.minimizer-a)<0.1*a || abs(res.minimizer)<0.1)
        rs = res.minimizer
        p1 = plot(0:(rb/300):rb,[q(i) for i in 0:(rb/300):rb],title="q profile",label=false)
        p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="q = $(m)/$(n)",line=:dash, xlabel="r (m)")
        p1 = vline!([res.minimizer],label="resonant surface location",line=:dash)
    else    
        verbose && @warn "Your chosen rational surface lies outside your minor radius!"
        p1 = plot(0:(rb/300):rb,[q(i) for i in 0:(rb/300):rb],title="q profile",label=false)
        p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="proposed rational q-value",line=:dash, xlabel="r (m)")
        p1 = vline!([res.minimizer],label="proposed resonant surface location",line=:dash)
    end

    verbose && display(p1)

    return res.minimizer[1] #, p1
end

function find_rs(q,m,n,rb;verbose=true, qtest=nothing, useOptim=false, return_plot=false)
    useOptim && (return find_rs_Optim(q,m,n,rb;verbose=verbose))
    
    if qtest isa Nothing
        qtest = m/n
    end

    f1 = x -> abs(q(first(x))-qtest)

    zeros = find_zeros(f1,0.0,rb)

    if length(zeros)>1
        verbose && print("Multiple rs values detected: \n") 
        verbose && print("rs = {$(zeros[1])")
        for i in 2:length(zeros)
            verbose && print(",$(zeros[i])")
        end
        print("}\n")

        return -2
    elseif length(zeros)==0
        verbose && print("No zeros detected.\n")
        verbose && @warn "Your chosen rational surface lies outside your minor radius!"

        #return find_rs_Optim(q,m,n,rb;verbose=verbose)
        return -1
    end

    rs = zeros[1]

    p1 = plot(0:(rb/300):rb,[q(i) for i in 0:(rb/300):rb],title="q profile",label=false)
    qtest isa Nothing && (p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="q = $(m)/$(n)",line=:dash, xlabel="r (m)"))
    !(qtest isa Nothing) && (p1 = plot!(0:(rb/300):rb,qtest.*ones(length(0:(rb/300):rb)),label="q = $(qtest)",line=:dash, xlabel="r (m)"))
    p1 = vline!([rs],label="resonant surface location",line=:dash)

    verbose && display(p1)

    return_plot && (return rs, p1)

    return rs
end
