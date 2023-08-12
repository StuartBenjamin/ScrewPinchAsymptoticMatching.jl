function call_deltac(e::Float64, f::Float64, g::Float64, h::Float64, k::Float64, m::Float64, taua::Float64, taur::Float64, v1::Float64, q_real::Float64, q_im::Float64; 
deltac_tol_in::Float64=1e-14, pfac_in::Float64=0.1 ,inps_xfac_in::Float64=1.0,
nx_in::Int=256, nq_in::Int=5, cutoff_in::Int=5, kmax_in::Int=8, interp_np_in::Int=10,
hermite_in_res_cells::Bool=false,noexp_in::Bool=true,restore_uh_in::Bool=true,
restore_us_in::Bool=true,restore_ul_in::Bool=true,diagnose_res_in::Bool=false,
diagnose_params_in::Bool=false,grid_diagnose_in::Bool=false,
deltac_bin_sol_in::Bool=false,deltac_out_sol_in::Bool=false)

    deltac_run_julia_out = zeros(Float64, 3, 2, 2)

    x = ccall((:__deltac_mod_MOD_deltac_run_julia, string(ScrewPinchAsymptoticMatching_path,"/src/deltac.so")),
                Cvoid,
                (Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
                Ref{Float64},Ref{Float64},Ref{Float64},
                Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
                Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
                Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
                Ptr{Array{Float64,3}}),
                e,f,g,h,k,m,taua,taur,v1,q_real,q_im,
                deltac_tol_in,pfac_in,inps_xfac_in,
                Int32(nx_in),Int32(nq_in),Int32(cutoff_in),Int32(kmax_in),Int32(interp_np_in),
                Int32(hermite_in_res_cells),Int32(noexp_in),Int32(restore_uh_in),
                Int32(restore_us_in),Int32(restore_ul_in),Int32(diagnose_res_in),
                Int32(diagnose_params_in),Int32(grid_diagnose_in),
                Int32(deltac_bin_sol_in),Int32(deltac_out_sol_in),
                deltac_run_julia_out)

    false && display(deltac_run_julia_out)

    Δo = deltac_run_julia_out[1,1,1]+im*deltac_run_julia_out[1,1,2]
    Δe = deltac_run_julia_out[1,2,1]+im*deltac_run_julia_out[1,2,2]
    return Δe,Δo
end

#mod_path=pathof(ScrewPinchAsymptoticMatching)
#call_deltac(0.1, 0., 0., 0., 0., 0., 1., 1., 1.,0.01,0.0)
#call_deltac(0.1, 0.005, 40., 0.05, 500., 0., 1., 1., 1.,3.0,0.0)