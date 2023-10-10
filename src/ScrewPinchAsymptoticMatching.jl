__precompile__()

module ScrewPinchAsymptoticMatching

using TaylorSeries
using MatrixEquations
using ForwardDiff
using Distances
using Plots
using Statistics
using NaNStatistics
using LinearAlgebra
using RootsAndPoles
using SpecialFunctions
using JLD2
using CSV
using DataFrames
#using CubicSplines

using QuadGK
using DelimitedFiles
using HypergeometricFunctions_Mod

using Interpolations
using StaticArrays
using DifferentialEquations
using OrdinaryDiffEq
using Optim
using Roots
using Random
using Distributions
using LaTeXStrings

const mu0 = 4*pi*1e-7
const H2_mass = 2*1.6605e-27
global mod_path = "" #once ScrewPinchAsymptoticMatching compiled run:# ScrewPinchAsymptoticMatching.mod_path=chop(pathof(ScrewPinchAsymptoticMatching);tail=31)

include("CubicSplines_Mod.jl")
export CubicSpline

include("InnerAsymptotics.jl")
export findXmax, checkXmax, generate_Umatrix

include("InnerIntegrator.jl")
export alpha3_on_alpha4

include("InnerOuterMatching.jl")
export contour_inner, generate_inner_ratios, generate_inner_ratios_Xvar, generate_D_Q 
export countour_Q, GRPF_Q, D_Q_Xvarying, countour_Q_Xvarying, cylinder_root_start_Glass75

include("ChandraScrewPinchEquilibrium.jl")
export Furth_Equil, Chandra_Equil, Scaffidi_Equil, print_equil_data, q_, local_beta, pm

include("ModularCurrentProfile.jl")
export integral, rJt_integral, initialise_internalInt_Spline, Bp_Spline, internalInt_Spline, Bt_Spline, Spline_Equil, diff_spln, forwarddiff_spln, find_rs

include("OuterIntegrator.jl")
export raw_delta_prime, Psi_w_scales, plot_full_Psis, k_, g_, find_rs

include("OuterAsymptotics.jl")
export plot_full_Psis_w_frobenius, Δl_Δr_calculator, Δl_Δr_calculator_zeroPressure, Δl_Δr_xmin, Δl_Δr_xmin_visualiser, ξ_plus_minus, f_

include("IdealStability.jl")
export test_ideal_stability, test_Suydam, ideal_stability_solin, ideal_stability_out, ideal_nowall_external_stability

#include("DataBaseGen.jl")

include("GeneratingInnerTerms.jl")
export generateInnerTerms_ScrewPinch_Scaffidi, generateInnerTerms_ScrewPinch_Glass75

include("deltac_interface.jl")
export call_deltac

end