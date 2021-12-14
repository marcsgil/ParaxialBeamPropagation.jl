module ParaxialBeamPropagation

using FFTW, Interpolations,
ClassicalOrthogonalPolynomials,SpecialFunctions, Plots, 
BenchmarkTools, PlutoUI,Printf,LaTeXStrings

export get_fourier_transform,display_fourier_transform,
propagate_beam,display_beam,animate_beam,HG,LG,airy_beam,vertical_obstacle,horizontal_obstacle,
circular_obstacle,rectangular_obstacle

include("Propagation.jl")
include("InitialProfiles.jl")

end
