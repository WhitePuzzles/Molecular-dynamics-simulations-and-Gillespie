module Mol_Gillespie

using Plots
using LinearAlgebra
using Random

include("constant.jl")
include("force.jl")
include("step_update.jl")
include("min_energy.jl")
include("mix.jl")
include("make_gif.jl")
include("gillespie.jl")



end
