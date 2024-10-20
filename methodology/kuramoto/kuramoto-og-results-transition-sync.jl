using DrWatson
using Revise
using OrdinaryDiffEq, JLD2, Attractors, PreallocationTools
using GLMakie, JLD2
# using CairoMakie, JLD2

include("kuramoto-and-attractors.jl")

using Distributions, LinearAlgebra, Random


set_theme!(theme_latexfonts())
update_theme!(fontsize=16, 
    Axis = (
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible=true,
        ygridvisible=true,
        spinewidth=1.5,
        xtickwidth=2.0,
        ytickwidth=2.0,
        xtickalign=1,
        ytickalign=1,
        xticklabelsize=18,
        yticklabelsize=18,
        xlabelsize=20,
        ylabelsize=20,
    ))


"""
Calculates R(t), given matrix of phases of Nxtimes
"""
function order_parameter(θs::AbstractMatrix{T}) where T
    NN, TT = size(θs)
    if NN > TT @warn "Size of phases maybe wrong as N = $NN > T = $TT" end
    mapslices(order_parameter, θs, dims=1)[1,:]
end

function order_parameter(θs::Vector{T}) where T
    if isempty(θs) return [-1] end
    θ_sum  = zero(ComplexF64)
    for θ in θs
        θ_sum += exp(im * θ)
    end
    return abs(θ_sum) / length(θs)
end

function stats_q(vec)
    return [mean(vec), std(vec), maximum(vec), minimum(vec), quantile(vec, 0.25), quantile(vec, 0.75), quantile(vec, 0.5)]
end

"""
Calculates Ri(t) and then the stats_q of Ri for i being the index of each community. Defines community as a chunk of 1:N.
"""
function order_parameter_community(θs, number_communities)
    N, T = size(θs)
    idxs_in_communities = chunk_array(1:N, number_communities)
    R_stats_communities = Dict{Int, Vector{Float64}}()
    for idx_community in eachindex(idxs_in_communities)
        Rs_community = mapslices(order_parameter, θs[idxs_in_communities[idx_community], :], dims=1 )[1,:]
        R_stats_communities[idx_community] = stats_q(Rs_community)
    end
    return idxs_in_communities, R_stats_communities
end

function kuramoto_featurizer_comms(A, ts)
    # Amat = Matrix(A)'[:, 1:10:end] #transforms into NXT matrix 
    Amat = Matrix(A)'[:, 1:1:end] #transforms into NXT matrix changedthisasofjun14!
    
    Rs = order_parameter(Amat)
    R_mean = mean(Rs)
    R_std = std(Rs)
    
    number_communities = 4
    _, R_stats_communities = order_parameter_community(Amat, number_communities)
    R_mean_communities = [R_stats[1] for (k, R_stats) in R_stats_communities]
    
    # Ωs = speeds(ts, A, params)
    Ωs = mapslices(diff, Matrix(A), dims=1) ./ (ts[2] - ts[1])
    spatial_deviation_frequencies = std(Ωs, dims=2)[:,1] #T els
    
    return [R_mean, R_std, R_mean_communities..., mean(spatial_deviation_frequencies)]
end


"""
Util function for calculating order parameter dividnig the network into local communities.
"""
function chunk_array(array, num_chunks)
    size_chunks = Int64(floor(length(array) / num_chunks)) #size of all (except, potentially, last chunk)
    num_regular_chunks = Int64(floor(length(array) / size_chunks))
    array_div = []
    for idx=1:num_regular_chunks-1
        idx_min = Int64((idx-1)*size_chunks+1 )
        idx_max = Int64(idx*size_chunks)
        push!(array_div, array[idx_min:idx_max])
    end

    push!(array_div, array[(num_regular_chunks-1)*size_chunks+1:end])
end

const DATADIR = "/dss/work/lufa9108/multistability-review/data"
# const DATADIR = datadir()
DrWatson.datadir() = DATADIR

pvals_kuramoto_og = dict_list(Dict(
    # :N => 500, 
    :N => 1000, 
    :ϵ => 1.0, 
    :topm => "wattsstrogatz", :k => 40, 
    :rewiring_prob => 0.0, 
    :topseed => 1,
    :freq_mean => 0.0, :freq_std => 1.0, :freq_seed => 1,
    :icmin => 0, :icmax => 2π, :icseed => 1, :ictype => "uniform",
    :samples_per_parameter => 30, 
    :diffeq => (alg = Tsit5(), reltol = 1e-8, abstol = 1e-8, maxiters = Inf),
    :Δt => 0.1, :T => 1000.0, :Ttr => 1000.0,
    :ss_xg => range(0, +2π, length=10), #for basins in featurizing 
    :ss_grid => Derived([:ss_xg, :N], (ss_xg, N)->ntuple(x->ss_xg, N)),
    :ic_grid => Derived([:ss_xg, :N], (ss_xg, N)->ntuple(x->ss_xg, N)),
    :featurizer => kuramoto_featurizer_comms,
    :grouping_kwargs => (threshold = 0.1, rescale_features = false), 
    :frequencies_description => "identical",
))[1]; 

mutable struct params_global
    ωs :: Array{Float64, 1}
    ϵ :: Float64
end

function kuramotonetwork_global!(du, u, p, t)
    @unpack ωs, ϵ = p
    N = length(ωs)
    z = mean(exp.(im .* u))
    r = abs(z)
    θ = angle(z)
    @inbounds for i in 1:N
        du[i] = ωs[i] + ϵ*r*sin(θ - u[i])
    end
    return
end


pvals = pvals_kuramoto_og
@unpack N, ϵ, diffeq = pvals
Icoup = PreallocationTools.dualcache(zeros(N));
ωs = randn(MersenneTwister(1), N)
ps = params_global(ωs, ϵ)
u0s = zeros(Float64, N)
prob = ODEProblem(kuramotonetwork_global!, u0s, (0.0, 2000.0), ps)
ds = CoupledODEs(prob, diffeq)

@unpack grouping_kwargs, T, Δt, Ttr, ic_grid, samples_per_parameter = pvals
prange = collect(range(0, 3.0; step=0.1)); push!(prange, [1.55, 1.65, 1.75]...); sort!(prange)
pidx = :ϵ # index of the parameter

# grouping_config = GroupViaPairwiseComparison(; grouping_kwargs...)
# mapper = AttractorsViaFeaturizing(ds, pvals[:featurizer], grouping_config; Δt, T, Ttr)
# 
# rng = MersenneTwister(1)
# sampler, = statespace_sampler_og(rng; min_bounds = minimum.(ic_grid), max_bounds = maximum.(ic_grid))
# ics = Dataset([sampler() for _ in 1:samples_per_parameter])
# 
# ascm = AttractorSeedContinueMatch(mapper)
# fractions_cont, attractors_cont = global_continuation(ascm, prange, pidx, sampler; samples_per_parameter = 100)
# Rmean_all = [mean(order_parameter(Matrix(atts[1])')) for atts in attractors_cont]


ic = rand(Float64, N)
attractors_cont = []
for (idx, p) in enumerate(prange)
    set_parameter!(ds, pidx, p)
    tr, ts = trajectory(ds, T, ic; Ttr, Δt)
    push!(attractors_cont, tr)
end

ts = Ttr:Δt:Ttr+T
# pvals[:featurizer](attractors_cont[1][1], ts)
Rmean_all = [mean(order_parameter(Matrix(atts)')) for atts in attractors_cont]

# res = @strdict Rmean_all prange fractions_cont attractors_cont
# res = @strdict Rmean_all prange 
filename = "kuramoto-og-continuation-singleic-eps_$prange.jld2"
jldsave(filename; res)
res = load(filename)["res"]
@unpack Rmean_all, prange = res
res = @strdict Rmean_all prange 
jldsave(filename; res)



size_fig = (400, 350)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)

fig = Figure(; size=size_fig)
ax = Axis(fig[1,1])
scatterlines!(ax, prange, Rmean_all; color=:black)
ax.ylabel = L"\langle R \rangle"
ax.xlabel = L"K"
ϵc = 1.596
vlines!(ax, ϵc; linestyle=:dash)
ax.xticks = ([0, 1, ϵc, 2, 3], [L"0", L"1", L"\epsilon_c", L"2", L"3"])
ax.yticks = 0:0.2:1
ylims!(ax, 0, 1)
filename = "kuramoto-og-continuation-singleic.png"
save(filename, fig; px_per_unit)
# 
# idx = 15 #chaotic looking
# idx = 20
# idx = 25
# att = attractors_cont[idx]
# Rs = order_parameter(Matrix(att)')
# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, ts, Rs)
# 
# heatmap(mod.(Matrix(att), 2pi))
 

# ϵs1 = range(0, ϵc, length=10)
# Rs1 = zeros(length(ϵs1))
# 
# ϵs2 = range(ϵc, 5,length=100)
# Rs2 = @. sqrt(1 - ϵc/ϵs2)
# 
# ϵs = [ϵs1; ϵs2]
# Rs = [Rs1; Rs2]
# fig = Figure()
# ax = Axis(fig[1,1])
# lines!(ax, ϵs, Rs; color=:black)
# ax.ylabel = L"\langle R \rangle"
# ax.xlabel = L"\epsilon"
# ϵc = 1.596
# vlines!(ax, ϵc; linestyle=:dash)
# ax.xticks = ([0, 1, ϵc, 2, 3, 4, 5], [L"0", L"1", L"\epsilon_c", L"2", L"3", L"4", L"5"])
# ax.yticks = 0:0.2:1
# ylims!(ax, 0, 1)