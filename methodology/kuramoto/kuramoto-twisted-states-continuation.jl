using DrWatson
using Revise
using OrdinaryDiffEq, JLD2, Attractors, PreallocationTools
using CairoMakie
using GLMakie

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

const DATADIR = "/dss/work/lufa9108/multistability-review/data"
# const DATADIR = datadir()
DrWatson.datadir() = DATADIR

pvals_kuramoto_twistedstates_cont = Dict(
    :N => 80, 
    :ϵ => 1.0, 
    :topm => "wattsstrogatz", :k => collect(range(1, 25; step=2)), 
    # :topm => "wattsstrogatz", :k => collect(range(1, 25; step=10)), 
    :rewiring_prob => 0.0, 
    :topseed => 1,
    :freq_mean => 0.0, :freq_std => 1.0, :freq_seed => 1,
    :icmin => 0, :icmax => 2π, :icseed => 1, :ictype => "uniform",
    :samples_per_parameter => 4000, 
    # :samples_per_parameter => 100, 
    # :samples_per_parameter => 20000,
    :diffeq => (alg = Tsit5(), reltol = 1e-8, abstol = 1e-8, maxiters = Inf),
    :Δt => 0.1, :T => 50.0, :Ttr => 500.0,
    :ss_xg => range(0, +2π, length=10), #for basins in featurizing 
    :ss_grid => Derived([:ss_xg, :N], (ss_xg, N)->ntuple(x->ss_xg, N)),
    :ic_grid => Derived([:ss_xg, :N], (ss_xg, N)->ntuple(x->ss_xg, N)),
    :featurizer => kuramoto_featurizer_twisted_states,
    :grouping_kwargs => (threshold = 0.01, rescale_features = false), 
    :frequencies_description => "identical",
); 
pvals = pvals_kuramoto_twistedstates_cont

atts_all = Vector{Dict}(undef, length(dict_list(pvals_kuramoto_twistedstates_cont)))
fs_all = Vector{Dict}(undef, length(dict_list(pvals_kuramoto_twistedstates_cont)))
# for (idx, pvals) in enumerate(dict_list(pvals_kuramoto_twistedstates_cont))
#     @unpack N, ϵ, diffeq,k,rewiring_prob, topseed = pvals
#     g = watts_strogatz(N, 2*k, rewiring_prob; is_directed=false, rng=MersenneTwister(topseed))
#     adjl = g.fadjlist; #out neighbors == receivers of each node
#     Icoup = PreallocationTools.dualcache(zeros(N));
#     ωs = zeros(Float64, N)
#     u0s = zeros(Float64, N)
#     ps = params_kuramoto(ωs, ϵ, adjl, Icoup)
#     prob = ODEProblem(kuramoto_network!, u0s, (0.0, 2000.0), ps)
#     ds = CoupledODEs(prob, diffeq)
# 
#     @unpack grouping_kwargs, T, Δt, Ttr, ic_grid, samples_per_parameter = pvals
#     grouping_config = GroupViaPairwiseComparison(; grouping_kwargs...)
#     mapper = AttractorsViaFeaturizing(ds, pvals[:featurizer], grouping_config; Δt, T, Ttr)
# 
#     rng = MersenneTwister(1)
#     sampler, = statespace_sampler_og(rng; min_bounds = minimum.(ic_grid), max_bounds = maximum.(ic_grid))
#     ics = Dataset([sampler() for _ in 1:samples_per_parameter])
#     fs, labels = basins_fractions(mapper, ics; show_progress=false)
#     attractors = extract_attractors(mapper)
#     atts_all[idx] = attractors
#     fs_all[idx] = fs
# end
# 
# fs_all_q0 = zeros(Float64, length(fs_all))
# for idx = 1:length(fs_all)
#     qs = Dict(kuramoto_featurizer_winding_number(att, nothing)=>k for (k, att) in atts_all[idx])
#     fs_all_q0[idx] = fs_all[idx][qs[0]]
# end
# 
param_vals = pvals_kuramoto_twistedstates_cont[:k]
filename = "twisted-states-continuation-k_$param_vals-$(pvals[:samples_per_parameter]).jld2"
# res = @strdict fs_all_q0 param_vals
# jldsave(filename; res)


res = load(filename)["res"]
@unpack fs_all_q0, param_vals = res


size_fig = (400, 300)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)

fig = Figure(; size=size_fig)
ax = Axis(fig[1,1])
# params_normalized = param_vals ./ pvals[:N]
scatterlines!(ax, param_vals, fs_all_q0; color=:black)
ax.ylabel = "fs(q=0)"
ax.xlabel = "k"
ax.yticks = 0:0.1:1.0 
# ax.xticks = 1:2:param_vals[end]
ax.xticks = param_vals
hlines!(ax, 1; color=:red, linestyle=:dash)
filename = "twisted-states-continuation.png"
save(filename, fig; px_per_unit)
