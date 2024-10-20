using DrWatson
using Revise
using OrdinaryDiffEq, JLD2, Attractors, PreallocationTools
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

pvals_kuramoto_twistedstates = dict_list(Dict(
    :N => 80, 
    :ϵ => 1.0, 
    :topm => "wattsstrogatz", :k => 2, 
    :rewiring_prob => 0.0, 
    :topseed => 1,
    :freq_mean => 0.0, :freq_std => 0.0, :freq_seed => 1,
    :icmin => 0, :icmax => 2π, :icseed => 1, :ictype => "uniform",
    :samples_per_parameter => 100_000, 
    :diffeq => (alg = Tsit5(), reltol = 1e-8, abstol = 1e-8, maxiters = Inf),
    :Δt => 0.1, :T => 50.0, :Ttr => 100.0,
    :ss_xg => range(0, +2π, length=10), #for basins in featurizing 
    :ss_grid => Derived([:ss_xg, :N], (ss_xg, N)->ntuple(x->ss_xg, N)),
    :ic_grid => Derived([:ss_xg, :N], (ss_xg, N)->ntuple(x->ss_xg, N)),
    # :featurizer => kuramoto_featurizer_comms, #TODO DOES THIS WORK??
    :featurizer => kuramoto_featurizer_twisted_states,
    :grouping_kwargs => (threshold = 0.01, rescale_features = false), 
    :frequencies_description => "identical",
))[1]; 


pvals = pvals_kuramoto_twistedstates


@unpack N, ϵ, diffeq,k,rewiring_prob, topseed = pvals
g = watts_strogatz(N, k, rewiring_prob; is_directed=false, rng=MersenneTwister(topseed))
adjl = g.fadjlist; #out neighbors == receivers of each node
Icoup = PreallocationTools.dualcache(zeros(N));
ωs = zeros(Float64, N)
u0s = zeros(Float64, N)
ps = params_kuramoto(ωs, ϵ, adjl, Icoup)
prob = ODEProblem(kuramoto_network!, u0s, (0.0, 2000.0), ps)
ds = CoupledODEs(prob, diffeq)

@unpack grouping_kwargs, T, Δt, Ttr, ic_grid, samples_per_parameter = pvals
grouping_config = GroupViaPairwiseComparison(; grouping_kwargs...)
mapper = AttractorsViaFeaturizing(ds, pvals[:featurizer], grouping_config; Δt, T, Ttr)

rng = MersenneTwister(1)
sampler, = statespace_sampler_og(rng; min_bounds = minimum.(ic_grid), max_bounds = maximum.(ic_grid))
ics = Dataset([sampler() for _ in 1:samples_per_parameter])
# fs, labels = basins_fractions(mapper, ics)
# attractors = extract_attractors(mapper)

# res = @strdict fs attractors 
filename = "twisted-states-distribution-fractions-k_$k-$samples_per_parameter.jld2"
# jldsave(filename; res)
res = load(filename)["res"]
@unpack fs, attractors = res

qs = Dict(k=>kuramoto_featurizer_winding_number(att, nothing) for (k, att) in attractors)
fs_with_qs = Dict(q=>fs[k] for (k,q) in qs)

size_fig = (400, 400)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)

xs = collect(keys(fs_with_qs))
ys = collect(values(fs_with_qs))
idxperm = sortperm(xs)
xs = xs[idxperm]
ys = ys[idxperm]

fig = Figure(; size=size_fig)
ax = Axis(fig[1,1])
scatterlines!(ax, xs, ys, color=:black)
ax.xticks = -7:1:7
ax.yticks = 0:0.05:0.25
ax.ylabel = "fs(q)"
ax.xlabel = "q"
save("twisted-states-distribution-fractions-k_$k.png", fig; px_per_unit)

# _adjl = results["details"]["mapper"].ds.p0.adjl
# 
# prange = 4.5:0.01:6
# pidx = 1 # index of the parameter
# ascm = AttractorSeedContinueMatch(mapper)
# fractions_cont, attractors_cont = global_continuation( ascm, prange, pidx, sampler; samples_per_parameter = 1_000)