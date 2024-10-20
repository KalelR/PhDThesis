using DrWatson
using Revise
using OrdinaryDiffEq, JLD2, Attractors, PreallocationTools
using GLMakie, JLD2

"""
    label_axes!(axs::Array{Axis};
        valign = :top, halign = :right, pad = 5, kwargs...
    )

Add labels (like a,b,c,...) to all axes.
Keywords customly adjust location, and `kwargs` are propagated to `Label`.
"""
function label_axes!(axs;
        labels = range('A'; step = 1, length = length(axs)),
        # transformation = x -> "("*string(x)*")",
        transformation = x -> ""*string(x)*" ",
        valign = :top, halign = :right,
        pad = 5, box_color="white", kwargs...,
    )

    lbs = @. string(transformation(labels))
    # Create padding from alignment options
    padding = [0,0,0,0]
    if halign == :right
        padding[2] = pad
    elseif halign == :left
        padding[1] = pad
    end
    if valign == :top
        padding[3] = pad
    elseif valign == :bottom
        padding[4] = pad
    end

    for (i, ax) in enumerate(axs)
        # @assert ax isa Axis
        gc = ax.layoutobservables.gridcontent[]
        x = gc.parent[gc.span.rows, gc.span.cols]
        # Currently `Label` has no way of having a box around it
        lab = Label(x, lbs[i];
            tellwidth=false, tellheight=false,
            valign, halign, padding, font = :bold, justification = :center,
            kwargs...
        )
        # but we can access the internals and get the box of the label,
        # and then make an actual box around it
        bx = Box(first(axs).parent; bbox = lab.layoutobservables.computedbbox, color = box_color)
        Makie.translate!(bx.blockscene, 0, 0, +1)
        Makie.translate!(lab.blockscene, 0, 0, +2)
    end
    return
end


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

fig = Figure(); axs =[]

# -------------------------------- PANELS A-B -------------------------------- #

using Colors, ColorSchemes
"""
Function to return a vector of colors to be used as the argument for color in Makie. Uses either COlors for distinguishable_colors or ColorSchemes for colorschemes
"""
function get_colors(number_colors, colormap=:viridis; rangemin=0, rangemax=1)
    if colormap == :distinguishable_colors
        cols = distinguishable_colors(number_colors, [RGB(1,1,1), RGB(0,0,0)], dropseed=true) #vector of RGB structures 
        return cols
    else
        colorscheme = colorschemes[colormap] 
        return get(colorscheme, range(rangemin, rangemax, length=number_colors))
    end
end

function complex_order_parameter_t(v_phase)
    θ  = 0 + 0im
    for i in v_phase
        θ += exp(im * i)
    end
    θ /= length(v_phase)
    return abs(θ), angle(θ)
end

R = 1.0
θs = range(0, 2pi, length=1000)
xs = R .* cos.(θs)
ys = R .* sin.(θs)

θs_units_all = [range(pi/1.5-pi/1.5, pi/1.5+pi/1.5, length=20), range(pi/2-pi/6, pi/2+pi/6, length=20)]
titles = ["weak phase sync", "strong phase sync"]

for idx=1:2
    θs_units = θs_units_all[idx]
    ax = Axis(fig[1,idx]; aspect=1); push!(axs, ax);

    xs_units = R .* cos.(θs_units)
    ys_units = R .* sin.(θs_units)
    colors_units = get_colors(length(θs_units), :distinguishable_colors)

    lines!(ax, xs, ys; color=:black)
    scatter!(ax, xs_units, ys_units; color=colors_units, markersize=15)

    r, ϕ =  complex_order_parameter_t(θs_units)
    arrow_head = [r*cos(ϕ), r*sin(ϕ)]

    arrows!(ax, [0], [0], [arrow_head[1]], [arrow_head[2]])#; linewidth=3, arrowsize=2.0)
    hidedecorations!(ax)
    hidespines!(ax)    
    ax.title = titles[idx]
end


# -------------------- Panel C: Og kuramoto with loernzian ------------------- #
ax = Axis(fig[2,1]); push!(axs, ax);
ϵc = 1.596
ϵs1 = range(0, ϵc, length=10)
Rs1 = zeros(length(ϵs1))

ϵs2 = range(ϵc, 5,length=100)
Rs2 = @. sqrt(1 - ϵc/ϵs2)

ϵs = [ϵs1; ϵs2]
Rs = [Rs1; Rs2]
lines!(ax, ϵs, Rs; color=:black)
ax.ylabel = L"r"
ax.xlabel = L"K"
vlines!(ax, ϵc; linestyle=:dash)
ax.xticks = ([0, 1, ϵc, 2, 3, 4, 5], [L"0", L"1", L"\epsilon_c", L"2", L"3", L"4", L"5"])
ax.yticks = 0:0.2:1
ylims!(ax, 0, 1)
ax.title = "Lorentzian distribution (theoretical)"


# -------------- Panel D: OG kuramoto with gaussian, simulations ------------- #
ax = Axis(fig[2,2]); push!(axs, ax);
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
    # :featurizer => kuramoto_featurizer_comms,
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

# ic = rand(Float64, N)
# attractors_cont = []
# for (idx, p) in enumerate(prange)
#     set_parameter!(ds, pidx, p)
#     tr, ts = trajectory(ds, T, ic; Ttr, Δt)
#     push!(attractors_cont, tr)
# end

# ts = Ttr:Δt:Ttr+T
# Rmean_all = [mean(order_parameter(Matrix(atts)')) for atts in attractors_cont]

# res = @strdict Rmean_all prange attractors_cont
# filename = "kuramoto-og-continuation-singleic-eps_$prange.jld2"
# jldsave(filename; res)
filename = "kuramoto-og-continuation-singleic-eps_$prange.jld2"
res = load(filename)["res"]
@unpack Rmean_all, prange = res

scatterlines!(ax, prange, Rmean_all; color=:black)
ax.ylabel = L"\langle r(t) \rangle"
ax.xlabel = L"K"
ϵc = 1.596
vlines!(ax, ϵc; linestyle=:dash)
ax.xticks = ([0, 1, ϵc, 2, 3], [L"0", L"1", L"\epsilon_c", L"2", L"3"])
ax.yticks = 0:0.2:1
ylims!(ax, 0, 1)

ax.title = "Gaussian distribution (N=$N)"


# ------------ panel E: twisted state - distribution of basin sies ----------- #
ax = Axis(fig[3,1]); push!(axs, ax);
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
@unpack k, samples_per_parameter = pvals

filename = "twisted-states-distribution-fractions-k_$k-$samples_per_parameter.jld2"
res = load(filename)["res"]
@unpack fs, attractors = res

qs = Dict(k=>kuramoto_featurizer_winding_number(att, nothing) for (k, att) in attractors)
fs_with_qs = Dict(q=>fs[k] for (k,q) in qs)

xs = collect(keys(fs_with_qs))
ys = collect(values(fs_with_qs))
idxperm = sortperm(xs)
xs = xs[idxperm]
ys = ys[idxperm]

scatterlines!(ax, xs, ys, color=:black)
ax.xticks = -7:2:7
ax.yticks = 0:0.05:0.25
ax.ylabel = "fs(q)"
ax.xlabel = "q"
ax.title = "Homogeneous frequencies, k = $k"

# ---------------------------- panel F -- fs(q=0) ---------------------------- #
ax = Axis(fig[3,2]); push!(axs, ax)
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
param_vals = pvals_kuramoto_twistedstates_cont[:k]
filename = "twisted-states-continuation-k_$param_vals-$(pvals[:samples_per_parameter]).jld2"
res = load(filename)["res"]
@unpack fs_all_q0, param_vals = res

scatterlines!(ax, param_vals, fs_all_q0; color=:black)
ax.ylabel = "fs(q=0)"
ax.xlabel = "k"
ax.yticks = 0:0.1:1.0 
ax.xticks = param_vals[1]:2:param_vals[end]
hlines!(ax, 1; color=:red, linestyle=:dash)
ax.title = "Homogeneous frequencies, q = 0"

# ------------------------------- beuatify fig! ------------------------------ #
size_fig = (675, 600)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)
resize!(fig, size_fig...)


label_axes!(axs; halign=:left)
filename = "kuramoto-thesis-figure.png"
save(filename, fig; px_per_unit)

fig