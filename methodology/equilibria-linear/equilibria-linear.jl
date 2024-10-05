using GLMakie, Attractors, DrWatson

using OrdinaryDiffEq, LinearAlgebra

set_theme!(theme_latexfonts())
update_theme!(fontsize=16, 
    Axis = (
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible=false,
        ygridvisible=false,
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

# include("$(srcdir())/phase-portrait/fixed_point_analysis.jl")

#unstable: integrate forward
#stable: interate backwards

linear_stable_manifold(fp, dist_to_fp, stable_eigvec) = fp .+ (dist_to_fp .* stable_eigvec)

function stable_manifold(fp, eigvec, max_dist_to_fp_negative, max_dist_to_fp_positive, ds, T, Δt, ds_back; num_ics=5)
    dists_to_fp = filter!(x->x!=0.0, collect(range(-max_dist_to_fp_negative, max_dist_to_fp_positive, length=num_ics)))
    pts_lsm = map(dist->linear_stable_manifold.(fp, dist, eigvec), dists_to_fp)
    backward_trajs_all = [trajectory(ds_back, T, u0; Δt=Δt)[1] for u0 in pts_lsm]
    return backward_trajs_all
end 

function unstable_manifold(fp, eigvec, max_dist_to_fp_negative, max_dist_to_fp_positive, ds, T, Δt, ds_back; num_ics=5)
    dists_to_fp = filter!(x->x!=0.0, collect(range(-max_dist_to_fp_negative, max_dist_to_fp_positive, length=num_ics)))
    pts_lsm = map(dist->linear_stable_manifold.(fp, dist, eigvec), dists_to_fp)
    forward_trajs_all = [trajectory(ds, T, u0; Δt)[1] for u0 in pts_lsm]
    return forward_trajs_all 
end 

function plot_manifold(args...; kwargs...)
    f = Figure()
    a = Axis(f[1,1])
    plot_manifold!(f, a, args...; kwargs...)
end

function plot_manifold!(f, a, rt, trajs_all; idx1 = 1, idx2 = 2, color=:white, color_fp = :black, markersize_fp=8, kwargs...)
    scatter!(a, rt...; color=color_fp, markersize=markersize_fp)
    if !isnothing(trajs_all)
    for idx=1:length(trajs_all) lines!(a, trajs_all[idx][:, idx1], trajs_all[idx][:, idx2]; color, kwargs...) end
    end
    return f, a
end

#TODO: handle eigvals with im part; handle case of no stable manfiold
function invariant_manifolds(pvals; which_manifolds="stable")
    @unpack ds_func, p, min_bounds, max_bounds, numics, stable_manifold_params, unstable_manifold_params, ds_func_back = pvals
    u0 = SVector{2}([0.0, 1.0])
    integ_kwargs = get(pvals, :integ_kwargs, NamedTuple())
    
    ds = CoupledODEs(ds_func, u0, p; diffeq = integ_kwargs)
    ds_back = CoupledODEs(ds_func_back, u0, p; diffeq = integ_kwargs)
    rts_all, eigs_all, stability, eig_vecs_all = equilibria(ds, min_bounds, max_bounds, numics)

    rt = rts_all[1]
    λs = eigs_all[1]
    vs = eig_vecs_all[1]
    @info "Using root $rt, with λs $λs"
    
    if "stable" ∈ which_manifolds
        @unpack T, Δt, max_dist_to_fp_negative, max_dist_to_fp_positive = stable_manifold_params
        if any(imag.(λs) .!= 0)
            trajs_all_stable = nothing
        else
            idx_negative_eigvals = findall(x->x < 0.0, real.(λs))
            if !isnothing(idx_negative_eigvals)
                trajs_all_stable = []
                for (idx, idx_negative_eigval) in enumerate(idx_negative_eigvals)
                    stable_eigvec = real.(vs[:,idx_negative_eigval])
                    trajs_all_stable_ = stable_manifold(rt, stable_eigvec, max_dist_to_fp_negative, max_dist_to_fp_positive, ds, T, Δt, ds_back)
                    push!(trajs_all_stable, trajs_all_stable_...)
                end
            else
                trajs_all_stable = nothing
            end
        end
    else
        trajs_all_stable = nothing
    end 

    if "unstable" ∈ which_manifolds
        @unpack T, Δt, max_dist_to_fp_negative, max_dist_to_fp_positive = unstable_manifold_params
        if any(imag.(λs) .!= 0)
            trajs_all_unstable = nothing
        else
            idx_positive_eigvals = findall(x->x > 0.0, real.(λs))
            if !isnothing(idx_positive_eigvals)
                trajs_all_unstable = []
                for (idx, idx_positive_eigval) in enumerate(idx_positive_eigvals)
                    unstable_eigvec = real.(vs[:, idx_positive_eigval])
                    trajs_all_unstable_ = unstable_manifold(rt, unstable_eigvec, max_dist_to_fp_negative, max_dist_to_fp_positive, ds, T, Δt, ds_back)
                    push!(trajs_all_unstable, trajs_all_unstable_...)
                end
            else
                trajs_all_unstable = nothing
            end
        end

    else 
        trajs_all_unstable = nothing
    end 
    
    extra = @strdict ds
    return trajs_all_stable, trajs_all_unstable, rt, extra
end

using NLsolve, ForwardDiff

function get_equilibria(pvals)
    @unpack N = pvals
    min_bounds = N == 2 ? (-110.0, -110, -0.5, -0.5) : (-110.0, -0.5)
    max_bounds = N == 2 ? (10., 10, 1.5, 1.5) : (10.0, 1.5)
    numics = 50
    ds = get_ds(pvals)
    # rts, eigs, stability = equilibria(ds, min_bounds, max_bounds, numics; ftol=1e-15, method=:newton, autodiff=:forward)
    rts, eigs, stability = equilibria(ds, min_bounds, max_bounds, numics; ftol=1e-8, method=:trust_region, iterations=100_000, num_digits_roots=5)
    # rts, eigs, stability = equilibria(ds, min_bounds, max_bounds, numics; ftol=1e-5, method=:anderson, beta=0.1, xtol=1e-5)
    check_convergence_roots(rts, ds)
    return rts, eigs, stability
end

function equilibria(ds::DynamicalSystem, args...; kwargs...)
    system_function = dynamic_rule(ds)
    p = ds.integ.sol.prob.p
    equilibria(system_function, p, args...; kwargs...)
end


"""
Currently for iip, but it works also for oop with a few small changes!
For a system function of size N
    * `min_bounds`: NTuple of size N, with the minimum values for the search box along each dimension
    * `min_bounds`: NTuple of size N, with the maximum values for the search box along each dimension
    * `numics`: number of randomly generated ics in the box
"""
function equilibria(system_function::Function, params, min_bounds, max_bounds, numics; kwargs...)
    sys_function(du, u) = system_function(du, u, params, 0.0)
    sys_jacobian(du, u) = ForwardDiff.jacobian(sys_function, du, u)
    rts = find_roots(sys_function, min_bounds, max_bounds, numics; kwargs...)
    eigs, stability, eig_vecs = stability_roots(rts, sys_jacobian)
    return rts, eigs, stability, eig_vecs
end


# include("$(srcdir())/dynamical-systems-analysis/attractors/sampler.jl")
function find_roots(system_function, min_bounds, max_bounds, numics; kwargs...)
    sampler, = statespace_sampler_og(Random.MersenneTwister(1); min_bounds, max_bounds)
    u0s = [sampler() for i=1:numics]
    find_roots(system_function, u0s; kwargs...)
end

all_converged(res) = res.x_converged && res.f_converged

function find_roots(system_function, u0s::Vector{Vector{A}}; num_digits_roots=10, kwargs...) where {A}
    roots = Vector{Vector{A}}()
    for (idx, u0) in enumerate(u0s)
        res = nlsolve(system_function, u0; kwargs...)
        rt = res.zero
        rt_rounded = round.(rt; digits=num_digits_roots)
        # if !all_converged(res) continue end
        if !res.f_converged || res.residual_norm >= 1e-5 continue end
        # @show res.x_converged, res.f_converged, res.residual_norm
        
        # if ! (rt_rounded ∈ round.(roots; digits=num_digits_roots)) push!(roots, rt) end
        roots_rounded = map(x->round.(x; digits=num_digits_roots), roots)
        if idx >= 2 if !(rt_rounded ∈ roots_rounded) push!(roots, rt) end end
    end
    return roots
end

function stability_roots(rts, jacobian; ds_type="continuous")
    eigs = Vector{Vector{Complex{Float64}}}(undef, length(rts))
    eig_vecs = Vector{Matrix{Complex{Float64}}}(undef, length(rts))
    
    for (idx, root) in enumerate(rts)
        du = zeros(Float64, length(root))
        J = jacobian(du, root)
        eigs[idx] = LinearAlgebra.eigvals(J)
        eig_vecs[idx] = LinearAlgebra.eigvecs(J)
    end
    
    stable = Bool[isstable_continuous(e) for e in eigs]
    
    return eigs, stable, eig_vecs
end

isstable_continuous(e) = maximum(real(x) for x in e) < 0

function check_convergence_roots(roots, ds)
    system_function = dynamic_rule(ds)
    p = ds.integ.sol.prob.p
    for (idx, root) in enumerate(roots)
        iteration = similar(root)
        system_function(iteration, root, p, 0.0)
        # @show iteration, root
    end
    nothing 
end

export statespace_sampler
using Distributions, LinearAlgebra, Random

"""
    statespace_sampler(rng = Random.GLOBAL_RNG; kwargs...) → sampler, isinside
Convenience function that creates two functions. `sampler` is a 0-argument function
that generates random points inside a state space region defined by the keywords.
`isinside` is a 1-argument function that decides returns `true` if the given 
state space point is inside that region.

The regions can be:
* **Rectangular box**, with edges `min_bounds` and `max_bounds`.
  The sampling of the points inside the box is decided by the keyword `method` which can
  be either `"uniform"` or `"multgauss"`.
* **Sphere**, of `spheredims` dimensions, radius `radius` and centered on `center`.
"""
function statespace_sampler_og(rng = Random.GLOBAL_RNG; 
        min_bounds=[], max_bounds=[], method="uniform", 
        radius::Number=-1,
        spheredims::Int=0, center=zeros(spheredims),
    )

    if min_bounds ≠ [] && max_bounds != []
        if method == "uniform" gen, isinside = boxregion(min_bounds, max_bounds, rng)
        elseif method == "multgauss" gen, isinside = boxregion_multgauss(min_bounds, max_bounds, rng)
        else @error("Unsupported boxregion sampling method")
        end
    elseif radius ≥ 0 && spheredims ≥ 1
        gen, isinside = sphereregion(radius, spheredims, center, rng)
    else
        @error("Incorrect keyword specification.")
    end
    return gen, isinside
end


function boxregion_multgauss(as, bs, rng)
    @assert length(as) == length(bs) > 0
    center = mean(hcat(as,bs), dims=2)
    gen() = [rand(rng, truncated(Normal(center[i]), as[i], bs[i])) for i=1:length(as)]
    isinside(x) = all(as .< x .< bs)
    return gen, isinside
end

# this has a docstring only because it was part of the expansionentropy api.
# It won't be exported in future versions
"""
    boxregion(as, bs, rng = Random.GLOBAL_RNG) -> sampler, isinside

Define a box in ``\\mathbb{R}^d`` with edges the `as` and `bs` and then
return two functions: `sampler`, which generates a random initial condition in that box
and `isinside` that returns `true` if a given state is in the box.
"""
function boxregion(as, bs, rng = Random.GLOBAL_RNG)
    @assert length(as) == length(bs) > 0
    gen() = [rand(rng)*(bs[i]-as[i]) + as[i] for i in 1:length(as)]
    isinside(x) = all(as .< x .< bs)
    return gen, isinside
end

# Specialized 1-d version
function boxregion(a::Real, b::Real, rng = Random.GLOBAL_RNG)
    a, b = extrema((a, b))
    gen() = rand(rng)*(b-a) + a
    isinside = x -> a < x < b
    return gen, isinside
end

#Algorithm is taken from https://math.stackexchange.com/questions/1585975/how-to-generate-random-points-on-a-sphere.
#It follows from the fact that a multivariate normal distribution is spherically symmetric.
function sphereregion(r, dim, center, rng)
    @assert r ≥ 0 
    gen() = normalize([( 2*randn(rng) - 1 ) for j=1:dim]) .* r .+ center
    isinside(x) = norm(x .- center) < r
    return gen, isinside
end

        

function linear!(du, u, p, t)
    A = p[1]
    # mul!(du, A, u)
    du .= A*u
    nothing
end

function linear_back!(du, u, p, t)
    A = p[1]
    # mul!(du, -A, u)
    du .= -A*u
    nothing
end

pvals_linear = Dict(
    :ds_func => linear!,
    :ds_func_back => linear_back!,
    :p => [[-1 0; 0 -1]],
    :min_bounds => [-0.5, -0.5],
    :max_bounds => [0.5, 0.5],
    :min_bounds_all => [-70, -0.05], :max_bounds_all => [-10, 1.0], #for all
    :numics => 100,
    :T => 20.0 ,
    :Δt => 0.01,
    :max_dist_to_fp => 0.001,
    :integ_kwargs => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
    # :stable_manifold_params => Dict( :max_dist_to_fp_negative => 0.00001, :max_dist_to_fp_positive => 0.1, :T => 5.0 , :Δt => 0.01),
    # :stable_manifold_params => Dict( :max_dist_to_fp_negative => 0.00001, :max_dist_to_fp_positive => 0.1, :T => 10.0 , :Δt => 0.01),
    :stable_manifold_params => Dict( :max_dist_to_fp_negative => 0.1, :max_dist_to_fp_positive => 0.1, :T => 10.0 , :Δt => 0.01),
    :unstable_manifold_params => Dict(:max_dist_to_fp_negative => 0.001, :max_dist_to_fp_positive => 0.001, :T => 15.0 , :Δt => 0.01),
    :xg => range(-3, 3; length = 400),
    :yg => range(-3, 3; length = 400),
)

pvals_linear_all = Dict(
    :ds_func => linear!,
    :ds_func_back => linear_back!,
    :p => [[[-1 0; 0 -1]], [[-1+1im 0; 0 -1-1im]], [[+1 0; 0 +1]], [[+1+1im 0; 0 +1-1im]], [[-1 0; 0 +1]]],
    :min_bounds => [[-0.5, -0.5]],
    :max_bounds => [[0.5, 0.5]],
    :min_bounds_all => [[-70, -0.05]], :max_bounds_all => [[-10, 1.0]], #for all
    :numics => 100,
    :T => 20.0 ,
    :Δt => 0.01,
    :max_dist_to_fp => 0.001,
    :integ_kwargs => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
    # :stable_manifold_params => Dict( :max_dist_to_fp_negative => 0.00001, :max_dist_to_fp_positive => 0.1, :T => 5.0 , :Δt => 0.01),
    # :stable_manifold_params => Dict( :max_dist_to_fp_negative => 0.00001, :max_dist_to_fp_positive => 0.1, :T => 10.0 , :Δt => 0.01),
    :stable_manifold_params => Dict( :max_dist_to_fp_negative => 0.1, :max_dist_to_fp_positive => 0.1, :T => 10.0 , :Δt => 0.01),
    :unstable_manifold_params => Dict(:max_dist_to_fp_negative => 0.001, :max_dist_to_fp_positive => 0.001, :T => 15.0 , :Δt => 0.01),
    :xg => range(-3, 3; length = 400),
    :yg => range(-3, 3; length = 400),
)

# pvals = pvals_linear
# stable_trajs_all, unstable_trajs_all, fp, extra = invariant_manifolds(pvals; which_manifolds = ["stable", "unstable"])
# 
# 
# fig = Figure()
# ax = Axis(fig[1,1])
# plot_manifold!(fig, ax, fp, stable_trajs_all; color=:green, markersize_fp=15)
# plot_manifold!(fig, ax, fp, unstable_trajs_all; color=:red, markersize_fp=15)
# xlims!(ax, -2.5, 2.5)
# ylims!(ax, -2.5, 2.5)
# elems = Any[LineElement(color = c, linestyle = :solid) for c in [:green, :red] ]
# push!(elems, MarkerElement(color=:black, marker=:circle))
# elems_str = ["Stable manifold", "Unstable manifold", "Saddle-point"]
# Legend(fig[1, 1], elems, elems_str, patchsize = (35, 35), rowgap = 10, tellwidth=false, tellheight=false, halign=:left, valign=:top)
# ax.xlabel = "x"
# ax.ylabel = "y"
# 
# makiesave("$(plotsdir())/thesis/duffing-bistability.png", fig; px_per_unit=4)

size_fig = (600, 400)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)

fig = Figure(; size=size_fig)
As = [[[-1.5 0; 0 -1]], [[-1 1; -1 -1]], [[-2 0; 0 +1]], [[+1.5 0; 0 +1]], [[+1 1; -1 +1]], ]
idxs_axis = [[1,1], [2,1], [1,2], [1,3], [2,3]]
axs = []
# ics = [[1.0, 1.0], [1.0, -1], [-1, 1], [-1, -1]]
ics = [
        [[2.5, 2.5], [2.5, -2.5], [-2.5, 2.5], [-2.5, -2.5]],
        [[2.5, 2.5], [2.5, -2.5], [-2.5, 2.5], [-2.5, -2.5]],
        # [[2.0, 2.0], [2.0, -2], [-2, 2], [-2, -2]],
        [[2.5, 0.3], [2.5, -0.3], [-2.5, 0.3], [-2.5, -0.3]],
        [[0.02, 0.1], [0.02, -0.1], [-0.02, 0.1], [-0.02, -0.1]],
        [[0.002, 0.002], [0.002, -0.002], [-0.002, 0.002], [-0.002, -0.002]],
    ]
for (idx, A) in enumerate(As)
    pvals = deepcopy(pvals_linear)
    pvals[:p] = A
    @info "idx = $idx, pvals = $pvals"
    
    stable_trajs_all, unstable_trajs_all, fp, extra = invariant_manifolds(pvals; which_manifolds = ["stable", "unstable"])
    ax = Axis(fig[idxs_axis[idx]...]); push!(axs, ax)
    
    plot_manifold!(fig, ax, fp, stable_trajs_all; color=:green, markersize_fp=15)
    plot_manifold!(fig, ax, fp, unstable_trajs_all; color=:red, markersize_fp=15)
    
    for ic in ics[idx]
        @unpack ds = extra 
        tr, ts = trajectory(ds, 10.0, ic)
        lines!(ax, tr[:,1], tr[:,2], color=:black)
        # scatter!(ax, tr[1:10:end-2][:,1], tr[1:10:end-2][:,2], color=:black, marker=:utriangle)
    end
    
    xlims!(ax, -2.5, 2.5)
    ylims!(ax, -2.5, 2.5)
    # ax.xlabel = "x"
    # ax.ylabel = "y"
end

for ax in [axs[1], axs[2], axs[5]]
    ax.xlabel = "x"
    ax.ylabel = "y"
end
axs[1].title = "r₋ = 2 \n stable"
axs[2].title = "c₋ = 2 \n stable"
axs[3].title = "r₊ = 1, r₋ = 1 \n unstable"
axs[4].title = "r₊ = 2 \n unstable"
axs[5].title = "c₊ = 2 \n unstable"

elems = Any[LineElement(color = c, linestyle = :solid) for c in [:green, :red, :black] ]
push!(elems, MarkerElement(color=:black, marker=:circle))
elems_str = ["Stable manifold", "Unstable manifold", "Typical trajectories", "Equilibrium"]
Legend(fig[2,2], elems, elems_str,  tellwidth=false, tellheight=false, halign=:center, valign=:center, labelsize=14)


save("hyperbolic-eq-2d.png", fig; px_per_unit=4)