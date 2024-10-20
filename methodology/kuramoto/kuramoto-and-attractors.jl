
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



circ_dist(x,y) =  angle(exp(1*im*x)./exp(1*im*y))
function circ_dist(vec)
    diffs = zeros(Float64, length(vec)-1)
    for idx = 1:length(vec)-1
        diffs[idx] = circ_dist(vec[idx+1], vec[idx])
    end
    return diffs 
end

function kuramoto_featurizer_twisted_states(A, t)
    Δθ = mean( circ_dist(A[end]) ) #mb do circ dist
    return Δθ
end

winding_number(Δθ, N) = round(Int, (N*Δθ)/(2π))

function kuramoto_featurizer_winding_number(A, t)
    # Δθ = mean( mod.( abs.( diff(A[end])) , 2pi) ) #mb do circ dist
    Δθ = mean( circ_dist(A[end]) ) #mb do circ dist
    q = winding_number(Δθ, length(A[end]))
    # @show Δθ,q,length(A[end])
    return q
end


function run_featurizing(pvals)
    prob, ps = get_problem_kuramoto(pvals)
    @unpack diffeq = pvals
    @show diffeq
    ds = CoupledODEs(prob, diffeq)
    fs, labels, attractors, details = run_featurizing(ds, pvals)
    res_featurizing = @strdict fs labels attractors details
end


function run_featurizing(ds, pvals)
    
    @unpack N, featurizer, ttrans, tend, samples_per_parameter, Δt, ic_grid, grouping_kwargs = pvals
    T = tend - ttrans; Ttr = ttrans;
    
    GroupingType = haskey(grouping_kwargs, :min_neighbors) ? GroupViaClustering : GroupViaPairwiseComparison
    @info "Finding att with grouping $GroupingType"
    grouping_config = GroupingType(; grouping_kwargs...)
    mapper = AttractorsViaFeaturizing(ds, featurizer, grouping_config; Δt, T, Ttr)
    
    rng = MersenneTwister(1)
    sampler, = statespace_sampler_og(rng; min_bounds = minimum.(ic_grid), max_bounds = maximum.(ic_grid))
    ics = Dataset([sampler() for _ in 1:samples_per_parameter])
    fs, labels = basins_fractions(mapper, ics)
    attractors = extract_attractors(mapper)
     
    details = @strdict mapper ics
    return fs, labels, attractors, details
end


function get_problem_kuramoto(d_params; _u0s=nothing, _tend=nothing, adjl=nothing, ωs=nothing)
    @unpack N, tend, freq_mean, freq_std, freq_seed, ϵ = d_params
    tend = isnothing(_tend) ? tend : _tend
    ωs = isnothing(ωs) ? get_frequencies_gauss(N, freq_mean, freq_std, freq_seed, 0) : ωs
    frequencies_description = get(d_params, :frequencies_description, "nonzero_mean")
    if frequencies_description == "zero_mean" ωs .-= mean(ωs) end
    if frequencies_description == "identical" ωs .= 0 end
    @show mean(ωs)
    adjl = isnothing(adjl) ? get_adjl(d_params) : adjl
    u0s = isnothing(_u0s) ? get_initial_conditions_uniform(d_params) : _u0s
    Icoup = PreallocationTools.dualcache(zeros(N));

    noise_amplitude = get(d_params, :noise_amplitude, nothing)
    # @show noise_amplitude
    if isnothing(noise_amplitude)
        idxs_decrease = get(d_params, :idxs_decrease, nothing)
        if isnothing(idxs_decrease)
            ps = params_kuramoto(ωs, ϵ, adjl, Icoup)
            prob = ODEProblem(kuramoto_network!, u0s, (0.0, tend), ps)
        else
            @unpack idxs_increase, weight = d_params
            @info "Running adj mat mode, with weight=$weight, idxs_decrease = $idxs_decrease, idxs_increase = $idxs_increase."
            Δθ = PreallocationTools.dualcache(zeros(N));
            A = adjlist_to_adjmat(adjl)
            ps = params_kuramoto_adjacency_weight_change(ωs, ϵ, A, Icoup, Δθ, idxs_decrease, idxs_increase, weight)
            prob = ODEProblem(kuramoto_network!, u0s, (0.0, tend), ps)
        end
    else
        ps = params_kuramoto_noise(ωs, ϵ, adjl, Icoup, noise_amplitude)
        # W = WienerProcess(0.0, zeros(N))
        # prob = SDEProblem(kuramoto_network!, kuramoto_noise!, u0s, (0.0, tend), ps, noise=W)
        prob = SDEProblem(kuramoto_network!, kuramoto_noise!, u0s, (0.0, tend), ps)
    end

    return prob, ps
end

mutable struct params_kuramoto
    ωs :: Array{Float64, 1}
    ϵ :: Float64
    adjl :: Array{Array{Int64,1}, 1}
    Icoup :: PreallocationTools.DiffCache{Vector{Float64}, Vector{Float64}}  #cache
end

function coupling_kuramoto!(u, p, t)
    @unpack adjl, ϵ = p; Icoup = get_tmp(p.Icoup, u)
    fill!(Icoup, 0.0)

    for i in eachindex(u)
        @inbounds inedges = adjl[i] #get neighborhood of i
        for j in inedges
            @inbounds Icoup[i] += sin(u[j] - u[i])
        end
    end
    Icoup .*= ϵ 
    return Icoup
end

function kuramoto_network!(du, u, p,t)
    coupling_kuramoto!(u, p, t)
    Icoup = get_tmp(p.Icoup, u)
    @unpack ωs = p
    @. du = ωs + Icoup
end

using Distributions, Graphs, Random, DelimitedFiles

function get_adjl(pvals)
    @unpack topm = pvals
    if topm isa Vector 
        adjl = topm
        return adjl
    else
        if topm == "ER"
            @unpack N, k, topseed = pvals
            filename = "$(datadir())/sims/inputs/N_$N/graph-$topm-N_$N-k_$k-seed_$topseed.jld2"
            if isfile(filename)
                g = load(filename)["g"]
                adjl = g.fadjlist;
            else
                @info("File $filename was not found. Generating it and saving now.")
                Ne = k*N
                g = erdos_renyi(N, Ne; is_directed=true, rng=MersenneTwister(topseed))
                safesave(filename, Dict("g"=>g))
                adjl = g.fadjlist; #out neighbors == receivers of each node
                return adjl
            end
        elseif topm == "121"
            adjl = [[2,3], [], []]
        elseif topm == "bidirectional"
            adjl = [[2], [1]]
        elseif topm == "wattsstrogatz"
            @unpack N, k, rewiring_prob, topseed = pvals
            filename = "$(datadir())/sims/inputs/N_$N/graph-$topm-N_$N-k_$k-rewiringprob_$rewiring_prob-seed_$topseed.jld2"
            if isfile(filename)
                g = load(filename)["g"]
                adjl = g.fadjlist;
            else
                g = watts_strogatz(N, k, rewiring_prob; is_directed=false, rng=MersenneTwister(topseed))
                safesave(filename, Dict("g"=>g))
                adjl = g.fadjlist; #out neighbors == receivers of each node
                return adjl
            end
        elseif topm == "wattsstrogatz_consecutive"
            @unpack N, k, topseed = pvals
            rewiring_prob = get(pvals, :rewiring_prob, nothing)
            if !isnothing(rewiring_prob) && !(rewiring_prob isa Int) 
                @unpack rewiring_prob = pvals
                filename = "$(datadir())/sims/inputs/N_$N/graph-$topm-N_$N-k_$k-rewiringprob_$rewiring_prob-seed_$topseed.jld2"
                number_rewired_conns = round(Int,N*k*rewiring_prob) 
            else
                @unpack number_rewired_conns = pvals
                filename = "$(datadir())/sims/inputs/N_$N/graph-$topm-N_$N-k_$k-numrewirings_$number_rewired_conns-seed_$topseed.jld2"
            end
            @info "Generating or loading $topm with N=$N, k = $k, with $number_rewired_conns rewired connections."
            if isfile(filename)
                g = load(filename)["g"]
                adjl = g.fadjlist;
            else
                g = watts_strogatz_consecutive(N, k, number_rewired_conns; is_directed=false, rng=MersenneTwister(topseed))
                safesave(filename, Dict("g"=>g))
                adjl = g.fadjlist; #out neighbors == receivers of each node
                return adjl
            end
        elseif topm[end-3:end] == ".dat"
            @unpack N = pvals
            adjv = readdlm(topm, Int64)[:,1]
            adjl = adjvec_to_adjlist(adjv, N)
            return adjl
        else
            error("No other topm found")
        end
    end
end

