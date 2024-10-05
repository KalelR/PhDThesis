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

function supertitle(fig, title)
    Label(fig[0, :], title, valign = :bottom,
        # padding = (0, 0, 5, 0), tellheight = true, tellwidth = false,
        tellheight = true, tellwidth = false,
        font = "TeX Gyre Heros Bold", # same font as Axis titles
    )
    return
end


activationfunc(V, Vh, k) = 1 / (1 + exp((Vh - V)/k))

@inbounds function inapk!(du, u, p, t)
    I = p
    El = -80; gl = 8; gna = 20; gk = 10; Vhm = -20; km =15; Vhn = -25; kn = 5; Ena = 60; Ek = -90; C = 1; N = 1; 
    τ = 0.16
    Vs = u[1:N]
    ns = u[N+1:2N]
    
    minf = activationfunc.(Vs, Vhm, km)
    ninf = activationfunc.(Vs, Vhn, kn)
    
    dV1 = (I - gl*(Vs[1]- El) - gna*minf[1]*(Vs[1] - Ena) - gk * ns[1] * (Vs[1] - Ek) ) / C
    dn1 = (ninf[1] - ns[1]) / τ 
    du[1] = dV1 
    du[2] = dn1 
    
    return nothing
end

@inbounds function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

function kuramoto_two!(du, u, p, t)
    ϵ, ω2, ω1 = p 
    du[1] = ω1 + ϵ*sin(u[2]-u[1])
    du[2] = ω2 + ϵ*sin(u[1]-u[2])
    nothing 
end

function vanderpol!(du,u,p,t)
    # μ = 0.1; β = 0.0; α = 1.0; g = 1/2; wf=2 
    μ = 0.1; β = 0.0; α = 1.0; g = 1/2; wf=sqrt(3)
    x, v = u
    du[1] = v 
    du[2] = μ*(1-x^2)*v - α*x - β*x^3 + g*cos(wf*t)
    nothing 
end
    

pvals_eq = Dict(
    :ds_func => inapk!, 
    :p => 2.0,
    :T => 10, :Δt => 0.001, :Ttr => 0,
    :diffeq => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
    :u0 => [-20, -0.1],
)

pvals_lc = Dict(
    :ds_func => inapk!, 
    :p => 6.0,
    :T => 10, :Δt => 0.001, :Ttr => 0,
    :diffeq => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
    :u0 => [-20, -0.1],
)

# pvals_torus = Dict(
#     :ds_func => kuramoto_two!, 
#     :p => [0.1, 1.0, 2.0],
#     :T => 50, :Δt => 0.001, :Ttr => 0,
#     :diffeq => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
#     :u0 => [0.0, 0.1],
# )

pvals_torus = Dict(
    :ds_func => vanderpol!, 
    :p => [],
    # :T => 100, :Δt => 0.001, :Ttr => 100,
    :T => 100, :Δt => 0.001, :Ttr => 0,
    :diffeq => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
    :u0 => [2.0,2],
)

pvals_chaos = Dict(
    :ds_func => lorenz!, 
    :p => [10, 28, 8/3],
    :T => 100, :Δt => 0.001, :Ttr => 0,
    :diffeq => (abstol=1e-10, reltol=1e-10, alg=Vern9()),
    :u0 => [0.1, 0.1, 0.1],
    :azimuth => 5.41,
    :elevation =>0.182,
)


size_fig = (1000, 500)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)

pvals_all = [pvals_eq, pvals_lc, pvals_torus, pvals_chaos]
fig = Figure(size=size_fig) 
axs = []
for (idx_sys, pvals) in enumerate(pvals_all)
# idx_sys = 4
# pvals = pvals_all[idx_sys]
@unpack ds_func, u0, p, T, diffeq, Ttr, Δt = pvals
ds = CoupledODEs(ds_func, u0, p; diffeq)

tr, ts = trajectory(ds, T, u0; Ttr, Δt)

if dimension(ds) == 2
    ax = Axis(fig[1,idx_sys]); push!(axs, ax)
    lines!(ax, tr[:,1], tr[:,2]; color=:black)
    if pvals == pvals_eq 
        scatter!(ax, tr[end]...; markersize=15, color=:black)
    end
else 
    @unpack elevation, azimuth = pvals
    ax = Axis3(fig[1,idx_sys]; azimuth, elevation); push!(axs, ax)
    lines!(ax, tr[:,1], tr[:,2], tr[:,3]; color=:black)
end
ax = Axis(fig[2,idx_sys]); push!(axs, ax)
lines!(ax, ts, tr[:,1]; color=:black)
end

axs[1].title = "Equilibrium"
axs[3].title = "Limit cycle"
axs[5].title = "Quasiperiodic"
axs[7].title = "Chaotic"

axs[1].ylabel = "y"
axs[1].xlabel = "x"
axs[3].xlabel = "x"
axs[5].xlabel = "x"
axs[2].ylabel = "x(t)"
for ax in [axs[2], axs[4], axs[6], axs[8]]
    ax.xlabel = "t"
end

axs[1].xticks = -60:20:-0
axs[3].xticks = -60:20:-0
axs[2].yticks = -60:20:-0
axs[4].yticks = -60:20:-0

supertitle(fig, "Basic types of attractors in nonlinear systems")

save("attractors-nonlinear.png", fig; px_per_unit=4)