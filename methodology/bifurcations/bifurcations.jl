using GLMakie, Attractors, DrWatson

using OrdinaryDiffEq, LinearAlgebra

include("../vis-code.jl")


size_fig = (600, 800)
pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)
fig = Figure(size=size_fig); axs=[]
# -------------------------------- Saddle-node ------------------------------- #
f_sn(x, α) = x^2 + α
fp_sn(α) = α <= 0 ? sqrt(-α) : NaN
xs = range(-1.2, 1.2, step=0.01)
αs = [-0.2, 0.0, 0.1]
idx_bif = 1
for (idx_col, α) in enumerate(αs)
    ys = f_sn.(xs, α)
    ax = Axis(fig[idx_bif, idx_col]); push!(axs, ax)
    lines!(xs, ys; color=:black)
    x_fp = fp_sn(α)
    if α != 0
        scatter!(ax, x_fp, 0; color=:red)
        scatter!(ax, -x_fp, 0; color=:green)
    else 
        scatter!(ax, 0,0;color=:gray)
    end
    hlines!(ax, 0; color=:black, linestyle=:dash)
    ylims!(-0.4, 0.4)
    ax.title = "α = $α"
end
ax = Axis(fig[idx_bif, 4]); push!(axs, ax)
αs = range(-1, 0, step=0.001)
xfp_all = fp_sn.(αs)
lines!(ax, αs, xfp_all; color=:green)
lines!(ax, αs, -xfp_all; color=:red)

for ax_idx in 1:3 axs[ax_idx].xlabel = "x"; axs[ax_idx].xticks=[-1.0, 0.0, 1.0] end 
axs[1].ylabel = L"\dot{x}"
axs[4].ylabel = L"x^\star" 
axs[4].xlabel = L"\alpha"


# ----------------------------------- Hopf ----------------------------------- #
fp_hb(α) = α <= 0 ? sqrt(-α) : NaN
xs = range(-1.2, 1.2, step=0.01)
αs = [-0.2, 0.0, 0.2]
idx_bif = 2
for (idx_col, α) in enumerate(αs)
    ax = Axis(fig[idx_bif, idx_col]); push!(axs, ax)
    x_fp = fp_sn(α)
    color_origin = α <= 0 ? :green : :red
    scatter!(0, 0; color=color_origin)
    if α > 0 
        radius = √α 
        θs = range(0, 2pi, step=0.01)
        xs = radius * cos.(θs)
        ys = radius * sin.(θs)
        lines!(xs, ys; color=:green)
    end 
    ylims!(-0.5, 0.5)
    xlims!(-0.5, 0.5)
    ax.title = "α = $α"
end
ax = Axis(fig[idx_bif, 4]); push!(axs, ax)
αs = range(-0.2, 0.2, step=0.0001)
radius_all = map(α->√α, αs[αs .> 0]) 
color_origin = map(α->(α <= 0 ? :green : :red), αs)
lines!(ax, αs, zeros(length(αs)); color=color_origin)
lines!(ax, αs[αs .> 0], radius_all; color=:green)
lines!(ax, αs[αs .> 0], -radius_all; color=:green)

for ax_idx in 5:7 axs[ax_idx].xlabel = "x" end 
axs[1+4].ylabel = L"y"
axs[4+4].ylabel = L"x_\mathrm{max}, x_\mathrm{min}" 
axs[4+4].xlabel = L"\alpha"
axs[4+4].xticks = [-0.2, 0.0, 0.2]


# ------------------------------------ HOM ----------------------------------- #
function sandstede!(du, u, p, t)
    α = p 
    x,y=u
    du[1] = -x + 2y + x^2 
    du[2] = (2-α)*x -y -3*x^2 + (3/2)*x*y
    nothing 
end


α = -0.2
T = 100.0; Ttr = 30.0; Δt=0.01
diffeq = (abstol=1e-10, reltol=1e-10, alg=Vern9())
αs = [-0.01, 0.00001, 0.1]
idx_bif = 3
for (idx_col, α) in enumerate(αs)
    u0 = [0.8, 0.1]
    ax = Axis(fig[idx_bif, idx_col]); push!(axs, ax)
    ds = CoupledODEs(sandstede!, u0, α; diffeq)
    tr, ts = trajectory(ds, T; Ttr, Δt)
    lines!(ax, tr[:,1], tr[:,2], color=:green)
    scatter!(ax, 0, 0; color=:red)
    xlims!(ax, -0.05, 1.1)
    ylims!(ax, -0.4, 0.4)
    ax.title = "α = $α"
end


#Please George forgive me for I have sinned
ax = Axis(fig[idx_bif, 4]); push!(axs, ax)
αs = range(0.0001, 0.1, step=0.001)
max_x = zeros(Float64, length(αs))
min_x = zeros(Float64, length(αs))
for (idx, α) in enumerate(αs)
    u0 = [0.8, 0.1]
    ds = CoupledODEs(sandstede!, u0, α; diffeq)
    tr, ts = trajectory(ds, T; Ttr, Δt)
    max_x[idx] = maximum(tr[:,1])
    min_x[idx] = minimum(tr[:,1])
end 
lines!(ax, αs, max_x; color=:green)
lines!(ax, αs, min_x; color=:green)
αs = range(-0.1, 0.1, step=0.001)
lines!(ax, αs, zeros(length(αs)); color=:red)

for ax_idx in 9:11 axs[ax_idx].xlabel = "x" end 
axs[1+8].ylabel = L"y"
axs[4+8].ylabel = L"x_\mathrm{max}, x_\mathrm{min}" 
axs[4+8].xlabel = L"\alpha"
axs[4+8].xticks = [-0.1, 0.0, 0.1]
# α = -0.2
# ds = CoupledODEs(sandstede!, [1.0, 1.0], α; diffeq)
# tr, ts = trajectory(ds, T; Ttr, Δt)
    
# supertitle(fig, "Basic types of attractors in nonlinear systems")
label_axes!(axs; halign=:left)


Label(fig[1, 0], "Saddle-\n node", valign = :center, tellwidth=false, tellheight=false)
Label(fig[2, 0], "Hopf",        valign = :center, tellwidth=false, tellheight=false)
Label(fig[3, 0], "Homo-\n clinic",  valign = :center, tellwidth=false, tellheight=false)

colsize!(fig.layout, 0, Relative(0.1))

size_fig = (600, 800)
resize!(fig, 800, 800)
save("bifurcations.png", fig; px_per_unit=4)

fig