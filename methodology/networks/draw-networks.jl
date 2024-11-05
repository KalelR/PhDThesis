using GLMakie
include("../vis-code.jl")

function adjmat_to_adjlist(A)
    N = size(A,1)
    adjlist = [Int64[] for i=1:N]
    for i=1:N 
        for j=1:N
            if A[i,j] == 1 
                push!(adjlist[i], j)
            end
        end
    end
    return adjlist
end

function plot_network(adjlist)
    fig = Figure()
    ax = Axis(fig[1,1]; aspect=1)
    plot_network!(ax, adjlist)
end

function plot_network!(ax, adjlist)
    N = length(adjlist)
    positions_units = plot_nodes!(ax, N)
    plot_connections!(ax, adjlist, positions_units)
    hidedecorations!(ax)
    hidespines!(ax)
    
    return fig, ax 
end

function plot_nodes!(ax, N; R=1, colors=:black, markersize=20)
    θs = range(0, 2π, length=N+1)[1:N]
    xs = @. R*cos(θs)
    ys = @. R*sin(θs)
    positions_units = [[xs[i], ys[i]] for i=1:N]
    scatter!(ax, xs, ys; color=colors, markersize)
    return positions_units
end

function plot_connections!(ax, adjlist, positions_units)
    for (i, neighbors) in enumerate(adjlist)
        pos_unit_i = positions_units[i]
        for j in neighbors 
            pos_unit_j = positions_units[j]
            draw_connection!(ax, pos_unit_i, pos_unit_j)
        end
    end
end

function draw_connection!(ax, pos_i, pos_j; color=:black)
    x_i, y_i = pos_i 
    x_j, y_j = pos_j 
    lines!(ax, [x_i, x_j], [y_i, y_j]; color)
end
        
# A = [0 1 1; 1 0 1; 1 0 0]
# adjlist = adjmat_to_adjlist(A)
# adjlist = [[2,4], [1,3,5], [2,4], [1,3], [2]]
# fig, ax = plot_network(adjlist)
# fig

al_regular = [[2,3, 9, 10], [1,3,4,10], [1,2,4,5], [2,3,5,6], [3,4,6,7], [4,5,7,8], [5,6,8,9], [6,7,9,10],[1,7,8,10], [1,2,8,9]]
al_sw = [[2,3,9], [1,3,4,10], [2,4,5], [2,5,6,9], [3,4,6,7], [4,5,7,8,10], [5,6,9], [6,9,10],[4,7,8,10], [2,6,8,9]]
al_random = [[6,5,9], [8,7,10], [4,5,6,9], [5,6,9], [1,3,4,7,8], [1,3,4,8,10], [2,5,9], [2,6,5,9,10],[1,3,4,7,8], [2,6,8]]

al_all = [al_regular, al_sw, al_random]
size_fig = (600,250)
fig = Figure(size=size_fig); axs=[]
for (idx, adjlist) in enumerate(al_all)
    ax = Axis(fig[1,idx]; aspect=1); push!(axs, ax)
    plot_network!(ax, adjlist)
end
fig     

label_axes!(axs; halign=:left)
titles = ["2-nearest-neighbors", "small-world", "random"]
for (ax,title) in zip(axs,titles) ax.title=title end

pixel_per_unit(out_size, dpi, res) = (out_size.*dpi./res)[1]
dpi = 600
out_width = 3.54 #inches 
px_per_unit = pixel_per_unit(out_width, dpi, size_fig)
save("wattsstrogatz.png", fig; px_per_unit)