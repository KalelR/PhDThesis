using GLMakie

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

r = 1.0
θs = range(0, 2pi, length=1000)
xs = r .* cos.(θs)
ys = r .* sin.(θs)

θs_units = range(pi/2-pi/3, pi/2+pi/3, length=20)
xs_units = r .* cos.(θs_units)
ys_units = r .* sin.(θs_units)
colors_units = get_colors(length(θs_units), :distinguishable_colors)

fig = Figure()
ax = Axis(fig[1,1]; aspect=1)
lines!(ax, xs, ys; color=:black)
scatter!(ax, xs_units, ys_units; color=colors_units, markersize=15)

r, ϕ =    complex_order_parameter_t(θs_units)
arrow_head = [r*cos(ϕ), r*sin(ϕ)]

arrows!(ax, [0], [0], [arrow_head[1]], [arrow_head[2]])#; linewidth=3, arrowsize=2.0)
hidedecorations!(ax)
hidespines!(ax)    

fig