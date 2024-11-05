using GLMakie, Random, ColorSchemes, Colors


"""
    lighten(c, f = 1.2)

Lighten given color `c` by multiplying its luminance with `f`.
If `f` is less than 1, the color is darkened.
"""
function lighten(c, f = 1.2)
    c = to_color(c)
    hsl = Makie.HSLA(c)
    neg = Makie.RGBAf(Makie.HSLA(hsl.h, hsl.s, clamp(hsl.l*f, 0.0, 1.0), hsl.alpha))
    neg = Makie.RGBf(Makie.HSL(hsl.h, hsl.s, clamp(hsl.l*f, 0.0, 1.0)))
    return neg
end

# Create a grid of x, y values
x = range(-8, 8, length=100)
y = range(-8, 8, length=100)
function meshgrid(x::AbstractVector, y::AbstractVector)
    X = [i for i in x, j in y]  # Matrix with x repeated in rows
    Y = [j for i in x, j in y]  # Matrix with y repeated in columns
    return X, Y
end
X, Y = meshgrid(x, y)
# X = x' .* ones(length(y))
# Y = y' .* ones(length(x))

function gaussian2D(A, x0, y0, sigma_x, sigma_y, X, Y)
    return A * exp.(-((X .- x0).^2 ./ (2 * sigma_x^2) .+ (Y .- y0).^2 ./ (2 * sigma_y^2)))
end

# using FastNoise2


# function generate_perlin_noise(x, y; scale, amplitude)
#     noise_model = FastNoise2.Noise(cellular_return_type=FastNoise2.ReturnType.Distance, 
#                                    noise_type=FastNoise2.NoiseType.Perlin)
#     Z_noise = [amplitude * noise_model(xi * scale, yi * scale) for xi in x, yi in y]  # Generate noise for each point
#     return Z_noise
# end

function generate_perlin_noise(x, y; scale, amplitude)
    noise_model = FastNoise2.Noise(style=:perlin, frequency=scale)  # Create a Perlin noise model
    Z_noise = [amplitude * noise_model(xi, yi) for xi in x, yi in y]  # Generate noise for each point
    return Z_noise
end

function generate_perlin_noise(x, y; scale, amplitude)
    noise_model = FastNoise2.perlin()  # Create a Perlin noise model
    Z_noise = [amplitude * noise_model(xi * scale, yi * scale) for xi in x, yi in y]  # Generate noise for each point
    return Z_noise
end

function generate_random_noise(X, Y, amplitude; randseed=1)
    Z_noise = amplitude * randn(MersenneTwister(randseed), size(X))  # Random noise with specified amplitude
    return Z_noise
end

# Define parameters for the central peak
A0, x0, y0, sigma0x, sigma0y = 5.0, 0.0, 0.0, 1.0, 1.0

function landscape(X,Y)
# Z .+= gaussian2D(A1, x1, y1, sigma1x, sigma1y, X, Y)  # First smaller peak
Z = gaussian2D(4.0, 0.0,0.0,1.0,1.0, X, Y)  # Central peak
Z .+= gaussian2D(3.0, 5.0, -2.5, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(2.2, 5.2, 1.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(1.0, 4.0, 2.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(2.0, -5.0, -2.5, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(2.0, -5.0, 0.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(2.0, -5.0, 2.5, 1.0, 1.0, X, Y)  # First smaller peak

Z .+= gaussian2D(2.0, -2.0, 3.0, 1.3, 1.3, X, Y)  # First smaller peak
Z .+= gaussian2D(1.0, 2.0, 5.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(0.8, 0.0, 5.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(1.5, 2.5, 2.5, 1.4, 1.4, X, Y)  # First smaller peak
Z .+= gaussian2D(2.0, -2.5, -5.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(2.0, 0, -5.0, 1.0, 1.0, X, Y)  # First smaller peak
Z .+= gaussian2D(3.5, 2.0, -3.0, 1.0, 1.0, X, Y)  # First smaller peak
return Z 
end

Z = landscape(X, Y)

# Z_noise = generate_perlin_noise(collect(x), collect(y), scale=1.0, amplitude=0.5)
Z_noise = generate_random_noise(X, Y, 0.1)
# Z_noise = generate_random_noise(X, Y, 0.0)
Z .+= Z_noise

cmap = ColorSchemes.darkterrain.colors
# Function to brighten the colormap by a factor
function brighten_colormap(colormap, factor)
    # Increase brightness by scaling the RGB components of each color
    return [RGB(clamp(c.r * factor, 0, 1), clamp(c.g * factor, 0, 1), clamp(c.b * factor, 0, 1)) for c in colormap]
end

using ColorSchemes
# Use a factor > 1 to brighten (e.g., 1.2)
brighter_colormap = brighten_colormap(cmap, 1.2)
# brighter_colormap = brighten_colormap(cmap, 1.3)


# Plot the surface
fig = Figure(;size=(600,500), margins=0); axs=[]
ax = Axis3(fig[1, 1], azimuth=31.22, elevation=1.13, viewmode=:stretch, protrusions=-5); push!(axs,ax)  # Set initial azimuth and elevation
# surface!(ax, x, y, Z; colormap=:darkterrain)
surface!(ax, x, y, Z; colormap=brighter_colormap)
hidedecorations!(ax)
hidespines!(ax)


# fig, ax, plt =surface(x, y, Z; colormap=:darkterrain)
# fig,ax,plt=surface(x, y, Z; colormap=:darkterrain, axis=(type=Axis3, xspinesvisible = false, yspinesvisible = false, zspinesvisible = false, xgridvisible=false, ygridvisible=false, zgridvisible=false, xticksvisible=false, yticksvisible=false, zticksvisible=false, xlabelvisible=false, ylabelvisible=false, zlabelvisible=false, xticks=[1000], yticks=[1000], zticks=[1000]))

function parametrize_trajectory(x0, y0, x1, y1, t)
    # Linear interpolation between (x0, y0) and (x1, y1)
    x_t = (1 - t) * x0 + t * x1
    y_t = (1 - t) * y0 + t * y1
    return x_t, y_t
end   

function plot_trajectory!(x0, y0, x1, y1; t_values = range(0, 1, length=100), plot_ball=false, color_traj=:red, color_ball=:purple)
    trajectory_x = []; trajectory_y = []; trajectory_z = []
    for t in t_values
        x_t, y_t = parametrize_trajectory(x0, y0, x1, y1, t)
        push!(trajectory_x, x_t)
        push!(trajectory_y, y_t)
        push!(trajectory_z, landscape([x_t], [y_t])[1])
    end
    # lift=0.2
    lift=0.25
    lines!(trajectory_x, trajectory_y, trajectory_z.+lift, color=color_traj, linewidth=3)
    if plot_ball
    scatter!(x0, y0, landscape([x0], [y0])[1].+lift; color=color_ball, markersize=15)
    end
end

plot_trajectory!(2.5, 2.5, 2.5, 0.0; plot_ball=true, color_traj=:red)
plot_trajectory!(2.5, 2.5, 0.5, 3.5; plot_ball=false, color_traj=:green)


plot_trajectory!(1.5, 1.5, 2.5, 0.0; plot_ball=true, color_traj=:red)
# plot_trajectory!(1.5, 1.5, 0.5, 3.5; plot_ball=false, color_traj=:green)
plot_trajectory!(1.5, 1.5, 0.5, 3.5; plot_ball=false, color_traj=:green)

plot_trajectory!(5.0, -0.5, 2.5, 0.0; plot_ball=true, color_traj=:red)
# plot_trajectory!(5.0, -0.5, 6.3, -0.5; plot_ball=false, color_traj=:green)
plot_trajectory!(5.0, -0.5, 6.3, -0.5; plot_ball=false, color_traj=:magenta)

plot_trajectory!(3.5, -3, 2.5, 0.0; plot_ball=true, color_traj=:red)
# plot_trajectory!(3.5, -3, 3.5, -5; plot_ball=false, color_traj=:green)
plot_trajectory!(3.5, -3, 3.5, -5; plot_ball=false, color_traj=:yellow)

plot_trajectory!(1.5, -1.5, 2.5, 0.0; plot_ball=true, color_traj=:red)
# plot_trajectory!(1.5, -1.5, -1.0, -3.0; plot_ball=true, color_traj=:green)
plot_trajectory!(1.5, -1.5, -1.0, -3.0; plot_ball=true, color_traj=:cyan)



function plot_trajectory_2d!(x0, x1; t_values = range(0, 1, length=100), plot_ball=false, color_traj=:red, color_ball=:purple)
    trajectory_x = []; trajectory_y = []; 
    xs = range(x0, x1; length=100)
    ys = potential.(xs)
    lines!(xs, ys, color=color_traj, linewidth=3)
    if plot_ball
    scatter!(x0, ys[1]; color=color_ball, markersize=15)
    end
end


# ax = Axis(fig[1,2]); push!(axs, ax)
ax = Axis(fig[2,1]); push!(axs, ax)
xs = range(-1.13, 1.13, length=100)
potential(x) = x^4 - x^2
ys = potential.(xs)
lines!(ax, xs, ys; color=:black)
# colsize!(fig.layout, 2, Relative(0.3))
rowsize!(fig.layout, 2, Relative(0.2))
plot_trajectory_2d!(0, -0.707; plot_ball=true)
plot_trajectory_2d!(0, +0.707; color_traj=:green)
hidedecorations!(ax)
hidespines!(ax)


include("../methodology/vis-code.jl")
label_axes!(axs; halign=:left)


save("landscape.png", fig; px_per_unit=4)