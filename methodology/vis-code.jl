
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