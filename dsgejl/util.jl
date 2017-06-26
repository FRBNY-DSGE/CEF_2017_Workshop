function plot_actual!(p::Plots.Subplot, v::Float64)
    # Get y-axis limits
    (y0, y1) = p.attr[:yaxis].d[:lims]
    plot!(p, [v, v], [y0, y1], linewidth = 3, label = "Actual")
end