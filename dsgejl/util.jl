function plot_actual!(p::Plots.Subplot{Plots.GRBackend}, v::Float64)
    # Get y-axis limits
    ylim = p.attr[:yaxis].d[:lims]
    
    # Plot actual value
    plot!(p, [v, v], [ylim[1], ylim[2]], linewidth = 3, label = "", color = :black)
end

function plot_actual!(p::Plots.Subplot{Plots.PlotlyBackend}, v::Float64)
    # Get y-axis limits
    ylim = p.attr[:yaxis].d[:extrema]
    
    # Plot actual value
    plot!(p, [v, v], [ylim.emin, ylim.emax], linewidth = 3, label = "", color = :black)
end