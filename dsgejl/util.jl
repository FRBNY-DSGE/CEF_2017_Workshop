function quarter_date_to_number(date::Date)
    y = Dates.year(date)
    m = Dates.month(date)
    if m == 3
        return y
    elseif m == 6
        return y + 0.25
    elseif m == 9
        return y + 0.5
    elseif m == 12
        return y + 0.75
    end
end

function get_date_ticks(start_date::Date, end_date::Date;
                        tick_size::Int = 5)
    dates = DSGE.quarter_range(start_date, end_date)
    get_date_ticks(dates, tick_size = tick_size)
end

function get_date_ticks(dates::AbstractArray{Date, 1};
                        tick_size::Int = 5)
    datenums = map(quarter_date_to_number, dates)
    t0 = ceil(datenums[1] / tick_size) * tick_size
    t1 = datenums[end]
    ticks = t0:tick_size:t1
    return ticks
end

function get_date_limits(start_date::Nullable{Date}, end_date::Nullable{Date},
                         dates::AbstractArray{Date, 1})
    if isnull(start_date)
        start_date = Nullable(dates[1])
    end
    if isnull(end_date)
        end_date = Nullable(dates[end])
    end

    return start_date, end_date
end

function get_date_limit_indices(start_date::Nullable{Date}, end_date::Nullable{Date},
                                dates::AbstractArray{Date, 1})
    start_ind = if isnull(start_date)
        1
    else
        if dates[1] <= get(start_date) <= dates[end]
            findfirst(dates, get(start_date))
        elseif get(start_date) < dates[1]
            1
        else
            error("start_date $(get(start_date)) cannot be after last forecast period $(dates[end])")
        end
    end
    end_ind = if isnull(end_date)
        length(dates)
    else
        if dates[1] <= get(end_date) <= dates[end]
            findfirst(dates, get(end_date))
        elseif get(end_date) > dates[end]
            length(dates)
        else
            error("end_date $(get(end_date)) cannot be before first historical period $(dates[1])")
        end
    end
    return start_ind, end_ind
end

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
