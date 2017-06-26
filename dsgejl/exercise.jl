using DSGE, Plots

do_setup    = true
do_estimate = true
do_plot     = true

if do_setup
    include("ma1.jl")
    include("generate_data.jl")
    include("util.jl")

    m = MA1()

    # Settings
    m <= Setting(:saveroot, dirname(@__FILE__))
    m <= Setting(:date_presample_start, DSGE.quartertodate("2001-Q1"))
    m <= Setting(:date_mainsample_start, DSGE.quartertodate("2001-Q1"))
    m <= Setting(:date_forecast_start, DSGE.quartertodate("2026-Q1"))
    m <= Setting(:n_mh_simulations, 500)
    m <= Setting(:n_mh_blocks, 10)
    m <= Setting(:n_mh_burn, 2)
end

if do_estimate
    # Sample from posterior distribution
    data = df_to_matrix(m, df)
    DSGE.estimate(m, data)
end

if do_plot
    # Initialize GR backend
    gr()

    # Plot prior distributions
    prior_draws = zeros(4500, 3)
    for i = 1:4000
        prior_draws[i, 1] = rand(get(m[:μ].prior))
        prior_draws[i, 2] = rand(get(m[:β].prior))
        prior_draws[i, 3] = rand(get(m[:σ].prior))
    end
    prior_plot = plot(prior_draws, t = [:histogram :histogram :histogram], layout = @layout([a b c]),
                      label = ["\\mu" "\\beta" "\\sigma"], title = ["" "Prior" ""])
    plot_actual!(prior_plot.subplots[1], μ)
    plot_actual!(prior_plot.subplots[2], β)
    plot_actual!(prior_plot.subplots[3], σ)

    # Plot posterior distributions
    post_draws = load_draws(m, :full)
    post_plot = plot(post_draws, t = [:histogram :histogram :histogram], layout = @layout([a b c]),
                     label = ["\\mu" "\\beta" "\\sigma"], title = ["" "Posterior" ""])
    plot_actual!(post_plot.subplots[1], μ)
    plot_actual!(post_plot.subplots[2], β)
    plot_actual!(post_plot.subplots[3], σ)

    # Plot both together
    all_plots = plot(prior_plot, post_plot, layout = @layout([a; b]))
end