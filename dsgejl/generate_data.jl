using DataFrames, Distributions

# True values of parameters
μ = 0.75
β = 0.9
σ = 0.25

dist = Normal(0, σ)

# Initialize states
u_t  = 0.0
u_t1 = 0.0

# Initialize DataFrame with 100 periods (2001-Q1 to 2025-Q4)
df = DataFrame(date = DSGE.quarter_range(DSGE.quartertodate("2001-Q1"), DSGE.quartertodate("2025-Q4")))
df[:x_t] = NaN

for t = 1:100
    # Set last period's u_t value to this period's u_{t-1}
    u_t1 = u_t

    # Draw new value of u_t
    u_t = rand(dist)

    # Apply measurement equation to get x_t
    x_t = μ + u_t + β*u_t1

    # Record in DataFrame
    df[t, :x_t] = x_t
end