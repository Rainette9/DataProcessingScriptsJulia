#!/usr/bin/env julia --project=.
"""
plot_raw_vs_processed.jl

For each day of SFC data, plot raw vs processed (despiked) time series
side by side for Ux, Uy, Uz, Ts, H2O. Saves one png per day.
"""

using Peddy, Dates, Statistics, CSV, DataFrames
using DimensionalData
using Plots

include("src/read_data,jl")

# === Configuration ===
sensor         = "SFC"
input_base     = "/home/engbers/Documents/PhD/EC_data_convert/2025/converted"
processed_base = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_HF/$sensor"
plot_dir       = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_HF/$sensor/plots"
mkpath(plot_dir)

# Variables to plot: (raw_column, processed_column, label)
plot_vars = [
    (:Ux,      :Ux,  "Ux [m/s]"),
    (:Uy,      :Uy,  "Uy [m/s]"),
    (:Uz,      :Uz,  "Uz [m/s]"),
    (:Ts,      :Ts,  "Ts [°C]"),
    (:LI_H2Om, :H2O, "H2O [mmol/m³]"),
]

# CSV options for reading processed .dat files
fo = Peddy.FileOptions(
    header = 1,
    delimiter = ",",
    comment = "#",
    timestamp_column = :timestamp,
    time_format = DateFormat("yyyy-mm-ddTHH:MM:SS.s"),
)

function read_processed_day(files::Vector{String}, opts::Peddy.FileOptions)
    dfs = DataFrame[]
    for f in files
        source = IOBuffer(replace(read(f, String), '"' => ""))
        df = CSV.read(source, DataFrame; header=opts.header, delim=opts.delimiter,
                      comment=opts.comment, dateformat=opts.time_format,
                      types=Dict(opts.timestamp_column => DateTime))
        push!(dfs, df)
    end
    combined = vcat(dfs...)
    sort!(combined, opts.timestamp_column)
    return combined
end

# === Discover processed files and group by date ===
date_pattern = r"_Fastdata_proc_(\d{4}-\d{2}-\d{2})_\d{4}\.dat$"
proc_files_by_date = Dict{Date, Vector{String}}()

for (root, _, files) in walkdir(processed_base)
    for f in files
        m = match(date_pattern, f)
        m === nothing && continue
        d = Date(m.captures[1], "yyyy-mm-dd")
        push!(get!(proc_files_by_date, d, String[]), joinpath(root, f))
    end
end

dates = sort(collect(keys(proc_files_by_date)))

if isempty(dates)
    error("No processed .dat files found under $processed_base")
end

@info "Found processed data for $(length(dates)) days"

# === Plot each day ===
for (i, date) in enumerate(dates)
    date_str = Dates.format(date, "yyyy-mm-dd")
    png_file = joinpath(plot_dir, "$(sensor)_raw_vs_proc_$(date_str).png")

    if isfile(png_file)
        println("  [$i/$(length(dates))] $date — plot exists, skipping")
        continue
    end

    println("--- [$i/$(length(dates))] Plotting $date ---")

    try
        # Read raw data for this day
        raw_data = read_fast_data(base_path=input_base, sensor=sensor, dates=[date])
        if nrow(raw_data) == 0
            @warn "No raw data for $date — skipping"
            continue
        end
        unique!(raw_data, :TIMESTAMP)
        sort!(raw_data, :TIMESTAMP)

        # Read processed data for this day
        day_files = sort(proc_files_by_date[date])
        proc_data = read_processed_day(day_files, fo)

        # Build subplots
        n_vars = length(plot_vars)
        subplots = []

        for (raw_col, proc_col, ylabel) in plot_vars
            # Skip if column doesn't exist in either dataset
            if !(string(raw_col) in names(raw_data)) || !(string(proc_col) in names(proc_data))
                continue
            end

            p = plot(raw_data.TIMESTAMP, Float64.(raw_data[!, raw_col]);
                label="raw", ylabel=ylabel, linewidth=0.3, alpha=0.5,
                color=:grey, legend=:topright, tickfontsize=6, guidefontsize=8,
                legendfontsize=6)
            plot!(p, proc_data.timestamp, Float64.(proc_data[!, proc_col]);
                label="processed", linewidth=0.3, color=:blue)

            push!(subplots, p)
        end

        if isempty(subplots)
            @warn "No matching variables for $date — skipping"
            continue
        end

        fig = plot(subplots...; layout=(length(subplots), 1),
                   size=(1400, 250 * length(subplots)),
                   plot_title="$sensor  $date  —  raw vs processed",
                   left_margin=10Plots.mm, bottom_margin=4Plots.mm)

        savefig(fig, png_file)
        println("  Saved $png_file")

    catch e
        @warn "Failed to plot $date" exception=(e, catch_backtrace())
    end

    GC.gc()
end

println("\nDone — $(length(dates)) days")
