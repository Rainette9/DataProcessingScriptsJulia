#!/usr/bin/env julia --project=.
"""
run_mrd.jl

For each day of processed SFC fast data:
1. Read all hourly .dat files for that day.
2. Concatenate into a single DimArray.
3. Run Orthogonal MRD (M=16, normalize=true).
4. Write per-scale summary statistics (median/q25/q75) -> MRD_SFC_<date>.dat
5. Save a summary plot -> MRD_SFC_<date>.png
"""

using Peddy, Dates, Statistics, CSV, DataFrames
using DimensionalData
using Plots
using LaTeXStrings
using Glob

# === Configuration ===
sensor       = "LOWER"
input_base   = "/home/engbers/Documents/PhD/EC_data_convert/2025/processed_HF/$sensor"
output_path  = "/home/engbers/Documents/PhD/EC_data_convert/2025/MRD/$sensor"
mkpath(output_path)

fo = Peddy.FileOptions(
    header = 1,
    delimiter = ",",
    comment = "#",
    timestamp_column = :timestamp,
    time_format = DateFormat("yyyy-mm-ddTHH:MM:SS.s"),
)

# === Helpers ===

function read_processed_dimarray(path::AbstractString, opts::Peddy.FileOptions;
                                 N::Type{T}=Float64, strip_quotes::Bool=true) where {T<:Real}
    tscol = opts.timestamp_column
    types_map = Dict(tscol => DateTime)
    source = strip_quotes ? IOBuffer(replace(read(path, String), '"' => "")) : path
    f = CSV.File(source; header=opts.header, delim=opts.delimiter, comment=opts.comment,
                 types=types_map, dateformat=opts.time_format)
    timestamps = f[tscol]
    vars = Symbol[x for x in f.names if x != tscol]
    data = Matrix{T}(undef, length(timestamps), length(vars))
    for (j, v) in pairs(vars)
        data[:, j] .= f[v]
    end
    DimArray(data, (Peddy.Ti(timestamps), Peddy.Var(vars)))
end

function summarize_mrd(res)
    scales = collect(res.scales)
    A = res.mrd
    M, _ = size(A)
    med = Vector{Float64}(undef, M)
    q25 = Vector{Float64}(undef, M)
    q75 = Vector{Float64}(undef, M)
    for i in 1:M
        vals = filter(!isnan, collect(@view A[i, :]))
        if isempty(vals)
            med[i] = NaN; q25[i] = NaN; q75[i] = NaN
        else
            med[i] = median(vals)
            q25[i] = quantile(vals, 0.25)
            q75[i] = quantile(vals, 0.75)
        end
    end
    return (; scales, median=med, q25, q75)
end

function write_mrd_dat(summary, filepath; delim=",")
    open(filepath, "w") do io
        println(io, join(["scale_s", "median", "q25", "q75"], delim))
        for i in eachindex(summary.scales)
            println(io, summary.scales[i], delim,
                    summary.median[i], delim,
                    summary.q25[i], delim,
                    summary.q75[i])
        end
    end
end

function plot_mrd_summary(summary, date::Date, output_file::String)
    median_vals = summary.median .* 1000.0
    q25_vals    = summary.q25 .* 1000.0
    q75_vals    = summary.q75 .* 1000.0

    finite_scales = filter(x -> isfinite(x) && x > 0, summary.scales)
    lo = floor(Int, log10(minimum(finite_scales)))
    hi = ceil(Int, log10(maximum(finite_scales)))
    decade_positions = [10.0^k for k in lo:hi]
    decade_labels = [LaTeXString("10^{$k}") for k in lo:hi]

    plt = plot(summary.scales, median_vals;
        title  = "MRD $sensor $date",
        xlabel = L"\tau [\mathrm{s}]",
        ylabel = L"C_{T_s w} [10^{-3} \mathrm{K m s}^{-1}]",
        xscale = :log10,
        xticks = (decade_positions, decade_labels),
        xminorgrid = true,
        minorgrid  = true,
        legend = :topright,
        label  = "median",
    )
    plot!(plt, summary.scales, q75_vals;
        fillrange = q25_vals,
        label     = "quartile range",
        fillalpha = 0.25,
        linealpha = 0.0,
        linecolor = :transparent,
    )
    savefig(plt, output_file)
end

# === Discover all processed .dat files and group by date ===
all_files = String[]
for (root, _, files) in walkdir(input_base)
    for f in files
        if occursin(r"_Fastdata_proc_\d{4}-\d{2}-\d{2}_\d{4}\.dat$", f)
            push!(all_files, joinpath(root, f))
        end
    end
end
sort!(all_files)

if isempty(all_files)
    error("No processed .dat files found under $input_base")
end

# Extract date from filename and group
date_pattern = r"_Fastdata_proc_(\d{4}-\d{2}-\d{2})_\d{4}\.dat$"
files_by_date = Dict{Date, Vector{String}}()
for f in all_files
    m = match(date_pattern, basename(f))
    m === nothing && continue
    d = Date(m.captures[1], "yyyy-mm-dd")
    push!(get!(files_by_date, d, String[]), f)
end

dates = sort(collect(keys(files_by_date)))
@info "Found $(length(all_files)) files across $(length(dates)) days"

# === Process each day ===
for (i, date) in enumerate(dates)
    date_str = Dates.format(date, "yyyy-mm-dd")
    dat_file = joinpath(output_path, "MRD_$(sensor)_$(date_str).dat")
    png_file = joinpath(output_path, "MRD_$(sensor)_$(date_str).png")

    # Skip if both outputs already exist
    if isfile(dat_file) && isfile(png_file)
        println("  [$i/$(length(dates))] $date — already done, skipping")
        continue
    end

    println("--- [$i/$(length(dates))] MRD for $date ---")

    try
        # Read and concatenate all hourly files for this day
        day_files = sort(files_by_date[date])
        arrays = [read_processed_dimarray(f, fo) for f in day_files]
        day_data = length(arrays) == 1 ? arrays[1] : cat(arrays...; dims=Ti)

        n_times = size(day_data, Ti)
        println("  $(length(day_files)) files, $n_times samples")

        # Run MRD
        shift = round(Int, 0.1 * 2^17)
        mrd = Peddy.OrthogonalMRD(M=17, shift=shift, normalize=true, regular_grid=true)
        Peddy.decompose!(mrd, day_data, nothing)
        res = Peddy.get_mrd_results(mrd)

        if res === nothing
            @warn "No MRD results for $date — skipping"
            continue
        end

        # Summarize and write
        summary = summarize_mrd(res)
        write_mrd_dat(summary, dat_file)
        plot_mrd_summary(summary, date, png_file)
        println("  Wrote $dat_file")
        println("  Wrote $png_file")

    catch e
        @warn "Failed MRD for $date" exception=(e, catch_backtrace())
    end

    GC.gc()
end

println("\nDone — $(length(dates)) days processed")
