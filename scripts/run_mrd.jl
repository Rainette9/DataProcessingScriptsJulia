#!/usr/bin/env julia --project=@.
"""
run_mrd.jl

For each day of processed fast data:
1. Read all hourly .dat files for that day.
2. Concatenate into a single DimArray.
3. Run Orthogonal MRD (M=17, normalize=true) for wT (Uz × Ts).
4. If LiCor H2O data is available and valid: also run MRD for wq (Uz × H2O).
5. Write per-scale summary statistics (median/q25/q75):
     MRD_<sensor>_<date>.dat       (wT)
     MRD_<sensor>_<date>_wq.dat    (wq, when LiCor available)
6. Save summary plots (.png) for each.
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
                                 number_type::Type{T}=Float64, strip_quotes::Bool=true) where {T<:Real}
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

function plot_mrd_summary(summary, date::Date, output_file::String, pair::String)
    median_vals = summary.median .* 1000.0
    q25_vals    = summary.q25 .* 1000.0
    q75_vals    = summary.q75 .* 1000.0

    finite_scales = filter(x -> isfinite(x) && x > 0, summary.scales)
    lo = floor(Int, log10(minimum(finite_scales)))
    hi = ceil(Int, log10(maximum(finite_scales)))
    decade_positions = [10.0^k for k in lo:hi]
    decade_labels = [LaTeXString("10^{$k}") for k in lo:hi]

    ylabel_str = pair == "wT" ?
        L"C_{T_s w}\ [10^{-3}\ \mathrm{K\,m\,s^{-1}}]" :
        L"C_{qw}\ [10^{-3}\ \mathrm{mmol\,m^{-2}\,s^{-1}}]"

    plt = plot(summary.scales, median_vals;
        title  = "MRD $sensor $date ($pair)",
        xlabel = L"\tau\ [\mathrm{s}]",
        ylabel = ylabel_str,
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

"""Return true if the :H2O column exists in `data` and has >10% finite positive values."""
function licor_available(data)
    :H2O ∉ dims(data, Var) && return false
    h2o = data[Var=At(:H2O)][:]
    valid = count(x -> isfinite(x) && x > 0, h2o)
    return valid / length(h2o) > 0.10
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

    wT_dat  = joinpath(output_path, "MRD_$(sensor)_$(date_str).dat")
    wT_png  = joinpath(output_path, "MRD_$(sensor)_$(date_str).png")
    wq_dat  = joinpath(output_path, "MRD_$(sensor)_$(date_str)_wq.dat")
    wq_png  = joinpath(output_path, "MRD_$(sensor)_$(date_str)_wq.png")
    wq_skip = joinpath(output_path, "MRD_$(sensor)_$(date_str)_wq.nolicor")

    wT_done = isfile(wT_dat) && isfile(wT_png)
    wq_done = (isfile(wq_dat) && isfile(wq_png)) || isfile(wq_skip)

    if wT_done && wq_done
        println("  [$i/$(length(dates))] $date — already done, skipping")
        continue
    end

    println("--- [$i/$(length(dates))] MRD for $date ---")

    try
        day_files = sort(files_by_date[date])
        arrays = [read_processed_dimarray(f, fo) for f in day_files]
        day_data = length(arrays) == 1 ? arrays[1] : cat(arrays...; dims=Ti)

        n_times = size(day_data, Ti)
        println("  $(length(day_files)) files, $n_times samples")

        shift = round(Int, 0.1 * 2^17)

        # --- wT MRD ---
        if !wT_done
            mrd_wT = Peddy.OrthogonalMRD(M=17, shift=shift, normalize=true, regular_grid=true,
                                          a=:Uz, b=:Ts)
            Peddy.decompose!(mrd_wT, day_data, nothing)
            res_wT = Peddy.get_mrd_results(mrd_wT)
            if res_wT === nothing
                @warn "No wT MRD results for $date — skipping"
            else
                summary_wT = summarize_mrd(res_wT)
                write_mrd_dat(summary_wT, wT_dat)
                plot_mrd_summary(summary_wT, date, wT_png, "wT")
                println("  Wrote $wT_dat")
                println("  Wrote $wT_png")
            end
        end

        # --- wq MRD ---
        if !wq_done
            if licor_available(day_data)
                mrd_wq = Peddy.OrthogonalMRD(M=17, shift=shift, normalize=true, regular_grid=true,
                                              a=:Uz, b=:H2O)
                Peddy.decompose!(mrd_wq, day_data, nothing)
                res_wq = Peddy.get_mrd_results(mrd_wq)
                if res_wq === nothing
                    @warn "No wq MRD results for $date — skipping"
                else
                    summary_wq = summarize_mrd(res_wq)
                    write_mrd_dat(summary_wq, wq_dat)
                    plot_mrd_summary(summary_wq, date, wq_png, "wq")
                    println("  Wrote $wq_dat")
                    println("  Wrote $wq_png")
                end
            else
                println("  No valid LiCor H2O data for $date — skipping wq")
                touch(wq_skip)
            end
        end

    catch e
        @warn "Failed MRD for $date" exception=(e, catch_backtrace())
    end

    GC.gc()
end

println("\nDone — $(length(dates)) days processed")
