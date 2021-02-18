#!/usr/bin/env julia

module compute_statistics_for_bwigs

using ArgParse
using BED
using BigWig
using Formatting
using Statistics
using GenomicFeatures
using Base.Threads

function parse_cmd_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--chr-subranges-bed"
        help = "Path to a BED file containing the genomic regions to be processed."
        arg_type = String
        "--extrusion-barriers-bed"
        help = "Path to a BED file containing the genomic coordinates of the extrusion barriers simulated (e.g. CTCF motifs)."
        arg_type = String
        "bigWigFiles"
        help = "bigWigFiles to be processed. Also accepts the path to one or more text files containing the list of bigWig files to be processed."
        arg_type = String
        nargs = '+'
        required = true
    end

    return parse_args(ARGS, s)
end

function parse_subranges(path_to_bed::Union{String,Nothing})::IntervalCollection
    if path_to_bed === nothing
        return IntervalCollection{Nothing}()
    end

    return IntervalCollection(
        [
            Interval(BED.chrom(b), BED.chromstart(b), BED.chromend(b)) for
            b in open(BED.Reader, path_to_bed)
        ],
        true,
    )
end

function parse_file_list(files::Array{String,1})::Vector{String}
    file_list = Vector{String}()

    for f in files
        if endswith(f, ".bw")
            if !ispath(f)
                throw(SystemError(f, 2))
            end
            push!(file_list, f)

        else
            open(f, "r") do fp
                for line in eachline(fp)
                    if !ispath(f)
                        throw(SystemError(f, 2))
                    end
                    push!(file_list, strip(line))
                end
            end
        end
    end
    return file_list
end

function parse_extr_barriers(
    fpath::Union{String,Nothing},
    ranges::IntervalCollection,
)::IntervalCollection
    if fpath === nothing
        return IntervalCollection{Nothing}()
    end

    barriers = IntervalCollection(
        [
            Interval(BED.chrom(b), BED.chromstart(b), BED.chromend(b)) for
            b in open(BED.Reader, fpath)
        ],
        true,
    )
    if ranges === IntervalCollection{Nothing}()
        return barriers
    end
    return IntervalCollection([b for (b, _) in eachoverlap(barriers, ranges)])
end

function compute_stats(
    fpath::String,
    ranges::IntervalCollection,
    barriers::IntervalCollection,
    chan::Base.Channel,
    fexpr,
)
    open(BigWig.Reader, fpath) do bw
        chroms =
            IntervalCollection([Interval(name, 0, size) for (name, size) in BigWig.chromlist(bw)])
        if ranges !== IntervalCollection{Nothing}()
            chroms = IntervalCollection([record for (_, record) in eachoverlap(chroms, ranges)])
        end
        for chrom in chroms
            cname = seqname(chrom)
            cstart = leftposition(chrom)
            cend = rightposition(chrom)
            if isempty(barriers)
                _mean = BigWig.mean(bw, cname, cstart, cend)
                _median =
                    Statistics.median([Float32(record.value) for record in eachoverlap(bw, chrom)])
                _std = BigWig.std(bw, cname, cstart, cend)
                _min = BigWig.minimum(bw, cname, cstart, cend)
                _max = BigWig.maximum(bw, cname, cstart, cend)
            else
                values = Vector{Float32}()
                for interval in
                    IntervalCollection([record for record in eachoverlap(barriers, chrom)])
                    push!(
                        values,
                        [Float32(record.value) for record in eachoverlap(bw, interval)][1],
                    )
                end
                _mean = Statistics.mean(values)
                _median = Statistics.median(values)
                _std = Statistics.std(values)
                _min = Statistics.minimum(values)
                _max = Statistics.maximum(values)
            end
            params = join(get_params_from_fname(fpath), '\t')
            put!(
                chan,
                format(fexpr, fpath, cname, cstart, cend, params, _mean, _std, _median, _min, _max),
            )
        end
    end
end

function get_param_labels_from_fname(fname::String, delim::Union{String,Char} = '~')::Vector{String}
    # Tokenize file path and remove tokens that do not refer to a parameter
    toks = split(fname, '/', keepempty = false)
    toks = toks[[delim in t for t in toks]]
    # Return param names
    return [getindex(split(t, delim), 1) for t in toks]
end

function get_params_from_fname(fname::String, delim::Union{String,Char} = '~')::Vector{String}
    # Tokenize file path and remove tokens that do not refer to a parameter
    toks = split(fname, '/', keepempty = false)
    toks = toks[[delim in t for t in toks]]
    # Return param names
    return [getindex(split(t, delim), 2) for t in toks]
end


function main()::Cint
    args = parse_cmd_args()
    files = parse_file_list(args["bigWigFiles"])
    ranges = parse_subranges(args["chr-subranges-bed"])
    barriers = parse_extr_barriers(args["extrusion-barriers-bed"], ranges)
    printfmtln(stderr, "# of files to be procesed: {:d}", length(files))
    ch = Channel{String}(Threads.nthreads() * 2)
    @async begin
        Threads.@threads for f in files
            compute_stats(f, ranges, barriers, ch, fexpr)
        end
        close(ch)
    end
    fexpr = FormatExpr("{:s}\t{:s}\t{:d}\t{:d}\t{:s}\t{:f}\t{:f}\t{:f}\t{:f}\t{:f}")

    params = join(get_param_labels_from_fname(files[1]), '\t')
    printfmtln("file\tchr\tstart\tend\t{:s}\tmean\tstd\tmedian\tmin\tmax", params)
    for msg in ch
        println(msg)
    end
    return 0
end

function julia_main()
    try
        main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module
