#!/usr/bin/env julia

module combine_bwig_stats

using ArgParse
using CodecZlib
using CSV
using DataFrames
using Formatting
using Mmap

function parse_cmd_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--column-labels-to-use-for-joining"
        help = "CSV of column labels used to join datasets."
        arg_type = String
        required = true
        "--column-suffixes"
        help = "CSV of suffixes to append to column labels."
        arg_type = String
        required = true
        "--sort-by"
        help = "List of column names to use for sorting purposes."
        arg_type = String
        required = true
        "--sort-descending"
        help = "Sort in descending order."
        default = true
        "TSVFiles"
        help = "TSV Files to be processed."
        arg_type = String
        nargs = '+'
        required = true
    end

    return parse_args(ARGS, s)
end

function read_df(
    fname::String,
    labels::Union{Vector{String},Set{String},Nothing} = nothing,
)::DataFrame
    df =
        endswith(fname, ".gz") ?
        DataFrame(CSV.File(transcode(GzipDecompressor, Mmap.mmap(fname)))) :
        DataFrame(CSV.File(fname))
    if labels === nothing
        return df
    end

    cnames = Set(names(df))
    missing_labels = Vector{String}()
    for label in labels
        if !(label in cnames)
            push!(missing_labels, label)
        end
    end

    if isempty(missing_labels)
        return df
    end
    error(
        format(
            "Unable to find the following columns in file \"{:s}\": {:s}\n",
            fname,
            join(missing_labels, ", "),
        ),
    )
end

function rename_dupl_columns(df::DataFrame, suffixes::Vector{String})
    labels = ["mean", "median", "std", "min", "max"]
    mappings = [Pair(Symbol(l), Symbol(join([l, strip(suffixes[1], '_')], '_'))) for l in labels]
    for i = 2:length(suffixes)
        append!(
            mappings,
            [
                Pair(
                    Symbol(format("{:s}_{:d}", l, i - 1)),
                    Symbol(join([l, strip(suffixes[i], '_')], '_')),
                ) for l in labels
            ],
        )
    end
    rename!(df, mappings)
end

function main()::Cint
    args = parse_cmd_args()
    files = args["TSVFiles"]
    suffixes = [String(strip(s)) for s in split(args["column-suffixes"], ',')]
    labels = [String(strip(s)) for s in split(args["column-labels-to-use-for-joining"], ',')]
    sort_by = [Symbol(l) for l in split(args["sort-by"], ',')]

    if length(files) != length(suffixes)
        error(
            format(
                "The number of files and suffixes should be the same. Got {:d} files and {:d} suffixes\n",
                length(files),
                length(suffixes),
            ),
        )
    end

    dfs = Vector{DataFrame}([read_df(files[1], labels)])

    for i = 2:length(files)
        push!(dfs, select(read_df(files[i], labels), Not(names(dfs[1])[1])))
    end

    df = innerjoin(dfs..., on = labels, makeunique = true)
    rename_dupl_columns(df, suffixes)
    sort!(df, sort_by, rev = [args["sort-descending"] for _ in sort_by])
    CSV.write(stdout, df, delim = '\t')
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
