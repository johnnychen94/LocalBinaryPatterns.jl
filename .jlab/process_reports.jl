function aggregate_reports(outfile; rootdir=".", pattern=raw"report.md$")
    report_files = filter(x->x!=outfile && !isnothing(match(Regex(pattern), x)), readdir(rootdir))

    for report in report_files
        report_name = replace(strip(split(lowercase(report), "report")[1], '_'), '_'=>' ')
        @info "Processing \"$report\"..."
        report_name = titlecase(report_name) * " Report"
        open(outfile, "a") do io
            println(io, "\n<details>\n<summary> $report_name </summary>\n")
            open(report, "r") do src_io
                write(io, read(src_io))
            end
            println(io, "\n</details>\n")
        end
    end
    return report_files
end

outfile = "cml_report.md"
report_files = aggregate_reports(outfile)

if !isempty(report_files)
    run(`cml-send-comment $outfile`)
end
