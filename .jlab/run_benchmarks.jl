using PkgBenchmark

# https://docs.gitlab.com/ee/ci/variables/predefined_variables.html
if haskey(ENV, "CI_MERGE_REQUEST_TARGET_BRANCH_NAME")
    # merge requests
    baseline = String(ENV["CI_MERGE_REQUEST_TARGET_BRANCH_NAME"])
    # This works around shallow git clone errors: GitError(Code:ENOTFOUND, Class:Reference, revspec 'master' not found)
    run(`git checkout $baseline`)
    run(`git checkout -`)
    results = PkgBenchmark.judge(".", baseline)
    export_markdown("benchmark_judge_report.md", results)
    export_markdown("benchmark_baseline_report.md", results.baseline_results)
    export_markdown("benchmark_target_report.md", results.target_results)
else
    # branch pushes
    results = PkgBenchmark.benchmarkpkg(".")
    export_markdown("benchmark_report.md", results)
end
