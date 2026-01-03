using MultiFlowSolver
using Documenter

DocMeta.setdocmeta!(MultiFlowSolver, :DocTestSetup, :(using MultiFlowSolver); recursive=true)

makedocs(;
    modules=[MultiFlowSolver],
    authors="Daniel Araujo <dc_junior@id.uff.br>",
    sitename="MultiFlowSolver.jl",
    format=Documenter.HTML(;
        canonical="https://dcunha95.github.io/MultiFlowSolver.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dcunha95/MultiFlowSolver.jl",
    devbranch="master",
)
