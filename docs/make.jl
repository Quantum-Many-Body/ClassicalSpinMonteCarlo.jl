using ClassicalSpinMonteCarlo
using Documenter

DocMeta.setdocmeta!(ClassicalSpinMonteCarlo, :DocTestSetup, :(using ClassicalSpinMonteCarlo); recursive=true)

makedocs(;
    modules=[ClassicalSpinMonteCarlo],
    authors="wwangnju <wwangnju@163.com> and contributors",
    repo="https://github.com/Quantum-Many-Body/ClassicalSpinMonteCarlo.jl/blob/{commit}{path}#{line}",
    sitename="ClassicalSpinMonteCarlo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Quantum-Many-Body.github.io/ClassicalSpinMonteCarlo.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Quantum-Many-Body/ClassicalSpinMonteCarlo.jl",
    devbranch="master",
)
