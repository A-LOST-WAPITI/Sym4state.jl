using Documenter
import Sym4state

DocMeta.setdocmeta!(Sym4state, :DocTestSetup, :(using Sym4state); recursive=true)

makedocs(;
    modules=[Sym4state],
    authors="Guolin Wan <glwan@iphy.ac.cn> and contributors",
    repo="https://github.com/a-lost-wapiti/Sym4state.jl/blob/{commit}{path}#{line}",
    sitename="Sym4state.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://a-lost-wapiti.github.io/Sym4state.jl",
        edit_link="main",
        assets=String[],
        mathengine = MathJax3(Dict(
            :loader => Dict("load" => ["[tex]/mhchem"]),
        ))
    ),
    pages=[
        "Home" => "home.md",
        "Manual" => "manual.md",
        "Types" => "types.md",
        "Utilities" => "utils.md",
        "API" => "index.md"
    ],
)

deploydocs(;
    repo="github.com/A-LOST-WAPITI/Sym4state.jl",
    devbranch="main",
)
