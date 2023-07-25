using Sym4state
using Documenter

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
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/a-lost-wapiti/Sym4state.jl",
    devbranch="main",
)
