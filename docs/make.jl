using CavitySpectroscopy
using Documenter

DocMeta.setdocmeta!(CavitySpectroscopy, :DocTestSetup, :(using CavitySpectroscopy); recursive=true)

makedocs(;
    modules=[CavitySpectroscopy],
    authors="Garrek Stemo <8449000+garrekstemo@users.noreply.github.com>",
    repo=Remotes.GitHub("garrekstemo", "CavitySpectroscopy.jl"),
    sitename="CavitySpectroscopy.jl",
    checkdocs=:exports,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://garrekstemo.github.io/CavitySpectroscopy.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Library" => "lib/public.md",
    ],
)

deploydocs(;
    repo="github.com/garrekstemo/CavitySpectroscopy.jl",
    devbranch="main",
)
