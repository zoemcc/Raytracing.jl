using Raytracing
using Documenter

makedocs(;
    modules=[Raytracing],
    authors="Zoe McCarthy <zoemccarthy12@gmail.com> and contributors",
    repo="https://github.com/zobot/Raytracing.jl/blob/{commit}{path}#L{line}",
    sitename="Raytracing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://zobot.github.io/Raytracing.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zobot/Raytracing.jl",
)
