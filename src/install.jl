using Pkg: add

for na in [
    "PyCall",
    "DataFrames",
    "Pandas",
]

    add(na)

end

add(url="https://github.com/KwatME/Kwat.jl")
