using Pkg: add

add(url="https://github.com/KwatME/Kwat.jl")

for na in [
    "PyCall",
    "DataFrames",
    "Pandas",
]

    add(na)

end
