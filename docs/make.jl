using Pkg

Pkg.develop(path = "..")

using Publish
using Hyperelastics

p = Publish.Project(BitSAD)

function build_and_deploy(label)
    rm(label; recursive = true, force = true)
    deploy(Hyperelastics; root="/Hyperelastics.jl", label=label)
end