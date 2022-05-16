
using Pkg

Pkg.develop(path = "..")

using Publish
using Hyperelastics
using Artifacts, LazyArtifacts

# override default theme
cp(artifact"flux-theme", "../_flux-theme"; force = true)

function build_and_deploy(label)
    rm(label; recursive = true, force = true)
    deploy(Hyperelastics; root = "/Hyperelastics.jl", label = label)
end