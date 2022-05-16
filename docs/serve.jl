using Pkg

Pkg.develop(path = "..")

using Revise
using Publish
using Hyperelastics
using Artifacts, LazyArtifacts

# override default theme
cp(artifact"flux-theme", "../_flux-theme"; force = true)

p = Publish.Project(Hyperelastics)

# serve documentation
serve(Hyperelastics)