push!(LOAD_PATH, ".")
using MathematicsToolsWYPIWYG
using WYPIWYGabTools

using Plots
using LaTeXStrings


pyplot()
default(dpi=1000) # Only for PyPlot

###################### Plot Configuration #################################

# font specification
fntsm = Plots.font("Computer Moddern", 18)
fntlg = Plots.font("Computer Modern", 18)

default(titlefont=fntlg,
        guidefont=fntlg,
        tickfont=fntsm,
        legendfont=fntsm)

# Size of the figure frame
dim1 = 600
dim2 = 600

# Lines and markers specification
mycolors = [:blue :red :black]
mylinewidth = 2.2
mymarkersize = 9.
mylynestyles = [:solid :dash :dot]
mymarkershapes = [:circle :rect :diamond ]
mymarkerstrokewidth = 2.
# Pictures save route

Pict_route = "./Pictures/arruda_wypiwyg_treloar"    # default  directory

###########################################################################
# Specify data route
Data_route ="./TestData/"   # default directory to find the experimental data
# Loading data from Treloar et al.
λuni, Puni = readCurve(Data_route*"treloar-uni.txt")
λbiax, Pbiax = readCurve(Data_route*"treloar-equibi.txt")
λshear, Pshear = readCurve(Data_route*"treloar-shear.txt")
np = 100
λuni_dense = [λuni[1]:(λuni[end]-λuni[1])/(convert(Float64,np)-1.):λuni[end] ...]
λbiax_dense = [λbiax[1]:(λbiax[end]-λbiax[1])/(convert(Float64,np)-1.):λbiax[end] ...]
λshear_dense = [λshear[1]:(λshear[end]-λshear[1])/(convert(Float64,np)-1.):λshear[end] ...]

################################################################################
#################### WYPIWYG using Bsplines and uniaxial data ##################
################################################################################
# Computing f_λch from uniaxial test, λu_block = 9.0
f_λch_uniax, df_λch_uniax = ffromtests(λuni, Puni, 9.0, 1, 2)
# Computing σ_w_uniax using f_λch from uniaxial test
σ_w_uniax = sigma_wypiwyg(λuni_dense,f_λch_uniax)

λch_plotf_uni = [1.:3.62/599.:4.62 ...]

############################## FIGURE #########################################
# plotting f_λch obtained from uniaxial test
p = plot(λch_plotf_uni, f_λch_uniax.(λch_plotf_uni),
        label = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from uniaxial"),
        linecolor = mycolors[1],
        linewidth = mylinewidth,
        linestyle = mylynestyles[1],
        xlabel= latexstring("\$ \\lambda_{ch}\\;[-]\$"),
        ylabel = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;[MPa]\$"),
        size = (dim1,dim2)
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from uniaxial test"))
display(p)
savefig(Pict_route*"/fchBSuni.eps")
savefig(Pict_route*"/fchBSuni.png")
savefig(Pict_route*"/fchBSuni.pdf")
###############################################################################

############################## FIGURE #########################################
# plotting treloar data for uniaxial test
p = scatter(λuni, Puni,
           label = latexstring("Test \\; data \\; uniaxial" ),
           markercolor = mycolors[1],
           markerstrokecolor = mycolors[1],
           markersize = mymarkersize,
           markerstrokewidth = mymarkerstrokewidth,
           markeralpha = 0.,
           xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
           ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
           size = (dim1,dim2)
           )
# plotting Nominal stress in uniaxial test computed with WYPIWYG from
#    f_λch obtained from uniaxial test
p = plot!(λuni_dense, sigma2P(λuni_dense,σ_w_uniax),
          label=latexstring("WYPiWYG \\; Spline \\; uniaxial" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1],
          )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from uniaxial test"))
display(p)
savefig(Pict_route*"/PLbsUnifUni.eps")
savefig(Pict_route*"/PLbsUnifUni.png")
savefig(Pict_route*"/PLbsUnifUni.pdf")
###############################################################################
# Computing σ_w_biax using f_λch from uniaxial test
σ_w_biax = sigma_wypiwyg(λbiax_dense,f_λch_uniax,2)

############################## FIGURE #########################################
# plotting treloar data for biaxial test
p = scatter(λbiax, Pbiax,
            label = latexstring("Test \\; data \\; biaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
            )
# plotting Nominal stress in biaxial test computed with WYPIWYG from
#    f_λch obtained from uniaxial test
p = plot!(λbiax_dense,sigma2P(λbiax_dense,σ_w_biax),
        label=latexstring("WYPiWYG \\; Spline \\; biaxial" ),
        linecolor = mycolors[2],
        linewidth = mylinewidth,
        linestyle = mylynestyles[1]
        )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from uniaxial test"))
display(p)
savefig(Pict_route*"/PLbsBiaxfUni.eps")
savefig(Pict_route*"/PLbsBiaxfUni.png")
savefig(Pict_route*"/PLbsBiaxfUni.pdf")
###############################################################################

# Computing σ_w_shear using f_λch from uniaxial test
σ_w_shear = sigma_wypiwyg(λshear_dense,f_λch_uniax,3)

############################## FIGURE #########################################
# plotting treloar data for shear test
p = scatter(λshear, Pshear,
            label = latexstring("Test \\; data \\; p.shear" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in shear test computed with WYPIWYG from
#    f_λch obtained from uniaxial test
p = plot!(λshear_dense,sigma2P(λshear_dense,σ_w_shear),
          label=latexstring("WYPiWYG \\; Spline \\; p.shear" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from uniaxial test"))
display(p)
savefig(Pict_route*"/PLbsPsfUni.eps")
savefig(Pict_route*"/PLbsPsfUni.png")
savefig(Pict_route*"/PLbsPsfUni.pdf")
###############################################################################


################################################################################
################### WYPIWYG using Bsplines and biaxial data ####################
################################################################################

# Computing f_λch from biaxial test , λu_block = 6.4
f_λch_biax, df_λch_uniax = ffromtests(λbiax, Pbiax, 6.4, 2, 2)
# Computing σ_w_biax using f_λch from biaxial test
σ_w_biax = sigma_wypiwyg(λbiax_dense,f_λch_biax,2)

λch_plotf_biax = [1.:3.62/599.:4.6 ...]

############################## FIGURE #########################################
# plotting f_λch obtained from biaxial test
p = plot(λch_plotf_biax,f_λch_biax.(λch_plotf_biax),
        label = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right)\$ from biaxial"),
        linecolor = mycolors[1],
        linewidth = mylinewidth,
        linestyle = mylynestyles[1],
        xlabel= latexstring("\$ \\lambda_{ch}\\;[-]\$"),
        ylabel = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;[MPa]\$"),
        size = (dim1,dim2)
        )

title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from biaxial test"))
display(p)
savefig(Pict_route*"/fchBSbiax.eps")
savefig(Pict_route*"/fchBSbiax.png")
savefig(Pict_route*"/fchBSbiax.pdf")
###############################################################################

############################## FIGURE #########################################
# plotting treloar data for biaxial test
p = scatter(λbiax, Pbiax,
            label = latexstring("Test \\; data  \\; biaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in biaxial test computed with WYPIWYG from
#    f_λch obtained from biaxial test
p = plot!(λbiax_dense, sigma2P(λbiax_dense,σ_w_biax),
          label=latexstring("WYPiWYG \\; Spline \\; biaxial" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from biaxial test"))
display(p)
savefig(Pict_route*"/PLbsBiaxfBiax.eps")
savefig(Pict_route*"/PLbsBiaxfBiax.png")
savefig(Pict_route*"/PLbsBiaxfBiax.pdf")
###############################################################################
# Computing σ_w_uniax using f_λch from biaxial test
σ_w_uniax = sigma_wypiwyg(λuni_dense,f_λch_biax,1)

############################## FIGURE #########################################
# plotting treloar data for uniaxial test
p = scatter(λuni, Puni,
            label = latexstring("Test \\; data  \\; uniaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in uniaxial test computed with WYPIWYG from
#    f_λch obtained from biaxial test
p = plot!(λuni_dense, sigma2P(λuni_dense,σ_w_uniax),
          label=latexstring("WYPiWYG \\; Spline \\; uniaxial" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from biaxial test"))
display(p)
savefig(Pict_route*"/PLbsUnifBiax.eps")
savefig(Pict_route*"/PLbsUnifBiax.png")
savefig(Pict_route*"/PLbsUnifBiax.pdf")
###############################################################################

# Computing σ_w_shear using f_λch from biaxial test
σ_w_shear = sigma_wypiwyg(λshear_dense,f_λch_biax,3)

############################## FIGURE #########################################
# plotting treloar data for shear test
p = scatter(λshear, Pshear,
            label = latexstring("Test \\; data  \\; p.shear" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in shear test computed with WYPIWYG from
#    f_λch obtained from biaxial test
p = plot!(λshear_dense, sigma2P(λshear_dense,σ_w_shear),
          label=latexstring("WYPiWYG \\; Spline \\; p.shear" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from biaxial test"))
display(p)
savefig(Pict_route*"/PLbsPsfBiax.eps")
savefig(Pict_route*"/PLbsPsfBiax.png")
savefig(Pict_route*"/PLbsPsfBiax.pdf")
###############################################################################

################################################################################
################# WYPIWYG using Bsplines and Pure Shear data ###################
################################################################################
# Computing f_λch from pure shear test, λu_block = 8.6
f_λch_shear, df_λch_shear = ffromtests(λshear, Pshear, 8.6, 3, 2)
# Computing σ_w_shear using f_λch from shear test
σ_w_shear = sigma_wypiwyg(λshear_dense,f_λch_shear,3)

λch_plotf_shear = [1.:3.62/599.:4.6 ...]

############################## FIGURE #########################################
# plotting f_λch obtained from shear test
p = plot(λch_plotf_shear, f_λch_shear.(λch_plotf_shear),
         label = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from p.shear"),
         linecolor = mycolors[1],
         linewidth = mylinewidth,
         linestyle = mylynestyles[1],
         xlabel= latexstring("\$ \\lambda_{ch}\\;[-]\$"),
         ylabel = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;[MPa]\$"),
         size = (dim1,dim2)
        )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from p.shear test"))
display(p)
savefig(Pict_route*"/fchBSsh.eps")
savefig(Pict_route*"/fchBSsh.png")
savefig(Pict_route*"/fchBSsh.pdf")
###############################################################################

############################## FIGURE #########################################
# plotting treloar data for shear test
 p = scatter(λshear, Pshear,
             label = latexstring("Test \\; data  \\; p.shear" ),
             markercolor = mycolors[1],
             markerstrokecolor = mycolors[1],
             markersize = mymarkersize,
             markerstrokewidth = mymarkerstrokewidth,
             markeralpha = 0.,
             xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
             ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
             size = (dim1,dim2)
            )
# plotting Nominal stress in shear test computed with WYPIWYG from
#    f_λch obtained from shear test
 p = plot!(λshear_dense, sigma2P(λshear_dense,σ_w_shear),
           label=latexstring("WYPiWYG \\; Spline \\; p.shear" ),
           linecolor = mycolors[2],
           linewidth = mylinewidth,
           linestyle = mylynestyles[1]
           )
 title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from p.shear test"))
 display(p)
 savefig(Pict_route*"/PLbsPsfSh.eps")
 savefig(Pict_route*"/PLbsPsfSh.png")
 savefig(Pict_route*"/PLbsPsfSh.pdf")
###############################################################################
# Computing σ_w_uniax using f_λch from shear test
σ_w_uniax = sigma_wypiwyg(λuni_dense,f_λch_shear,1)

############################## FIGURE #########################################
# plotting treloar data for uniaxial test
p = scatter(λuni, Puni,
            label = latexstring("Test \\; data  \\; uniaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in uniaxial test computed with WYPIWYG from
#    f_λch obtained from shear test
p = plot!(λuni_dense, sigma2P(λuni_dense,σ_w_uniax),
          label=latexstring("WYPiWYG \\; Spline \\; Uniaxial" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from p.shear test"))
display(p)
savefig(Pict_route*"/PLbsUnifSh.eps")
savefig(Pict_route*"/PLbsUnifSh.png")
savefig(Pict_route*"/PLbsUnifSh.pdf")
###############################################################################

# Computing σ_w_biax using f_λch from shear test
σ_w_biax = sigma_wypiwyg(λbiax_dense,f_λch_shear,2)

############################## FIGURE #########################################
# plotting treloar data for biaxial test
p = scatter(λbiax, Pbiax,
            label = latexstring("Test \\; data  \\; biaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in biaxial test computed with WYPIWYG from
#    f_λch obtained from shear test
p = plot!(λbiax_dense, sigma2P(λbiax_dense,σ_w_biax),
          label=latexstring("WYPiWYG \\; Spline \\; Uniaxial" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ from p.shear test"))
display(p)
savefig(Pict_route*"/PLbsBiaxfSh.eps")
savefig(Pict_route*"/PLbsBiaxfSh.png")
savefig(Pict_route*"/PLbsBiaxfSh.pdf")
###############################################################################

################################################################################
############# Comparing the Energy Functions from different tests ##############
################################################################################

############################## FIGURE #########################################
# plotting f_λch obtained from uniaxial test
p = plot(λch_plotf_uni,f_λch_uniax.(λch_plotf_uni),
         label = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;\$ from uniaxial"),
         linecolor = mycolors[1],
         linewidth = mylinewidth,
         linestyle = mylynestyles[1],
         xlabel= latexstring("\$\\lambda_{ch}\\;[-]\$"),
         ylabel = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;[MPa]\$"),
         size = (dim1,dim2)
         )
# plotting f_λch obtained from biaxial test
p = plot!(λch_plotf_biax,f_λch_biax.(λch_plotf_biax),
          label = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;\$ from biaxial"),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[2]
          )
# plotting f_λch obtained from shear test
p = plot!(λch_plotf_shear,f_λch_shear.(λch_plotf_shear),
          label = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;\$ from p.shear"),
          linecolor = mycolors[3],
          linewidth = mylinewidth,
          linestyle = mylynestyles[3]
         )

display(p)
savefig(Pict_route*"/allfchBS.eps")
savefig(Pict_route*"/allfchBS.png")
savefig(Pict_route*"/allfchBS.pdf")
###############################################################################

################################################################################
########################### Computin the avg f(λch) ############################
################################################################################

# Read treloar data
λuni, Puni = readCurve(Data_route*"treloar-uni.txt")
λbiax, Pbiax = readCurve(Data_route*"treloar-equibi.txt")
λshear, Pshear = readCurve(Data_route*"treloar-shear.txt")
np = 100
λuni_dense = [λuni[1]:(λuni[end]-λuni[1])/(convert(Float64,np)-1.):λuni[end] ...]
λbiax_dense = [λbiax[1]:(λbiax[end]-λbiax[1])/(convert(Float64,np)-1.):λbiax[end] ...]
λshear_dense = [λshear[1]:(λshear[end]-λshear[1])/(convert(Float64,np)-1.):λshear[end] ...]

# specify the λch_block for the average function
λch_block = 5.2

# Compute the point of the f(λch) from experimental data
λch_uni, fch_uni = ffromtests_points(λuni, Puni, 1 , 2)
λch_biax, fch_biax = ffromtests_points(λbiax, Pbiax, 2 , 2)
λch_shear, fch_shear = ffromtests_points(λshear, Pshear, 3 , 2)


λch_glb = vcat(λch_uni, λch_biax, λch_shear)
fch_glb = vcat(fch_uni, fch_biax, fch_shear)


p = sortperm(λch_glb)


λch_glb_ordered  = λch_glb[p]
fch_glb_ordered = fch_glb[p]

# Compute the average f(λch) from the points cloud (points of the f(λch) computed from diff tests)
f_avg, df_avg = rational_joint(λch_glb_ordered, fch_glb_ordered , λch_glb_ordered[end] , λch_block , 3, 1)

λch_plot = [1.:3.62/999.:5...]
############################## FIGURE #########################################
# plotting f(λch) avg
p = plot(λch_plot, f_avg.(λch_plot),
        label = latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ avg"),
        linecolor = :green,
        size = (dim1,dim2),
        linewidth = mylinewidth,
        linestyle = :solid,
        xlabel= latexstring("\$\\lambda_{ch}\\;[-]\$"),
        ylabel = latexstring("\$ f\\;\\left( { \\lambda  }_{ ch } \\right) \\;[MPa]\$")
        )
# plotting points for the cloud f(λch) from uniaxial test
scatter!(λch_uni, fch_uni,
        label = latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ points from uniaxial"),
        markershape = mymarkershapes[1],
        markercolor = mycolors[1],
        markerstrokecolor = mycolors[1],
        markerstrokewidth = mymarkerstrokewidth,
        markersize = mymarkersize,
        markeralpha = 0.
        )
# plotting points for the cloud f(λch) from biaxial test
scatter!(λch_biax, fch_biax,
        label = latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ points from biaxial"),
        markershape = mymarkershapes[2],
        markercolor = mycolors[2],
        markerstrokecolor = mycolors[2],
        markerstrokewidth = mymarkerstrokewidth,
        markersize = mymarkersize,
        markeralpha = 0.
        )
# plotting points for the cloud f(λch) from shear test
scatter!(λch_shear, fch_shear,
        label = latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ points from p.shear"),
        markershape = mymarkershapes[3],
        markercolor = mycolors[3],
        markerstrokecolor = mycolors[3],
        markerstrokewidth = mymarkerstrokewidth,
        markersize = mymarkersize,
        markeralpha = 0.
        )

title!(latexstring("\$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ points cloud"))
display(p)
savefig(Pict_route*"/fchfrom3tests.eps")
savefig(Pict_route*"/fchfrom3tests.png")
savefig(Pict_route*"/fchfrom3tests.pdf")
###############################################################################

σ_w_uniax_avg = sigma_wypiwyg(λuni_dense,f_avg,1)
σ_w_biax_avg = sigma_wypiwyg(λbiax_dense,f_avg,2)
σ_w_shear_avg = sigma_wypiwyg(λshear_dense,f_avg,3)

############################## FIGURE #########################################
# plotting treloar data for uniaxial test
p = scatter(λuni, Puni,
            label = latexstring("Test \\; data  \\; uniaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in uniaxial test computed with WYPIWYG from
#    f_λch obtained with avg
p = plot!(λuni_dense,sigma2P(λuni_dense,σ_w_uniax_avg),
         label=latexstring("WYPiWYG \\; Spline \\; uniaxial" ),
         linecolor = mycolors[2],
         linewidth = mylinewidth,
         linestyle = mylynestyles[1]
        )
title!(latexstring(" Curves obtained from \$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ avg"))
display(p)
savefig(Pict_route*"/PLuniavg.eps")
savefig(Pict_route*"/PLuniavg.png")
savefig(Pict_route*"/PLuniavg.pdf")
###############################################################################

############################## FIGURE ########################################
# plotting treloar data for biaxial test#
p = scatter(λbiax, Pbiax,
            label = latexstring("Test \\; data  \\; biaxial" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in biaxial test computed with WYPIWYG from
#    f_λch obtained with avg
p = plot!(λbiax_dense, sigma2P(λbiax_dense,σ_w_biax_avg),
          label=latexstring("WYPiWYG \\; Spline \\; biaxial" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
         )
title!(latexstring(" Curves obtained from \$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ avg"))
display(p)
savefig(Pict_route*"/PLbiaxavg.eps")
savefig(Pict_route*"/PLbiaxavg.png")
savefig(Pict_route*"/PLbiaxavg.pdf")
###############################################################################

############################## FIGURE #########################################
# plotting treloar data for shear test
p = scatter(λshear, Pshear,
            label = latexstring("Test \\; data  \\; p.shear" ),
            markercolor = mycolors[1],
            markerstrokecolor = mycolors[1],
            markersize = mymarkersize,
            markerstrokewidth = mymarkerstrokewidth,
            markeralpha = 0.,
            xlabel= latexstring("Stretch: \$\\lambda\\;[-]\$"),
            ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
            size = (dim1,dim2)
           )
# plotting Nominal stress in shear test computed with WYPIWYG from
#    f_λch obtained with avg
p = plot!(λshear_dense, sigma2P(λshear_dense,σ_w_shear_avg),
          label=latexstring("WYPiWYG \\; Spline \\; p.shear" ),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[1]
)
title!(latexstring(" Curves obtained from \$f\\;\\left( { \\lambda  }_{ ch } \\right) \$ avg"))
display(p)
savefig(Pict_route*"/PLshavg.eps")
savefig(Pict_route*"/PLshavg.png")
savefig(Pict_route*"/PLshavg.pdf")
###############################################################################
