push!(LOAD_PATH, ".")
using MathematicsToolsWYPIWYG
using WYPIWYGabTools

using Plots
using LaTeXStrings
using LinearAlgebra


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
mymarkersize1 = 11.
mymarkersize2 = 7.
mymarkerstrokewidth = 4.5
mylynestyles = [:solid :dash :dot]
mymarkershapes = [:circle :rect :diamond ]

# Pictures save route

Pict_route = "./Pictures/arruda_wypiwyg"   # default directory

###########################################################################

## Inverse Langevin functions

# Create Inverse Langevin function
ILF = createIL()


# Checks for Langevin and ILF calculations


# Arruda-Boyce parameters  N and G=nkT
Nch = 26.5;  G = 0.27

########################### Energy from uniaxial test ##########################

    # 8 Chain Model for σ computing
        λmax_uniax = 8.0
        λu_uniax = [1.:(λmax_uniax-1.)/99.:λmax_uniax ...]
        σ_ab_uniax = arrBoy(G, Nch, ILF, λu_uniax) # :uniaxial test default

    # WYPIWYG using previously computed σ_ab_uniax and λmax_uni
        f_λch_uniax = ffromModel(λu_uniax, σ_ab_uniax, 1, 1) # 1:uniaxial test, 1: cauchy stress
        σ_w_uniax = sigma_wypiwyg(λu_uniax,f_λch_uniax) # :uniaxial test default


########################### Figure ####################################
    plot_gap = 3
    p = scatter(λu_uniax[1:plot_gap:end],
                sigma2P(λu_uniax[1:plot_gap:end], σ_ab_uniax[1:plot_gap:end]),
                linecolor = mycolors[1],
                markersize = mymarkersize1,
                markerstrokecolor = mycolors[1],
                markerstrokewidth = mymarkerstrokewidth,
                markeralpha = 0.,
                label = latexstring("\\mathrm{8-Chain \\; model}"),
                xlabel = latexstring("Stretch: \$\\lambda\\;[-]\$"),
                ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
                size = (dim1,dim2)
               )

    p = scatter!(λu_uniax[1:plot_gap:end],
                 sigma2P(λu_uniax[1:plot_gap:end],σ_w_uniax[1:plot_gap:end]),
                 marker=:circle,
                 markercolor = mycolors[2],
                 markerstrokecolor = mycolors[2],
                 markersize = mymarkersize2,
                 markeralpha = 1.,
                 label=latexstring("WYPiWYG")
                )

    title!(latexstring("Uniaxial \\; test}"))
    display(p)
savefig(Pict_route*"/uniax1.eps")
savefig(Pict_route*"/uniax1.pdf")
savefig(Pict_route*"/uniax1.png")
#####################################################################

######################### Energy from equibiaxial test #########################

    # 8 Chain Model for σ computing
        λmax_biax = 6.0
        λu_biax = [1.:(λmax_biax-1.)/99.:λmax_biax...]; # Sampling uniaxial stretches
        σ_ab_biax = arrBoy(G, Nch, ILF, λu_biax, 2)  # 2 indicates biaxial test

    # WYPIWYG using previously computed σ_ab_biax and λmax_biax
        f_λch_biax = ffromModel(λu_biax, σ_ab_biax, 2, 1) # 2:biaxial test, 1: cauchy stress
        σ_w_biax = sigma_wypiwyg(λu_biax, f_λch_biax, 2) # 2:biaxial

########################### Figure ####################################
    p = scatter(λu_biax[1:plot_gap:end],
                sigma2P(λu_biax[1:plot_gap:end],σ_ab_biax[1:plot_gap:end]),
                linecolor = mycolors[1],
                markersize = mymarkersize1,
                markerstrokecolor = mycolors[1],
                markerstrokewidth = mymarkerstrokewidth,
                markeralpha = 0.,
                label = latexstring("\\mathrm{8-Chain \\; model}"),
                xlabel = latexstring("Stretch: \$\\lambda\\;[-]\$"),
                ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
                size = (dim1,dim2)
             )

    p = scatter!(λu_biax[1:plot_gap:end],
                 sigma2P(λu_biax[1:plot_gap:end],σ_w_biax[1:plot_gap:end]),
                 marker=:circle,
                 markercolor = mycolors[2],
                 markerstrokecolor = mycolors[2],
                 markersize = mymarkersize2,
                 markeralpha = 1.,
                 label=latexstring("WYPiWYG")
                 )

    title!(latexstring("Equibiaxial \\; test}"))
    display(p)
savefig(Pict_route*"/equib1.eps")
savefig(Pict_route*"/equib1.pdf")
savefig(Pict_route*"/equib1.png")
#####################################################################


###################### Energy from pure shear test #############################

    # 8 Chain Model for σ computing
        λmax_ps = 8.0
        λu_ps = [1.:(λmax_ps-1.)/99.:λmax_ps...];      # Sampling uniaxial stretches
        σ_ab_ps = arrBoy(G, Nch, ILF, λu_ps, 3) # 3 indicates pure shear test
    # WYPIWYG using previously computed σ_ab_biax and λmax_biax
        f_λch_ps = ffromModel(λu_ps, σ_ab_ps, 3, 1) # 3: ps test, 1: cauchy stress
        σ_w_ps = sigma_wypiwyg(λu_ps, f_λch_ps, 3) # 3: ps test

    p = scatter(λu_ps[1:plot_gap:end],
             sigma2P(λu_ps[1:plot_gap:end],σ_ab_ps[1:plot_gap:end]),
             linecolor = mycolors[1],
             markersize = mymarkersize1,
             markerstrokecolor = mycolors[1],
             markerstrokewidth = mymarkerstrokewidth,
             markeralpha = 0.,
             label = latexstring("\\mathrm{8-Chain \\; model}"),
             xlabel = latexstring("Stretch: \$\\lambda\\;[-]\$"),
             ylabel = latexstring("Nominal stress: \$ P \\;[MPa]\$ "),
             size = (dim1,dim2)
             )

    p = scatter!(λu_ps[1:plot_gap:end],
                 sigma2P(λu_ps[1:plot_gap:end],σ_w_ps[1:plot_gap:end]),
                 marker=:circle,
                 markercolor = mycolors[2],
                 markerstrokecolor = mycolors[2],
                 markersize = mymarkersize2,
                 markeralpha = 1.,
                 label=latexstring("WYPiWYG")
                 )

    title!(latexstring("Pure \\; Shear \\; test}"))
    display(p)

savefig(Pict_route*"/push1.pdf")
savefig(Pict_route*"/push1.eps")
savefig(Pict_route*"/push1.png")
######################## fλch obtained from different test #####################



λch_plotf = [1.:0.1:4.6 ...]


########################### Figure ####################################
plot_gap = 2
p = scatter(λch_plotf[1:plot_gap:end],
         f_λch_uniax.(λch_plotf[1:plot_gap:end]),
         linecolor = mycolors[1],
         linewidth = mylinewidth,
         linestyle = mylynestyles[1],
         markershape = mymarkershapes[1],
         markersize = 10.5,
         markerstrokecolor = mycolors[1],
         markeralpha = 0.,
         markerstrokewidth = mymarkerstrokewidth,
         label = latexstring("Uniaxial"),
         xlabel= latexstring("\$ \\lambda_{ch}\\;[-]\$"),
         ylabel= latexstring("\$ f\\;(λ_{ch})\\;[MPa]\$"),
         size = (dim1,dim2)
         )

p = scatter!(λch_plotf[1:plot_gap:end],
          f_λch_biax.(λch_plotf[1:plot_gap:end]),
          linecolor = mycolors[2],
          linewidth = mylinewidth,
          linestyle = mylynestyles[2],
          markerstrokecolor = mycolors[2],
          markersize = 6.4,
          markershape = mymarkershapes[2],
          markeralpha = 0.,
          markerstrokewidth = mymarkerstrokewidth,
          label = latexstring("Biaxial"),
          )


p = scatter!(λch_plotf[1:plot_gap:end],
          f_λch_ps.(λch_plotf[1:plot_gap:end]),
          linecolor = mycolors[3],
          linewidth = mylinewidth,
          linestyle = mylynestyles[3],
          markershape = mymarkershapes[3],
          markersize = 4.5,
          markercolor = mycolors[3],
          markerstrokecolor = mycolors[3],
          markeralpha = 0.,
          markerstrokewidth = mymarkerstrokewidth,
          label = latexstring("Pure \\; shear"),
          )

display(p)
savefig(Pict_route*"/fch1i.eps")
savefig(Pict_route*"/fch1i.png")
savefig(Pict_route*"/fch1i.pdf")
#########################################################################


########################### Figure ####################################
p = scatter(λch_plotf[1:plot_gap:end],
         f_λch_uniax.(λch_plotf[1:plot_gap:end]),
         linecolor = mycolors[1],
         linestyle = mylynestyles[1],
         markershape = mymarkershapes[1],
         markersize = 9.,
         markerstrokewidth = 2.,
         markerstrokecolor = mycolors[1],
         markeralpha = 0.,
         label = latexstring("Discrete \\; data"),
         xlabel= latexstring("\$ \\lambda_{ch}\\;[-]\$"),
         ylabel=latexstring("\$ f\\;(λ_{ch})\\;[MPa]\$"),
         size = (dim1,dim2)
         )

display(p)

p = plot!(λch_plotf,f_λch_uniax.(λch_plotf),
         linecolor = :green,
         linewidth = mylinewidth,
         linestyle = mylynestyles[1],
         label = latexstring("Spline \\; curve"),
         xlabel= latexstring("\$ \\lambda_{ch}\\;[-]\$"),
         size = (dim1,dim2)
         )

         display(p)
savefig(Pict_route*"/fch1ii.eps")
savefig(Pict_route*"/fch1ii.pdf")
savefig(Pict_route*"/fch1ii.png")
#########################################################################

################################ General Load Case ############################

# Generic deformation gradient

#=

F = [
            (1+γ)^n         γ^n          0.
               0.         (1+γ)^n        0.
               0.            0.        (1+γ)^(-2n)
    ]


=#

γ0 = -0.5
γf = 1.5
γ_vector = [prevfloat(γ0):(γf-γ0)/(19.):nextfloat(γf)...]

        σ11_w = zeros(Float64,length(γ_vector))
        σ12_w = zeros(Float64,length(γ_vector))
        σ13_w = zeros(Float64,length(γ_vector))
        σ21_w = zeros(Float64,length(γ_vector))
        σ22_w = zeros(Float64,length(γ_vector))
        σ23_w = zeros(Float64,length(γ_vector))
        σ31_w = zeros(Float64,length(γ_vector))
        σ32_w = zeros(Float64,length(γ_vector))
        σ33_w = zeros(Float64,length(γ_vector))

        σ11_ab = zeros(Float64,length(γ_vector))
        σ12_ab = zeros(Float64,length(γ_vector))
        σ13_ab = zeros(Float64,length(γ_vector))
        σ21_ab = zeros(Float64,length(γ_vector))
        σ22_ab = zeros(Float64,length(γ_vector))
        σ23_ab = zeros(Float64,length(γ_vector))
        σ31_ab = zeros(Float64,length(γ_vector))
        σ32_ab = zeros(Float64,length(γ_vector))
        σ33_ab = zeros(Float64,length(γ_vector))

        n = 1


for i = 1:length(γ_vector)
        γ = γ_vector[i]

        # Computing the eigenvalues

        λ2_1 = (γ+1.)^(2*n) +
               ((1.)/(2.))*γ^(2*n) +
               (((1.)/(2.))*γ^n)*sqrt((4.)*((γ+1.)^(2*n))+γ^(2*n))


        λ2_2 = (γ+1.)^(2*n) +
              ((1.)/(2.))*γ^(2*n) -
              (((1.)/(2.))*γ^n)*sqrt((4.)*((γ+1.)^(2*n))+γ^(2*n))


        λ2_3 = (γ+1.)^(-4*n)

        λi2 = [λ2_1, λ2_2, λ2_3]

        num1 = γ^n+sqrt((4.)*(γ+1.)^(2*n)+γ^(2*n))
        num2 = γ^n-sqrt((4.)*(γ+1.)^(2*n)+γ^(2*n))
        den = (2.)*(γ+1.)^n

        t1 = num1/den
        t2 = num2/den

        # Computing the eigenvectors
        v1 = [t1, 1. , 0.]
        v2 = [t2, 1. , 0.]
        v3 = [0., 0., 1. ]

        # Normalized eigenvectors

        u1 = v1/norm(v1)
        u2 = v2/norm(v2)
        u3 = v3

        R = hcat(u1, u2, u3)



        λch_γ = sqrt((λi2[1]+λi2[2]+λi2[3])/(3.))
        λch_γN = λch_γ/sqrt(Nch)


##################### Cauchy Stress computed with WYPIWYG ######################
        σ1 = G/(3.)*sqrt(Nch)/λch_γ*ILF(λch_γN)*(λi2[1]-λi2[3])
        σ2 = G/(3.)*sqrt(Nch)/λch_γ*ILF(λch_γN)*(λi2[2]-λi2[3])



        σvect = [σ1, σ2, 0. ]
        σ_ppal = Diagonal(σvect)
        σnoppal = R*σ_ppal*R'

        σ11_ab[i] = σnoppal[1,1]
        σ12_ab[i] = σnoppal[1,2]
        σ13_ab[i] = σnoppal[1,3]
        σ21_ab[i] = σnoppal[2,1]
        σ22_ab[i] = σnoppal[2,2]
        σ23_ab[i] = σnoppal[2,3]
        σ31_ab[i] = σnoppal[3,1]
        σ32_ab[i] = σnoppal[3,2]
        σ33_ab[i] = σnoppal[3,3]

##################### Cauchy Stress computed with WYPIWYG ######################
        σ1 = (λi2[1]-λi2[3])*f_λch_uniax(λch_γ)
        σ2 = (λi2[2]-λi2[3])*f_λch_uniax(λch_γ)

        σvect = [σ1, σ2, 0. ]
        σ_ppal = Diagonal(σvect)
        σnoppal = R*σ_ppal*R'

        σ11_w[i] = σnoppal[1,1]
        σ12_w[i] = σnoppal[1,2]
        σ13_w[i] = σnoppal[1,3]
        σ21_w[i] = σnoppal[2,1]
        σ22_w[i] = σnoppal[2,2]
        σ23_w[i] = σnoppal[2,3]
        σ31_w[i] = σnoppal[3,1]
        σ32_w[i] = σnoppal[3,2]
        σ33_w[i] = σnoppal[3,3]

end # end for


########################### Figure ####################################


p = scatter(γ_vector, σ11_ab,
         linecolor = mycolors[1],
         linewidth = mylinewidth,
         markersize = mymarkersize1,
         markerstrokecolor = mycolors[1],
         markerstrokewidth = mymarkerstrokewidth,
         markeralpha = 0.,
         label = latexstring("\\mathrm{8-Chain \\; model}"),
         xlabel = latexstring("Stretch: \$\\gamma\\;[-]\$"),
         ylabel = latexstring("Cauchy stress: \$ \\sigma_{11} \\;[MPa]\$ "),
         size = (dim1,dim2)
         )

p = scatter!(γ_vector, σ11_w,
             marker=:circle,
             markercolor = mycolors[2],
             markerstrokecolor = mycolors[2],
             markersize = mymarkersize2,
             markeralpha = 1.,
             label=latexstring("WYPiWYG")
             )

title!(latexstring("General \\; load \\; case}"))
display(p)
savefig(Pict_route*"/cauchy11.png")
savefig(Pict_route*"/cauchy11.pdf")
savefig(Pict_route*"/cauchy11.eps")
#######################################################################


########################### Figure ####################################


p = scatter(γ_vector, σ22_ab,
         linecolor = mycolors[1],
         linewidth = mylinewidth,
         markersize = mymarkersize1,
         markerstrokecolor = mycolors[1],
         markerstrokewidth = mymarkerstrokewidth,
         markeralpha = 0.,
         label = latexstring("\\mathrm{8-Chain \\; model}"),
         xlabel = latexstring("Stretch: \$\\gamma\\;[-]\$"),
         ylabel = latexstring("Cauchy stress: \$ \\sigma_{22} \\;[MPa]\$ "),
         size = (dim1,dim2)
         )

p = scatter!(γ_vector, σ22_w,
             marker=:circle,
             markercolor = mycolors[2],
             markerstrokecolor = mycolors[2],
             markersize = mymarkersize2,
             markeralpha = 1.,
             label=latexstring("WYPiWYG")
             )

title!(latexstring("General \\; load \\; case}"))
display(p)

savefig(Pict_route*"/cauchy22.png")
savefig(Pict_route*"/cauchy22.pdf")
savefig(Pict_route*"/cauchy22.eps")
#######################################################################


########################### Figure ####################################


p = scatter(γ_vector, σ12_ab,
         linecolor = mycolors[1],
         linewidth = mylinewidth,
         markersize = mymarkersize1,
         markerstrokecolor = mycolors[1],
         markerstrokewidth = mymarkerstrokewidth,
         markeralpha = 0.,
         label = latexstring("\\mathrm{8-Chain \\; model}"),
         xlabel = latexstring("Stretch: \$\\gamma\\;[-]\$"),
         ylabel = latexstring("Cauchy stress: \$ \\sigma_{12} \\;[MPa]\$ "),
         size = (dim1,dim2)
         )

p = scatter!(γ_vector, σ12_w,
             marker=:circle,
             markercolor = mycolors[2],
             markerstrokecolor = mycolors[2],
             markersize = mymarkersize2,
             markeralpha = 1.,
             label=latexstring("WYPiWYG")
             )

title!(latexstring("General \\; load \\; case}"))
display(p)

savefig(Pict_route*"/cauchy12.pdf")
savefig(Pict_route*"/cauchy12.png")
savefig(Pict_route*"/cauchy12.eps")
#######################################################################
