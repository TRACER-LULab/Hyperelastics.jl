"""
###################################
#      Module WYPIWYGabTools      #
###################################

This module contains the following functions:

1. ffromtests_points
2. ffromtests
3. ffromModel
4. arrboy
5. sigma2P
6. sigma_wypiwyg
"""

module WYPIWYGabTools  #  begin of the module WYPIWYGabTools

include("./MathematicsToolsWYPIWYG.jl")
using .MathematicsToolsWYPIWYG
######################################################################
##################### Functions that are exported ####################
######################################################################
export ffromtests
export ffromtests_points
export ffromModel
export arrBoy
export sigma2P
export sigma_wypiwyg
######################################################################
######################## Functions definition ########################
######################################################################

##############
#     1      #
#########

######################### begin: ffromtests_points ###################

"""
ffromtests_points(λu, σP, flagTest, flagStress = 1)-> λch, fλch

Computes the points of f(λch) given experimental the stretches(λu) and its
corresponding stresses(σP).

Input:
    λu::Array{Float64,1} = stretches test points
    σP::Array{Float64,1} = stresses test points corresponding to λu
    flagTest::Int64 = specifies which test are the points obtained from
        flagTest == 1 Uniaxial test
        flagTest == 2 Biaxial test
        flagTest == 3 Pure Shear test

Outputs:
    λch::Array{Float64,1} : λch corresponding to λu
    fch::Array{Float64,1} : fch points computed from σP and λu
Optionals:

    flagStress::Int64 = it specifies if the stresses are:
        flagStress == 1 Cauchy stress(σ) (Default)
        flagStress == 2 Nomial stress(P)

# Examples
```
julia>  λu = [1. 1.2 1.4 1.8]; σP = [0.0 0.19 0.30 0.5]
julia>  λch, fch = ffromtests_points(λu, σP, flagTest, flagStress = 1)
```
"""
function ffromtests_points(λu::Array{Float64,1}, σP::Array{Float64,1},
                           flagTest::Int64, flagStress::Int64 = 1)
     λu2 = λu .^2
     iλu = 1. ./ λu

if flagTest == 1 # Uniaxial test

     λch = sqrt.((λu2 + 2. *iλu)/3.)

      # Computes fλch deending on the  flagStress

     if flagStress == 1 # σ_Pu Cauchy stress                              DEFAULT
      fλch = σP ./ (λu2 - iλu);  quitNaN!(fλch)
     elseif flagStress == 2 # σ_Pu Nominal (1st Piola kirchoff) stress
      fλch = σP ./ (λu - iλu.^2);  quitNaN!(fλch)
     else
      error("flagStress no válido en ffromtest in uniax")
     end

  # THAT'S IT!!
elseif flagTest == 2 # biaxial test

    λch = sqrt.((2. * λu2 + iλu.^4)/3.)


    if flagStress == 1 # The input is cauchy stress
        fλch =  σP ./ (λu2 - iλu.^4)  ;  quitNaN!(fλch)

    elseif flagStress == 2 # The input is nominal stress
        fλch =  σP ./ (λu - iλu.^5)  ;  quitNaN!(fλch)
    else
        error("flagStress no válido en ffrombiax")
    end

elseif flagTest == 3 # Pure shear

    λch = sqrt.((λu2 + iλu.^2 .+ 1.0)/3.)

    if flagStress == 1 # The input is cauchy stress
        fλch =  σP ./ (λu2 - iλu.^2)  ;  quitNaN!(fλch)

    elseif flagStress == 2 # The input is nominal stress
        fλch =  σP ./ (λu - iλu.^3)  ;  quitNaN!(fλch)
    else
     error("flagStress no válido en ffromtest in shear")
    end

else
    error("Incorrect flag in ffromtests")
end



return λch::Array{Float64,1}, fλch::Array{Float64,1}
end


######################### end:  ffromtests_points ####################

##############
#     2      #
##############

############################ begin: ffromtests #######################
"""
ffromtests_points(λu, σP, λu_block, flagTest, flagStress = 1)
                                ->f::Function, df::Function
Computes a function and its derivative which are a union between
    a B-spline and a rational function with asymptote at λu_block.


Input:
    λu::Array{Float64,1} = stretches test points
    λu_block::Array{Float64,1} = test stretch at which the user consider
        the chains are blocked (this value is obtained analyzing the test curve)
    σP::Array{Float64,1} = stresses test points corresponding to λu
    flagTest::Int64 = specifies which test are the points obtained from
        flagTest == 1 Uniaxial test
        flagTest == 2 Biaxial test
        flagTest == 3 Pure Shear test

Outputs:
    f::Array{Float64,1} : function returned
    df::Array{Float64,1} : function derivative returned
Optionals:

    flagStress::Int64 = it specifies if the stresses are:
        flagStress == 1 Cauchy stress(σ) (Default)
        flagStress == 2 Nomial stress(P)

# Examples
```
julia>  λu = [1. 1.2 1.4 1.8]; σP = [0.0 0.19 0.30 0.5]
julia>  λch, fch = ffromtests_points(λu, σP, flagTest, flagStress = 1)
```
"""

function ffromtests(λu::Array{Float64,1}, σP::Array{Float64,1},
                    λu_block::Float64, flagTest::Int64, flagStress::Int64 = 1)


if flagTest == 1
    λch_block = sqrt( (λu_block^2 + 2. * (1. / λu_block ) ) / 3.)
elseif flagTest == 2
    λch_block = sqrt((2. * λu_block^2 + (1. / λu_block )^4)/3.)
elseif flagTest == 3
    λch_block = sqrt.((λu_block^2 + (1. / λu_block)^2 + 1.)/3.)
else
    error("flagTest incorrect in fromtests")
end

    λch, fλch = ffromtests_points(λu, σP, flagTest, flagStress)


    f, df = rational_joint(λch, fλch, λch[end], λch_block, 3, 1)
    return f::Function, df::Function
end


############################# end:  ffromtests #######################


##############
#     3      #
##############

############################ begin: ffromModel #######################
"""
ffromModel(λu, σP, flagTest, flagStress = 1)->f::Function

Computes the cubic spline fλch(λch) for that, it first computes λch, fλch
points from λu and σP from a particular test.


Input:
    λu::Array{Float64,1} = stretches test points
    σP::Array{Float64,1} = stresses test points corresponding to λu
    flagTest::Int64 = specifies which test are the points obtained from
        flagTest == 1 Uniaxial test
        flagTest == 2 Biaxial test
        flagTest == 3 Pure Shear test

Outputs:
    f::Array{Float64,1} : cubic spline function returned

Optionals:

    flagStress::Int64 = it specifies if the stresses are:
        flagStress == 1 Cauchy stress(σ) (Default)
        flagStress == 2 Nomial stress(P)

# Examples
```
julia>  λu = [1. 1.2 1.4 1.8]; σP = [0.0 0.19 0.30 0.5]
julia>   f = ffromModel(λu, σP, flagTest, flagStress = 1)
```
"""

function ffromModel(λu::Array{Float64,1}, σP::Array{Float64,1},
                    flagTest::Int64, flagStress::Int64 = 1)

    λch, fλch = ffromtests_points(λu, σP, flagTest, flagStress)

    f = csplinef(λch, fλch)
    return f::Function
end


############################ end:  ffromModel ###################


##############
#     4      #
##############

############################## begin: arrBoy #########################

"""
arrBoy(G, Nch, ILF, λu, flagTest = 1)->σuAB
Uses the 8-chain Model to compute cauchy stress correcponding to
 a vector of stretches.

 Input:
    G, Nch::Float64 : Material parameters
    ILF::Function : Inverse Langevin function
    λu::Array{Float64,1} : streteches vector which define where the
     function is goint to determine the stress
    flagTest::Int64 : Indicates which test we want to compute with
        the 8-model

        flagTest == 1 Uniaxial test
        flagTest == 2 Biaxial test
        flagTest == 3 Pure Shear test

    Output:

        σuAB::Array{Float64,1} : computed cauchy stresses with 8-chain model
```
"""

function arrBoy(G::Float64, Nch::Float64, ILF::Function,
               λu::Array{Float64,1}, flagTest::Int64 = 1)

 λu2 = λu.^2               # Squared stretches (used repeatedly below)
 iλu = (1.)./ λu

# Compute λch

 # Computes fλch deending on the  flagStress

 if flagTest == 1 # Uniaxial test                             DEFAULT
   λch = sqrt.((λu2 + (2.)*iλu)/(3.));  λchN = λch./sqrt(Nch)
   σuAB = G/(3.)*sqrt(Nch) ./ λch .* ILF.(λchN) .* (λu2 - iλu)

 elseif flagTest == 2 # Biaxial test
  λch = sqrt.(((2.)*λu2 + iλu.^4)/(3.));  λchN = λch./sqrt(Nch)
  σuAB = G/(3.)*sqrt(Nch) ./ λch .* ILF.(λchN) .* (λu2 - iλu.^4)


 elseif flagTest == 3 # Pure shear test

  λch = sqrt.((λu2 + iλu.^2 .+ 1.0)/(3.));  λchN = λch./sqrt(Nch)
  σuAB = G/(3.)*sqrt(Nch) ./ λch .* ILF.(λchN) .* (λu2 - iλu.^2)

 else
  error("flagTest no válido en ffromuniax")
 end


return σuAB::Array{Float64,1}
end
################################# end: arrBoy ########################


##############
#     5      #
##############

############################# begin: sigma2P #########################
"""
sigma2P(λu::Array{Float64,1}, σu::Array{Float64,1})->P::{Float64,1}

Computes nominal stress from the stretches and the corresponding
cauchy stress.

# Examples
```
julia>  λu = [1. 1.2 1.4 1.8]; σP = [0.0 0.19 0.30 0.5]
julia>  P = sigma2P(λu, σu)
```
"""
function sigma2P(λu::Array{Float64,1}, σu::Array{Float64,1})

P = σu ./ λu

return P::Array{Float64,1}
end

######################### end: sigma2P ###############################


##############
#     6      #
##############

########################## begin: sigma_wypiwyg ######################
"""
sigma_wypiwyg(λu::Array{Float64,1}, f::Function, flagTest::Int64 = 1)->σw
Given the test stretches in which we want to compute the cauchy stress
and the function f(λch) it computes the cauchy stresses for a certain test

 Input:
     λu::Array{Float64,1} : streteches vector which define where the
     function is goint to determine the stress
     f::Function : f(λch) which is unique for a certain material
    flagTest::INt64 : indicates the test in which we want to compute the
        the σ.
        flagTest == 1 Uniaxial test
        flagTest == 2 Biaxial test
        flagTest == 3 Pure Shear test

    Output:

        σw::Array{Float64,1} : computed cauchy stresses with WYPIWYG model
```
"""

function sigma_wypiwyg(λu::Array{Float64,1}, f::Function, flagTest::Int64 = 1)

  λu2 = λu  .^2
  iλu = 1. ./ λu

 if (flagTest == 1)  # Uniaxial stress

  λch_uniax = sqrt.((λu2 + (2.)*iλu)/(3.))
  σ_w = f.(λch_uniax).*(λu2 - iλu)

 elseif (flagTest == 2)  # Biaxial stress

  λch_biax = sqrt.(((2.)*λu2 + iλu.^4)/(3.))
  σ_w = f.(λch_biax).*(λu2 - iλu.^4)

 elseif (flagTest == 3)  # shear stress
  λch_ps = sqrt.((λu2 + iλu.^2 .+ 1.0)/(3.))
  σ_w = f.(λch_ps).*(λu2 - iλu.^2)

 else
  error("flagTest no válido en sigma_wypiwyg")
 end

return σ_w::Array{Float64,1}
end

###################### end:  sigma_wypiwyg ###########################


end  #  end of the module WYPIWYGabTools
