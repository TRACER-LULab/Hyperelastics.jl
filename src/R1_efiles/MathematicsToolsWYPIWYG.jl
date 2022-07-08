"""
###################################
# Module MathematicsToolsWYPIWYG #
###################################

This module contains the following functions:

1. csplinef
2. createIL()
3. quitNaN!
4. readCurve
5. bsplinefdf
6. rational
7. rational_joint
"""

module MathematicsToolsWYPIWYG  #  begin of the module MathematicsToolsWYPIWYG

################################################################################
########################## Functions that are exported #########################
################################################################################

export csplinef
export createIL
export quitNaN!
export readCurve
export bsplinefdf
export rational_joint
export rational

################################################################################
############################# Functions definition #############################
################################################################################

##############
#     1      #
##############

################################### begin: csplinef ##############################
"""
    csplinef(x, y, EndCond1, EC1, EndCond2, EC2) -> csplinefval::Function

Creates an interpolating cubic spline function with given end conditions

(N-1) Polynomials are a + b x  + c x^2 + d x^3  with N = length(x)

Input:
    x::Array{Float64,1} = Abscisae (increasing), length N
    y::Array{Float64,1} = Ordinates corresponding to x

Optionals (defaults to natural splines)
    Endcond1::Int64     = order of derivative prescribed (1,2)
    EC1::Float64        = Value of end condition at start
    Endcond2::Int64     = order of derivative prescribed (1,2)
    EC2::Float64        = Value of end condition at end

Returns the function name which evaluates the spline

# Examples
```
julia> x = [1. 2. 3. 4. 5.]; y = [-1. 1. 2. -1. 0.5]
julia> f = csplinef(x, y, EndCond1, EC1, EndCond2, EC2)
julia> xvalue = [1.,0.1,5. ...]
julia> plot(x,y); plot!(xvalue,f.(xvalue))
```

"""
function csplinef(x::Array{Float64,1}, y::Array{Float64,1},
                  EndCond1::Int64 = 2, EC1::Float64 = 0.,
                  EndCond2::Int64 = 2,EC2::Float64 = 0.)

    # Check that data has correct dimensions
    N = length(x)
    Ny= length(y)
    (N != Ny) && error("Leng of y (=$Ny) no equal to length of x (=$N)")
    (N < 3)   && error("Length of data (=$N) must be > 3")
    for i=2:N
        (x[i] <= x[i-1]) && error("x values must be monotonically increasing")
    end

    # Dimensioning of arrays
    D=zeros(Float64,N)            # h.r.s. lengths
    h=zeros(Float64,N-1)          # Interval lengths
    l=zeros(Float64,N)           # Lower terms of the tridiagonal matrix
    u=zeros(Float64,N)            # Upper terms of the tridiagonal matrix
    di=zeros(Float64,N)           # Diagonal terms
    Y =zeros(Float64,N)           # Second derivatives
                                     # This is y(x) = a + bx + cx^2 + dx^3
    a = zeros(Float64,N-1)        # Independent coeff of splines
    b = zeros(Float64,N-1)        # Linear coeff of splines
    c = zeros(Float64,N-1)       # Quadratic coeff of splines
    d = zeros(Float64,N-1)        # Cubic coeff of splines

    # Compute h coefficients
    for i=N-1:-1:1
        h[i] = x[i+1] - x[i]        # interval length
    end
    for i=N-1:-1:1
        d[i] = (y[i+1]-y[i])/(x[i+1]-x[i])  # Differences
    end

    # End conditions
    if EndCond1==1        # First derivative prescribed at start
        D[1] = 6.0 * ( d[1] - EC1 )
        u[1] = h[1]
        di[1]= 2.0*h[1]
    elseif EndCond1==2    # Second derivative prescribed at start
        D[1] = EC1
        u[1] = 0.0
        di[1]= 1.0
    else
        error("EndCond1 cannot be $EndCond1")
    end
    if EndCond2==1        # End derivative prescribed at end
        D[N] = 6.0 * ( EC2 - d[N-1] )
        l[N] = h[N-1]
        di[N]= 2.0*h[N-1]
    elseif EndCond2==2    # Second derivative prescribed at end
        D[N] = EC2
        l[N] = 0.0
        di[N]= 1.0
    else
        error("EndCond2 cannot be $EndCond2")
    end

    # r.h.s. equal to d coefficients; calculate also r, s terms of matrix
    for i=N-1:-1:2
        D[i] = 6.0 * ( d[i] - d[i-1] )
        l[i] = h[i-1]
        u[i] = h[i]
        di[i]= 2.0 * ( h[i-1] + h[i] )
    end

    # Solve tridiagonal system (forward reduction + backsubstitution)
    u[1] /= di[1]
    D[1] /= di[1]
    for i=2:N-1
        z    = di[i] - l[i] * u[i-1]
        u[i] /= z
        D[i] = ( D[i] - l[i] * D[i-1] ) / z
    end
    D[N] = ( D[N] - l[N] * D[N-1] ) / ( di[N] - l[N] * u[N-1] )
    Y[N] = D[N]
    for i=N-1:-1:1
        Y[i] = D[i] - u[i] * Y[i+1]
    end

    # Obtain the rest of the coefficients; these are the coeffs of
    #   y(x) = a + b (x-x0) + c (x-x0)^2 + d (x-x0)^3 with x0 origin of interval
    for i=1:N-1
        a[i]  = y[i]
        b[i]  = d[i] - h[i] / 6.0 * (2.0 * Y[i] + Y[i+1])
        c[i]  = Y[i] / 2.0
        d[i]  = (Y[i+1] - Y[i]) / (6.0 * h[i])
    end

    # This is a flag to tell that it is a uniform spline
    aN = 0.0; if all(h[1] .== h); aN = 1.0; end

    # This is saved in a single variable
    C = [[a;aN] [b;y[1]] [c;[y[N]]] [d;0.0] x]

    # THIS IS THE FUNCTION RETURNED

    """
        csplinefval(x::Float64)

    Returns the value of a cubic spline given the coefficients C
    N is the number of break points, N-1 the number of polynomia
    C has dimension of [N,5], the last column are the break points
    if C[N,1] = -1

    See function csplinef() generating the coefficients

    """
    function csplinefval(xx)

        # Get dimensions of problem
        N = size(C)[1]        # Number of points of the spline (N-1) are the pieces
        xx = convert(Float64, xx)         # Guarantee float

        # Check if uniform spline or not
        if C[N,1] == 1.

            # This is a uniform spline, just take uniform h and compute indexing
            x0 = C[1,5]
            h = C[2,5] - x0
            #
            # For each value in x
                # Compute applicable interval s
                if xx <= x0
                    s = 1
                elseif xx >= x0 + h * (N-1)
                    s = N-1
                else
                    s = floor(Int,(xx - x0)/h) + 1
                end
                # Compute solution
                z = xx - C[s,5]; z2 = z*z
                yy = C[s,1] + C[s,2] * z + C[s,3] * z2 + C[s,4] * z2 * z
        else

            # This is a non-uniform spline, must check applicable interval
            #
            # For each value in x

                # Search the applicable interval s by bisection method
                s0 = 1
                s1 = N-1
                if xx <= C[s0,5]
                    s = s0
                elseif xx >= C[s1,5]
                    s = s1
                else
                    s0 = 1; s1 = N-1
                    while s1-s0 > 1
                        si = s0 + floor(Int,(s1-s0)/2)
                        (xx == C[si,5]) && (s0 = si; s1 = si)
                        (xx > C[si,5])  && (s0 = si)
                        (xx < C[si,5])  && (s1 = si)
                    end
                    s = s0
                end
                # Compute solution
                z = xx - C[s,5]; z2 = z*z
                yy = C[s,1] + C[s,2] * z + C[s,3] * z2 + C[s,4] * z2 * z
        end
        return yy
    end # of csplineval
    #
    # Return the name of the function that evaluates
    return csplinefval::Function
end # of csplinef

##################################### end: csplinef ############################

##############
#     2      #
##############
################################# begin: createIL() ############################
"""
    createIL() -> InvLangevin(x)::Function

Creates the inverse langevin function
    returns a function InvLangevin(x)

**Examples**
```
julia> ILf = createIL()
julia> x = 0.
julia> plot(x, ILf(x))
```
"""
function createIL()
    # STEP 0: Selection of some parameters
    #np  = 1001. ; xir = 0.943; yir = 17.543859649122449  # Points & cut-off
    np  = 100001.; xir = 0.980; yir = 50.      # This option is almost emach
    dyir= yir^2/(1+yir^2-yir^2*(coth(yir))^2)  # Derivative @ cutoff

    # STEP 1: Evaluate the Langevin Function for some y values
    y1 = collect(1e-14:(2yir-1e-14)/(2np-2):2yir)  # 2np-1 points between 1e-14 and 2yir
    x1 = coth.(y1)-1.0./y1          # this is the Langevin function for those points
    x1[1] = 0.0; y1[1] = 0.0        # (the inverse Langevin function, but not good for efficient evaluation)

    # STEP 2 and 3: Initial spline of the Inverse Langevin function
    sp1 = csplinef(x1,y1) # Creates a nonuniform spline function (uniform in ordinates)

    # STEP 4: Optimized spline for fast evaluation
    x2 = collect(0:xir/(np+1):xir); # Uniform spacing for efficiency in eval.
    dx = x2[2] - x2[1]       # This is the uniform abscisae spacing
    y2 = sp1.(x2)             # Evaluates values for uniform spacing
    sp2= csplinef(x2,y2,1,3.0,1,dyir) # Uniform spline w/ prescribed end der.

    # STEP 5: Rational function to accomodate asymptote at 1 [ax+b]/(1-x^2)
    btemp = yir*(1.0 -xir^2)
    aILrf = -2.0*btemp*xir/(1.0-xir^2) + dyir*(1.0-xir^2)
    bILrf = btemp - aILrf*xir

    # STEP 6: Create the Inverse Langevin function to be returned
    """
        InvLangevin(xval) -> yval::Float64
    Returns the inv of Langevin at a certain number xval

    **Examples**
    ```
    julia> y = InvLangevin(xval)
    ```
    """
    function InvLangevin(xval)
        xval = convert(Float64, xval)
            if xval <= 0.
                yval = 0.
            elseif xval >=1.
                yval = Inf64
            elseif xval < xir
                yval = sp2(xval)
            else
                yval = (aILrf*xval + bILrf) / (1.0 -xval^2)
            end
        return yval::Float64
    end
    return InvLangevin::Function
end

################################### end: createIL() ############################

##############
#     3      #
##############
################################# begin: quitNaN!() ############################


"""
    quitNaN!(x::Array{Float64,1})

This function detects NaN in the vector x and replaces them with the value in its
    right



# Examples
```
julia> x = [NaN 1. 1. 0. ]
julia> quitNaN!(x)
julia> x
[1. 1. 1. 0.]

```

"""

function quitNaN!(x::Array{Float64,1})

    for i=1:length(x)

        if (isnan(x[i]))

            x[i] = x[i+1]

        else
        end

    end
end

################################### end: quitNaN!() ############################

##############
#     4      #
##############
################################ begin: readCurve() ############################
"""
    readCurve(Nombre::String)-> x::Array{Float64,1}, y::Array{Float64,1}

Read a files with 2 columns data and saves each column in a vector.

Input:
    Nombre::String = it includes the name and the extension of the file

Output:
    x::Array{Float64,1} : contains the first column in the file
    Y::Array{Float64,1} : contains the first column in the file

# Examples
```
julia> x, y = readCurve("Nombre_Curva.extension")

```
"""

function readCurve(Nombre::String)

f = readlines(open(Nombre,"r"));
x = zeros(Float64,length(f))
y = zeros(Float64,length(f))

for i = 1:length(f)
     s1, s2= split(f[i], "\t");
     s1 = parse(Float64, s1);
     s2 = parse(Float64, s2);
     x[i] = s1;
     y[i] = s2;
end
 return x::Array{Float64,1}, y::Array{Float64,1}
end
################################# end: readCurve() #############################

##############
#     5     #
##############
########################### begin: bsplinefdf() ################################
"""
bsplinefdf(x_data, y_data,flag1 = 1, flag2 = 1, nvertices = 11)->f::Function, df::Function

Builds b-spline with penalization of D2 (second derivative) using periodic splines
as base functions. Options can be specified to increase the penalization of D2
in those vertices which does not satisfy the condition D1 (first derivative)>0

Input:
    x_data::Array{Float64,1} : x values from which b-spline is computed
    y_data::Array{Float64,1} : y values corresponding to those in x_data.

Optional:
    flag1::Int64 : Flag to modify the penalization on the second derivative

            flag1 == 1  (Default)
                - Penalizes all the vertices with a value.
                - Increases Penalization at the beginning and at the end.
            flag1 == 2
                - No penalization
            flag1 == 3
            - Penalizes all the vertices with a value.
            - Increases Penalization at the beginning and at the end.
            - Penalizes more those vertices in which 1st derivative is negative

    flag2::Int64 : Flag to modify the penalization applied when
    (yB_spline(x_data)-y_data)^2 = 0 is not satisfied.

            flag2 == 1
                - Increases penalization at the end
                - same penalization for the others
            flag2 == 2
                - Increases penalization at the beginning
                - Same penalization for the others
            flag2 == 3
                - Increases penalization at the beginning and at the end
                - Same penalization for all the vertices
            flag2 ==4
                - Same penalization for all the vertices

    nvertices : Number of Vertices to build the B-spline (11 vertices by default)

Output:
    Bsplinefval::Function : B-spline function obtained from data

    Bsplinedffval::Function : derivative of the B-spline function obtained from data

# Examples
```
julia> f, df = bsplinefdf(x_data,y_data,flag1,flag2,11)
julia> y = f(x) # x and y are numbers
julia> dy = df(x) # x and dy are numbers
julia> y = f.(x) # x and y are arrays
julia> dy = df.(x) # x and dy are arrays

```
"""
function bsplinefdf(x_data::Array{Float64,1}, y_data::Array{Float64,1},
                    flag1::Int64 = 1, flag2::Int64 = 1, nvertices::Int64 = 11)
    #=
        - Knots are considered to be at the x axis
        - Periodic Base functions are used
        - The coordinates of the control polygon are given by (KnotV[j],B[j])
            where B[j] are to be determined.
    =#

    m = nvertices-3   # segments in the evaluation range
    x0 = x_data[1]    # beginning of the evaluation range
    xf = x_data[end]  # ending of the evaluation range
    lsubd = (xf-x0)/(convert(Float64,nvertices)-3.) # segment length(distance between adjacent knots)
    knotVin = [x0:lsubd:nextfloat(xf) ...] # Knots within the evaluation range
    knotV = [x0-lsubd; knotVin; xf+lsubd ] # All the Knots (x coord of knots)


    nodesaffected = 2 # Dj(2) requires Bj+2 and Bj+1 apart from Bj

# Matrix that weights the importance of the error in the y coordinate for every vertex
    W = zeros(Float64, length(x_data),length(x_data))

    for i = 1:length(x_data)
        W[i,i] = 1.
    end
    vert1 = Int64(floor(nvertices/5))  # Auxiliar variable to introduce smooth penalty variation
    vert2 = Int64(floor(nvertices/2.5)) # Auxiliar variable to introduce smooth penalty variation


if flag2 == 1

    # Increase penalization at the end
    W[end,end] = 80.
    W[end-vert1,end-vert1] = 60.
    W[end-vert2,end-vert2] = 40.

elseif flag2 == 2
    # Increase penalization at the beginning
    W[1,1] = 100.
    W[1+vert1,2+vert1] = 40.
    W[1+vert2,1+vert2] = 20.

elseif flag2 == 3
    # Increase penalization at the beginning
    W[1,1] = 100.
    W[1+vert1,2+vert1] = 40.
    W[1+vert2,1+vert2] = 20.
    # Increase penalization at the end
    W[end,end] = 80.
    W[end-vert1,end-vert1] = 60.
    W[end-vert2,end-vert2] = 40.

elseif flag2 == 4
     # Do nothing
else

    error("flag2 in bsplinefdf has a wrong value")

end



# Matrix that weights the importance of having D(2) different from zero for every vertex
    Ω = zeros(Float64, nvertices-nodesaffected, nvertices-nodesaffected)

    if flag1 == 1 # begin if flag1
    # Uniform penalization
        for i = 1:nvertices-nodesaffected
            Ω[i,i] = 2.0
        end # end for

    # Increase penalization at the beginning
        Ω[1,1] = 6.
        Ω[1+vert1,2+vert1] = 4.
        Ω[1+vert2,1+vert2] = 2.

    # Increase penalization at the end
        Ω[end,end] = 6.
        Ω[end-vert1,end-vert1] = 4.
        Ω[end-vert2,end-vert2] = 3.

    elseif flag1 == 2 # else if flag1

    # No penalization
        for i = 1:nvertices-nodesaffected
            Ω[i,i] = 0.
        end # end for

    elseif flag1 == 3 # else if flag1

    # Uniform penalization
        for i = 1:nvertices-nodesaffected
            Ω[i,i] = 2.0
        end # end for
    # Increase penalization at the beginning
        Ω[1,1] = 6.
        Ω[1+vert1,2+vert1] = 4.
        Ω[1+vert2,1+vert2] = 2.
    # Increase penalization at the end
        Ω[end,end] = 6.
        Ω[end-vert1,end-vert1] = 4.
        Ω[end-vert2,end-vert2] = 3.

        # flag1 == 3 also adds additional penalization see Iteration "Iteration to reduce D(1)<0"

    else # else  flag1
        error(" flag1  wrong value in bsplinefdf ")
    end # end if flag1



################################################################################
### functions to build the matrices needed to solve  the system of equations ###


############################## begin: buildNmatrix() ###########################
function buildNmatrix(xdata,knotVin)
# Compute the length of data set

l = length(xdata)

# Basis functions definition
Transf_matrix = (1. / 6.) * [ -1. 3. -3. 1.;
                               3. -6. 3. 0.;
                               -3. 0. 3. 0.;
                                1. 4. 1. 0.]

# Number of knots to be less than number of data
nvertices = length(knotVin) + 2 # 2 used nodes has to be added
# number of used vertices will be m + 1 + 2

lsubd = knotVin[2] - knotVin[1]
Nmatrix = zeros(Float64,l,nvertices)

    for i = 1:l
        # Compute applicable interval s
        if nextfloat(xdata[i]) < knotVin[1]
            error("The value is out of the domain (left)")
        elseif xdata[i] > nextfloat(knotVin[end])
            error("The value is out of the domain (right)")
        else
            s = floor(Int,(xdata[i]-knotVin[1])/lsubd) + 1

            if s > nvertices-3
                s = nvertices-3
            end
        end

        # Once s has been determined compute local coord and evaluate base functions

        psi = (xdata[i]-knotVin[s])/lsubd
        xi_local = [psi^3, psi^2, psi, 1.]
        Nmatrix[i,s:s+3] = xi_local'*Transf_matrix
    end

    return Nmatrix::Array{Float64,2}
end
############################## end: buildNmatrix() ###########################


############################## begin: buildD3matrix() ##########################
#=
buildD3matrix(nvertices) -> Dmatrix::Array{Float64,2}
Computes a matrix of coeff that multiplied by te vertices y coord gives an approx
of the 3rd derivative with the polygon vertices.
 Dmatrix*B = 2nd deriv approx
=#
function buildD3matrix(nvertices)

Dmatrix = zeros(Float64,nvertices-3,nvertices)

    for i=1:nvertices-3
        Dmatrix[i,i] = -1.
        Dmatrix[i,i+1] = 3.
        Dmatrix[i,i+2] = -3.
        Dmatrix[i,i+3] = 1.
    end

    return Dmatrix::Array{Float64,2}

end
############################### end: buildD3matrix() ###########################

############################## begin: buildD2matrix() ##########################
#=
buildD2matrix(nvertices) -> Dmatrix::Array{Float64,2}

Computes a matrix of coeff that multiplied by te vertices y coord gives an approx
of the 2nd derivative with the polygon vertices.
 Dmatrix*B = 2nd deriv approx
=#
function buildD2matrix(nvertices::Int64)

    Dmatrix = zeros(Float64,nvertices-2,nvertices)

    for i=1:nvertices-2
        Dmatrix[i,i] = 1.
        Dmatrix[i,i+1] = -2.
        Dmatrix[i,i+2] = 1.
    end

    return Dmatrix::Array{Float64,2}

end


########################€###### end: buildD2matrix() ###########################

############################## begin: buildD1matrix ############################
#=
buildD1matrix(nvertices) -> Dmatrix::Array{Float64,2}

Computes a matrix of coeff that multiplied by te vertices y coord gives an approx
of the 1st derivative with the polygon vertices.
 Dmatrix*B = 1st deriv approx
=#

function buildD1matrix(nvertices::Int64)

Dmatrix = zeros(Float64,nvertices-1,nvertices)

    for i=1:nvertices-1
        Dmatrix[i,i] = -1.
        Dmatrix[i,i+1] = 1.
    end

    return Dmatrix::Array{Float64,2}

end

########################€###### end: buildD1matrix() ###########################

################################################################################

    D1 = buildD1matrix(nvertices)
    D2 = buildD2matrix(nvertices)
    D3 = buildD3matrix(nvertices)
    N = buildNmatrix(x_data,knotVin)
################################################################################
####### functions to construct the matrix A and solve the system Ax = b #######

########################### begin: compute_vertices() ##########################
#=
compute_vertices(x,y,N,D,W,Ω) -> B::Array{Float64,1}

The y coord of the vertices given the matrices N,D,W and Ω and the experimental points
given by the vectors x and y.
=#
function compute_vertices(x::Array{Float64,1}, y::Array{Float64,1},
                          N::Array{Float64,2}, D::Array{Float64,2},
                          W::Array{Float64,2}, Ω::Array{Float64,2})


A = N'*W*N + D'*Ω*D
# independent terms vector
b = N'*W*y

# Solving the linear system
B = A\b

return B::Array{Float64,1}
end

############################# end: compute_vertices() ##########################
################################################################################

# Iteration to reduce D(1)<0

# It increases the penalization on the 2nd derivative on those vertices where
# D(1)<0

if flag1 == 3  # if flag1 == 3

    Δ = 0.2
# Try to eliminate D(1)<0 incrementing up to 60 times the penalization Ω
    for j=1:60  # for j 1 to 60
    # with initial Ω compute vertices
    B =  compute_vertices(x_data, y_data, N, D2, W, Ω)
    # with the computed vertices compute approximations to 1st, 2nd and 3rd deriv
    D1B = D1*B
    D2B = D2*B
    D3B = D3*B
    # find vertices where 1st deriv < 0
    actinvertices = findall(D1B.<0)

    # where 1st deriv < 0 modify the penalty
        if(isempty(actinvertices)) # if isempty(actinvertices)
            break
        else # else isempty(actinvertices)
            for i =1:length(actinvertices) # for i

                if (actinvertices[i] == 1) # if actinvertices[i]
                    Ω[actinvertices[i],actinvertices[i]] = Ω[actinvertices[i],actinvertices[i]] +  Δ
                    Ω[actinvertices[i]+1,actinvertices[i]+1] = Ω[actinvertices[i]+1,actinvertices[i]+1] +  Δ/2.
                elseif (actinvertices[i] >= nvertices-nodesaffected) # else if actinvertices[i]
                    Ω[nvertices-nodesaffected-1,nvertices-nodesaffected-1] = Ω[nvertices-nodesaffected-1,nvertices-nodesaffected-1] +  Δ/2.
                    Ω[nvertices-nodesaffected,nvertices-nodesaffected] = Ω[nvertices-nodesaffected,nvertices-nodesaffected] +  Δ
                else # else actinvertices[i]
                    Ω[actinvertices[i]-1,actinvertices[i]-1] = Ω[actinvertices[i]-1,actinvertices[i]-1] +  Δ/2.
                    Ω[actinvertices[i],actinvertices[i]] = Ω[actinvertices[i],actinvertices[i]] +  Δ
                    Ω[actinvertices[i]+1,actinvertices[i]+1] = Ω[actinvertices[i]+1,actinvertices[i]+1] +  Δ/2.
                end # end actinvertices[i]
            end # end for i
        end # end if isempty(actinvertices)
    end # end for j 1 to 60

end # end if flag1 == 3

# Once the iteration process finishes, compute the final vertices y coord
    B =  compute_vertices(x_data, y_data, N, D2, W, Ω)
    D1B = D1*B


############################## begin: Bsplinefval() ############################

"""
Bsplinefval(xeval::Float64)->yval::Float64

Evaluates the value of the Bspline computed in the main program bsplinefdf
in a certain value.

# Examples
```
julia> yeval = Bsplinefval(xeval) # xeval and yeval are numbers
julia> yeval = Bsplinefval.(xeval) # xeval and yeval are arrays

```
"""
    function Bsplinefval(xeval::Float64)

        Transf_matrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                       3. -6. 3. 0.;
                                       -3. 0. 3. 0.;
                                        1. 4. 1. 0.]

        nvertices = length(B)
        lsubd = knotVin[2]-knotVin[1]

       if nextfloat(xeval) < knotVin[1]
           error("The value is out of the domain (left)")
       elseif xeval > nextfloat(knotVin[end])
           error("The value is out of the domain (right)")
       else
           s = floor(Int,(xeval - knotVin[1])/lsubd) + 1

           if s > nvertices-3
               s = nvertices-3
           end
       end
       # Once s has been determined

       psi = (xeval - knotVin[s])/lsubd

       #print(s)
        B_local = [ B[s], B[s+1], B[s+2], B[s+3] ]
        xi_local = [ psi^3, psi^2, psi, 1.]

       yfinal = xi_local'*Transf_matrix*B_local


        return yfinal::Float64
    end


################################################################################

############################# begin: Bsplinedffval() ###########################
"""
Bsplinedffval(xeval::Float64)->yval::Float64

Evaluates the derivative  of the Bspline computed in the main program bsplinefdf
in a certain value.

# Examples
```
julia> dyeval = Bsplinedffval(xeval) # xeval and dyeval are numbers
julia> dyeval = Bsplinedffval.(xeval) # xeval and dyeval are arrays

```
"""


    function Bsplinedffval(xeval::Float64)

        nvertices = length(B)
        lsubd = knotVin[2]-knotVin[1]
        Transf_matrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                       3. -6. 3. 0.;
                                       -3. 0. 3. 0.;
                                        1. 4. 1. 0.]


       if nextfloat(xeval) < knotVin[1]
           error("The value is out of the domain (left)")
       elseif xeval > nextfloat(knotVin[end])
           error("The value is out of the domain (right)")
       else
           s = floor(Int,(xeval - knotVin[1])/lsubd) + 1

           if s > nvertices-3
               s = nvertices-3
           end
       end
       # Once s has been determined
       psi = (xeval - knotVin[s])/lsubd

       #print(s)

       B_local = [ B[s], B[s+1], B[s+2], B[s+3] ]
       xi_local = [3. * psi^2, 2. * psi, 1., 0.]


       dydxi = xi_local'*Transf_matrix*B_local

       dyfinal = dydxi / lsubd

    return dyfinal::Float64
    end
################################################################################



    return Bsplinefval::Function, Bsplinedffval::Function
    #, D1B::Array{Float64,1}
end

################################ end: bsplinefdf() #############################


##############
#     6     #
##############
############################### begin: rational() ##############################
"""
    rational(xl::Float64, yl::Float64, dyl::Float64, λch_b::Float64) -> f::Function, df::Function

Returns a function which is a rational passing through the point (xl,yl),
    with derivative (xl,dyl) and with asymptote at λch_b.

Input:
    xl::Float64 = abscisae the connection point
    yl::Float64 = function value at the connection point
    dyl::Float64 = function derivative at the connection point
    λch_b::Float64 = value in which the asymptote is placed

Output:
    f::Function : rational function obtained
    df::Function : derivative of the rational function obtained

# Examples
```
julia> f, df = rational(xl,yl,dyl,λch_b)

```

"""
function rational(xl::Float64, yl::Float64, dyl::Float64, λch_b::Float64)

    a = dyl*(λch_b^2 - xl^2) - 2. * yl * xl
    b = yl*(λch_b^2 - xl^2) - a*xl

    f(x) = (a*x + b)/(λch_b^2 - x^2)
    df(x) = (a*x^2 + 2. * b * x + a*λch_b^2)/(λch_b^2 - x^2)^2

    return f::Function, df::Function

end


################################# end: rational() ##############################


##############
#     7     #
##############
############################ begin: rational_joint() ###########################
"""
    rational_joint(x, y, xl, xlim) -> f::Function, df::Function

Given 2 vectors x and y it returns a piecewise function which is the B_spline with
    computed with x and y from x(1) to xl and a rational curve which has C1 continuity
    and vertical asymptote at xlim.

Input:
    x::Array{Float64,1} = abcisae vector to construct B-spline
    y::Float64{Float64,1} = y values corresponding to x used to construct the B-spline
    xl::Float64 = point at which continutity C1 is imosed between the B-spline
     and the rational
    xlim::Float64 = value in which the asymptote is placed
Output:

    f::Function : piecewise function (Bspline and rational function)
    df::Function : piecewise function ( derivative of the Bspline and of the rational function)


# Examples
```
julia> f, df = rational_joint(x,y,xl,xlim)

```

"""

function rational_joint(x::Array{Float64,1}, y::Array{Float64,1},
                        xl::Float64 , xlim::Float64,
                        cond1::Int64 = 3, cond2::Int64 = 1)

if length(x)>11
    f, df= bsplinefdf(x, y, cond1, cond2)
else
    f, df= bsplinefdf(x, y, cond1, cond2, length(x)-3)
end

    yl = f(xl)
    dyl = df(xl)

    f2, df2 = rational(xl, yl, dyl, xlim)


    function raional_joint_eval(xeval)


            if xeval > xl #if 1
                if xeval > xlim
                    error(" Trying to evaluate the function after blocking")
                end
                yeval = f2(xeval)
            elseif xeval <= xl  #elseif 1
                yeval = f(xeval)
            else
            end # end if 1

        return yeval
    end # end of the rational_joint_eval

    function raional_joint_eval_df(xeval)


            if xeval > xl #if 1
                if xeval > xlim
                    error(" Trying to evaluate the function after blocking")
                end
                dyeval = df2(xeval)
            elseif xeval <= xl  #elseif 1
                dyeval = df(xeval)
            else
            end # end if 1

        return dyeval
    end # end of the rational_joint_eval


    return raional_joint_eval::Function, raional_joint_eval_df::Function
end # end of the rational_joint


############################## end: rational_joint() ###########################

end  #  end of the module MathematicsToolsWYPIWYG
