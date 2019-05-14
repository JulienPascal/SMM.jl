#-------------------------------------------------------------------------------
# WORK IN PROGRESS
#-------------------------------------------------------------------------------
"""
    NW_kernel(j::Int64,m::Int64)

Newey-West kernel function.
source: https://www.federalreserve.gov/pubs/ifdp/2012/1060/ifdp1060.pdf
"""
function NW_kernel(j::Int64,m::Int64)

    output=0.::Float64

    if j<= m
        output = 1. - (j/(m+1.0))
    end

    return output

end

"""

Function to calculate a Heteroskedasticity and Autocorrelation Consistent (HAC)
covariance matrix. The input X is a matrix with each column corresponding to
a variable and a row corresponding to a given period. m is the number of autocovariances
used to calculate the HAC covariance matrix. corrected is a Boolean. If true,
divide by T-1 instead of T.
"""
function hac_cov(X::Array{Float64,2}; m::Int64=4, kernel::Symbol=:NW, corrected::Bool=false)

    # Initialization:
    HAC_matrix = zeros(size(X,2),size(X,2))

    # Safety checks:
    if m < 0
        error("m has be greater than 0.")
    end

    # Newey-West formula
    #-------------------
    if kernel ==:NW
        # Calculate the first term (variance-covariance)
        HAC_matrix += cov(X)

        # Loop over different time lags
        #------------------------------
        for j=1:m
            Omega_j=Omega(X, j, corrected=corrected)
            HAC_matrix += NW_kernel(j,m)*(Omega_j + transpose(Omega_j))
        end
    else
        error("Only kernel=:NW is supported for the moment.")
    end

    return HAC_matrix

end

function hac_cov(df::DataFrame; m::Int64=4, kernel::Symbol=:NW, corrected::Bool=false)
    hac_cov(convert(Array, df), m=m, kernel=kernel, corrected=corrected)
end

"""
    Omega(X::Array{Float64,2}, j::Int64; corrected::Bool=false)

Function to calculate the "Omaga j" in Whitney K. Newey and Kenneth D. West (1987)
paper. Source: https://www.jstor.org/stable/1913610?seq=1#metadata_info_tab_contents
"""
function Omega(X::Array{Float64,2}, j::Int64; corrected::Bool=false)

    T, p = size(X)
    O = zeros(p,p)
    mean_X = zeros(p)
    for i=1:p
        mean_X[i] = mean(X[:, i])
    end
    for h=1:p, s = 1:p
        for t = j+1:T
            @inbounds O[h, s] += (X[t, s] - mean_X[s])*(X[t-j, h] - mean_X[h])
        end
    end
    if corrected==true
        O=(1/(T-1)).*O
    else
        O=(1/T).*O
    end
    return O
end


function Omega(X::Array{Int64,2}, j::Int64; corrected::Bool=false)
    Omega(convert(Array{Float64,2},X), j::Int64, corrected=corrected)
end
