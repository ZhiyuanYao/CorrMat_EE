#=============================================================================#
#   This program use local correlation matrix with range-k and basis
#   Li = σ_1 ⊗ σ_2 ⊗ σ_3 ⋯ σ_k to study the PXP model.
#
#   Currently, the program is limitted to PBC case since we are calling Fortran
#   program to feed us with the eigenstates in the (k,I)=(0,+) sector.
#
#   Version:  1.0
#   Created:  2022-04-26 11:09
#    Author:  Zhiyuan Yao, zhiyuan.yao@me.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#=============================================================================#
# Change Log:
#   *
#==============================================================================#
using LinearAlgebra
using SparseArrays
using OffsetArrays
using DelimitedFiles

#------------------------------------------------------------------------------
# system parameters
#------------------------------------------------------------------------------
const L = 20
const range = 5
const dimH = 2^L
const dimLH = 4^(range) - 1
const ⊗ = kron

#-----------------------------------------------------------------------------#
# Utility functions
#-----------------------------------------------------------------------------#
function speye(n::Integer, ::Type{T}=Float64) where T <: Union{Int32, Int64, Float64, ComplexF64}
    if n < 1
        error("speye(): argument error")
    end
    return sparse(T(1)*I, n, n)
end
#------------------------------------------------------------------------------
function getThermalEntropy(beta::Real, eigvals::Vector{Float64})
    if beta < 0
        println("!!!==== getThermalEntropy(): input beta < 0 ====!!!")
    end
    Elist = (issorted(eigvals) ? eigvals : sort(eigvals))
    if beta < 0
        # to make sure in the for sum below p is decreasing from 1
        Elist = reverse(Elist)
    end
    E_0 = Elist[1]
    Z, sum = 0.0, 0.0
    for E in Elist
        p = exp(-beta*(E-E_0))
        Z += p
        # sum += -p*log(p)     # unstable for large beta, gives NaN
        sum += p*beta*(E-E_0)
    end
    return sum/Z + log(Z)
end
#------------------------------------------------------------------------------
function getThermalBeta(E::Real, eigvals::Vector{Float64}, beta0=1.0::Real; std=1E-6)
    function func(beta::Real)
        Elist = (issorted(eigvals) ? eigvals : sort(eigvals))
        if beta < 0
            # to make sure in the for sum below p is decreasing from 1
            Elist = reverse(Elist)
        end
        E_0 = Elist[1]
        Z, Esum = 0.0, 0.0
        for E in Elist
            p = exp(-beta*(E-E_0))
            Z += p
            Esum += E*p
        end
        E = Esum/Z
        return E
    end
    beta, std_x, std_y = doNewtonBisect(beta0, E, func; step=0.1, std_x=std, std_y=std)
    if beta < 0
        println("!!!==== getThermalBeta(): beta < 0 ====!!!")
    end
    return beta
end
#------------------------------------------------------------------------------
function doNewtonBisect(x0::T, ystar::T, f::Function; step=0.1, std_x=1E-6, std_y=Inf, Nmax=1000) where T <: Real
    if step < 0
        error("doNewtonBisect(): step must be bigger than zero!")
    end
    y0 = f(x0)
    x1 = x0 + step
    y1 = f(x1)
    slope = (y1 > y0 ? 1 : -1)
    if slope == 1
        if ystar > y0 && ystar > y1
            while ystar > y1
                step *= 2
                x0, y0 = x1, y1
                x1 = x1 + step
                y1 = f(x1)
            end
        elseif ystar < y0 && ystar < y1
            while ystar < y0
                step *= 2
                x1, y1 = x0, y0
                x0 = x0 - step
                y0 = f(x0)
            end
        end
        xa, ya, xb, yb = x0, y0, x1, y1
    else
        if ystar > y0 && ystar > y1
            while ystar > y0
                step *= 2
                x1, y1 = x0, y0
                x0 = x0 - step
                y0 = f(x0)
            end
        elseif ystar < y0 && ystar < y1
            while ystar < y1
                step *= 2
                x0, y0 = x1, y1
                x1 = x1 + step
                y1 = f(x1)
            end
        end
        xa, ya, xb, yb = x1, y1, x0, y0
    end
    itime = 0
    #------------------------------------------------------------------------------
    # make sure ya < yb to facillate robust programming; although this should
    # be guranteed from the above code
    if ya > yb
        xa, xb = xb, xa
        ya, yb = yb, ya
    end
    #------------------------------------------------------------------------------
    while true
        if itime > Nmax
            error("doNewtonBisect(): bisection performed $(Nmax) times without convergence!")
        end
        x  = (xa + xb)/2
        dx = abs(xb - xa)/2
        y  = f(x)
        dy = abs(y - ystar)
        if dx < std_x && dy < std_y
            return x, dx, dy
        elseif y > ystar
            xb = x
        else
            xa = x
        end
        itime += 1
    end
end
#------------------------------------------------------------------------------
function checkDensityMatrix(rho::Matrix{T}; tol=1E-10) where T <: Union{Float64, ComplexF64}
    dimRA = size(rho, 1)
    if dimRA != size(rho, 2)
        error("checkDensityMatrix(): argument size error")
    end
    if eltype(rho) == Float64
        rho = Symmetric(rho)
        for i in 1:dimRA
            for j in 1:i-1
                if abs(rho[i, j] - rho[j, i]) > tol
                    error("checkDensityMatrix(): real rho not symmetric")
                end
            end
        end
    else
        for i in 1:dimRA
            for j in 1:i-1
                if abs(rho[i, j] - conj(rho[j, i])) > tol
                    error("checkDensityMatrix(): complex rho off-diagonal not Hermitian")
                end
            end
            if abs(imag(rho[i, i])) > tol
                error("checkDensityMatrix(): complex rho diagonal not real")
            end
        end
    end
end


#=============================================================================#
# call Fortran program
#=============================================================================#
const momentum, inversion = 0, 1
const mass = 0.0

run(`make clean`); run(`make`)
ccall((:__symmetryed_MOD_initialize, "./symmetryed.so"), Cvoid, (Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Float64}), Ref(L), Ref(momentum), Ref(inversion), Ref(mass))
println("range = $(range)")

#------------------------------------------------------------------------------
# The cglobal access global variable and returns a pointer to its address.
#------------------------------------------------------------------------------
const basisN = unsafe_load(cglobal((:__symmetryed_MOD_basisn, "./symmetryed.so"), Int32))
const stateN = unsafe_load(cglobal((:__symmetryed_MOD_staten, "./symmetryed.so"), Int64))
const state_istate = unsafe_wrap(Vector{Int64}, ccall((:get_statetable, "./symmetryed.so"), Ptr{Int64}, ()), stateN)
const index_istate = (state_istate .+ 1)

function geteigvect(n::Int32)
    psi = unsafe_wrap(Vector{Float64}, ccall((:get_psi, "./symmetryed.so"), Ptr{Float64}, (Ptr{Int32}, ), Ref(n)), stateN)
    return sparse(psi)
end

#=============================================================================#
# Filter State and prepare state_istate
#=============================================================================#
const NOTHING = -1
const FibonacciTable = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610,
987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393]
function get_index_istate(L::Int, hasPBC::Bool)
    stateN = (hasPBC ? (FibonacciTable[L-1] + FibonacciTable[L+1]) : FibonacciTable[L+2])
    state_istate = Vector{Int64}(undef, stateN)
    istate_state = OffsetVector(fill(NOTHING, 2^L), 0:2^L-1)
    invalid_stateN = 0
    istate = 0
    for state in 0:2^L-1
        bits = (digits(state, base = 2, pad=L) |> reverse)
        #---------------------------------------------------------------------
        # site :   1 2 3    ... L-1 L
        #          _________________
        #         |_|_|_|_| ... |_|_|
        #
        # bit :  high weight     low weight
        #---------------------------------------------------------------------
        valid_state = true
        for i in 1:L-1
            if bits[i]*bits[i+1] == 1
                valid_state = false
                break
            end
        end
        if hasPBC && bits[1]*bits[L] == 1
            valid_state = false
        end
        if valid_state
            istate += 1
            state_istate[istate] = state
            istate_state[state]  = istate
        else
            invalid_stateN += 1
        end
    end
    if stateN + invalid_stateN != 2^L
        error("stateN + invalid_stateN != 2^L")
    end
    return (state_istate .+ 1)
end

#=============================================================================#
# use the SU(2) generators to construct the general Hamiltonian
# the Fortran convention for basis index is used |↓> = [1, 0] & |↑> = [0, 1]
#=============================================================================#
const σ_I = spzeros(2, 2); σ_I[1, 1] = σ_I[2, 2] = 1
const σ_x = spzeros(2, 2); σ_x[1, 2] = σ_x[2, 1] = 1
const σ_y = spzeros(Complex{Float64}, 2, 2); σ_y[1, 2] = im; σ_y[2, 1] = -im
const σ_z = spzeros(2, 2); σ_z[1, 1] = -1; σ_z[2, 2] = 1

const M1 = [σ_I, σ_x, σ_y, σ_z]

#------------------------------------------------------------------------------
# i takes from 1 to dimLH, i.e. from I ⊗ I ⊗... σ_x to σ_z ⊗ σ_z ⊗ ... σ_z
#------------------------------------------------------------------------------
function getL(i::Int)
    abcs = reverse(digits(i, base=4, pad=range) .+ 1)
    Li = foldl(⊗, [M1[n] for n in abcs[1:end]]) ⊗ speye(2^(L-range))
    Li = Li[:, index_istate]
    return Li
end

#=============================================================================#
# Correlation Spectrum for restricted case; use sparse properly to speed up
#=============================================================================#
function getCorrelationSpectrum(psi::SparseVector{T}) where T <: Union{Float64, ComplexF64}
    Ls = zeros(ComplexF64, dimLH)
    Mat = zeros(ComplexF64, dimLH, dimLH)
    # store Li|ψ⟩ as φ_i
    phiList = Vector{SparseVector{ComplexF64, Int64}}(undef, dimH)
    # calculate Ls
    for i in 1:dimLH
        phiList[i] = getL(i)*psi
        Ls[i] = psi'*(phiList[i][index_istate])  # ⟨ψ|L_i|ψ⟩
    end
    # calculate Mat
    for i in 1:dimLH
        for j in 1:i-1
            Mat[i, j] = (phiList[i])'*phiList[j] - conj(Ls[i])*Ls[j]  # general form is used
            Mat[j, i] = conj(Mat[i, j])
        end
        # diagonal term
        Mat[i, i] = norm(phiList[i])^2 - abs(Ls[i])^2
    end
    vals, vects = eigen(Hermitian(Mat))
    nzero = count(x -> abs(x) < 1E-10, vals)
    # println("Correlation Spectrum N₀ = $(nzero)")
    return  nzero, vals, vects
end

#=============================================================================#
# construct local
#=============================================================================#
function getlocalL(i::Int)
    abcs = reverse(digits(i, base=4, pad=range) .+ 1)
    Li = foldl(⊗, [M1[n] for n in abcs[1:end]])
    return Li
end

function constructLocalHam(psi::SparseVector{Float64}, nzero::Int, vects::Matrix{ComplexF64}, nterm::Int)
    Ham = zeros(ComplexF64, 2^range, 2^range)
    for n in nzero+1:min(nzero+nterm, dimLH)
        A  = zeros(ComplexF64, 2^range, 2^range)
        for i in 1:dimLH
            A .+= vects[i, n]*getlocalL(i)
        end
        O  = sparse(A) ⊗ speye(2^(L-range))
        # O |ψ⟩ = ξ |ψ⟩ + |ψ⟂⟩
        ξ = psi'*O[index_istate, index_istate]*psi
        Ham .+= (A - ξ*speye(2^range))'*(A - ξ*speye(2^range))
    end
    # make local Ham not too big so that the resulting spectrum of range sites is order 1
    return Ham/(dimLH*nterm)
end

const LA = unsafe_load(cglobal((:__symmetryed_MOD_la, "./symmetryed.so"), Int32))
const dimRA = unsafe_load(cglobal((:__symmetryed_MOD_dimra, "./symmetryed.so"), Int32))
const index_istate_A = get_index_istate(Int(LA), false)
function constructH_A(psi::SparseVector{Float64}, nterm::Int)
    nzero, vals, vects = getCorrelationSpectrum(psi)
    H = constructLocalHam(psi, nzero, vects, nterm)
    if LA < range
        error("constructH_A(): LA < range, increase L")
    end
    Ham = spzeros(ComplexF64, 2^LA, 2^LA)
    for n in 1:LA-range+1
        Ham .+= speye(2^(n-1)) ⊗ H ⊗ speye(2^(LA-n-range+1))
    end
    # select index_istate to make sure the constaints of underlying are
    # satisfied, this should place more constaints on the low energy
    # eigenvectors, thus gives a better bound on subsequent Sv
    return Ham[index_istate_A, index_istate_A]
end

#=============================================================================#
# get the entropy bound of given state
#=============================================================================#
function getEntropyBound(n::Int, nterm::Int)
    n = Int32(n)
    psi = geteigvect(n)
    H_A = constructH_A(psi, nterm)
    rhoA = zeros(dimRA, dimRA)
    ccall((:get_rhoa, "./symmetryed.so"), Cvoid, (Ptr{Int32}, Ptr{Int32}, Ptr{Matrix{Float64}}), Ref(n), Ref(dimRA), rhoA)
    checkDensityMatrix(rhoA)
    E = real(tr(rhoA*H_A))
    eigvals, = eigen(Hermitian(Matrix(H_A)))
    beta = getThermalBeta(E, eigvals)
    if beta < 0
        error("beta < 0!")
    end
    return getThermalEntropy(beta, eigvals)
end

function getMinimumEntropyBound(n::Int, nterm::Int)
    ntermList = [nterm + i for i in -10:10]
    ntermList = ntermList[findfirst(x -> x > 0, ntermList):end]
    entropyList = zeros(length(ntermList))
    for (in, nterm) in enumerate(ntermList)
        entropyList[in] = getEntropyBound(n, nterm)
    end
    perm = sortperm(entropyList)
    sort!(entropyList)
    println(ntermList[perm][1:5])
    println(entropyList[1:5])
    return entropyList[1], ntermList[perm][1]
end

#------------------------------------------------------------------------------
# L = 20, k= 3
# nList = [1, 2, 7, 30, 102]; nterm = 15    # N_op - N_0 = 39 - 24
# L = 20, k= 4
# nList = [1, 2, 7, 30, 102]; nterm = 64    # N_op - N_0 = 207/191 - 128
# L = 20, k= 5
# nList = [1, 2, 7, 30, 102]; nterm = 240   # N_op - N_0 = 895 - 608
#------------------------------------------------------------------------------
function getnListEntropyBound(nList::Vector{Int}, ntermList::T) where T <:Union{Int, Vector{Int}}
    entropyList  = zeros(length(nList))
    NoptimalList = zeros(Int, length(nList))
    if typeof(ntermList) == Int
        ntermList = fill(ntermList, length(nList))
    end
    for (in, n) in enumerate(nList)
        nterm = ntermList[in]
        entropyList[in], NoptimalList[in] = getMinimumEntropyBound(n, nterm)
    end
    run(`mkdir -p ./data/`)
    open("./data/L$(L)_k$(range)_EE_bound.dat", "a") do io
        writedlm(io, [nList entropyList NoptimalList])
    end
end
nList = [1, 2, 7, 30, 102];
ntermList = [297, 329, 203, 336, 315];  # k = 5
# ntermList = [81, 93, 76, 86, 69];  # k = 4
# ntermList = [10, 24, 12, 16, 19];  # k = 3
getnListEntropyBound(nList, ntermList)
