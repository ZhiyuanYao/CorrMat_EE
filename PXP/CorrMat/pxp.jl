#=============================================================================#
#   This program calculate the local correlation matrix spectrum of the the PXP
#   model using range-k basis (excluding the trivial identity operator)
#
#                       Li = σ_1 ⊗ σ_2 ⊗ σ_3 ⋯ σ_k
#
#   Currently, the program is limitted to PBC case since we are calling Fortran
#   program to feed us with the eigenstates in the (k,I)=(0,+) sector.
#
#   Version:  1.0
#   Created:  2022-04-26 10:38
#    Author:  Zhiyuan Yao, zhiyuan.yao@me.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#=============================================================================#
using LinearAlgebra
using SparseArrays
using OffsetArrays
using DelimitedFiles

#-----------------------------------------------------------------------------#
# System parameters
#-----------------------------------------------------------------------------#
const ⊗ = kron
const range = 5
const dimLH = 4^(range) - 1

#-----------------------------------------------------------------------------#
# Utility functions
#-----------------------------------------------------------------------------#
function speye(n::Integer, ::Type{T}=Float64) where T <: Union{Int32, Int64, Float64, ComplexF64}
    if n < 1
        error("speye(): argument error")
    end
    return sparse(T(1)*I, n, n)
end

#=============================================================================#
# call Fortran program
#=============================================================================#
run(`make clean`); run(`make`)
const L = ccall((:get_l, "./symmetryed.so"), Cint, ())
const mass = ccall((:get_mass, "./symmetryed.so"), Cdouble, ())
const dimH = 2^L
ccall((:__symmetryed_MOD_initialize, "./symmetryed.so"), Cvoid, ())
println("range = $(range),  mass = $(mass)")

#------------------------------------------------------------------------------
# The cglobal access global variable and returns a pointer to its address.
#------------------------------------------------------------------------------
const basisN = unsafe_load(cglobal((:__symmetryed_MOD_basisn, "./symmetryed.so"), Int32))
const stateN = unsafe_load(cglobal((:__symmetryed_MOD_staten, "./symmetryed.so"), Int64))
const state_istate = unsafe_wrap(Vector{Int64}, ccall((:get_statetable, "./symmetryed.so"), Ptr{Int64}, ()), stateN)
const index_istate = (state_istate .+ 1)

function geteigvect(n::Int)
    psi = unsafe_wrap(Vector{Float64}, ccall((:get_psi, "./symmetryed.so"), Ptr{Float64}, (Ptr{Int64}, ), Ref(n)), stateN)
    return sparse(psi)
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
    return  nzero, vals
end

#=============================================================================#
# eigenList: calculate the correlation spectrum for n in eigenList
#=============================================================================#
function getCorrSpectrum(; eigenList=collect(1:basisN))
    nzeroList = zeros(Int, basisN)
    lamvalMat = zeros(dimLH, basisN)
    for n in eigenList
        nzero, lamvals = getCorrelationSpectrum(geteigvect(n))
        nzeroList[n] = nzero
        lamvalMat[:, n] = lamvals
        println("n = $n, nzero = $nzero")
    end
    run(`mkdir -p ./data`)
    writedlm("./data/L$(L)_k$(range)_nzeros.dat", nzeroList)
    writedlm("./data/L$(L)_k$(range)_lamMat.dat", lamvalMat)
    writedlm("./data/L$(L)_k$(range)_ibasis.dat", eigenList)
end

# getCorrSpectrum(eigenList=collect(30:102))
getCorrSpectrum()
