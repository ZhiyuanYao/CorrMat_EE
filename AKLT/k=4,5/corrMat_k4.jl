#==============================================================================#
# Description
# -----------
#   This is program use MPO × MPS method to calculate the η-pairing states in
#   the AKLT model, |Ψ_s⟩ = (Q^†)^N |G⟩, and then calculate the correlation
#   matrix of range_-4 operators.
#
#   We also check to make sure that for N=1 scar different system size produce
#   the same eigen-operator space, Ker(M).
#
#   Version:  1.0
#   Created:  2022-06-13 15:30
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
# Change Log:
#   * This version is based on version 1.0 of k=3 corrMat.jl
#==============================================================================#
using Printf
using OffsetArrays
using SparseArrays
using LinearAlgebra
using DelimitedFiles

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
# Gram–Schmidt column vectors stored in matrix A
#------------------------------------------------------------------------------
function gram_schmidt(vecList::Vector{Vector{T}}; tol = 1E-10, info=true) where T <: Union{Float64, ComplexF64}
    B = fill(Vector{T}(undef, length(vecList[1])), length(vecList))
    dim = 0
    for i = 1:length(vecList)
        vec = vecList[i]
        if length(vec) != length(vecList[1])
            println("gram_schmidt(): Input vec in vecList size NOT the same")
            return
        end
        for j = 1:dim
            vec -= (B[j]'*vecList[i]) * B[j]
        end
        if norm(vec) > tol
            dim += 1
            B[dim] = vec/norm(vec)
        end
    end
    if dim < length(vecList) && info
        println("gram_schmidt(): Input $nc vectors are linearly dependent with dim = $dim")
    end
    return B[1:dim]
end
#------------------------------------------------------------------------------
function gram_schmidt(A::Matrix{T}; tol = 1E-10, info=true) where T <: Union{Float64, ComplexF64}
    nr, nc = size(A)
    B = similar(A)
    dim = 0
    for i = 1:nc
        vec = A[:, i]
        for j = 1:dim
            vec -= (B[:, j]'*A[:, i]) * B[:, j]
        end
        if norm(vec) > tol
            dim += 1
            B[:, dim] = vec/norm(vec)
        end
    end
    if dim < nc && info
        println("gram_schmidt(): Input $nc vectors are linearly dependent with dim = $dim")
    end
    return B[:, 1:dim]
end
#------------------------------------------------------------------------------
# calculate W - V
#------------------------------------------------------------------------------
function complement_span(W::Matrix{T}, V::Matrix{T}) where T <: Union{Float64, ComplexF64}
    if size(W, 1) != size(V, 1)
        error("common_span(): argument dimension error")
    end
    dimW = rank(W); dimV = rank(V)
    A = gram_schmidt(V, info=false)
    B = deepcopy(W)
    # project out the component of orthonormal basis in V
    for i in 1:size(B, 2)
        for k in 1:size(A, 2)
            B[:, i] -= (A[:, k]'*B[:, i])*A[:, k]
        end
    end
    W_V = gram_schmidt(B, info=false)
    if dimW - dimV != size(W_V, 2)
        error("complement_span(): dimension error")
    end
    return W_V
end
#------------------------------------------------------------------------------
function complement_span(W::Matrix{T}, V::Vector{T}) where T <: Union{Float64, ComplexF64}
    if size(W, 1) != length(V)
        error("common_span(): argument dimension error")
    end
    dimW = rank(W)
    B = deepcopy(W)
    A = normalize(V)
    # project out the component of orthonormal basis in V
    for i in 1:size(B, 2)
        B[:, i] -= (A'*B[:, i])*A
    end
    W_V = gram_schmidt(B, info=false)
    if dimW - 1 != size(W_V, 2)
        error("complement_span(): dimension error")
    end
    return W_V
end
#------------------------------------------------------------------------------
function pruneMat!(M::AbstractMatrix{T}, tol=1E-14) where T <: Union{Float64, ComplexF64}
    if eltype(M) == Float64
        for i in 1:size(M, 1)
            for j in 1:size(M, 2)
                if abs(M[i, j]) < tol
                    M[i, j] = 0.0
                end
            end
        end
    else
        for i in 1:size(M, 1)
            for j in 1:size(M, 2)
                realpart = real(M[i, j])
                imagpart = imag(M[i, j])
                realpart = (abs(realpart) < tol ? 0.0 : realpart)
                imagpart = (abs(imagpart) < tol ? 0.0 : imagpart)
                M[i, j] = realpart + im*imagpart
            end
        end
    end
end
function pruneMat!(M::SparseMatrixCSC{T, Int64}, tol=1E-14) where T <: Union{Float64, ComplexF64}
    droptol!(M, tol)
end

#------------------------------------------------------------------------------
# settings
#------------------------------------------------------------------------------
const ⊗ = kron
const range_=4

#-----------------------------------------------------------------------------#
# use the SU(3) generators to construct the general Hamiltonian
# λ_0 has been modified so that tr(λ_a λ_b) = 2 δ_ab valid for it
#-----------------------------------------------------------------------------#
const λ_0 = spzeros(ComplexF64, 3, 3); λ_0[1, 1] =λ_0[2, 2] = λ_0[3, 3] = sqrt(2.0/3)
const λ_1 = spzeros(ComplexF64, 3, 3); λ_1[1, 2], λ_1[2, 1] = 1, 1
const λ_2 = spzeros(ComplexF64, 3, 3); λ_2[1, 2], λ_2[2, 1] = -im, im
const λ_3 = spzeros(ComplexF64, 3, 3); λ_3[1, 1], λ_3[2, 2] = 1, -1
const λ_4 = spzeros(ComplexF64, 3, 3); λ_4[1, 3], λ_4[3, 1] = 1, 1
const λ_5 = spzeros(ComplexF64, 3, 3); λ_5[1, 3], λ_5[3, 1] = -im, im
const λ_6 = spzeros(ComplexF64, 3, 3); λ_6[2, 3], λ_6[3, 2] = 1, 1
const λ_7 = spzeros(ComplexF64, 3, 3); λ_7[2, 3], λ_7[3, 2] = -im, im
const λ_8 = spzeros(ComplexF64, 3, 3); λ_8[1, 1], λ_8[2, 2], λ_8[3, 3] = [1, 1, -2]/sqrt(3)
const  M1 = [λ_0, λ_1, λ_2, λ_3, λ_4, λ_5, λ_6, λ_7, λ_8]

#------------------------------------------------------------------------------
# define 4-body operator index list
#------------------------------------------------------------------------------
function getabcdList()
    abcdList = Vector{Int}[]
    #--------------------------------------------------------------------------
    # get all single, two-body and three-body operators
    #--------------------------------------------------------------------------
    # for a in 1:9
    #     for b in 1:9
    #         for c in 1:9
    #             for d in 1:9
    #                 push!(abcdList, [a, b, c, d])
    #             end
    #         end
    #     end
    # end
    # abcdList = abcdList[2:end]

    # --------------------------------------------------------------------------
    # get all pure 3-body operators
    #--------------------------------------------------------------------------
    for a in 2:9
        for b in 2:9
            for c in 2:9
                for d in 2:9
                    push!(abcdList, [a, b, c, d])
                end
            end
        end
    end

    return abcdList
end
const abcdList = getabcdList()
const dimLH = length(abcdList)

#------------------------------------------------------------------------------
# get MPS representation of the (N+1)-th scar state |Ψ_s⟩ = (Q^†)^N |G⟩
#------------------------------------------------------------------------------
function getScarMPS(L::Int, N::Int)
    if N < 1 || N > Int(floor(L/2))
        error("getScarMPS(): argument invalid")
    end
    #------------------------------------------------------------------------------
    # prepare A matrix for the AKLT ground state |G⟩
    #------------------------------------------------------------------------------
    # define standard A matrix for AKLT ground state: As[1],As[2],As[3] for m = -1,0,1
    As = [spzeros(Float64, 2, 2) for _ in 1:3]
    As[1][1, 2] = sqrt(2.0/3.0)
    As[2][1, 1], As[2][2, 2] = -1/sqrt(3.0), 1/sqrt(3.0)
    As[3][2, 1] = -sqrt(2.0/3.0)
    #------------------------------------------------------------------------------
    # MPO × MPS representation of (Q^†)^N |G⟩
    #------------------------------------------------------------------------------
    # MPS representation of Q^†, refer to (57) in Phys. Rev. B 98, 235156 (2018).
    M = Matrix{SparseMatrixCSC{Float64,Int64}}(undef, 3, 3)
    for n in 1:3
        for m in 1:3
            M[n, m] = spzeros(Float64, N+1, N+1)
            for i in 1:N
                M[n, m][i, i]   = (n == m ? (-1)^i : 0)
                M[n, m][i, i+1] = ((n == 3 && m == 1) ? 2*(-1)^i : 0)
            end
            M[n,m][N+1, N+1]  = (n == m ? (-1)^(N+1) : 0)
        end
    end
    # new MPS matrix from by MPO × MPS, here N for one time acting of (Q^†)^N on |G⟩
    AList   = [spzeros(Float64, 2*(N+1), 2*(N+1)) for _ in 1:3]
    A_LList = [spzeros(Float64, 2, 2*(N+1)) for _ in 1:3]
    A_RList = [spzeros(Float64, 2*(N+1), 2) for _ in 1:3]
    b_L = prepend!(zeros(N), [1.0])
    b_R = append!(zeros(N), [1.0])
    for n in 1:3
        for m in 1:3
            AList[n] .+= M[n, m] ⊗ As[m]
            A_LList[n] .+= (transpose(b_L)*M[n, m]) ⊗ As[m]
            A_RList[n] .+= (M[n, m]*b_R) ⊗ As[m]
        end
    end
    return AList, A_LList, A_RList
end

#------------------------------------------------------------------------------
# Calculate the correlation matrix of (Q^†)^N |G⟩ and its eigen-operators
#------------------------------------------------------------------------------
function getAarray(L::Int; N::Int)
    AList, A_LList, A_RList = getScarMPS(L, N)
    #  ----o----       depending on the position of the contracted A[m]
    #      |           we have EL, EA and ER
    #  ----o----
    # determine normalization
    EA = sum((A ⊗ A for A in AList))
    EL = sum((A ⊗ A for A in A_LList))
    ER = sum((A ⊗ A for A in A_RList))
    # EA_L_2 = EA^(L-2)
    normalization = sum(Diagonal(EL*EA^(L-2)*ER))   # normalization ⟨ψ | ψ⟩
    if normalization < 1E-4
        println("The state is a null state with ⟨ψ|ψ⟩=$normalization")
        return nothing
    end
    # store for later use
    EL6 = EA^(L-6)
    #--------------------------------------------------------------------------
    # determine opeartor average ⟨ψ | O | ψ⟩
    #--------------------------------------------------------------------------
    χ = size(AList[1], 1); χ2 = χ^2
    # AlambdasList[i][n] = A^{[n]}_{λᵢ}. The index sequence is defined in the
    # way that the outmost index, n, is in the outmost position for performance
    AlambdasList = [[spzeros(ComplexF64, χ, χ) for n in 1:3] for i in 1:9]
    for m in 1:3
        for i in 1:length(M1)
            for n in 1:3
                AlambdasList[i][m] .+= AList[n]*M1[i][m, n]
            end
        end
    end
    # transfer matrix from contraction of complex(A^{[n]}[λi]) and A^{[n]}[λj]
    Eij = Matrix{SparseMatrixCSC{ComplexF64,Int64}}(undef, 9, 9)
    for i in 1:9
        for j in 1:9
            Eij[i, j] = spzeros(ComplexF64, χ2, χ2)
            for n in 1:3
                Eij[i, j] .+= conj(AlambdasList[i][n]) ⊗ AlambdasList[j][n]
            end
        end
    end
    Ls = zeros(ComplexF64, dimLH)
    Mat = zeros(ComplexF64, dimLH, dimLH)
    # calculate Mat
    for j in 1:dimLH
        aj, bj, cj, dj = abcdList[j]
        Ls[j] = sum(Diagonal(EL*Eij[1, aj]*Eij[1, bj]*Eij[1, cj]*Eij[1, dj]*EL6*ER))/normalization
        for i in 1:j-1
            ai, bi, ci, di = abcdList[i]
            Mat[i, j]  = sum(Diagonal(EL*Eij[ai, aj]*Eij[bi, bj]*Eij[ci, cj]*Eij[di, dj]*EL6*ER))/normalization
            Mat[i, j] -= conj(Ls[i])*Ls[j]  # general form is used
            Mat[j, i]  = conj(Mat[i, j])
        end
        # diagonal term
        Mat[j, j]  = sum(Diagonal(EL*Eij[aj, aj]*Eij[bj, bj]*Eij[cj, cj]*Eij[dj, dj]*EL6*ER))/normalization
        Mat[j, j] -= abs(Ls[j])^2
    end
    vals, vects = eigen(Hermitian(Mat))
    nzero = count(x -> abs(x) < 1E-10, vals)
    println("Correlation Spectrum N₀ = $(nzero)")
    # return  nzero, vals, vects
    #----------------------------------------------------------------------
    # get λ of the eigen-operator A ψ = λ ψ
    #----------------------------------------------------------------------
    Aarray = vects[:, 1:nzero]
    lamvals = ComplexF64[]
    for n in 1:nzero
        lamval = 0.0
        for j in 1:dimLH
            lamval += Aarray[j, n]*Ls[j]
        end
        if abs(imag(lamval)) > 1E-12
            error("getAarray(): imag(λ) != 0")
        end
        # if we want to double check that  A |ψ⟩ = ξ |ψ⟩, one can calculate
        # ⟨ψ|A^† A|ψ⟩ and compare with |ξ|^2 ⟨ψ|ψ⟩
        push!(lamvals, real(lamval))
    end
    return Aarray, lamvals
end

#------------------------------------------------------------------------------
# use all the local operators to construct the semi-definite Hamiltonian
#------------------------------------------------------------------------------
function getProj(L::Int; N::Int)
    Aarray, lamvals = getAarray(L, N=N)
    Proj = zeros(ComplexF64, 3^range_, 3^range_)
    for (i, λ) in enumerate(lamvals)
        vec = Aarray[:, i]
        A = spzeros(ComplexF64, 3^range_, 3^range_)
        for i in 1:dimLH
            a, b, c, d = abcdList[i]
            A .+= vec[i]*(M1[a] ⊗ M1[b] ⊗ M1[c] ⊗ M1[d])
        end
        A .-= λ*speye(3^range_)
        Proj .+= A'*A
    end
    # rescale to make the biggest element be order 1
    Proj ./= maximum(abs.(Proj))
    pruneMat!(Proj, 1E-14)
    return Proj
end

L, N = 256, 58
L, N = 16, 3
L, N = 24, 5
Proj = getProj(L, N=N)
run(`mkdir -p data`)
writedlm("./data/Proj_k$(range_)_L$(L)_N$(N).txt", Proj)
