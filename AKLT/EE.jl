#==============================================================================#
# Description
# -----------
#   This is program use MPO × MPS method to calculate the η-pairing states in
#   the AKLT model, |Ψ_s⟩ = (Q^†)^N |G⟩, and then calculate its energy and
#   entanglement entropy.
#
#   Version:  1.0
#   Created:  2021-11-10 09:01
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
using OffsetArrays
using SparseArrays
using LinearAlgebra
using DelimitedFiles
using ITensors

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
# settings
#------------------------------------------------------------------------------
const ⊗ = kron
const Sx = spzeros(Complex{Float64}, 3, 3); Sx[1,2] = Sx[2,1] = Sx[2,3] = Sx[3,2] = 1/sqrt(2.0)
const Sy = spzeros(Complex{Float64}, 3, 3); Sy[1,2] = Sy[2,3] =  im/sqrt(2.0)
                                            Sy[2,1] = Sy[3,2] = -im/sqrt(2.0)
const Sz = spzeros(Int, 3, 3); Sz[1,1], Sz[3, 3] = -1, 1
const Ss = [Sx, Sy, Sz]

#------------------------------------------------------------------------------
# get the MPS representation of the scar state |Ψ_s⟩ = (Q^†)^N |G⟩
#------------------------------------------------------------------------------
function getScarMPS(L::Int, N::Int)
    if isodd(L) || N > Int(L/2)
        error("getScarMPS(): argument invalid")
    end
    #------------------------------------------------------------------------------
    # prepare A matrix for the AKLT ground state |G⟩
    #------------------------------------------------------------------------------
    # define standard A matrix for AKLT ground state: As[1],As[2],As[3] for m = -1,0,1
    As = [spzeros(Float64, 2, 2) for i in 1:3]
    As[1][1, 2] = sqrt(2.0/3.0)
    As[2][1, 1], As[2][2, 2] = -1/sqrt(3.0), 1/sqrt(3.0)
    As[3][2, 1] = -sqrt(2.0/3.0)
    if N == 0
        return As, As, As
    end
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
    AList   = [spzeros(Float64, 2*(N+1), 2*(N+1)) for i in 1:3]
    A_LList = [spzeros(Float64, 2, 2*(N+1)) for i in 1:3]
    A_RList = [spzeros(Float64, 2*(N+1), 2) for i in 1:3]
    b_L = prepend!(zeros(N), [1.0])
    b_R = append!(zeros(N), [1.0])
    for n in 1:3
        for m in 1:3
            AList[n] .+= M[n, m] ⊗ As[m]
            A_LList[n] .+= (transpose(b_L)*M[n, m]) ⊗ As[m]
            A_RList[n] .+= (M[n, m]*b_R) ⊗ As[m]
        end
    end
    return A_LList, AList, A_RList
end

#------------------------------------------------------------------------------
# Calculate the energy of the scar state (Q^†)^N |G⟩
#------------------------------------------------------------------------------
function getEnergy_MPS(L::Int, A_LList::T, AList::T, A_RList::T) where T <: Vector{SparseMatrixCSC{Float64,Int64}}
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
    #--------------------------------------------------------------------------
    # determine opeartor average ⟨ψ | O | ψ⟩
    #--------------------------------------------------------------------------
    χ = size(AList[1], 1); χ2 = χ^2
    EO = spzeros(ComplexF64, χ2, χ2)
    # (S * S)/2  term
    #  -------o-------
    #         |
    #         □   E is the generalized transfer matrix for operator O
    #         |
    #  -------o-------
    for i in 1:3
        E = spzeros(ComplexF64, χ2, χ2)
        for m = 1:3
            for n = 1:3
                E .+= (AList[m] ⊗ AList[n]) * Ss[i][m, n]
            end
        end
        EO .+= E*E/2
    end
    # (S * S)^2/6  term
    for i in 1:3
        for j in 1:3
            E = spzeros(ComplexF64, χ2, χ2)
            for m = 1:3
                for n = 1:3
                    E .+= (AList[m] ⊗ AList[n]) * (Ss[i]*Ss[j])[m, n]
                end
            end
            EO .+= E*E/6
        end
    end
    # times L since only one term of the Hamiltonian is calculated
    return L*(sum(Diagonal(EL*EO*EA^(L-4)*ER))/normalization + 1/3.0)
end

#-----------------------------------------------------------------------------#
# convert PBC MPS to OBC MPS format
#-----------------------------------------------------------------------------#
function convertMPS_PBC2OBC(A_LList::T, AList::T, A_RList::T) where T <: Vector{SparseMatrixCSC{Float64,Int64}}
    q = length(A_LList); D = size(A_LList[1], 1); χ = size(A_LList[1], 2)
    if length(AList) != q || length(A_RList) != q
        error("convertMPS_PBC2OBC(): physical dimension of argument error")
    elseif size(A_LList[1]) != reverse(size(A_RList[1]))
        error("convertMPS_PBC2OBC(): boundary array dimension does not match")
    elseif χ != size(AList[1], 1) || χ != size(AList[1], 2)
        error("convertMPS_PBC2OBC(): argument A array dimension does not match")
    end
    A_LList_OBC = [transpose(A_LList[n])[:] for n in 1:q]
    AList_OBC   = [speye(D) ⊗ AList[n] for n in 1:q]
    A_RList_OBC = [A_RList[n][:] for n in 1:q]
    #--------------------------------------------------------------------------
    # normalize A so that the weight of ψ does not blow up
    #--------------------------------------------------------------------------
    # method 1: Σₘ A^m† A^m = λ I
    # λ = sum([sum(Diagonal(A'*A)) for A in AList_OBC])/size(AList[1], 1)
    # AList_OBC = [A/sqrt(abs(λ)) for A in AList_OBC]
    #--------------------------------------------------------------------------
    # method 2: E = Σₘ A^m† ⊗ A^m  then find the largest eigenvalue
    #--------------------------------------------------------------------------
    # E = sum([(conj(A) ⊗ A) for A in AList_OBC])
    # eigvals, = eigen(Matrix(E))
    # λ = maximum(abs.(eigvals))
    # AList_OBC_new = AList_OBC/sqrt(abs(λ))  # AList_OBC_new = [A/sqrt(abs(λ)) for A in AList_OBC]
    # F = sum([(conj(A) ⊗ A) for A in AList_OBC_new])
    # eigvals, = eigen(Matrix(F))
    # λ = maximum(abs.(eigvals))
    # println("---", λ, "---")
    return A_LList_OBC, AList_OBC, A_RList_OBC
end

#------------------------------------------------------------------------------
# convert A_LList, AList, A_RList representation to ITensor MPS
#------------------------------------------------------------------------------
function getITensorMPS(N::Int, A_LList::T1, AList::T, A_RList::T1) where {T <: Vector{SparseMatrixCSC{Float64,Int64}}, T1 <: Vector{SparseVector{Float64,Int64}}}
    q = size(A_LList,  1)
    χ = size(AList[1], 1)
    sites = siteinds(q, N)
    psi = MPS(sites, linkdims=χ)
    #--------------------------------------------------------------------------
    # set value of psi
    #--------------------------------------------------------------------------
    TL = spzeros(dims(psi[1])...)
    Tm = Array{Float64}(undef, dims(psi[2])...)
    TR = spzeros(dims(psi[N])...)
    for m in 1:q
        TL[:, m] = A_LList[m]
        TR[:, m] = A_RList[m]
        Tm[:, m, :] = AList[m]
    end
    psi[1] = ITensor(TL, inds(psi[1]))    # Array(psi[1], inds(psi[1])...)
    psi[N] = ITensor(TR, inds(psi[N]))    # Array(psi[N], inds(psi[N])...)
    for i in 2:N-1
        psi[i] =  ITensor(Tm, inds(psi[i]))
    end
    return psi
end

#-----------------------------------------------------------------------------#
# get the entanglement of subsystem of size LA for OBC MPS with left/right
# boundary vector A_L/A_R and middle site-independent matrix A
#-----------------------------------------------------------------------------#
function getEntropy(L::Int, LA::Int, A_LList::T1, AList::T, A_RList::T1) where {T <: Vector{SparseMatrixCSC{Float64,Int64}}, T1 <: Vector{SparseVector{Float64,Int64}}}
    psi = getITensorMPS(L, A_LList, AList, A_RList)
    # normalization
    norm2 = inner(psi, psi)
    if norm2 < 1E-4
        println("The state is a null state with ⟨ψ|ψ⟩=$norm2")
        return nothing
    end
    b = Int(L/2)
    orthogonalize!(psi, b)
    # println("orthogonalize -> ", inner(psi, psi))
    U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
    S_von::Float64 = 0.0
    for n in 1:dim(S, 1)
        p = S[n, n]^2/norm2
        if (p < 1E-17)
            break
        end
        S_von += -p*log(p)
    end
    return S_von
end

function getEntropyMaxList()
    LList = [6, 8, 10, 12, 14, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    entropyMaxList = zeros(length(LList))
    indexMaxList   = zeros(Int, length(LList))
    for (iL, L) in enumerate(LList)
        println("In calculating entropy for L=$L")
        Nmax = iseven(Int(L/2)) ? Int(L/2) : Int(L/2) - 1
        entropyList = OffsetVector(zeros(Nmax+1), 0:Nmax)
        for N in 0:Nmax
            A_LList, AList, A_RList = getScarMPS(L, N)
            # energy = getEnergy_MPS(L, A_LList, AList, A_RList)
            # println("energy = $(energy)")
            A_LList_OBC, AList_OBC, A_RList_OBC = convertMPS_PBC2OBC(A_LList, AList, A_RList)
            entropyList[N] = getEntropy(L, Int(L/2), A_LList_OBC, AList_OBC, A_RList_OBC)
            # println("N = $(N), entropy = $(entropyList[N])")
        end
        entropyMaxList[iL], indexMaxList[iL] = findmax(entropyList)
    end
    run(`mkdir -p ./data`)
    writedlm("./data/entropyMaxITensor.dat", [LList indexMaxList entropyMaxList])
end

getEntropyMaxList()
