#==============================================================================#
# Description
# -----------
#   This is an implementation of the recursive method in section III of the
#   supplementary material that determines the ground state dengeracy of a
#   projector Hamiltonian recursively.
#
#   Version:  1.0
#   Created:  2021-04-28 14:44
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
using Random
using Printf
using OffsetArrays
using SparseArrays
using LinearAlgebra
using DelimitedFiles

#-----------------------------------------------------------------------------#
# Utility functions
#-----------------------------------------------------------------------------#
function eye(n::Integer, ::Type{T}=Float64) where T <: Union{Int32, Int64, Float64, ComplexF64}
    if n < 1
        error("eye(): argument error")
    end
    return Matrix{T}(I, n, n)
end
function speye(n::Integer, ::Type{T}=Float64) where T <: Union{Int32, Int64, Float64, ComplexF64}
    if n < 1
        error("speye(): argument error")
    end
    return sparse(T(1)*I, n, n)
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
Random.seed!(1234)

#------------------------------------------------------------------------------
# global parameters
#------------------------------------------------------------------------------
const PBC = true
const q = 3       # Hilbert space dimension of each site
const range = 3   # operator range

#=============================================================================#
# use the SU(3) generators to construct the general Hamiltonian
# λ_0 has been modified so that tr(λ_a λ_b) = 2 δ_ab valid for it
#=============================================================================#
# The general form of range-2 Hamiltonian can be written as H = Σ c_ab^i λ_a^i λ_b^{i+1}
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
# define 3-body operator index list
#------------------------------------------------------------------------------
function getabcList()
    abcList = Vector{Int}[]
    for a in 2:9
        for b in 2:9
            for c in 2:9
                push!(abcList, [a, b, c])
            end
        end
    end
    return abcList
end
const abcList = getabcList()
const dimLH = length(abcList)

#------------------------------------------------------------------------------
# use all the local operators to construct the semi-definite Hamiltonian
#------------------------------------------------------------------------------
function construct3SiteA(vec::Vector{ComplexF64}, λ::ComplexF64)
    A = spzeros(ComplexF64, 3^3, 3^3)
    for i in 1:dimLH
        a, b, c = abcList[i]
        A .+= vec[i]*(M1[a] ⊗ M1[b] ⊗ M1[c])
    end
    A = A - λ*speye(3^3)
    return A
end

#------------------------------------------------------------------------------
# use all the local operators to construct the semi-definite Hamiltonian
#------------------------------------------------------------------------------
function get3SiteHamiltonian(Aarray::Matrix{ComplexF64}, Alamvals::Vector{ComplexF64}; L::Int, n::Int)
    Ham = spzeros(ComplexF64, 3^L, 3^L)
    nzero = size(Aarray, 2)
    if n <= L-2
        Ham3 = zeros(ComplexF64, 3^3, 3^3)
        for i in 1:nzero
            A = construct3SiteA(Aarray[:, i], Alamvals[i])
            Ham3 .+= rand()*A'*A
        end
        pruneMat!(Ham3, 1E-16)
        Ham = speye(3^(n-1)) ⊗ Ham3 ⊗ speye(3^(L-n-2))
    else
        #----------------------------------------------------------------------
        # add term for site (n-1)--n--1  &  n--1--2
        #----------------------------------------------------------------------
        Ham = get3SiteHamiltonian(Aarray, Alamvals, L=L, n=L-2)
        if n == L-1
            permutation = getOperatorPermutation(L, 1)
        else
            permutation = getOperatorPermutation(L, 2)
        end
        Ham = Ham[permutation, :]; Ham = Ham[:, permutation]
    end
    return Ham
end

#------------------------------------------------------------------------------
# smarter way of get value (1-based) of shifted qudits (3 here)
# shift > 0 : right shift operator --> H[iL i1 ... i(L-1), jL j1 ... j(L-1)]
#------------------------------------------------------------------------------
function nextBasis!(L::Int, basisvec::Vector{Int})
    for i in L:-1:1
        if (basisvec[i] < 2)
            basisvec[i] += 1
            basisvec[i+1:L] .= 0
            return basisvec
        end
    end
end
#------------------------------------------------------------------------------
function getOperatorPermutation(L::Int, shift::Int)
    sgn = sign(shift); shift = abs(shift)
    shiftweight= 3^shift
    weightList = [3^(L-i) for i in 1:L]
    basisvec = zeros(Int, L)
    permutation = zeros(Int, 3^L)
    for idx in 1:3^L
        index_shift = (idx - 1) ÷ shiftweight
        index_shift += sum([weightList[i]*basisvec[L-shift+i]  for i in 1:shift]) + 1
        permutation[idx] = index_shift
        nextBasis!(L, basisvec)
    end
    #--------------------------------------------------------------------------
    # Permutation[i1 i2 .... iL] =  iL i1 ... i(L-1)  for shift = 1
    # For operator, we need to inverse H[iL i1 ... i(L-1), jL j1 ... j(L-1)]
    #--------------------------------------------------------------------------
    return sgn > 0 ? sortperm(permutation) : permutation
end

function constructHam(Aarray::Matrix{ComplexF64}, Alamvals::Vector{ComplexF64}, L::Int; hasPBC::Bool)
    Ham = spzeros(ComplexF64, 3^L, 3^L)
    for n in 1:L-range+1
        Ham .+= get3SiteHamiltonian(Aarray, Alamvals, L=L, n=n)
    end
    if hasPBC
        for n in L-range+2:L
            Ham .+= get3SiteHamiltonian(Aarray, Alamvals, L=L, n=n)
        end
    end
    return Ham
end

#------------------------------------------------------------------------------
# Utility Program
#------------------------------------------------------------------------------
function getSwapIndices(indices::Vector{Int}, m::Int, n::Int)
    #----------------------------------------------------------------------
    # Let A and B be two Hilbert spaces of dimensions m and n,
    # respectively. The basis for A\otimes B or
    # B\otimes A are chosen according to the canonical Kronecker product
    # convention. For example, basis vectors of A\otimes B are ordered as
    # |1>|1>,|1>|2>,...,|1>|n>,|2>|1>,|2>|2>,...,|2>|n>,...,|m>|n>.
    #
    # the input indices is a column of basis vector indices in A\otimes B,
    # the output NewIndices is the corresponding basis vector indices in
    # B\otimes A.
    #----------------------------------------------------------------------
    iB, iA = Base._ind2sub((n,m), indices)
    NewIndices = Base._sub2ind((m,n), iA, iB)
    return NewIndices
end
function getSwapOperator(M, dimA::Int, dimB::Int)
    #----------------------------------------------------------------------
    # Let A and B be two Hilbert spaces. There is an isomorphism between
    # A\otimes B and B\otimes A. Let M be an operator acting on A\otimes B
    # under the kron product representation, Mnew is the corresponding
    # operator acting on B\otimes A.
    #----------------------------------------------------------------------
    dim = dimA*dimB
    ijList = findall(!iszero, M)
    values = M[ijList]
    iList = [x[1] for x in ijList]
    jList = [x[2] for x in ijList]
    iListnew  = getSwapIndices(iList, dimA, dimB)
    jListnew  = getSwapIndices(jList, dimA, dimB)
    return sparse(iListnew, jListnew, values, dim, dim)
end

#------------------------------------------------------------------------------
# get zero-energy manifold
#------------------------------------------------------------------------------
function getdspace(H)
    eigvals, eigvecs = eigen(Hermitian(Matrix(H)))
    nzero = count(x -> x < 1E-12, eigvals)
    return nzero, eigvecs[:, 1:nzero]
end

#=============================================================================#
# Determine the degeneracy of ground states of H with OBC and PBC if specified.
#=============================================================================#
function KerProjH(Aarray::Matrix{ComplexF64}, Alamvals::Vector{ComplexF64}, LMin::Int, LMax::Int)
    # first define the projection operator P that acts on 3 sites
    P = get3SiteHamiltonian(Aarray, Alamvals, L=3, n=1)

    nzero_OBC = zeros(Int, LMax)
    if PBC
        nzero_PBC = zeros(Int, LMax)
    end
    dspace = nothing   # to be used later

    rightProj = Vector{Any}(undef, LMax)
    leftProj  = Vector{Any}(undef, LMax)

    #--------------------------------------------------------------------------
    # Compute the kernel dimension.
    #--------------------------------------------------------------------------
    for L in LMin:LMax
        #----------------------------------------------------------------------
        # Minimal system size.
        #----------------------------------------------------------------------
        if L == LMin
            H_OBC = constructHam(Aarray, Alamvals, LMin, hasPBC=false)
            nzero_OBC[L], dspace = getdspace(H_OBC)
            rightProj[L] = dspace
            leftProj[L]  = dspace

            if PBC
                H_PBC = H_OBC
                for n = L-range+2:L
                    H_PBC .+= get3SiteHamiltonian(Aarray, Alamvals, L=L, n=n)
                end
                nzero_PBC[L], = getdspace(H_PBC)
            end
        else
            #------------------------------------------------------------------
            # Larger system sizes; OBC Left extension calculation.
            #------------------------------------------------------------------
            if L < LMin + range
                # the matrix 'space' here is a particular choice of orthonormal basis.
                Previousdspace = dspace
                # For small system sizes, construct the full Hamiltonian.
                H_OBC = P ⊗ eye(q^(L-range))
                # Tensor the previous dspace with one more site.
                subspace = eye(q) ⊗ Previousdspace
            else # L >= LMin + range
                # Construct the Hamiltonian in a small subspace.
                H_OBC = P ⊗ eye(nzero_OBC[L-range])
                subspace = 1
                for j = 1:(range-1)
                    subspace = subspace*leftProj[L-range+j]
                    subspace = eye(q) ⊗ subspace
                end
            end

            # Project the Hamiltonian to the smaller subspace.
            H_OBC_Sub = Hermitian(subspace'*H_OBC*subspace)
            nzero_OBC[L], dspace_1Step = getdspace(H_OBC_Sub)

            # Update dspace.
            if L < LMin + range # dspace is only needed for small system sizes.
                dspace = subspace*dspace_1Step; # Embed back to the full Hilbert space.
            end
            # Save the single-step embedding matrix.
            leftProj[L] = dspace_1Step

            if PBC
                #--------------------------------------------------------------
                # PBC calculation.
                #--------------------------------------------------------------
                if L < LMin + range
                    H_PBC = spzeros(ComplexF64, q^L, q^L)
                    for n = L-range+2:L
                        H_PBC .+= get3SiteHamiltonian(Aarray, Alamvals, L=L, n=n)
                    end
                    H_PBC_Sub = Hermitian(dspace'*H_PBC*dspace)
                    # Be careful not to rewrite dspace here.
                    nzero_PBC[L],  = getdspace(H_PBC_Sub)
                else # L >= LMin + range
                    H_PBC_Sub = zeros(nzero_OBC[L], nzero_OBC[L])
                    for j = 1:(range-1)
                        M = P ⊗ eye(nzero_OBC[L-range])
                        dimA = q^j
                        dimB = (q^(range-j))*nzero_OBC[L-range]
                        M = getSwapOperator(M, dimA, dimB)
                        subspace = rightProj[L-range+1]
                        for iR = 2:j
                            subspace = subspace ⊗ eye(q)
                            subspace = subspace*rightProj[L-range+iR]
                        end
                        for iL = (j+1):range
                            subspace = eye(q) ⊗ subspace
                            subspace = subspace*leftProj[L-range+iL]
                        end
                        H_PBC_Sub = H_PBC_Sub + subspace'*M*subspace
                    end
                    H_PBC_Sub = Hermitian(H_PBC_Sub)
                    # Be careful not to rewrite dspace here.
                    nzero_PBC[L], = getdspace(H_PBC_Sub)
                end
                #--------------------------------------------------------------
                # OBC Right extension calculation.
                #--------------------------------------------------------------
                if L < LMin + range
                    # Previousdspace has been defined, current dspace is also updated.
                    rightProj[L] = (Previousdspace ⊗ eye(q))'*dspace
                else # L> = LMin + range
                    rightProj[L] = (leftProj[L-1] ⊗ eye(q))'*((eye(q) ⊗ rightProj[L-1])*leftProj[L])
                end
            end
        end
    end
    return nzero_OBC, nzero_PBC
    # return nzero_OBC, nzero_PBC, leftProj, rightProj
end

function getDegeneracyList()
    LMin, LMax = 6, 512
    Aarray = readdlm("./data/Aarray.txt",  '\t', ComplexF64, '\n')
    Alamvals = readdlm("./data/Alamvals.txt",  '\t', ComplexF64, '\n')[:, 1]
    nzero_OBC, nzero_PBC = KerProjH(Aarray, Alamvals, LMin, LMax)
    writedlm("./data/degeneracy.dat", [LMin:LMax nzero_OBC[LMin:LMax] nzero_PBC[LMin:LMax]])
end

getDegeneracyList()
