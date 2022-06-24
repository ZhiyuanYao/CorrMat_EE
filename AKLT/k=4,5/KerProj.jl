#==============================================================================#
# Description
# -----------
#   This program is to construct the range_-k projection operators that
#   annihilates all the SGA scar states Q†^N |G⟩ semi-analytically. The key
#   observation is that for a projector Hamiltonian, the ground state |G⟩ can
#   be decomposed in the following way
#
#                         |G⟩ = ∑_η | Ψ_k^η⟩ ⊗ | Φ ⟩
#
#   where η=1, 2, 4, 4 labels four ground state of a k-site AKLT model with OBC.
#   Therefore, all we need to do is the find the complement space of space
#   formed by { |Ψ_k^η⟩, Q† |Ψ_k^η⟩, (Q†)^2 |Ψ_k^η⟩, ⋯ }.
#   Note my way of getting Proj differs from that of Shang Liu's a little bit.
#   But the final result should be the same.
#
#   After we identified the projection operators, we can then use the recursive
#   method to find the ground state degeneracy of the auxiliary Hamiltonian.
#
#   Version:  1.0
#   Created:  2022-06-14 13:40
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
# Change Log:
#   *
#==============================================================================#
# todo:
#   *
#==============================================================================#
using Printf
using OffsetArrays
using SparseArrays
using LinearAlgebra
using MoreLinearAlgebra
using DelimitedFiles

#------------------------------------------------------------------------------
# settings
#------------------------------------------------------------------------------
const ⊗ = kron
# const range_=4

#==============================================================================#
# blocks to get the Hamiltonian of the system
#==============================================================================#
# note the basis index follows my Fortran convention
# |-1⟩ = [1, 0, 0], |0⟩ = [0, 1, 0] and |1⟩ = [0, 0, 1]
const S_x = spzeros(ComplexF64, 3, 3); S_x[1,2] = S_x[2,1] = S_x[2,3] = S_x[3,2] = 1/sqrt(2.0)
const S_y = spzeros(ComplexF64, 3, 3); S_y[1,2] = S_y[2,3] =  im/sqrt(2.0)
                                       S_y[2,1] = S_y[3,2] = -im/sqrt(2.0)
const S_z = spzeros(ComplexF64, 3, 3); S_z[1,1], S_z[3, 3] = -1, 1
#-----------------------------------------------------------------------------#
# construct the pure AKLT model
#-----------------------------------------------------------------------------#
function constructAKLTHam(L::Int, hasPBC::Bool)
    Ham = spzeros(ComplexF64, 3^L, 3^L)
    #=========================================================================#
    # H =  Σ[1/2*S_i*S_j + 1/6*(S_i*S_j)^2 + 1/3] = Σ P_2(ij)
    #=========================================================================#
    SS = kron(S_x, S_x) + kron(S_y, S_y) + kron(S_z, S_z)
    S2 = SS/2 + SS*SS/6 + sparse(Matrix(I, 9, 9))/3
    for n = 1:L-1
        Ham .+= speye(3^(n-1)) ⊗ S2 ⊗ speye(3^(L-n-1))
    end
    if hasPBC
        H1L = spzeros(ComplexF64, 3^L, 3^L)
        H1L .+= S_x ⊗ speye(3^(L-2)) ⊗ S_x
        H1L .+= S_y ⊗ speye(3^(L-2)) ⊗ S_y
        H1L .+= S_z ⊗ speye(3^(L-2)) ⊗ S_z
        H1L   = H1L/2 + (H1L*H1L)/6 + sparse(I, 3^L, 3^L)/3
        Ham .+= H1L
    end
    return Ham
end

#------------------------------------------------------------------------------
# construct the space of  |Ψ_k^η⟩, Q†*|Ψ_k^η⟩
#------------------------------------------------------------------------------
function constructStateSpace(range_::Int)
    L, hasPBC = range_, false
    eigvals, eigvects = eigen(Hermitian(Matrix(constructAKLTHam(L, hasPBC))))
    dspace = eigvects[:, 1:4]
    #--------------------------------------------------------------------------
    Qdagger = spzeros(Float64, 3^range_, 3^range_)
    Splus   = real.(S_x + im*S_y)
    for n in 1:range_
        Qdagger .+= (-1)^n*speye(3^(n-1)) ⊗ (Splus*Splus) ⊗ speye(3^(range_-n))
    end
    #--------------------------------------------------------------------------
    N_max = floor(Int, (range_ + 1)/2)
    for N in 1:N_max
        dspace = [dspace (Qdagger)^N*eigvects[:, 1:4]]
    end
    return gram_schmidt(dspace)
end

#------------------------------------------------------------------------------
# get the complement of dspace, so that each vector is orthogonal to all scar
# states of arbitrary system size
#------------------------------------------------------------------------------
function getProj(range_::Int)
    dspace = constructStateSpace(range_)
    Kernal = complement_span(eye(3^range_, ComplexF64), dspace)
    Proj = zeros(ComplexF64, 3^range_, 3^range_)
    for i in 1:size(Kernal, 2)
        Proj .+= Kernal[:, i]*(Kernal[:,i])'
    end
    Proj ./= maximum(abs.(Proj))
    pruneMat!(Proj)
    run(`mkdir -p ./data`)
    writedlm("./data/Proj_k$(range_).txt", Proj)
    return Proj
end

#------------------------------------------------------------------------------
# check whether the last state |↑↑↑↑↑⟩ from (Q^†)^N_max |G_k,η⟩ belongs to the
# space spanned by {|G_k,η⟩, Q^† |G_k,η⟩, ... (Q^†)^(N_max-1) |G_k,η⟩}
#------------------------------------------------------------------------------
function check_span()
    for range_ in [3, 5, 7]
        L, hasPBC = range_, false
        eigvals, eigvects = eigen(Hermitian(Matrix(constructAKLTHam(L, hasPBC))))
        dspace = eigvects[:, 1:4]
        #--------------------------------------------------------------------------
        Qdagger = spzeros(Float64, 3^range_, 3^range_)
        Splus   = real.(S_x + im*S_y)
        for n in 1:range_
            Qdagger .+= (-1)^n*speye(3^(n-1)) ⊗ (Splus*Splus) ⊗ speye(3^(range_-n))
        end
        #--------------------------------------------------------------------------
        N_max = Int(range_ + 1)/2 - 1
        for N in 1:N_max
            dspace = [dspace (Qdagger)^N*eigvects[:, 1:4]]
        end
        N_max += 1
        if rank(dspace) != rank([dspace (Qdagger)^N*eigvects[:, 1:4]])
            println("k = $(range_) failed")
            return
        else
            println("k = $(range_) passed")
        end
   end
end

#------------------------------------------------------------------------------
# compare projector space from that obtained from CorrMat method
#------------------------------------------------------------------------------
function compare_Proj(L::Int, N::Int, range_::Int)
    Proj = readdlm("./data/Proj_k$(range_)_L$(L)_N$(N).txt", '\t', ComplexF64, '\n')
    eigvals, eigvects = eigen(Hermitian(Proj))
    nzero = count(x -> abs(x) < 1E-10, eigvals)
    Proj = eigvects[:, nzero+1:end]

    Proj_k = readdlm("./data/Proj_k$(range_).txt",  '\t', ComplexF64, '\n')
    eigvals, eigvects = eigen(Hermitian(Proj_k))
    nzero = count(x -> abs(x) < 1E-10, eigvals)
    Proj_C = eigvects[:, nzero+1:end]

    if size(Proj_C, 2) ==  size(Proj, 2)  && size(complement_span(Proj_C, Proj), 2) == 0 &&  size(complement_span(Proj, Proj_C), 2) == 0
        println("Projector space V_O are the same")
        return true
    else
        println("Projector space V_O are different")
        return false
    end
end
