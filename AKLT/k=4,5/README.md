# Description of Files

## Files
* `corrMat_k4.jl`: Use MPO × MPS method to calculate the correlation matrix of
  a given state. In choosing a proper scar state, the constructed projection
  operators will annihilate all scar states. Here `_k4` stands for range-4
  oeprators. 
* `KerProj.jl`: semi-analytical way to determine the projection operators of a
  general range that annihilate all scar states.
