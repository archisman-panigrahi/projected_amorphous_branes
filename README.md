### Julia Codes for Projected Amorphous Branes by Archisman Panigrahi

The file `PTB_chern_BHZ_random_sites_joint_LDOS.jl` contains the codes to generate the Hamiltonian for Projected Amorphous Branes and compute its eigenspectrum, Bott index, and local Chern marker. The script can be systematically run with `run_joint_LDOS.sh`.

The file `PTB_chern_BHZ_random_joint_LDOS_disorder_chern.jl` contains the analogous codes for on-sight disorder potential. The script can be systematically run with `run_joint_LDOS_disorder.sh`.

The files `angmom.jl`, `generate_matrices2D.jl`, `generate_matrices1D.jl` as well as `Hermitian_Check.jl` contain library functions to generate angular momentum matrices as well as lattice tight binding matrices. They are called by `PTB_chern_BHZ_random_sites_joint_LDOS.jl`.

All these codes are written in the [Julia](https://julialang.org/) programming language.
