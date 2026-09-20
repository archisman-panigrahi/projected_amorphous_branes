include("Hermitian_Check.jl")

function generate_matrices_1D(Lx::Int64, periodic_x::Int64)
    nSites = Lx
    CX1D = zeros(nSites, nSites)
    SX1D = zeros(ComplexF64,nSites, nSites)

    ### Constant matrix
    const_1D = Diagonal(ones(nSites))

    ## Bulk
    for ii = 1:(Lx-1)
        CX1D[ii+1,ii] = 0.5;
        CX1D[ii,ii+1] = 0.5;

        SX1D[ii+1,ii] = -0.5im;
        SX1D[ii,ii+1] = 0.5im;
    end

    ## boundary
    if(periodic_x == 1)
        CX1D[1,Lx] = 0.5;
        CX1D[Lx,1] = 0.5;

        SX1D[1,Lx] = -0.5im;
        SX1D[Lx,1] = 0.5im;
    end

    if (Hermitian_Check(CX1D)==1 && Hermitian_Check(SX1D)==1 && Hermitian_Check(const_1D)==1)
        println("Verified: Building Block Matrices are Hermitian")
        return const_1D, CX1D, SX1D
    else
        println("Error!!! Building Block Matrices are not Hermitian!")
    end
end