include("Hermitian_Check.jl")

function generate_matrices_2D(Lx::Int64, Ly::Int64, periodic_x::Int64, periodic_y::Int64)
    nSites = Lx*Ly
    CX2D = zeros(nSites, nSites)
    CY2D = zeros(nSites, nSites)
    SX2D = zeros(ComplexF64,nSites, nSites)
    SY2D = zeros(ComplexF64,nSites, nSites)
    CXCY2D = zeros(nSites, nSites)
    SXCY2D = zeros(ComplexF64,nSites, nSites)
    CXSY2D = zeros(ComplexF64,nSites, nSites)

    ### Constant matrix
    const_2D = Diagonal(ones(nSites))

    ### 2Dx
    for jj = 1:Ly
        for ii = 1:Lx-1
            CX2D[(jj-1)*Lx + ii + 1, (jj-1)*Lx + ii] = 0.5;
            CX2D[(jj-1)*Lx + ii, (jj-1)*Lx + ii + 1] = 0.5;

            ### Definition of sin function: Assuming <R + a| Exp(I G a) |R> = 1,
            ### the sin function has to be <R + a| sin(G a) | R> = 1/I = -I
            SX2D[(jj-1)*Lx + ii + 1, (jj-1)*Lx + ii] = -0.5im;
            SX2D[(jj-1)*Lx + ii, (jj-1)*Lx + ii + 1] = 0.5im;
        end
    ### For periodic systems
        if(periodic_x == 1)
            CX2D[(jj-1)*Lx + 1, (jj-1)*Lx + Lx] = 0.5;
            CX2D[(jj-1)*Lx + Lx, (jj-1)*Lx + 1] = 0.5;

            SX2D[(jj-1)*Lx + 1, (jj-1)*Lx + Lx] = -0.5im;
            SX2D[(jj-1)*Lx + Lx, (jj-1)*Lx + 1] = 0.5im;
        end
    end

    ### 2Dy
    for ii = 1:Lx
        for jj = 1:Ly-1
            CY2D[jj*Lx + ii, (jj-1)*Lx + ii] = 0.5;
            CY2D[(jj-1)*Lx + ii, jj*Lx + ii] = 0.5;

            SY2D[jj*Lx + ii, (jj-1)*Lx + ii] = -0.5im;
            SY2D[(jj-1)*Lx + ii, jj*Lx + ii] = 0.5im;
        end
    ### For periodic systems
        if(periodic_y == 1)
            CY2D[ii, (Ly-1)*Lx + ii] = 0.5;
            CY2D[(Ly-1)*Lx + ii, ii] = 0.5;

            SY2D[ii, (Ly-1)*Lx + ii] = -0.5im;
            SY2D[(Ly-1)*Lx + ii, ii] = 0.5im;
        end
    end

    ## 2nd neighbor hopping

    ####################################
    ### Couplings which are like
                       ###
                       ###
                       ###
                       ###
                       ###
            ##############
            ##############

    ### Bulk
    for jj = 1:(Ly-1)
        for ii = 1:(Lx-1)
            CXCY2D[(jj-1 +1)*Lx + ii + 1, (jj-1)*Lx + ii] = 0.25;
            CXCY2D[(jj-1)*Lx + ii, (jj-1 +1)*Lx + ii + 1] = 0.25;

            SXCY2D[(jj-1 +1)*Lx + ii + 1, (jj-1)*Lx + ii] = -0.25im;
            SXCY2D[(jj-1)*Lx + ii, (jj-1 +1)*Lx + ii + 1] = 0.25im;

            CXSY2D[(jj-1 +1)*Lx + ii + 1, (jj-1)*Lx + ii] = -0.25im;
            CXSY2D[(jj-1)*Lx + ii, (jj-1 +1)*Lx + ii + 1] = 0.25im;
        end
    end
    ## Periodic in x
    if(periodic_x == 1)
        for jj = 1:(Ly-1)
            CXCY2D[(jj-1 +1)*Lx + 1, (jj-1)*Lx + Lx] = 0.25;
            CXCY2D[(jj-1)*Lx + Lx, (jj-1+1)*Lx + 1] = 0.25;

            SXCY2D[(jj-1 +1)*Lx + 1, (jj-1)*Lx + Lx] = -0.25im;
            SXCY2D[(jj-1)*Lx + Lx, (jj-1+1)*Lx + 1] = 0.25im;

            CXSY2D[(jj-1 +1)*Lx + 1, (jj-1)*Lx + Lx] = -0.25im;
            CXSY2D[(jj-1)*Lx + Lx, (jj-1+1)*Lx + 1] = 0.25im;
        end
    end
    ## Periodic in y
    if(periodic_y == 1)
        for ii = 1:(Lx-1)
            CXCY2D[0*Lx + ii + 1, (Ly-1)*Lx + ii] = 0.25;
            CXCY2D[(Ly-1)*Lx + ii, 0*Lx + ii + 1] = 0.25;

            SXCY2D[0*Lx + ii + 1, (Ly-1)*Lx + ii] = -0.25im;
            SXCY2D[(Ly-1)*Lx + ii, 0*Lx + ii + 1] = 0.25im;

            CXSY2D[0*Lx + ii + 1, (Ly-1)*Lx + ii] = -0.25im;
            CXSY2D[(Ly-1)*Lx + ii, 0*Lx + ii + 1] = 0.25im;
        end
    end
    ## Periodic in both
    if(periodic_x == 1 && periodic_y == 1)
        CXCY2D[1, Lx*Ly] = 0.25;
        CXCY2D[Lx*Ly, 1] = 0.25;

        SXCY2D[1, Lx*Ly] = -0.25im;
        SXCY2D[Lx*Ly, 1] = 0.25im;

        CXSY2D[1, Lx*Ly] = -0.25im;
        CXSY2D[Lx*Ly, 1] = 0.25im;
    end

    ####################################
    ### Couplings which are like
        ###############
        ###############
                    ###
                    ###
                    ###
                    ###
                    ###

    ### Bulk
    for jj = 2:Ly
        for ii = 1:(Lx-1)
            CXCY2D[(jj-1 -1)*Lx + ii + 1, (jj-1)*Lx + ii] = 0.25;
            CXCY2D[(jj-1)*Lx + ii, (jj-1 -1)*Lx + ii + 1] = 0.25;

            SXCY2D[(jj-1 -1)*Lx + ii + 1, (jj-1)*Lx + ii] = -0.25im;
            SXCY2D[(jj-1)*Lx + ii, (jj-1 -1)*Lx + ii + 1] = 0.25im;

            CXSY2D[(jj-1 -1)*Lx + ii + 1, (jj-1)*Lx + ii] = 0.25im;
            CXSY2D[(jj-1)*Lx + ii, (jj-1 -1)*Lx + ii + 1] = -0.25im;
        end
    end
    ## Periodic in x
    if(periodic_x == 1)
        for jj = 2:Ly
            CXCY2D[(jj-1 -1)*Lx + 1, (jj-1)*Lx + Lx] = 0.25;
            CXCY2D[(jj-1)*Lx + Lx, (jj-1-1)*Lx + 1] = 0.25;

            SXCY2D[(jj-1 -1)*Lx + 1, (jj-1)*Lx + Lx] = -0.25im;
            SXCY2D[(jj-1)*Lx + Lx, (jj-1-1)*Lx + 1] = 0.25im;

            CXSY2D[(jj-1 -1)*Lx + 1, (jj-1)*Lx + Lx] = 0.25im;
            CXSY2D[(jj-1)*Lx + Lx, (jj-1-1)*Lx + 1] = -0.25im;
        end
    end
    ## Periodic in y
    if(periodic_y == 1)
        for ii = 1:(Lx-1)
            CXCY2D[(Ly-1)*Lx + ii + 1, 0*Lx + ii] = 0.25;
            CXCY2D[0*Lx + ii, (Ly-1)*Lx + ii + 1] = 0.25;

            SXCY2D[(Ly-1)*Lx + ii + 1, 0*Lx + ii] = -0.25im;
            SXCY2D[0*Lx + ii, (Ly-1)*Lx + ii + 1] = 0.25im;

            CXSY2D[(Ly-1)*Lx + ii + 1, 0*Lx + ii] = 0.25im;
            CXSY2D[0*Lx + ii, (Ly-1)*Lx + ii + 1] = -0.25im;
        end
    end
    ## Periodic in both
    if(periodic_x == 1 && periodic_y == 1)
        CXCY2D[(Ly-1)*Lx+1, Lx] = 0.25;
        CXCY2D[Lx, (Ly-1)*Lx+1] = 0.25;

        SXCY2D[(Ly-1)*Lx+1, Lx] = -0.25im;
        SXCY2D[Lx, (Ly-1)*Lx+1] = 0.25im;

        CXSY2D[(Ly-1)*Lx+1, Lx] = 0.25im;
        CXSY2D[Lx, (Ly-1)*Lx+1] = -0.25im;
    end

    if (Hermitian_Check(CX2D)==1 && Hermitian_Check(CY2D)==1 && Hermitian_Check(SX2D)==1 && Hermitian_Check(SY2D)==1 && Hermitian_Check(CXCY2D)==1 && Hermitian_Check(SXCY2D)==1 && Hermitian_Check(CXSY2D)==1 && Hermitian_Check(const_2D)==1)
        println("Verified: Building Block Matrices are Hermitian")
        return const_2D, CX2D, SX2D, CY2D, SY2D, CXCY2D, SXCY2D, CXSY2D
    else
        println("Error!!! Building Block Matrices are not Hermitian!")
    end
end