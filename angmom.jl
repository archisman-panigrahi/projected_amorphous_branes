### Generates spin matrices for arbitrary spin-j
eye2 = [1 0;
        0 1]

function generate_spin_matrices(j)
    if(isinteger(2*j) == false)
        println("Spin j must be half-integer or integer")
        println("Here, j was ",j)
        return
    end
    Jz = zeros(Int(2*j+1), Int(2*j+1))
    Jx = zeros(Int(2*j+1), Int(2*j+1))
    Jy = zeros(Int(2*j+1), Int(2*j+1))
    J_plus = zeros(Int(2*j+1), Int(2*j+1))
    J_minus = zeros(Int(2*j+1), Int(2*j+1))

    r = deepcopy(j)
    for a=1:Int(2*j)
        Jz[a,a] = r;
        J_minus[a+1,a] = sqrt(j*(j+1) - r*(r-1))
        r = r-1
    end
    Jz[Int(2*j+1), Int(2*j+1)] = -j

    J_plus = J_minus'
    
    Jx = (J_plus + J_minus)/2
    Jy = (J_plus - J_minus)/(2*im)

    return Jx,Jy,Jz
end