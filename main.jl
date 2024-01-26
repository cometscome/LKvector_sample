include("Chebyshev.jl")
include("LKvector.jl")
using .Chebyshev
using SparseArrays
using .LK

function calc_A(Nx, Ny, μ, Δ, aa)
    Ln = Nx * Ny * 2
    A = spzeros(Nx * Ny, Nx * Ny)

    for ix = 1:Nx
        for iy = 1:Ny
            #Diagonal element
            ii = (iy - 1) * Nx + ix
            jx = ix
            jy = iy
            jj = (jy - 1) * Nx + jx
            A[ii, jj] = -μ
            #+1 in x direction
            jx = ifelse(ix == Nx, 1, ix + 1)
            jy = iy
            jj = (jy - 1) * Nx + jx
            A[ii, jj] = -1.0
            #-1 in x direction
            jx = ifelse(ix == 1, Nx, ix - 1)
            jy = iy
            jj = (jy - 1) * Nx + jx
            A[ii, jj] = -1.0
            #+1 in y direction
            jx = ix
            jy = ifelse(iy == Ny, 1, iy + 1)
            jj = (jy - 1) * Nx + jx
            A[ii, jj] = -1.0
            #-1 in y direction
            jx = ix
            jy = ifelse(iy == 1, Ny, iy - 1)
            jj = (jy - 1) * Nx + jx
            A[ii, jj] = -1.0

        end
    end
    Abdg = [
        A Δ
        Δ' -conj.(A)
    ]

    #=
    for ii = 1:Nx*Ny
        for jj = 1:Nx*Ny
            A[ii+Nx*Ny, jj+Nx*Ny] = -conj(A[ii, jj])
            A[ii, jj+Nx*Ny] = Δ[ii, jj]
            A[ii+Nx*Ny, jj] = conj(Δ[jj, ii])
        end
    end
    =#

    return Abdg / aa

end


function main()
    μ = -1.5

    aa = 10
    Nx = 64 * 4
    Ny = 64 * 4
    Ln = 2 * Nx * Ny
    Δ0 = 0.5
    Δ = spzeros(Nx * Ny, Nx * Ny)
    for i = 1:Nx*Ny
        Δ[i, i] = Δ0
    end
    @time A = calc_A(Nx, Ny, μ, Δ, aa)

    left_i = 1
    right_j = 1
    nc = 200

    x = zeros(Float64, Ln)
    println("vector")
    vec_a = Chebyshev.calc_polynomials(nc, left_i, right_j, A)
    @time vec_a = Chebyshev.calc_polynomials(nc, left_i, right_j, A)
    println(vec_a[1:10])
    println(sum(abs.(vec_a)))
    println("fast vector")

    x = LKvector(Float64, Ln, 1e-6)

    #vec_a = Chebyshev.calc_polynomials(nc, left_i, right_j, A)
    #@time vec_a = Chebyshev.calc_polynomials(nc, left_i, right_j, A)
    #println(vec_a[1:10])
    vec_a = Chebyshev.calc_greenfunction(A, x, left_i, right_j; nc)
    @time vec_a = Chebyshev.calc_greenfunction(A, x, left_i, right_j; nc)
    println(vec_a[1:10])
    println(sum(abs.(vec_a)))

end
main()