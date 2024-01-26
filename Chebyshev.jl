module Chebyshev
using LinearAlgebra


function calc_greenfunction(A, x::T, left_i, right_j; nc=200) where {T}
    tempvector = Vector{T}(undef, 3)
    for i = 1:3
        tempvector[i] = zero(x)
    end
    vec_jnmm = tempvector[1]
    vec_jnm = tempvector[2]
    vec_jn = tempvector[3]
    vec_jn[right_j] = 1
    vec_ai = zeros(eltype(A), nc)

    @inbounds for nn = 0:nc-1
        if nn == 0
            vec_jnmm[right_j] = 1
        elseif nn == 1
            mul!(vec_jnmm, A, vec_jnm)
        else
            #mul!(C, A, B, α, β) -> C , A B α + C β
            mul!(vec_jnmm, A, vec_jnm, 2, -1)
        end

        vec_ai[nn+1] = vec_jnmm[left_i]
        #println(nn, "\t", vec_ai[nn+1], "\t", vec_jnmm.nonzeros)
        vec_jnm, vec_jnmm = vec_jnmm, vec_jnm
    end
    return vec_ai
end

function calc_polynomials(nc, left_i, right_j, A)
    Ln = size(A, 1)
    typeA = typeof(A[1, 1])
    vec_jnmm = zeros(typeA, Ln)
    vec_jnm = zeros(typeA, Ln)
    vec_jn = zeros(typeA, Ln)
    vec_jn[right_j] = 1.0
    vec_ai = zeros(typeA, nc)
    @inbounds for nn = 0:nc-1
        if nn == 0
            vec_jn[right_j] = 1.0
        elseif nn == 1
            mul!(vec_jn, A, vec_jnm)
            #A_mul_B!(vec_jn,A,vec_jnm)
        else
            mul!(vec_jn, A, vec_jnm)
            #A_mul_B!(vec_jn,A,vec_jnm)
            #                vec_jn *= 2
            for i = 1:Ln
                vec_jn[i] = vec_jn[i] * 2 - vec_jnmm[i]
            end
        end

        vec_ai[nn+1] = vec_jn[left_i]
        #println("vec_a, ",vec_ai[nn+1])
        for i = 1:Ln
            vec_jnmm[i] = vec_jnm[i]
            vec_jnm[i] = vec_jn[i]
        end


    end
    return vec_ai
end
end