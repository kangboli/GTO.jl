export make_r_tensor

const IntegralKey = Tuple{Float64,Vector{Float64}}

function batched_two_electron_integral(hermites::Vector{HermiteGaussian})

    r_tensor_sizes = Dict{IntegralKey,Vector{Int}}()

    for (i, p) in enumerate(hermites)
        for (j, q) in enumerate(hermites)
            i <= j || continue
            αc = α(p) * α(q) / (α(p) + α(q))
            pq = c(p) - c(q)
            required_size = [i(p) + i(q), j(p) + j(q), k(p) + k(q)]
            r_tensor_sizes[(αc, pq)] = map(max, get(r_tensor_sizes,
                    (αc, pq), [1, 1, 1]), required_size)
        end
    end

    batched_integrals = Dict{IntegralKey,Array}(
        k => make_r_tensor(k..., r_tensor_sizes[k]) for k in keys(r_tensor_sizes))

    return batched_integrals
end

"""
"""
function make_r_tensor(α::Float64, pq::Vector{Float64}, r_tensor_size::Vector{Int})
    i_max, j_max, k_max = r_tensor_size
    n_max = sum(r_tensor_size)
    r_tensor = zeros(Float64, n_max, i_max, j_max, k_max)
    a, b, c = pq
    T = α * (a^2 + b^2 + c^2)

    for n = 0:n_max-1
        r_tensor[n+1, 1, 1, 1] = (-2α)^n * F(n, T)
    end

    for k = 1:k_max-1
        n_lim = n_max - k
        r_tensor[1:n_lim, 1, 1, k+1] = c * r_tensor[2:n_lim+1, 1, 1, k]
        if k != 1
            r_tensor[1:n_lim, 1, 1, k+1] += (k - 1) * r_tensor[2:n_lim+1, 1, 1, k-1]
        end
    end

    for j = 1:j_max-1
        n_lim = n_max - k_max - j
        r_tensor[1:n_lim, 1, j+1, :] = b * r_tensor[2:n_lim+1, 1, j, :]
        if j != 1
            r_tensor[1:n_lim, 1, j+1, :] += (j - 1) * r_tensor[2:n_lim+1, 1, j-1, :]
        end
    end

    for i = 1:i_max-1
        n_lim = n_max - k_max - j_max - i
        r_tensor[1:n_lim, i+1, :, :] = a * r_tensor[2:n_lim+1, i, :, :]
        if i != 1
            r_tensor[1:n_lim, i+1, :, :] += (i - 1) * r_tensor[2:n_lim+1, i-1, :, :]
        end
    end
    return r_tensor
end
