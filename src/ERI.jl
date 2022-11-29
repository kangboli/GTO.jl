export eri_cpqr

function eri_cpqr(gps::Vector{CHG}, ϵ::Float64=1e-3)
    norm(g) = g' * g
    norms = norm.(gps)
    # sort!(gps, by=norm)
    norm_perm = sortperm(norms, rev=true)
    for i = 1:length(norm_perm)
        ni = norm_perm[i]
        norm(gps[ni]) > ϵ || return gps, i
        for j = i+1:length(norm_perm)
            nj = norm_perm[j]
            gps[nj] = gps[nj] - (gps[ni]' * gps[nj]) * gps[ni]
            norms[nj] = norm(gps[nj]) 
            norms[nj] > ϵ || break
        end
        norm_perm = sortperm(norms, rev=true)
    end
end

