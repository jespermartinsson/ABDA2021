using PyPlot
using Random

include((@__DIR__)*"/../src/stats.jl")

# data
y = [0, 1, 0, 1, 0, 0, 0, 0, 0, 1]

# log-likelihood bernoulli
function log_likelihood(θ::Float64,y::Union{Vector{Float64},Vector{Int64}})
    if (0<θ<1)
        n = length(y)
        llh = 0.0
        for i in 1:n
            llh += y[i]*log(θ) + (1-y[i])*log(1-θ)
        end
        return llh
    else
        return -Inf
    end
    
end

# log prior
function log_beta(θ,α,β)
    if 0<θ<1
        return (α-1)*log(θ) + (β-1)*log(1-θ)
    else
        return -Inf
    end
end

# log posterior
log_posterior(θ) = log_likelihood(θ[1],y) + log_beta(θ[1],1,1)



## sample the posterior
θ_init = [0.5]
N_burn_in  = 500
θ_samp, lps = batch_sample(copy(θ_init), ones(length(θ_init)), log_posterior, 100_000, N_burn_in; printing=true)

# remove "burn in" phase
θ_samp = θ_samp[:,N_burn_in+1:end]
lps = lps[N_burn_in+1:end] 



# plot mcmc chain
figure()
subplot(1,2,1), plot(θ_samp[1,:])
subplot(1,2,2), hist(θ_samp[1,:],100)
tight_layout()