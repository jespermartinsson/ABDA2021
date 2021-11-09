using PyPlot
using Statistics
using BenchmarkTools
include((@__DIR__)*"/../src/stats.jl")


close("all")

run_mcmc = true


script_name = split(@__FILE__, "/")[end]
pathdir = prod([k * "/" for k in split(@__FILE__, "/")[1:end - 1]])

y = [607.0, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395, 346, 361, 310, 300, 382, 294, 315, 250, 320, 335, 297, 315, 316, 319, 326, 280, 330, 272, 317, 311, 336, 308, 330, 311, 304, 327, 367, 414, 411, 391, 333, 425, 416, 379, 368, 284, 260, 294, 265, 275, 270, 270, 282, 281, 283, 307, 344, 318, 326, 308, 306, 294, 342, 323, 325, 402, 219, 285, 277, 283, 271, 280, 289, 288, 199, 267, 354, 234, 270, 320, 214, 252, 234, 223, 332, 268, 286, 292, 295, 292, 275, 300, 256, 282, 319, 288, 316, 265, 306, 347, 374, 328, 305, 327, 344, 275, 218, 263, 282, 386, 307, 267, 282, 314, 328, 332, 386, 462, 368, 354, 283, 335, 264, 304, 248, 239, 288, 239, 249, 272, 289, 274, 281, 232, 279, 308, 260, 309, 312, 307, 296, 280, 331, 298, 342, 297, 297, 305, 308, 510, 490, 458, 425, 522, 927, 555, 550, 516, 548, 560, 545, 633, 496, 498, 555, 387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327]
const ind = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23]
logy = log.(y)

const I = length(y)
const J = ind[end]
const zlogy = (logy .- mean(logy)) ./ std(logy)



function logprior_theta(theta::Array{Float64}, mu::Float64, tau::Float64)
    return sum(-log(tau) .- 0.5 * ((theta .- mu) ./ tau).^2)
end

function logprior_sigma(sigma::Float64)
    return log(sigma > 0)
end

function logprior_tau(tau::Float64)
    return log(tau > 0)
end

function logprior_mu(mu::Float64)
    return 0
end
function logll(beta::Array{Float64})
    theta = beta[1:J]
    mu = beta[J + 1]
    sigma = beta[J + 2]
    tau = beta[J + 3]

    lp = 0.0
    for i in 1:I
        j = ind[i] # here I obtain the individual j corresponding to the i:th measurement
        lp += -log(sigma) - 0.5 * ((zlogy[i] - theta[j]) / sigma)^2
    end
    return lp
end

function logpost(beta::Array{Float64})
    theta = beta[1:J]
    mu = beta[J + 1]
    sigma = beta[J + 2]
    tau = beta[J + 3]
    if (sigma > 0) && (tau > 0)
        return logll(beta) + logprior_theta(theta, mu, tau) + logprior_mu(mu) + logprior_sigma(sigma) + logprior_tau(tau)
    else
        return -Inf
    end
end



beta = ones(J + 3)
zeta_samp = 0
if run_mcmc
    #@btime logpost(beta)
    #zeta_samp, lp = slice_sample(beta, ones(length(beta)), logpost, 12_500; printing=true)
    @time zeta_samp, lp = batch_sample(beta, ones(length(beta)), logpost, 10_0000; printing=true)

end



figure()
N = Int(ceil(sqrt(length(beta))))
for i in 1:length(beta)
    subplot(N, N, i)
    myhist(zeta_samp[i,:], 100)
    xlabel(string("\$\\theta_{", i, "}\$"))
end
tight_layout()

theta_samp = zeta_samp[1:J,:]
mu_samp = zeta_samp[J + 1,:]
sigma_samp = zeta_samp[J + 2,:]
tau_samp = zeta_samp[J + 3,:]

# (ln(y)-my)/sy = theta + s*e
# ln(y) = sy*(theta + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*theta+my, sigma = sy*s
# y = exp(sy*theta + my + sy*s*e)
# y = exp(sy*theta+my)*exp(sy*s*e)
# y = exp(theta2)*exp(sigma2*e)

sigma2_samp = std(logy) .* sigma_samp
theta2_samp = std(logy) .* theta_samp .+ mean(logy)

# mean value of log-normal pdf
# https://en.wikipedia.org/wiki/Log-normal_distribution
mu_y_samp = exp.(theta2_samp .+ repeat(sigma2_samp', size(theta2_samp, 1), 1).^2 ./ 2)



figure(figsize = (8 * 2, 6 * 2), dpi = 80)
N = Int(ceil(sqrt(J)))
col = "b"
for i in 1:J
    subplot(N, N, i)
    x = mu_y_samp[i,:]
    myhist(x, color = col)
    xlabel(string("\$\\exp(\\theta_{", i, "}+\\sigma^2/2)\$"))
end
tight_layout()

figure()
col = "b"
i = 4
x = mu_y_samp[i,:]
myhist(x, color = col)
xlabel(string("\$\\exp(\\theta_{", i, "}+\\sigma^2/2)\$"))
title("The expected reaction time for the dude")
tight_layout()



mu2_samp = std(logy) .* mu_samp .+ mean(logy)
tau2_samp = std(logy) .* tau_samp


# Derivation for the expected reaction time for the group
# For a more analytical approach you may also let:
# theta = mu + tau*xi ~ N(mu,tau)
# x = theta + sigma*eta ~ N(theta,sigma)
# where eta~N(0,1) and xi~N(0,1) assuming independence. Then
# x = mu + sigma*eta + tau*xi ~ N(mu, sqrt(sigma^2 + tau^2))
# and sqrt(sigma^2 + tau^2) is the std.
# If y = exp(x) then E(y) = exp(mu + (sigma^2 + tau^2)/2), see also simulation here:
# https://github.com/jespermartinsson/ABDA.jl/blob/master/dev/test_var.jl

exp_mu2_samp = exp.(mu2_samp .+ 0.5 .* tau2_samp.^2 .+ 0.5 .* sigma2_samp.^2)



figure()
myhist(exp_mu2_samp, color = "b")
xlabel(string("\$\\exp(\\mu + \\tau^2/2 + \\sigma^2/2)\$"))
title("expected reaction time for the group")
tight_layout()

my = mean(logy)
sy = std(logy)
y_pred = zeros(size(zeta_samp, 2))
for i in 1:size(zeta_samp, 2)
    theta_i = mu_samp[i] + tau_samp[i] * randn()
    logzy_i = theta_i + sigma_samp[i] * randn()
    logy_i = sy * logzy_i + my
    y_pred[i] = exp(logy_i)
end

figure()
myhist(y_pred)
title("posterior predicted reaction time a random individual")

figure()
for j in 1:J
    plot(log.(y[ind .== j]))
end

mle = []
for j in 1:J
    push!(mle, mean(logy[ind .== j]))
end

function logprior_thetas(theta, mu, tau) # p(theta[1],...,theta[J]) ~ N(mu,tau)
    return (-log.(tau) .- 0.5 .* ((theta .- mu) ./ tau).^2)
end

offset = 2000
figure()
for j in 1:J
#    fr,bins = np.myhistogram(theta2_samp[j,:],Int(round(sqrt(length(theta2_samp[j,:])))),normed=true)
#    x,y = get_mystep(bins,fr)
#    plot(x,y*0.3+j)
    myhist(theta2_samp[j,:], baseline = j * offset, color = "b", scale = 100.0)
    plot(mle[j], j * offset, "ro")
end
plot(mle[1], 1 * offset, "ro", label = "MLE estimate \$\\hat{\\theta}_{\\mathrm{MLE}}\$")
xl = xlim()
thetas = range(xl[1], stop=xl[2], length=1000)
plot(thetas, 10*offset * exp.(logprior_thetas(thetas, mean(mu2_samp), mean(tau2_samp))), label = "\$p(\\hat{\\mu},\\hat{\\tau)}\$", color = "m")
legend()
yticks([])
xlabel(string("\$\\theta_j\$ expected log(reaction time)"))
tight_layout()



figure()
subplot(121)
myhist(tau_samp)
xlabel(string("\$\\tilde{\\tau}\$"))
tight_layout()

# (ln(y)-my)/sy = theta + phi*c + s*e
# ln(y) = sy*(theta + phi*c + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*theta + sy*phi*c + my, sigma = sy*s
# y = exp(sy*theta + sy*phi*c + my + sy*s*e)
# y = exp(sy*theta + sy*phi*c +my)*exp(sy*s*e)
# y = exp(theta2+phi2*c)*exp(sigma2*e)
# theta2 = sy*theta + my
# phi2 = sy*phi


subplot(122)
myhist(std(logy) * tau_samp)
xlabel(string("\$\\tau\$"))
tight_layout()



if false
    my = mean(logy)
    sy = std(logy)
    mean_reaction_time = zeros(size(zeta_samp, 2))
    for i in 1:size(zeta_samp, 2)
        theta_i = mu_samp[i] + tau_samp[i] * randn()
        theta2_i = sy * theta_i + my
        mean_reaction_time[i] = exp(theta2_i + 0.5 * sigma2_samp[i]^2)
    end

    figure()
    myhist(mean_reaction_time)
    title("posterior predicted reaction time for the group")

    figure()
    plot(tau2_samp, sigma2_samp, ".")
end


figure()
subplot(121)
myhist(sigma_samp)
xlabel(string("\$\\tilde{\\sigma}\$"))
tight_layout()
subplot(122)
myhist(std(logy) * sigma_samp)
xlabel(string("\$\\sigma\$"))
tight_layout()