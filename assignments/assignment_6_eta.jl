using PyPlot
using Statistics
using BenchmarkTools
include((@__DIR__) * "/../src/stats.jl")

run_mcmc = true


script_name = split(@__FILE__, "/")[end]
pathdir = prod([k * "/" for k in split(@__FILE__, "/")[1:end-1]])

y = [607.0, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395, 346, 361, 310, 300, 382, 294, 315, 250, 320, 335, 297, 315, 316, 319, 326, 280, 330, 272, 317, 311, 336, 308, 330, 311, 304, 327, 367, 414, 411, 391, 333, 425, 416, 379, 368, 284, 260, 294, 265, 275, 270, 270, 282, 281, 283, 307, 344, 318, 326, 308, 306, 294, 342, 323, 325, 402, 219, 285, 277, 283, 271, 280, 289, 288, 199, 267, 354, 234, 270, 320, 214, 252, 234, 223, 332, 268, 286, 292, 295, 292, 275, 300, 256, 282, 319, 288, 316, 265, 306, 347, 374, 328, 305, 327, 344, 275, 218, 263, 282, 386, 307, 267, 282, 314, 328, 332, 386, 462, 368, 354, 283, 335, 264, 304, 248, 239, 288, 239, 249, 272, 289, 274, 281, 232, 279, 308, 260, 309, 312, 307, 296, 280, 331, 298, 342, 297, 297, 305, 308, 510, 490, 458, 425, 522, 927, 555, 550, 516, 548, 560, 545, 633, 496, 498, 555, 387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327]
const ind = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23]
logy = log.(y)

const I = length(y)
const J = ind[end]
const zlogy = (logy .- mean(logy)) ./ std(logy)

const child_j = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
const child_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


script_name = split(@__FILE__, "/")[end]
pathdir = prod([k * "/" for k in split(@__FILE__, "/")[1:end-1]])






function logprior_theta(eta::Array{Float64}, mu::Float64, tau::Float64)
    return sum(-log.(tau) .- 0.5 * ((eta .- mu) ./ tau) .^ 2)
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
    eta = beta[1:J]
    b0 = beta[J+1]
    b1 = beta[J+2]
    sigma = beta[J+3]
    tau = beta[J+4]

    mu = b0 .+ b1 * child_j
    theta = mu .+ tau * eta

    lp = 0.0
    for i = 1:I
        j = ind[i]
        lp += -log(sigma) - 0.5 * ((zlogy[i] - theta[j]) / sigma)^2
    end
    return lp
end

function logpost(beta::Array{Float64})
    eta = beta[1:J]
    sigma = beta[J+3]
    tau = beta[J+4]

    if (sigma > 0) && (tau > 0)
        return logll(beta) + logprior_theta(eta, 0.0, 1.0) + logprior_sigma(sigma) + logprior_tau(tau)
    else
        return -Inf
    end

end



beta = ones(J + 4)
#beta[J+2]=0.1
#@enter logpost(beta)
#error()
N = 100_000
if run_mcmc
    @time zeta_samp, lp = batch_sample(beta, ones(length(beta)), logpost, N; printing = true)
end



# j = 1
# y_j = []
# y_i = []
# for i in 1:I
#     if ind[i] == j
#         push!(y_i, y[i])
#     else
#         push!(y_j, y_i)
#         y_i = []
#         push!(y_i, y[i])
#         j += 1
#     end
# end
# push!(y_j, y_i)

# figure()
# col = ""
# for j in 1:J
#     if child_j[j] == 1
#         col = "r"
#     else
#         col = "k"
#     end
#     subplot(121)
#     plot(1:length(y_j[j]), y_j[j], col * "-")
#     ylabel("reaction time")
#     xlabel("attempt nr")
#     xlim([1,10])
#     subplot(122)
#     plot(1:length(y_j[j]), log.(y_j[j]), col * "-")
#     ylabel("log(reaction time)")
#     xlabel("attempt nr")
#     xlim([1,10])
# end
# tight_layout()





eta_samp = zeta_samp[1:J, :]
b0_samp = zeta_samp[J+1, :]
b1_samp = zeta_samp[J+2, :]
sigma_samp = zeta_samp[J+3, :]
tau_samp = zeta_samp[J+4, :]
mu_samp = repeat(b0_samp, 1, J)' .+ repeat(b1_samp, 1, J)' .* repeat(child_j, 1, size(b1_samp, 1))
theta_samp = mu_samp .+ eta_samp .* repeat(tau_samp, 1, J)'


# (ln(y)-my)/sy = theta + phi*c + s*e
# ln(y) = sy*(theta + phi*c + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*theta + sy*phi*c + my, sigma = sy*s
# y = exp(sy*theta + sy*phi*c + my + sy*s*e)
# y = exp(sy*theta + sy*phi*c +my)*exp(sy*s*e)
# y = exp(theta2+phi2*c)*exp(sigma2*e)
# theta2 = sy*theta + my
# phi2 = sy*phi


phi_samp = b1_samp
phi2_samp = std(logy) * b1_samp

figure()
subplot(121)
myhist(phi_samp)
xlabel(raw"$\tilde{\varphi}$")
tight_layout()
subplot(122)
myhist(phi2_samp)
xlabel(raw"$\varphi$")
tight_layout()





# (ln(y)-my)/sy = theta + s*e
# ln(y) = sy*(theta + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*theta+my, sigma = sy*s
# y = exp(sy*theta + my + sy*s*e)
# y = exp(sy*theta+my)*exp(sy*s*e)
# y = exp(theta2)*exp(sigma2*e)

sigma2_samp = std(logy) * sigma_samp
theta2_samp = std(logy) * theta_samp .+ mean(logy)

mu_y_samp = exp.(theta2_samp + repeat(sigma2_samp', size(theta2_samp, 1), 1) .^ 2 / 2)


figure()
Np = Int(ceil(sqrt(J)))
for i = 1:J
    subplot(Np, Np, i)
    myhist(mu_y_samp[i, :])
    xlabel(string("\$\\theta_{", i, "}\$"))
end
tight_layout()

mu2_samp = std(logy) * mu_samp .+ mean(logy)
tau2_samp = std(logy) * tau_samp

exp_mu2_samp = exp.(mu2_samp .+ repeat(tau2_samp, 1, J)' .^ 2 / 2 .+ repeat(sigma2_samp, 1, J)' .^ 2 / 2)

figure()
Np = Int(ceil(sqrt(J)))
i = 0
for i = 1:J
    subplot(Np, Np, i)
    myhist(exp_mu2_samp[i, :])
    xlabel(string("\$\\theta_{", i, "}\$"))
end
tight_layout()


#eta = beta[1:J]
#b0 = beta[J+1]
#b1 = beta[J+2]
#sigma = beta[J+3]
#tau = beta[J+4]
#mu = b0 + b1*child_j
#theta = mu + tau*eta

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


using Distributions
a, b = 1, 1
if true
    z = sum(child_j)
else
    z = Int(round(J / 2)) # test half kids half adults
end
thetas = rand(Beta(z + a, J - z + b), size(zeta_samp, 2))
child_rand = thetas .> rand(size(zeta_samp, 2))



my = mean(logy)
sy = std(logy)

y_pred = zeros(size(zeta_samp, 2), 3)
for i = 1:size(zeta_samp, 2)
    child_i = 0.0
    mu_i = b0_samp[i] + b1_samp[i] * child_i
    theta_i = mu_i + tau_samp[i] * randn()
    logzy_i = theta_i + sigma_samp[i] * randn()
    logy_i = sy * logzy_i + my
    y_pred[i, 1] = exp(logy_i)

    child_i = 1.0
    mu_i = b0_samp[i] + b1_samp[i] * child_i
    theta_i = mu_i + tau_samp[i] * randn()
    logzy_i = theta_i + sigma_samp[i] * randn()
    logy_i = sy * logzy_i + my
    y_pred[i, 2] = exp(logy_i)

    child_i = child_rand[i]
    mu_i = b0_samp[i] + b1_samp[i] * child_i
    theta_i = mu_i + tau_samp[i] * randn()
    logzy_i = theta_i + sigma_samp[i] * randn()
    logy_i = sy * logzy_i + my
    y_pred[i, 3] = exp(logy_i)
end

figure()
myhist(y_pred[:, 1], color = "b")
myhist(y_pred[:, 2], color = "r")
myhist(y_pred[:, 3], color = "k")
title("posterior predicted reaction time a random individual")






mle = []
for j = 1:J
    push!(mle, mean(logy[ind.==j]))
end

function logprior_thetas(theta, mu, tau) # p(theta[1],...,theta[J]) ~ N(mu,tau)
    return (-log(tau) - 0.5 * ((theta - mu) / tau) .^ 2)
end


figure()
for j = 1:J
    #    fr,bins = np.histogram(theta2_samp[j,:],Int(round(sqrt(length(theta2_samp[j,:])))),normed=true)
    #    x,y = get_mystep(bins,fr)
    #    plot(x,y*0.3+j)
    if child_j[j] == 0
        myhist(theta2_samp[j, :], baseline = j * 10, color = "k")
        plot(mle[j], j * 10, "ko")
    else
        myhist(theta2_samp[j, :], baseline = j * 10, color = "r")
        plot(mle[j], j * 10, "ro")
    end
end
plot(mle[1], 1 * 10, "ro", label = "MLE estimate \$\\hat{\\theta}_{\\mathrm{MLE}}\$")
xl = plt["xlim"]()
thetas = range(xl[1], xl[2], 1000)
plot(thetas, 80 * exp.(logprior_thetas.(thetas, mean(mu2_samp, dims = 2)[3], mean(tau2_samp))), label = "\$p(\\hat{\\mu},\\hat{\\tau)}\$", color = "k")
plot(thetas, 80 * exp.(logprior_thetas.(thetas, mean(mu2_samp, dims = 2)[1], mean(tau2_samp))), label = "\$p(\\hat{\\mu}+\\hat{\\varphi},\\hat{\\tau)}\$", color = "r")
plt["legend"]()
plt["yticks"]([])
xlabel(string("\$\\theta_j\$ expected log(reaction time)"))
plt["tight_layout"]()







figure()
thetas = range(xl[1], xl[2], 1000)
plot(thetas, exp.(logprior_thetas.(thetas, mean(mu2_samp, dims = 2)[3], mean(tau2_samp))), label = "\$p(\\hat{\\mu},\\hat{\\tau)}\$", color = "k")
plot(thetas, exp.(logprior_thetas.(thetas, mean(mu2_samp, dims = 2)[1], mean(tau2_samp))), label = "\$p(\\hat{\\mu}+\\hat{\\varphi},\\hat{\\tau)}\$", color = "r")
mean_pdf_kids = zeros(1000)
mean_pdf_adults = zeros(1000)
for i = 1:1000
    mean_pdf_adults .+= exp.(logprior_thetas.(thetas, mu2_samp[3, i], tau2_samp[i]))
    mean_pdf_kids .+= exp.(logprior_thetas.(thetas, mu2_samp[1, i], tau2_samp[i]))
    plot(thetas, exp.(logprior_thetas.(thetas, mu2_samp[3, i], tau2_samp[i])), color = "k", alpha = 0.01)
    plot(thetas, exp.(logprior_thetas.(thetas, mu2_samp[1, i], tau2_samp[i])), color = "r", alpha = 0.01)
end
mean_pdf_adults .= mean_pdf_adults / 1000
mean_pdf_kids .= mean_pdf_kids / 1000
plot(thetas, mean_pdf_adults, "--", label = "\$\\frac{1}{N}\\sum_{i=1}^N p(\\mu_i,\\tau_i)\$", color = "k")
plot(thetas, mean_pdf_kids, "--", label = "\$\\frac{1}{N}\\sum_{i=1}^Np(\\mu_i+\\varphi_i,\\tau_i)\$", color = "r")



plt["legend"]()
plt["yticks"]([])
xlabel(string("\$\\theta_j\$ expected log(reaction time)"))
plt["tight_layout"]()




# (ln(y)-my)/sy = theta + phi*c + s*e
# ln(y) = sy*(theta + phi*c + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*theta + sy*phi*c + my, sigma = sy*s
# y = exp(sy*theta + sy*phi*c + my + sy*s*e)
# y = exp(sy*theta + sy*phi*c +my)*exp(sy*s*e)
# y = exp(theta2+phi2*c)*exp(sigma2*e)
# theta2 = sy*theta + my
# phi2 = sy*phi

figure()
subplot(121)
myhist(tau_samp)
xlabel(raw"$\tilde{\tau}$")
tight_layout()
subplot(122)
myhist(tau2_samp)
xlabel(raw"$\tau$")
tight_layout()








figure()
subplot(121)
myhist(sigma_samp)
xlabel(string("\$\\tilde{\\sigma}\$"))
plt["tight_layout"]()
subplot(122)
myhist(std(logy) * sigma_samp)
xlabel(string("\$\\sigma\$"))
plt["tight_layout"]()






## Fake data
a, b = 1, 1
if false
    z = sum(child_j)
else
    z = Int(round(J * 0.5)) # test half kids half adults
end
thetas = rand(Beta(z + a, J - z + b), N)
child_rand = thetas .> rand(N)

y_pred_fake = zeros(N, 3)
for i = 1:N
    k = 1
    for child_i in [0, 1, child_rand[i]]
        theta_i = (mu_samp[i] + phi_samp[i] * child_i) + tau_samp[i] * randn()
        logzy_i = theta_i + sigma_samp[i] * randn()
        logy_i = sy * logzy_i + my
        y_pred_fake[i, k] = exp(logy_i)
        k += 1
    end
end

figure()
myhist(y_pred_fake[:, 1], color = "k", label = "adults")
myhist(y_pred_fake[:, 2], color = "r", label = "kids")
myhist(y_pred_fake[:, 3], color = "b", label = "mixed")
title("posterior predicted reaction time a random individual (50 % chidren)")
legend()






figure()
pval = sum(exp.(phi2_samp) .> 1) / N
myhist(exp.(phi2_samp), label = (raw"$\mathrm{Pr}\{\mathrm{exp}(\varphi)>1\} = $" * string(pval)))
xlabel(raw"$\mathrm{exp}(\varphi)$")
title("kid's mutiplicative effect on average reaction time")
legend()
tight_layout()
