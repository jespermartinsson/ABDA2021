using PyPlot
using Statistics
using BenchmarkTools
include((@__DIR__) * "/../src/stats.jl")

run_mcmc = true


script_name = split(@__FILE__, "/")[end]
pathdir = prod([k * "/" for k in split(@__FILE__, "/")[1:end-1]])


const y = [607.0, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395, 346, 361, 310, 300, 382, 294, 315, 250, 320, 335, 297, 315, 316, 319, 326, 280, 330, 272, 317, 311, 336, 308, 330, 311, 304, 327, 367, 414, 411, 391, 333, 425, 416, 379, 368, 284, 260, 294, 265, 275, 270, 270, 282, 281, 283, 307, 344, 318, 326, 308, 306, 294, 342, 323, 325, 402, 219, 285, 277, 283, 271, 280, 289, 288, 199, 267, 354, 234, 270, 320, 214, 252, 234, 223, 332, 268, 286, 292, 295, 292, 275, 300, 256, 282, 319, 288, 316, 265, 306, 347, 374, 328, 305, 327, 344, 275, 218, 263, 282, 386, 307, 267, 282, 314, 328, 332, 386, 462, 368, 354, 283, 335, 264, 304, 248, 239, 288, 239, 249, 272, 289, 274, 281, 232, 279, 308, 260, 309, 312, 307, 296, 280, 331, 298, 342, 297, 297, 305, 308, 510, 490, 458, 425, 522, 927, 555, 550, 516, 548, 560, 545, 633, 496, 498, 555, 387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327]
const y_j = [[607.0, 583, 521, 494, 369], [782.0, 570, 678, 467, 620], [425.0, 395, 346, 361, 310, 300, 382, 294, 315], [250.0], [320.0, 335, 297, 315, 316, 319, 326, 280, 330, 272, 317, 311, 336, 308, 330], [311.0, 304, 327, 367, 414, 411, 391, 333, 425, 416, 379, 368], [284.0, 260, 294, 265, 275, 270, 270, 282, 281, 283], [307.0, 344, 318, 326, 308, 306, 294, 342, 323, 325, 402], [219.0, 285, 277, 283, 271, 280, 289, 288, 199, 267], [354.0, 234, 270, 320, 214, 252, 234, 223, 332, 268], [286.0, 292, 295, 292, 275, 300, 256, 282, 319, 288, 316, 265, 306], [347.0, 374, 328, 305, 327], [344.0, 275, 218, 263, 282], [386.0, 307, 267, 282, 314], [328.0, 332, 386, 462, 368], [354.0, 283, 335, 264, 304], [248.0, 239, 288, 239, 249], [272.0, 289, 274, 281, 232, 279, 308, 260, 309, 312, 307, 296, 280, 331, 298], [342.0, 297, 297, 305, 308], [510.0, 490, 458, 425, 522], [927.0, 555, 550, 516, 548], [560.0, 545, 633, 496, 498], [555.0, 387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327]]
const ind = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23]
logy = log.(y)

const I = length(y)
const J = ind[end]
const zlogy = (logy .- mean(logy)) ./ std(logy)

const child_j = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
const child_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

const x = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
const zx = (x .- mean(x)) ./ std(x)


#cd(dirname(PROGRAM_FILE))
#pathdir = "/media/jesper/ssd/LTU/research/olin/julia/"
using PyPlot
close("all")

function logprior_theta(theta::Array{Float64}, mu::Float64, sigma::Float64)
    #mu = b0 + b1*child_j
    return sum(-log(sigma) .- 0.5 * ((theta .- mu) ./ sigma) .^ 2)
end

function logprior_sigma(sigma::Float64)
    return log(sigma > 0)
end

function logprior_tau(tau::Float64)
    return log(tau > 0)
end

function logprior_mu(mu::Float64)
    return 0.0
end

function logll(beta::Array{Float64})
    eta0 = beta[1:J]

    b00 = beta[2*J+1]
    b01 = beta[2*J+2]

    eta1 = beta[J+1:2*J]
    b10 = beta[2*J+3]
    b11 = beta[2*J+4]

    sigma = beta[2*J+5]
    tau0 = beta[2*J+6]
    tau1 = beta[2*J+7]

    mu0 = b00 .+ b01 * child_j
    mu1 = b10 .+ b11 * child_j

    theta0 = mu0 .+ eta0 * tau0
    theta1 = mu1 .+ eta1 * tau1

    lp = 0.0
    for i = 1:I
        j = ind[i]
        mu = theta0[j] + theta1[j] * zx[i]
        lp += -log(sigma) - 0.5 * ((zlogy[i] - mu) / sigma)^2
    end
    return lp
end

function logpost(beta::Array{Float64})
    eta0 = beta[1:J]

    #b00 = beta[2*J+1]
    #b01 = beta[2*J+2]

    eta1 = beta[J+1:2*J]
    #b10 = beta[2*J+3]
    #b11 = beta[2*J+4]

    sigma = beta[2*J+5]
    tau0 = beta[2*J+6]
    tau1 = beta[2*J+7]

    #mu0 = b00 + b01*child_j
    #mu1 = b10 + b11*child_j

    #theta0 = mu0 + eta0*tau0
    #theta1 = mu1 + eta1*tau1

    if (sigma > 0) && (tau0 > 0) && (tau1 > 0)
        return logll(beta) + logprior_theta(eta0, 0.0, 1.0) + logprior_theta(eta1, 0.0, 1.0) + logprior_sigma(sigma) + logprior_tau(tau0) + logprior_tau(tau1)
    else
        return -Inf
    end

end


beta = ones(2 * J + 7)

N = 100_000
if run_mcmc
    @time zeta_samp, lp = batch_sample(beta, ones(length(beta)), logpost, N; printing = true)
end


eta0_samp = zeta_samp[1:J, :]

b00_samp = zeta_samp[2*J+1, :]
b01_samp = zeta_samp[2*J+2, :]

eta1_samp = zeta_samp[J+1:2*J, :]
b10_samp = zeta_samp[2*J+3, :]
b11_samp = zeta_samp[2*J+4, :]

sigma_samp = zeta_samp[2*J+5, :]
tau0_samp = zeta_samp[2*J+6, :]
tau1_samp = zeta_samp[2*J+7, :]

mu0_samp = repeat(b00_samp, 1, J)' + repeat(b01_samp, 1, J)' .* repeat(child_j, 1, size(b01_samp, 1))
mu1_samp = repeat(b10_samp, 1, J)' + repeat(b11_samp, 1, J)' .* repeat(child_j, 1, size(b11_samp, 1))

theta0_samp = mu0_samp + eta0_samp .* repeat(tau0_samp, 1, J)'
theta1_samp = mu1_samp + eta1_samp .* repeat(tau1_samp, 1, J)'





# (ln(y)-my)/sy = theta0 + theat1*zx + s*e
# ln(y) = sy*(theta0 + theat1*zx + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*(theta0 + theat1*zx) + my, sigma = sy*s
# y = exp(sy*(theta0 + theat1*zx) + sy*s*e + mu)
# y = exp( (sy*theta0 + my) + theat1*zx + sy*s*e)
# y = exp( (sy*theta0 + my) + theat1*zx) * exp(sy*s*e)
# y = exp(theta02 + theta12*zx) * exp(sigma2*e)
# y = exp(theta02 + theta12*(x-mx)/sx) * exp(sigma2*e)
# y = exp(theta02 + theta12*x/sx - theta12*mx/sx) * exp(sigma2*e)
# y = exp(theta02  - theta12*mx/sx + theta12*x/sx) * exp(sigma2*e)
# y = exp(theta03  + theta13*x) * exp(sigma2*e)

sigma2_samp = std(logy) .* sigma_samp
theta02_samp = std(logy) .* theta0_samp .+ mean(logy)
theta12_samp = std(logy) .* theta1_samp

theta03_samp = theta02_samp .- theta12_samp .* mean(x) ./ std(x)
theta13_samp = theta12_samp ./ std(x)


figure(figsize = (8 * 2, 6 * 2), dpi = 80)
N = Int(ceil(sqrt(J)))
for i = 1:J
    #subplot(N,N,i)
    subplot(9, 4, i)
    x = theta0_samp[i, :]
    myhist(x, color = "r")
    xlabel(string("\$\\theta_0[", i, "]\$"))
    yl = ylim()
    ylim([yl[1], yl[2] * 1.3])
end
tight_layout()


figure(figsize = (8 * 2, 6 * 2), dpi = 80)
N = Int(ceil(sqrt(J)))
for i = 1:J
    #subplot(N,N,i)
    subplot(9, 4, i)
    x = theta1_samp[i, :]
    myhist(x, color = "b")
    xlabel(string("\$\\theta_1[", i, "]\$"))
    yl = ylim()
    ylim([yl[1], yl[2] * 1.3])
end
plt["tight_layout"]()







figure(figsize = (8 * 2, 6 * 2), dpi = 80)
N = Int(ceil(sqrt(J)))
color_list = ["b", "r", "b", "g", "m"]
col = "k"
cntr = 0
for k in [1, 5]
    global cntr
    mu_y_samp = exp.(theta03_samp .+ k * theta13_samp .+ repeat(sigma2_samp', size(theta13_samp, 1), 1) .^ 2 ./ 2)
    for i = 1:J
        #subplot(N,N,i)
        subplot(9, 4, i)
        x = mu_y_samp[i, :]
        col = color_list[mod(cntr, length(color_list))+1]
        myhist(x, color = col)
        plot(-10, 0, color = col, label = string("\$x=", k, "\$"))
        xlabel(string("\$\\exp(\\theta_0[", i, "]+\\theta_1[{", i, "}]x + \\sigma^2/2)\$"))
        xlim([100, 800])
        yl = ylim()
        ylim([yl[1], yl[2] * 1.3])
    end
    cntr += 1
end
for i = 1:J
    #subplot(N,N,i)
    subplot(9, 4, i)
    legend(numpoints = 1, fancybox = true)
end
plt["tight_layout"]()




figure(figsize = (8 * 2, 6 * 2), dpi = 80)
N = Int(ceil(sqrt(J)))
color_list = ["k", "r"]
for j = 1:J
    subplot(4, 9, j)
    k = 0:0.1:21
    col = color_list[Int(mod(child_j[j], 2) + 1)]
    for i = 1:10:2000
        mu_y_samp = exp.(theta03_samp[j, i] .+ k .* theta13_samp[j, i] .+ sigma2_samp[i] .^ 2 ./ 2)
        plot(k, mu_y_samp, "b-", alpha = 0.1)
    end
    plot(1:length(y_j[j]), y_j[j], col * "o", label = string("j = ", j))
    ylim([100, 800])
    xlim([0, 21])
    if j in [1, 10, 19, 28]
        ylabel("reaction time")
    end
end
for j = 1:J
    #subplot(N,N,i)
    subplot(4, 9, j)
    legend(numpoints = 1, fancybox = true)
end
plt["tight_layout"]()




figure(figsize = (8 * 2, 6 * 2), dpi = 80)
N = Int(ceil(sqrt(J)))
color_list = ["b", "r", "b", "g", "m"]
col = "k"
cntr = 0
for k in [1, 5]
    global cntr, col, color_list
    mu_y_samp = exp.(theta03_samp .+ k * theta13_samp .+ repeat(sigma2_samp', size(theta13_samp, 1), 1) .^ 2 ./ 2)
    jcntr = 1
    for i in [1, 3, 4]
        #subplot(N,N,i)
        subplot(2, 3, jcntr)
        x = mu_y_samp[i, :]
        col = color_list[mod(cntr, length(color_list))+1]
        myhist(x, color = col)
        plot(-10, 0, color = col, label = string("\$x=", k, "\$"))
        xlabel(string("\$\\exp(\\theta_0[", i, "]+\\theta_1[{", i, "}]x + \\sigma^2/2)\$"))
        xlim([100, 800])
        yl = ylim()
        ylim([yl[1], yl[2] * 1.3])
        jcntr += 1
    end
    cntr += 1
end
jcntr = 1
for i in [1, 3, 4]
    global jcntr
    #subplot(N,N,i)
    subplot(2, 3, jcntr)
    legend(numpoints = 1, fancybox = true)
    jcntr += 1
end
plt["tight_layout"]()

#figure(figsize=(8*2, 6*2), dpi=80)
N = Int(ceil(sqrt(J)))
color_list = ["k", "r"]
jcntr = 1
for j in [1, 3, 4]
    global jcntr
    subplot(2, 3, jcntr + 3)
    k = 0:0.1:21
    col = color_list[Int(mod(child_j[j], 2) + 1)]
    for i = 1:10:2000
        mu_y_samp = exp.(theta03_samp[j, i] .+ k .* theta13_samp[j, i] .+ sigma2_samp[i] .^ 2 ./ 2)
        plot(k, mu_y_samp, "b-", alpha = 0.1)
    end
    plot(1:length(y_j[j]), y_j[j], col * "o", label = string("j = ", j))
    plot([1, 1], [100, 800], "b")
    plot([5, 5], [100, 800], "r")
    ylim([100, 800])
    xlim([0, 21])
    if j in [1]
        ylabel("reaction time")
    end
    xlabel("attempt nr")
    jcntr += 1
end
jcntr = 1
for i in [1, 3, 4]
    global jcntr
    #subplot(N,N,i)
    subplot(2, 3, jcntr)
    legend(numpoints = 1, fancybox = true)
    jcntr += 1
end
plt["tight_layout"]()



figure()
N = Int(ceil(sqrt(J)))
color_list = ["b", "r", "b", "g", "m"]
col = "k"
cntr = 0
mu_y_samp = exp.(theta13_samp)
jcntr = 1
for i in [1, 3, 4]
    global cntr, jcntr
    #subplot(N,N,i)
    subplot(1, 3, jcntr)
    x = mu_y_samp[i, :]
    col = color_list[mod(cntr, length(color_list))+1]
    myhist(x, color = col)
    xlabel(string("\$\\exp(\\theta_1[{", i, "}])\$"))
    yl = ylim()
    ylim([yl[1], yl[2] * 1.3])
    jcntr += 1
end
cntr += 1

jcntr = 1
for i in [1, 3, 4]
    global jcntr
    #subplot(N,N,i)
    subplot(1, 3, jcntr)
    legend(numpoints = 1, fancybox = true)
    jcntr += 1
end
plt["tight_layout"]()


figure()
subplot(121)
myhist(tau0_samp)
xlabel(string("\$\\tilde{\\tau}_0\$"))
plt["tight_layout"]()

# (ln(y)-my)/sy = (theta0 + phi0*c) + (theta1 + phi1*c)*zx + s*e
# ln(y) = sy*((theta0 + phi0*c) + (theta1 + phi1*c)*zx + s*e) + my
# ln(y) = beta + sigma*e, beta = sy*theta + sy*phi*c + my, sigma = sy*s
# y = exp(sy*theta0 + sy*phi0*c + my + sy*(theta1 + phi1*c)*zx + sy*s*e)
# y = exp(sy*theta0 + sy*phi0*c + my + sy*(theta1 + phi1*c)*(x-mx)/sx + sy*s*e)
# y = exp(sy*theta0 + sy*phi0*c + my + sy/sx*(theta1 + phi1*c)*x - sy/sx*mx*(theta1 + phi1*c) + sy*s*e)
# y = exp(sy*theta0 + sy*phi0*c + my + sy/sx*(theta1 + phi1*c)*x - sy/sx*mx*(theta1 + phi1*c))*exp(sy*s*e)
# y = exp(sy*theta0 + my - sy/sx*mx*theta1 + (sy*phi0 - sy/sx*mx*phi1)*c + (sy/sx*theta1 + sy/sx*phi1*c)*x )*exp(sy*s*e)
# y = exp( (theta02+phi02*c) + (theta12 +phi12*c)*x)*exp(sigma2*e)
# theta02 = sy*theta0  + my - sy/sx*mx*theta1 + (sy*phi0 - sy/sx*mx*phi1)*c
# theta12 = (sy/sx*theta1 + sy/sx*phi1*c)
# phi02 = sy*phi0 - sy/sx*mx*phi1
# phi12 = sy/sx*phi1

# theta02 = sy*theta0  + my - sy/sx*mx*theta1 + (sy*phi0 - sy/sx*mx*phi1)*c
# theta12 = (sy/sx*theta1 + sy/sx*phi1*c)
# phi02 = sy*phi0 - sy/sx*mx*phi1
# phi12 = sy/sx*phi1


mu0_samp = repeat(b00_samp, 1, J)' .+ repeat(b01_samp, 1, J)' .* repeat(child_j, 1, size(b01_samp, 1))
mu1_samp = repeat(b10_samp, 1, J)' .+ repeat(b11_samp, 1, J)' .* repeat(child_j, 1, size(b11_samp, 1))


phi0_samp = b01_samp
phi1_samp = b11_samp

phi02_samp = std(logy) .* phi0_samp .- std(logy) / std(x) * mean(x) .* phi1_samp
phi12_samp = std(logy) / std(x) .* phi1_samp




subplot(122)
myhist(std(logy) * tau0_samp)
xlabel(string("\$\\tau_0\$"))
plt["tight_layout"]()




figure()
plot(theta03_samp[1, :], theta13_samp[1, :], ".")
xlabel(string("\$\\theta_0\$"))
ylabel(string("\$\\theta_1\$"))
plt["tight_layout"]()

figure()
plot(theta02_samp[1, :], theta12_samp[1, :], ".")
xlabel(string("\$\\theta_0\$"))
ylabel(string("\$\\theta_1\$"))
plt["tight_layout"]()

figure()
plot(theta03_samp[1, :], theta13_samp[1, :], ".")
xlabel(string("\$\\theta_0\$"))
ylabel(string("\$\\theta_1\$"))
plt["tight_layout"]()


figure()
subplot(121)
myhist(sigma_samp)
xlabel(string("\$\\tilde{\\sigma}\$"))
plt["tight_layout"]()
subplot(122)
myhist(std(logy) * sigma_samp)
xlabel(string("\$\\sigma\$"))
plt["tight_layout"]()




figure()
subplot(121)
myhist(phi02_samp)
xlabel(string("\$\\varphi_0 = std(logy)\\tilde{\\varphi}_0 - \\tilde{\\varphi}_1 std(logy)mean(x)/std(x)\$"))
plt["tight_layout"]()
subplot(122)
myhist(phi12_samp)
xlabel(string("\$\\varphi_1 = std(logy)/std(x)\\tilde{\\varphi}_1\$"))
plt["tight_layout"]()

figure()
subplot(121)
myhist(exp.(phi02_samp))
xlabel(string("\$\\exp(\\varphi_0)\$"))
plt["tight_layout"]()
subplot(122)
myhist(exp.(phi12_samp))
xlabel(string("\$\\exp(\\varphi_1)\$"))
plt["tight_layout"]()


figure()
subplot(121)
myhist(phi0_samp)
xlabel(string("\$\\tilde{\\varphi}_0\$"))
plt["tight_layout"]()
subplot(122)
myhist(phi1_samp)
xlabel(string("\$\\tilde{\\varphi}_1\$"))
plt["tight_layout"]()



# mu = exp((theta0+phi0*c) + (theta1+phi1*c)*x + 0.5s^2) 
# mu = exp(phi0*c) * exp(phi1*c*x) * exp((theta0 + theta1*x + 0.5s^2) 

figure(),
subplot(211),
plot(Array(eta0_samp'))
ylabel(raw"$\eta_0$")
xlabel(raw"sample")
xlim([5000, 5200])
subplot(212),
plot(Array(eta1_samp'))
ylabel(raw"$\eta_1$")
xlabel(raw"sample")
xlim([5000, 5200])
tight_layout()

println("ess(eta0_samp): ", ess(eta0_samp))
println("ess(eta1_samp): ", ess(eta1_samp))