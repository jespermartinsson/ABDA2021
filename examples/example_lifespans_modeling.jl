

age = [18, 40, 60, 80, 95]
lifespans = [[82], [82, 84, 81], [83, 84], [90], [98, 100, 101]]

include("plot_lifespans.jl")


# let a = age, and l = lifespan
# p(l|a) = p(a|l)p(l)/p(a) 
# we can not be older than our lifespan
#   p(a|l) = 1 if a<=l, 0 otherwise 
#   p(a) = int_l=a^inf  p(a|l)p(l) dl = = int_l=a^inf p(l) dl = 1-int_l=0^a p(l) dl


# https://www.scb.se/hitta-statistik/sverige-i-siffror/manniskorna-i-sverige/doda-i-sverige/#doda-efter-alder
d = [
	122,
	12,
	4,
	9,
	8,
	2,
	4,
	6,
	5,
	5,
	2,
	4,
	3,
	8,
	9,
	7,
	13,
	20,
	22,
	36,
	32,
	34,
	45,
	38,
	33,
	39,
	43,
	42,
	50,
	45,
	50,
	56,
	42,
	51,
	47,
	65,
	46,
	59,
	62,
	54,
	69,
	74,
	75,
	73,
	73,
	83,
	61,
	98,
	95,
	115,
	122,
	147,
	138,
	141,
	219,
	233,
	248,
	280,
	281,
	308,
	338,
	357,
	402,
	420,
	499,
	547,
	592,
	615,
	683,
	764,
	813,
	935,
	1034,
	1192,
	1219,
	1355,
	1477,
	1687,
	1659,
	1705,
	1748,
	1595,
	1760,
	1783,
	1706,
	1763,
	1652,
	1634,
	1593,
	1503,
	1494,
	1457,
	1181,
	1007,
	870,
	683,
	501,
	392,
	258,
	169,
	293,
]
d = d ./ sum(d)

function p(l)
	if 0 <= l < length(d)
		return d[round(Int, l + 1)]
	else
		return 0.0
	end
end

function nan2zero(x)
	if isnan(x)
		return 0.0
	else
		return x
	end
end

function p(ls, a)
	return ([p(l) / (1 - sum(p.(0:a))) * (l .>= a) for l in ls])
end
l = 0:150
figure(),
plot(l, p.(l))

a = 70
figure(),
plot(l, p(l, a))


l2 = 0:120

y = []
for a ∈ l2
	push!(y, sum(l .* p(l, a)) / sum(p(l, a)))
end
figure(), plot(l2, y)
