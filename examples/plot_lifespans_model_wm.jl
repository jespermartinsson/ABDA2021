using PyPlot
using Stats
close("all")

figure()

for k ∈ 1:2
	subplot(212)
	dx = 0.25
	for i in 1:length(age)
		y = lifespans_m[i]
		x = age[i]
		if i == 1 && k == 1
			labels = ["obs. male", "mean", nothing]
		else
			labels = [nothing, nothing, nothing]
		end
		plot(age[i] * ones(length(y)) .- dx, y, "k."; label = labels[1], alpha = 0.5)
		if length(y) > 0
			plot(age[i] .- dx, mean(y), "dk"; label = labels[2])
			plot(age[i] * [1, 1] .- dx, [minimum(y), maximum(y)], "k-"; alpha = 0.5, label = labels[3])
		end
	end

	for i in 1:length(age)
		y = lifespans_w[i]
		x = age[i]
		if i == 1 && k == 1
			labels = ["obs. women", "mean", nothing]
		else
			labels = [nothing, nothing, nothing]
		end
		plot(age[i] * ones(length(y)) .+ dx, y, "g."; label = labels[1], alpha = 0.5)
		if length(y) > 0
			plot(age[i] .+ dx, mean(y), "dg"; label = labels[2])
			plot(age[i] * [1, 1] .+ dx, [minimum(y), maximum(y)], "g-"; alpha = 0.5, label = labels[3])
		end
	end
	axis([0, 120, 0, 140])
	grid("on")
	xlabel("age (years)")
	ylabel("lifespan (years)")





	# let a = age, and l = lifespan
	# p(l|a) = p(a|l)p(l)/p(a) 
	# we can not be older than our lifespan
	#   p(a|l) = 1 if a<=l, 0 otherwise 
	#   p(a) = int_l=a^inf  p(a|l)p(l) dl = = int_l=a^inf p(l) dl = 1-int_l=0^a p(l) dl


	# https://www.scb.se/hitta-statistik/sverige-i-siffror/manniskorna-i-sverige/doda-i-sverige/#doda-efter-alder
	dm = [
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

	dw = [
		107,
		12,
		3,
		9,
		2,
		3,
		4,
		4,
		0,
		4,
		2,
		4,
		1,
		7,
		11,
		11,
		4,
		9,
		11,
		11,
		8,
		5,
		12,
		17,
		5,
		17,
		16,
		16,
		12,
		19,
		27,
		20,
		24,
		30,
		22,
		22,
		20,
		22,
		35,
		31,
		29,
		41,
		34,
		30,
		36,
		38,
		51,
		61,
		73,
		79,
		91,
		112,
		113,
		115,
		127,
		128,
		182,
		194,
		200,
		220,
		227,
		233,
		258,
		279,
		293,
		366,
		385,
		420,
		459,
		537,
		597,
		644,
		744,
		819,
		911,
		1008,
		1187,
		1210,
		1318,
		1382,
		1321,
		1380,
		1517,
		1578,
		1638,
		1722,
		1751,
		1773,
		1841,
		1916,
		1918,
		1984,
		1828,
		1778,
		1570,
		1357,
		1140,
		1017,
		761,
		529,
		1086,
	]

	if k == 1
		d0 = dm
	else
		d0 = dw
	end
	d = d0 ./ sum(d0)

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
	# l = 0:150
	# figure(),
	# plot(l, p.(l))

	# a = 70
	# figure(),
	# plot(l, p(l, a))

	l = 0:120
	l2 = 0:120
	y = []
	for a ∈ l2
		push!(y, sum(l .* p(l, a)) / sum(p(l, a)))
	end
	if k == 1
		markercolor = "k-"
	else
		markercolor = "g-"
	end
	slabel = Dict(1 => "men", 2 => "women")
	plot(l2, y, markercolor, label = "BMMSE $(slabel[k])")
	legend()

	subplot(211)
	plot(0:length(d0)-1, d0, markercolor, label = "SCB $(slabel[k]) death by age (2022)")
	legend()
	grid("on")
	xlabel("age (years)")
	ylabel("nr deaths 2022")
	xlim([0, 120])

end
