using PyPlot
using Stats
close("all")

dx = 0.25
figure()
for i in 1:length(age)
	y = lifespans_m[i]
	x = age[i]
	if i == 1
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
	if i == 1
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
legend()
