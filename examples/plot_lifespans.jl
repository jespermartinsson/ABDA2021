using PyPlot
using Stats
close("all")

figure()
for i in 1:length(age)
	y = lifespans[i]
	x = age[i]
	if i == 1
		labels = ["observations", "mean", "range"]
	else
		labels = [nothing, nothing, nothing]
	end
	plot(age[i] * ones(length(lifespans[i])), lifespans[i], "k."; label = labels[1])
	plot(age[i], mean(lifespans[i]), "db"; alpha = 0.5, label = labels[2])
	plot(age[i] * [1, 1], [minimum(lifespans[i]), maximum(lifespans[i])], "b-"; alpha = 0.5, label = labels[3])
end
axis([0, 120, 0, 140])
grid("on")
xlabel("age (years)")
ylabel("lifespan (years)")
legend()
