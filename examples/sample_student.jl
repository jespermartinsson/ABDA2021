using PyPlot
using Statistics
using Random
using Dates
close("all")

function gen_seed()
    t = time()
    return round(Int,(t-floor(t))*1e6)
end

function sample_student(rows, n=1, rng=Xoshiro())
    row = (rows[1]-0.5) .+ (rows[end]-rows[1]+1)*rand(rng,n)
    return round.(Int,row)
end

function test_sampling(rows, n = 100_000, rng=Xoshiro())
    x = sample_student(rows,n,rng)
    bins = rows[1]:rows[end]
    nr = zeros(length(bins))
    for i in 1:length(bins)
        nr[i] = sum(x.==bins[i])
    end
    
    
    figure()
    subplot(211)
    plot(x[1:1000],1:1000,"o-",alpha=0.3)
    xlim([rows[1]-0.5,rows[end]+0.5])
    ylabel("sample index")
    subplot(212)
    bar(bins,nr)
    xlim([rows[1]-0.5,rows[end]+0.5])
    xlabel("row")
    ylabel("frequency")
end


rows = [4,23]
test_sampling(rows)

seed = gen_seed();
rng = Xoshiro(seed);
printstyled("seed: ", seed, "\n", color = :green)
printstyled("student at row: ", sample_student(rows,1,rng), "\n", color=:magenta)



