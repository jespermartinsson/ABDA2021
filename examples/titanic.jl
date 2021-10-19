using CSV
using DataFrames
df = CSV.read(DataFrame,(@__DIR__)*"/data/titanic.csv")
