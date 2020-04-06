push!(LOAD_PATH,"../Modules")
push!(LOAD_PATH,"Modules")

using Revise
using JuliRay

include("../Modules/JuliRay.jl")

Sphere([1,1,32],3)

Torus(3,3)

a=Torus([1,1,1],[1,3,5],[3,23,-3],0.3)

translate2pov(a)

f(x)=x

f(i for i in 1:5)

@less sum(i for i in 1:3)

@less sum(1)

sum2(f, a) = mapreduce(f, +, a)
sum2(a) = sum2(identity, a)
sum2(a::AbstractArray{Bool}) = count(a)

sum2(i for i in 1:3)

Base.add_sum
