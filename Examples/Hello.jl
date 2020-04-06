push!(LOAD_PATH,"../Modules")
push!(LOAD_PATH,"Modules")

using Revise
using JuliRay


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

sum2(f, a) = mapreduce(f, +, a)

hoge(x::Real)=3

hoge(i for i in 1:5)

hoge(i for i in 1:5)

identity(3)


(i for i in 1:5)

sum2(identity, i for i in 1:5)

sum2(identity, i for i in 1:5)


function fuga(bg::Base.Generator)
    return bg(1)
end

fuga(i for i in 1:5)

@less mapreduce(identity,+,1:3)

reduce(push!, [1,2,3]; init=[])

collect(i for i in 1:3)


csgUnion()

csgUnion(Sphere([0,0,0],1))

csgUnion(Sphere([0,0,0],i) for i in 1:3)
