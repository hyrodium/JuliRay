push!(LOAD_PATH,"../Modules")
push!(LOAD_PATH,"Modules")

using Revise
using JuliRay

Sphere([1,1,32],3)

Torus(3,3)

a=Torus([1,1,1],[1,3,5],[3,23,-3],0.3)

povray_script(a)
