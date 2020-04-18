module JuliRay

using LinearAlgebra
using Printf

include("_BasicFunctions.jl")

## Translation to POV-Ray code

"""
2.3 -> "2.3"
"""
function povray_script(x::Real)
    if isfinite(x)
        return repr(x)
        # return @sprintf "%1.24f" x
    else
        error("numerical value must be finite")
    end
end

"""
[2.3, -5.2] -> "<2.3,-5.2>"
"""
function povray_script(x::Array{T,1}) where T <: Real
    return "<"*join(povray_script.(x), ", ")*">"
end

"""
[[2.3, -5.2], [-1.1, 8.2]] -> "<2.3,-5.2>, <-1.1,8.2>"
"""
function povray_script(x::Array{Array{T,1},1}) where T <: Real
    return join(povray_script.(x), ", ")
end

# Types
abstract type JR end
abstract type Object <: JR end

include("_PrimitiveObject.jl")
include("_CSG.jl")
include("_Transformation.jl")
include("_Pigment.jl")
include("_CompoundObject.jl")
include("_ParametricObject.jl")
include("_Camera.jl")

BASE_DIR = "."

# Rendering
export render
function render(object::Object; name="new", index::Int=0, camera::Camera=LngLatCamera(), lights::Array{T,1}=Light[]) where T <: Light
    width=camera.width
    height=camera.height

    if 0 < index < 1000000
        Index = "_"*(@sprintf "%06d" index)
    elseif index == 0
        Index = ""
    else
        error("index must be non-negative and less than 1000000")
    end
    if endswith(name,".pov")
        Name = name[1:end-4]*Index
    else
        Name = name*Index
    end

    str = "#version 3.7;\nglobal_settings{assumed_gamma 1.0}\n"
    str = str*povray_script(camera)*"\n"
    str = str*(join(povray_script.(lights)))
    str = str*"background{rgb<1,1,1>}"*"\n"
    str = str*povray_script(object)*"\n"
    path_pov = BASE_DIR * "/" * Name * ".pov"
    io = open(path_pov,"w")
    write(io,str)
    close(io)

    str = "Width=$width\nHeight=$height\nAntialias=On\n"
    path_ini = BASE_DIR * "/" * "povray.ini"
    io = open(path_ini,"w")
    write(io,str)
    close(io)

    run(`povray $path_pov`)
end

end # module
