module JuliRay

using LinearAlgebra
using Printf

include("BasicFunctions.jl")

## Translation to POV-Ray code
function translate2pov(x::Real)
    if (isfinite(x))
        return @sprintf "%1.24f" x
    else
        error("numerical value must be finite")
    end
end

function translate2pov(x::Array{T,1}) where T <: Real
    y=repr(translate2pov.(x))
    y=replace(y,"\""=>"")
    return "<"*y[2:end-1]*">"
end

# Types
abstract type JR end
abstract type Object <: JR end

include("PrimitiveObject.jl")
include("CSG.jl")
include("Transformation.jl")
include("Pigment.jl")
include("CompoundObject.jl")
include("Camera.jl")

# Rendering
function render(object::Object;name="new", index::Int=0, camera::Camera=LngLatCamera(), lights::Array{T,1}=Light[]) where T <: Light
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
        Name = name[1:end-4]*Index*".pov"
    else
        Name = name*Index*".pov"
    end

    str = "#version 3.7;\nglobal_settings{assumed_gamma 1.0}\n"
    str = str*translate2pov(camera)*"\n"
    str = str*(join(translate2pov.(lights)))
    str = str*"background{rgb<1,1,1>}"*"\n"
    str = str*translate2pov(object)*"\n"
    io = open(Name,"w")
    write(io,str)
    close(io)

    str = "Width=$width\nHeight=$height\nAntialias=On\n"
    io = open("povray.ini","w")
    write(io,str)
    close(io)

    run(`povray $Name`)
end

end # module
