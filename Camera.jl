## Camera

abstract type Camera <: JuliRay end
abstract type Light <: JuliRay end

struct PerspectiveCamera <: Camera
    attitude :: Array{Float64,2}
    position :: Array{Float64,1}
    area::Float64
    width::Int
    height::Int
    color::Color
    PerspectiveCamera(attitude,position,area,width,height,color)=
    if (area ≤ 0)
        error("area must be positive")
    elseif (area<0)
        error("だめです")
    elseif (!isorthogonal(attitude))
        error("だめです")
    else
        new(attitude,position,area,width,height,color)
    end
end

function translate2pov(camera::PerspectiveCamera)
    e₁=camera.attitude[:,1]
    e₂=camera.attitude[:,2]
    e₃=camera.attitude[:,3]
    position=camera.position
    up=sqrt(camera.area*camera.height/camera.width)/2
    right=-sqrt(camera.area*camera.width/camera.height)/2
    color=camera.color
    light=PointLight(position,color)
    str="camera{perspective location "*translate2pov(position)*" right "*translate2pov(e₁*right)*" up "*translate2pov(e₂*up)*" direction "*translate2pov(e₃)*" sky "*translate2pov(e₂)*" look_at "*translate2pov(position-e₃)*"}\n"
    str=str*translate2pov(light)
    return str
end

struct ParallelCamera <: Camera
    attitude :: Array{Float64,2}
    position :: Array{Float64,1}
    area::Float64
    width::Int
    height::Int
    color::Color
    PerspectiveCamera(attitude,position,area,width,height,color)=
    if (area ≤ 0)
        error("area must be positive")
    elseif (area<0)
        error("だめです")
    elseif (!isorthogonal(attitude))
        error("だめです")
    else
        new(attitude,position,area,width,height,color)
    end
end



function LngLatCamera(;lng=-π/2,lat=π/2,tilt=0,pers=0.1,zoom=1.0,lookat=[0.0,0.0,0.0],width::Int=500,height::Int=500,color::Color=RGB(1,1,1))
    e₃=[cos(lat)*cos(lng),cos(lat)*sin(lng),sin(lat)]
    e₁=rotatematrix(e₃,tilt)*[-sin(lng),cos(lng),0.0]
    e₂=cross(e₃,e₁)
    attitude=hcat(e₁,e₂,e₃)
    position=lookat+e₃/(zoom*pers)
    area=4*pers*pers
    return PerspectiveCamera(attitude,position,area,width,height,color)
end

## Light
struct PointLight <: Light
    position :: Array{Float64,1}
    color :: Color
end

function translate2pov(light::PointLight)
    str="light_source{"*translate2pov(light.position)*" "*translate2pov(light.color)*"}"
    return str
end
