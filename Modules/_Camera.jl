# Camera

abstract type Camera <: JR end
abstract type Light <: JR end

export PerspectiveCamera
struct PerspectiveCamera <: Camera
    attitude :: Array{Float64,2}
    position :: Array{Float64,1}
    area::Float64
    width::Int
    height::Int
    color::Color
    function PerspectiveCamera(attitude,position,area,width,height,color)
        if area ≤ 0
            error("area must be positive")
        elseif !isorthogonal(attitude)
            error("だめです")
        else
            new(attitude,position,area,width,height,color)
        end
    end
end

function povray_script(camera::PerspectiveCamera)
    e₁=camera.attitude[:,1]
    e₂=camera.attitude[:,2]
    e₃=camera.attitude[:,3]
    position=camera.position
    up=sqrt(camera.area*camera.height/camera.width)/2
    right=-sqrt(camera.area*camera.width/camera.height)/2
    color=camera.color
    light=PointLight(position,color)
    str="camera{perspective location "*povray_script(position)*" right "*povray_script(e₁*right)*" up "*povray_script(e₂*up)*" direction "*povray_script(e₃)*" sky "*povray_script(e₂)*" look_at "*povray_script(position-e₃)*"}\n"
    str=str*povray_script(light)
    return str
end

export ParallelCamera
struct ParallelCamera <: Camera
    attitude :: Array{Float64,2}
    position :: Array{Float64,1}
    area::Float64
    width::Int
    height::Int
    color::Color
    function PerspectiveCamera(attitude,position,area,width,height,color)
        if area ≤ 0
            error("area must be positive")
        elseif !isorthogonal(attitude)
            error("だめです")
        else
            new(attitude,position,area,width,height,color)
        end
    end
end

export LngLatCamera
function LngLatCamera(;lng=-π/2,lat=π/2,tilt=0,pers=0.1,zoom=1.0,lookat=[0.0,0.0,0.0],width::Int=500,height::Int=500,color::Color=RGB(1,1,1))
    e₃=[cos(lat)*cos(lng),cos(lat)*sin(lng),sin(lat)]
    e₁=rotatematrix(e₃,tilt)*[-sin(lng),cos(lng),0.0]
    e₂=cross(e₃,e₁)
    attitude=hcat(e₁,e₂,e₃)
    position=lookat+e₃/(zoom*pers)
    area=4*pers*pers
    return PerspectiveCamera(attitude,position,area,width,height,color)
end

# Light
export PointLight
struct PointLight <: Light
    position :: Array{Float64,1}
    color :: Color
end

function povray_script(light::PointLight)
    str="light_source{"*povray_script(light.position)*" "*povray_script(light.color)*"}"
    return str
end
