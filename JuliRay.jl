module JuliRay

using LinearAlgebra
using Colors
using Printf

export


## Basic functions
function DeleteDuplicates(v::Array{T,1})::Array{T,1} where T <: Any
    w=Array{eltype(v),1}()
    if (length(v)>0)
        for e ∈ v
            if e ∉ w
                push!(w,e)
            end
        end
    end
    return w
end

function OrthogonalVector(v::Array{T,1}) where T <: Real
    if(norm(v)==0)
        error("zero vector")
    else
        v₁=normalize(v)
        i=argmin(abs.(v))
        v₂=zeros(length(v))
        v₂[i]=1.0
        v₂=normalize(v₂-dot(v₁,v₂)*v₁)
        return v₂
    end
end

function Circumcenter(point₁,point₂,point₃)
    l₁²=(norm(point₂-point₃))^2
    l₂²=(norm(point₃-point₁))^2
    l₃²=(norm(point₁-point₂))^2
    c₁=l₁²*(l₂²+l₃²-l₁²)
    c₂=l₂²*(l₃²+l₁²-l₂²)
    c₃=l₃²*(l₁²+l₂²-l₃²)
    center=(c₁*point₁+c₂*point₂+c₃*point₃)/(c₁+c₂+c₃)
    return center
end

function NormalVector(v₁::Array{Float64,1},v₂::Array{Float64,1},v₃::Array{Float64,1})
    return normalize(cross(v₁,v₂)+cross(v₂,v₃)+cross(v₃,v₁))
end


## Types
abstract type JuliRay end
abstract type Object <: JuliRay end
abstract type PrimitiveObject <: Object end
abstract type csgObject <: Object end
abstract type ColoredObject <: Object end
abstract type TransformedObject <: Object end

struct FT <: JuliRay
    filter::Float64
    transmit::Float64
    FT(filter,transmit) =
    if(filter<0 || transmit<0 || filter+transmit>1)
        error("だめです")
    else
        new(filter,transmit)
    end
end





## Primitive shapes
struct Empty <: PrimitiveObject end
struct Sphere <: PrimitiveObject
    center :: Array{Float64,1}
    radius :: Float64
    Sphere(center,radius) = (radius ≠ 0.0) ? new(center,abs(radius)) : Empty
end
struct Cylinder <: PrimitiveObject
    end1 :: Array{Float64,1}
    end2 :: Array{Float64,1}
    radius :: Float64
    Cylinder(end1,end2,radius) = (norm(end2-end1)*radius ≠ 0) ? new(end1,end2,abs(radius)) : Empty
end
struct Cone <: PrimitiveObject
    end1 :: Array{Float64,1}
    end2 :: Array{Float64,1}
    radius :: Float64
    Cone(end1,end2,radius) = (norm(end2-end1)*radius ≠ 0) ? new(end1,end2,abs(radius)) : Empty
end
struct Box <: PrimitiveObject
    vertex1 :: Array{Float64,1}
    vertex2 :: Array{Float64,1}
    Box(vertex1,vertex2) = (norm(vertex2-vertex1) ≠ 0) ? new(vertex1,vertex2) : Empty
end
struct Disc <: PrimitiveObject
    center :: Array{Float64,1}
    normal :: Array{Float64,1}
    radius :: Float64
    Disc(center,normal,radius) =
    if(norm(normal) == 0)
        error()
    elseif(radius == 0)
        Empty
    else
        new(center,normal,abs(radius))
    end
end
struct Torus <: PrimitiveObject
    radius1 :: Float64
    radius2 :: Float64
    Torus(radius1,radius2) =
    if(radius2 == 0)
        Empty
    elseif(radius1 == 0)
        Sphere([0,0,0],radius2)
    else
        new(radius1,radius2)
    end
end
struct Polygon <: PrimitiveObject
    vertices :: Array{Array{Float64,1},1}
    Polygon(vertices) =
    if(rank(hcat(vertices...)-repeat(+(vertices...)/length(vertices),1,length(vertices))) ≠ 2)
        Empty
    else
        new(vertices)
    end
end

struct csgUnion <: csgObject
    objects :: Array{Object,1}
    csgUnion(objects) =
    if(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 0)
        Empty
    elseif(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 1)
        objects[1]
    else
        new(DeleteDuplicates(deleteat!(objects, objects.== Empty)))
    end
end
struct csgIntersection <: csgObject
    objects :: Array{Object,1}
    csgIntersection(objects) =
    if(Empty ∈ objects)
        Empty
    elseif(length(DeleteDuplicates(objects)) == 0)
        Empty
    elseif(length(DeleteDuplicates(objects)) == 1)
        objects[1]
    else
        new(DeleteDuplicates(objects))
    end
end
struct csgMerge <: csgObject
    objects :: Array{Object,1}
    csgMerge(objects) =
    if(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 0)
        Empty
    elseif(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 1)
        objects[1]
    else
        new(DeleteDuplicates(deleteat!(objects, objects.== Empty)))
    end
end
struct csgDifference <: csgObject
    objects :: Array{Object,1}
    csgDifference(objects) =
    if(length(objects) ≠ 2)
        error("だめです")
    elseif(objects[1] == Empty)
        Empty
    elseif(objects[2] == Empty)
        objects[1]
    else
        new(objects)
    end
end

# Transform object to object
struct rgbColor <: ColoredObject
    object :: Object
    color :: Color
    rgbColor(object,color) = (object ≠ Empty) ? new(object,color) : Empty
end
struct rgbftColor <:ColoredObject
    object :: Object
    color :: Color
    transparence :: FT
    rgbftColor(object,color,ft) = (object ≠ Empty) ? new(object,color,ft) : Empty
end


struct AffineTransform <:TransformedObject
    object :: Object
    A :: Array{Float64,2}
    b :: Array{Float64,1}
    AffineTransform(object,A,b) = (object ≠ Empty && det(A) ≠ 0) ? new(object,A,b) : Empty
end





## Translation to POV-Ray code
function vector2pov(x::Array{T,1}) where T <: Real
    "<"*repr(x)[2:end-1]*">"
end

function object2pov(sphere :: Sphere)
    return "sphere{"*vector2pov(sphere.center)*","*repr(sphere.radius)*"}"
end
function object2pov(cylinder :: Cylinder)
    return "cylinder{"*vector2pov(cylinder.end1)*","*vector2pov(cylinder.end2)*","*repr(cylinder.radius)*"}"
end
function object2pov(cone :: Cone)
    return "cone{"*vector2pov(cone.end1)*","*repr(cone.radius)*","*vector2pov(cone.end2)*",0}"
end
function object2pov(box :: Box)
    return "box{"*vector2pov(box.vertex1)*","*vector2pov(box.vertex2)*"}"
end
function object2pov(disc :: Disc)
    return "disc{"*vector2pov(disc.center)*","*vector2pov(disc.normal)*","*repr(disc.radius)*"}"
end
function object2pov(torus :: Torus)
    return "torus{"*repr(torus.radius1)*","*repr(torus.radius2)*"}"
end
function object2pov(polygon :: Polygon)
    n=length(polygon.vertices)
    return "polygon{"*repr(n)* (*(reshape([(repeat([","],n),vector2pov.(polygon.vertices))[i][j] for i in 1:2, j in 1:n],2n)...))*"}"
end
# csg
function object2pov(csg :: csgUnion)
    return "union{"* *(object2pov.(csg.objects)...)*"}"
end
function object2pov(csg :: csgIntersection)
    return "intersection{"* *(object2pov.(csg.objects)...)*"}"
end
function object2pov(csg :: csgMerge)
    return "merge{"* *(object2pov.(csg.objects)...)*"}"
end
function object2pov(csg :: csgDifference)
    return "difference{"* *(object2pov.(csg.objects)...)*"}"
end
# Transform
function object2pov(rgbcolor :: rgbColor)
    r=string(Float64(rgbcolor.color.r))
    g=string(Float64(rgbcolor.color.g))
    b=string(Float64(rgbcolor.color.b))
    return "object{"*object2pov(rgbcolor.object)*" pigment{rgb<"*r*","*g*","*b*">}}"
end
function object2pov(rgbftcolor :: rgbftColor)
    r=string(Float64(rgbftcolor.color.r))
    g=string(Float64(rgbftcolor.color.g))
    b=string(Float64(rgbftcolor.color.b))
    f=string(Float64(rgbftcolor.transparence.filter))
    t=string(Float64(rgbftcolor.transparence.transmit))
    return "object{"*object2pov(rgbftcolor.object)*" pigment{rgbft<"*r*","*g*","*b*","*f*","*t*">}}"
end
function object2pov(affinetransform :: AffineTransform)
    return "object{"*object2pov(affinetransform.object)*" matrix"*vector2pov(vcat(reshape(affinetransform.A,9),affinetransform.b))*"}"
end

## Verargs objects
function Polygon(vertices...)
    return Polygon([vertices...])
end
function csgUnion(objects::Object...)
    return csgUnion([objects...])
end
function csgIntersection(objects::Object...)
    return csgIntersection([objects...])
end
function csgMerge(objects::Object...)
    return csgMerge([objects...])
end
function csgDifference(object1,object2)
    return csgDifference([object1,object2])
end


## Compound objects
function Arrow(end1, end2, radius)
    n=normalize(end2-end1)
    end3=end2-6*radius*n
    return csgUnion(Cylinder(end1,end3,radius),Cone(end3,end2,2*radius))
end
function Torus(center::Array{Float64,1}, normal::Array{Float64,1}, radius1::Float64, radius2::Float64)
    e₁=normalize(normal)
    e₂=OrthogonalVector(e₁)
    e₃=cross(e₁,e₂)
    A=hcat(e₃,e₁,e₂)
    return AffineTransform(Torus(radius1,radius2),A,center)
end
function Torus(point1::Array{Float64,1}, point2::Array{Float64,1}, point3::Array{Float64,1}, radius::Float64)
    center=Circumcenter(point1,point2,point3)
    normal=NormalVector(point1,point2,point3)
    R=norm(point1-center)
    return Torus(center,normal,R,radius)
end
function Disc(point1::Array{Float64,1}, point2::Array{Float64,1}, point3::Array{Float64,1})
    center=Circumcenter(point1,point2,point3)
    normal=NormalVector(point1,point2,point3)
    R=norm(point1-center)
    return Disc(center,normal,R)
end
function Blocks³(point1::Array{Float64,1}, point2::Array{Float64,1}, point3::Array{Float64,1}, thickness::Float64)
    center=Circumcenter(point1,point2,po,[2,3,4]int3)
    e₃=NormalVector(point1,point2,point3)
    R=norm(point1-center)
    e₁=normalize(point1-center)
    e₂=cross(e₃,e₁)
    box=Box([R+thickness,R+thickness,thickness],-[R+thickness,0,thickness])
    box1=AffineTransform(box,hcat(e₁,e₂,e₃),center)
    e₁=-normalize(point3-center)
    e₂=cross(e₃,e₁)
    box2=AffineTransform(box,hcat(e₁,e₂,e₃),center)
    # println(det(hcat(point1-center,point3-center,e₃)))
    if(det(hcat(point1-center,point3-center,e₃))>0)
        return csgIntersection(box1,box2)
    else
        return csgMerge(box1,box2)
    end
end
function Arc(point1::Array{Float64,1}, point2::Array{Float64,1}, point3::Array{Float64,1}, radius::Float64)
    return csgIntersection(Torus(point1,point2,point3,radius),Blocks³(point1,point2,point3,2radius))
end


function Transparent(object::PrimitiveObject,ft::FT)
    return object
end
function Transparent(csg::csgObject,ft::FT)
    csgtype=typeof(csg)
    return csgtype(map(obj->Transparent(obj,ft),csg.objects))
end

function Transparent(rgbcolor::rgbColor,ft::FT)
    return rgbftColor(rgbcolor.object,rgbcolor.color,ft)
end
function Transparent(rgbftcolor::rgbftColor,ft::FT)
    f0=rgbftcolor.transparence.filter
    t0=rgbftcolor.transparence.transmit
    f1=ft.filter
    t1=ft.transmit
    f2=1-(1-f0)*(1-f1)
    t2=1-(1-t0)*(1-t1)
    return rgbftColor(rgbftcolor.object,rgbftcolor.color,FT(f2,t2))
end
function Transparent(affinetransform::AffineTransform,ft::FT)
    return AffineTransform(Transparent(affinetransform.object,ft),affinetransform.A,affinetransform.b)
end


# Rendering
function render(obj::T;name="new", index::Int=0, width::Int=500, height::Int=500) where T <: Object
    if(0<index<1000000)
        Index="_"*(@sprintf "%06d" index)
    elseif(index==0)
        Index=""
    else
        error("index must be non-negative and less than 1000000")
    end
    if(endswith(name,".pov"))
        Name=name[1:end-4]*Index*".pov"
    else
        Name=name*Index*".pov"
    end

    str="#version 3.7;\nglobal_settings{assumed_gamma 1.0}\n\n#include \"Hy_constants.inc\"\n#include \"Hy_functions.inc\"\n#include \"Hy_colors.inc\"\n\n#declare Lng=30;\n#declare Lat=30;\n#declare Pers=0.1;\n#declare Zoom=0.9;\n#declare LookAt=<0,0,0>;\n#include \"Hy_camera.inc\"\n\n"
    str=str*object2pov(obj)
    io=open(Name,"w")
    write(io,str)
    close(io)

    str="Width=$width\nHeight=$height\nAntialias=On\n"
    io=open("povray.ini","w")
    write(io,str)
    close(io)

    run(`povray $Name`)
end

end
