## Primitive Object

abstract type PrimitiveObject <: Object end

export Empty
struct Empty <: PrimitiveObject end

export Sphere
struct Sphere <: PrimitiveObject
    center :: Array{Float64,1}
    radius :: Float64

    function Sphere(center,radius)
        if length(center) ≠ 3
            error("center of sphere must be an element of ℝ³")
        elseif radius == 0
            Empty()
        else
            new(center,abs(radius))
        end
    end
    function Sphere(point₁::RealVector,point₂::RealVector,point₃::RealVector,point₄::RealVector)
        center=Circumcenter(point₁,point₂,point₃,point₄)
        radius=Circumradius(point₁,point₂,point₃,point₄)
        return Sphere(center,radius)
    end
end

export Cylinder
struct Cylinder <: PrimitiveObject
    end1 :: Array{Float64,1}
    end2 :: Array{Float64,1}
    radius :: Float64

    function Cylinder(end1,end2,radius)
        if length(end1) ≠ 3
            error("end of cylinder must be an element of ℝ³")
        elseif length(end2) ≠ 3
            error("end of cylinder must be an element of ℝ³")
        elseif radius == 0
            Empty()
        elseif norm(end2-end1) == 0
            Empty()
        else
            new(end1,end2,abs(radius))
        end
    end
end

export Cone
struct Cone <: PrimitiveObject
    end1 :: Array{Float64,1}
    end2 :: Array{Float64,1}
    radius :: Float64

    function Cone(end1,end2,radius)
        if length(end1) ≠ 3
            error("end of cone must be an element of ℝ³")
        elseif length(end2) ≠ 3
            error("end of cone must be an element of ℝ³")
        elseif radius == 0
            Empty()
        elseif norm(end2-end1) == 0
            Empty()
        else
            new(end1,end2,abs(radius))
        end
    end
end

export Box
struct Box <: PrimitiveObject
    vertex1 :: Array{Float64,1}
    vertex2 :: Array{Float64,1}

    function Box(vertex1,vertex2)
        if length(vertex1) ≠ 3
            error("vertex of box must be an element of ℝ³")
        elseif length(vertex2) ≠ 3
            error("vertex of box must be an element of ℝ³")
        elseif norm(vertex2-vertex1) == 0
            Empty()
        else
            new(vertex1,vertex2)
        end
    end
end

export Disc
struct Disc <: PrimitiveObject
    center :: Array{Float64,1}
    normal :: Array{Float64,1}
    radius :: Float64

    function Disc(center,normal,radius)
        if length(center) ≠ 3
            error("center of disc must be an element of ℝ³")
        elseif length(normal) ≠ 3
            error("normal of disc must be an element of ℝ³")
        elseif norm(normal) == 0
            error("normal vector must be non-zero")
        elseif radius == 0
            Empty()
        else
            new(center,normal,abs(radius))
        end
    end
    function Disc(p₁::RealVector, p₂::RealVector, p₃::RealVector)
        center=Circumcenter(p₁,p₂,p₃)
        normal=NormalVector(p₁,p₂,p₃)
        R=norm(p₁-center)
        return Disc(center,normal,R)
    end
end

export Torus
struct Torus <: PrimitiveObject
    radius1 :: Float64
    radius2 :: Float64

    function Torus(radius1,radius2)
        if radius2 == 0
            Empty()
        elseif radius1 == 0
            Sphere([0,0,0],abs(radius2))
        else
            new(abs(radius1),abs(radius2))
        end
    end
    function Torus(center::RealVector, normal::RealVector, R::Float64, r::Float64)
        e₁=normalize(normal)
        e₂=OrthogonalVector(e₁)
        e₃=cross(e₁,e₂)
        A=hcat(e₃,e₁,e₂)
        return AffineTransform(Torus(R,r),A,center)
    end
    function Torus(p₁::RealVector, p₂::RealVector, p₃::RealVector, r::Float64; ε=10^(-4))
        α₁,α₂,α₃=Angles(p₁,p₂,p₃)
        if α₁>π-ε
            return Line(p₂,p₃,r)
        elseif α₂>π-ε
            return Line(p₃,p₁,r)
        elseif α₃>π-ε
            return Line(p₁,p₂,r)
        else
            center=Circumcenter(p₁,p₂,p₃)
            normal=NormalVector(p₁,p₂,p₃)
            R=(norm(p₁-center)+norm(p₂-center)+norm(p₃-center))/3
            return Torus(center,normal,R,r)
        end
    end

end

export Polygon
struct Polygon <: PrimitiveObject
    vertices :: Array{Array{Float64,1},1}

    function Polygon(vertices)
        if !all(e->e==3,length.(vertices))
            error("vertices of polygon must be an element of ℝ³")
        elseif rank(hcat(vertices...)-repeat(+(vertices...)/length(vertices),1,length(vertices)),atol=10^(-10)) ≠ 2
            Empty()
        else
            new(vertices)
        end
    end
    function Polygon(vertices::RealVector...)
        return Polygon([vertices...])
    end
end

function povray_script(sphere :: Sphere)
    return "sphere{"*povray_script(sphere.center)*","*povray_script(sphere.radius)*"}"
end
function povray_script(cylinder :: Cylinder)
    return "cylinder{"*povray_script(cylinder.end1)*","*povray_script(cylinder.end2)*","*povray_script(cylinder.radius)*"}"
end
function povray_script(cone :: Cone)
    return "cone{"*povray_script(cone.end1)*","*povray_script(cone.radius)*","*povray_script(cone.end2)*",0}"
end
function povray_script(box :: Box)
    return "box{"*povray_script(box.vertex1)*","*povray_script(box.vertex2)*"}"
end
function povray_script(disc :: Disc)
    return "disc{"*povray_script(disc.center)*","*povray_script(disc.normal)*","*povray_script(disc.radius)*"}"
end
function povray_script(torus :: Torus)
    return "torus{"*povray_script(torus.radius1)*","*povray_script(torus.radius2)*"}"
end
function povray_script(polygon :: Polygon)
    n=length(polygon.vertices)
    return "polygon{"*povray_script(n)* (*(reshape([(repeat([","],n),povray_script.(polygon.vertices))[i][j] for i in 1:2, j in 1:n],2n)...))*"}"
end

import Base: isempty

function isempty(object::Object)
    return repr(object)=="Empty()"
end
