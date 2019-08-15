## Primitive Object

abstract type PrimitiveObject <: Object end

struct Empty <: PrimitiveObject end
struct Sphere <: PrimitiveObject
    center :: Array{Float64,1}
    radius :: Float64
    Sphere(center,radius) =
    if (length(center) ≠ 3)
        error("center of sphere must be an element of ℝ³")
    elseif (radius == 0)
        Empty
    else
        new(center,abs(radius))
    end
end
struct Cylinder <: PrimitiveObject
    end1 :: Array{Float64,1}
    end2 :: Array{Float64,1}
    radius :: Float64
    Cylinder(end1,end2,radius) =
    if (length(end1) ≠ 3)
        error("end of cylinder must be an element of ℝ³")
    elseif (length(end2) ≠ 3)
        error("end of cylinder must be an element of ℝ³")
    elseif (radius == 0)
        Empty
    elseif (norm(end2-end1) == 0)
        Empty
    else
        new(end1,end2,abs(radius))
    end
end
struct Cone <: PrimitiveObject
    end1 :: Array{Float64,1}
    end2 :: Array{Float64,1}
    radius :: Float64
    Cone(end1,end2,radius) =
    if (length(end1) ≠ 3)
        error("end of cone must be an element of ℝ³")
    elseif (length(end2) ≠ 3)
        error("end of cone must be an element of ℝ³")
    elseif (radius == 0)
        Empty
    elseif (norm(end2-end1) == 0)
        Empty
    else
        new(end1,end2,abs(radius))
    end
end
struct Box <: PrimitiveObject
    vertex1 :: Array{Float64,1}
    vertex2 :: Array{Float64,1}
    Box(vertex1,vertex2) =
    if (length(vertex1) ≠ 3)
        error("vertex of box must be an element of ℝ³")
    elseif (length(vertex2) ≠ 3)
        error("vertex of box must be an element of ℝ³")
    elseif (norm(vertex2-vertex1) == 0)
        Empty
    else
        new(vertex1,vertex2)
    end
end
struct Disc <: PrimitiveObject
    center :: Array{Float64,1}
    normal :: Array{Float64,1}
    radius :: Float64
    Disc(center,normal,radius) =
    if (length(center) ≠ 3)
        error("center of disc must be an element of ℝ³")
    elseif (length(normal) ≠ 3)
        error("normal of disc must be an element of ℝ³")
    elseif (norm(normal) == 0)
        error("normal vector must be non-zero")
    elseif (radius == 0)
        Empty
    else
        new(center,normal,abs(radius))
    end
end
struct Torus <: PrimitiveObject
    radius1 :: Float64
    radius2 :: Float64
    Torus(radius1,radius2) =
    if (radius2 == 0)
        Empty
    elseif (radius1 == 0)
        Sphere([0,0,0],abs(radius2))
    else
        new(abs(radius1),abs(radius2))
    end
end
struct Polygon <: PrimitiveObject
    vertices :: Array{Array{Float64,1},1}
    Polygon(vertices) =
    if (!all(e->e==3,length.(vertices)))
        error("vertex of polygon must be an element of ℝ³")
    elseif (rank(hcat(vertices...)-repeat(+(vertices...)/length(vertices),1,length(vertices))) ≠ 2)
        Empty
    else
        new(vertices)
    end
end

function Polygon(vertices::RealVector...)
    return Polygon([vertices...])
end

function translate2pov(sphere :: Sphere)
    return "sphere{"*translate2pov(sphere.center)*","*translate2pov(sphere.radius)*"}"
end
function translate2pov(cylinder :: Cylinder)
    return "cylinder{"*translate2pov(cylinder.end1)*","*translate2pov(cylinder.end2)*","*translate2pov(cylinder.radius)*"}"
end
function translate2pov(cone :: Cone)
    return "cone{"*translate2pov(cone.end1)*","*translate2pov(cone.radius)*","*translate2pov(cone.end2)*",0}"
end
function translate2pov(box :: Box)
    return "box{"*translate2pov(box.vertex1)*","*translate2pov(box.vertex2)*"}"
end
function translate2pov(disc :: Disc)
    return "disc{"*translate2pov(disc.center)*","*translate2pov(disc.normal)*","*translate2pov(disc.radius)*"}"
end
function translate2pov(torus :: Torus)
    return "torus{"*translate2pov(torus.radius1)*","*translate2pov(torus.radius2)*"}"
end
function translate2pov(polygon :: Polygon)
    n=length(polygon.vertices)
    return "polygon{"*translate2pov(n)* (*(reshape([(repeat([","],n),translate2pov.(polygon.vertices))[i][j] for i in 1:2, j in 1:n],2n)...))*"}"
end
