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
    center=Circumcenter(point1,point2,point3)
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
