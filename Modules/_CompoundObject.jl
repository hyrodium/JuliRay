## Compound objects

export Line
function Line(p₁::RealVector, p₂::RealVector, r::Real; infty=10^5)
    v = normalize(p₂-p₁)
    return Cylinder(p₁-infty*v,p₂+infty*v,r)
end

export HalfLine
function HalfLine(p₁::RealVector, p₂::RealVector, r::Real; infty=10^5)
    v = normalize(p₂-p₁)
    return Cylinder(p₁,p₂+infty*v,r)
end

export Arrow
function Arrow(end1::RealVector, end2::RealVector, r::Real)
    n = normalize(end2-end1)
    end3 = end2-6*r*n
    return csgUnion(Cylinder(end1,end3,r),Cone(end3,end2,2*r))
end

function Blocks³(p₁::RealVector, p₂::RealVector, p₃::RealVector, thickness::Float64)
    center = Circumcenter(p₁,p₂,p₃)
    e₃ = NormalVector(p₁,p₂,p₃)
    R = norm(p₁-center)
    e₁ = normalize(p₁-center)
    e₂ = cross(e₃,e₁)
    box = Box([R+thickness,R+thickness,thickness],-[R+thickness,0,thickness])
    box1 = AffineTransform(box,hcat(e₁,e₂,e₃),center)
    e₁ = -normalize(p₃-center)
    e₂ = cross(e₃,e₁)
    box2 = AffineTransform(box,hcat(e₁,e₂,e₃),center)
    if det(hcat(p₁-center,p₃-center,e₃)) > 0
        return csgIntersection(box1,box2)
    else
        return csgMerge(box1,box2)
    end
end

export Arc
function Arc(p₁::RealVector, p₂::RealVector, p₃::RealVector, r::Float64; ε=10^(-4))
    α₁,α₂,α₃ = Angles(p₁,p₂,p₃)
    if (α₁>π-ε)
        return csgUnion(HalfLine(p₁,p₁+(p₂-p₃),r),HalfLine(p₃,p₃-(p₂-p₃),r))
    elseif (α₂>π-ε)
        return Cylinder(p₁,p₃,r)
    elseif (α₃>π-ε)
        return csgUnion(HalfLine(p₁,p₁+(p₁-p₂),r),HalfLine(p₃,p₃-(p₁-p₂),r))
    else
        return csgBound(csgIntersection(Torus(p₁,p₂,p₃,r),Blocks³(p₁,p₂,p₃,2r)),Torus(p₁,p₂,p₃,r))
    end
end

export RoundedBox
function RoundedBox(vertex1::RealVector, vertex2::RealVector, radius::Real)
    if radius == 0
        return Box(vertex1, vertex2)
    elseif radius < 0
        error("radius must be non-negative")
    end
    x=vertex1[1], vertex2[1]
    y=vertex1[2], vertex2[2]
    z=vertex1[3], vertex2[3]
    x_min, x_max=minimum(x), maximum(x)
    y_min, y_max=minimum(y), maximum(y)
    z_min, z_max=minimum(z), maximum(z)

    v_min=[x_min, y_min, z_min] .+ radius
    v_max=[x_max, y_max, z_max] .- radius

    x_r=v_min[1], v_max[1]
    y_r=v_min[2], v_max[2]
    z_r=v_min[3], v_max[3]

    spheres=[Sphere([x_r[i], y_r[j], z_r[k]], radius) for i in 1:2, j in 1:2, k in 1:2]

    cylinders_x=[Cylinder([x_r[1], y_r[j], z_r[k]], [x_r[2], y_r[j], z_r[k]], radius) for j in 1:2, k in 1:2]
    cylinders_y=[Cylinder([x_r[i], y_r[1], z_r[k]], [x_r[i], y_r[2], z_r[k]], radius) for i in 1:2, k in 1:2]
    cylinders_z=[Cylinder([x_r[i], y_r[j], z_r[1]], [x_r[i], y_r[j], z_r[2]], radius) for i in 1:2, j in 1:2]

    box_x=Box([x_min, y_r[1], z_r[1]], [x_max, y_r[2], z_r[2]])
    box_y=Box([x_r[1], y_min, z_r[1]], [x_r[2], y_max, z_r[2]])
    box_z=Box([x_r[1], y_r[1], z_min], [x_r[2], y_r[2], z_max])

    return csgMerge(spheres..., cylinders_x..., cylinders_y..., cylinders_z..., box_x, box_y, box_z)
end
