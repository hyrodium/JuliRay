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
        return csgIntersection(Torus(p₁,p₂,p₃,r),Blocks³(p₁,p₂,p₃,2r))
    end
end
