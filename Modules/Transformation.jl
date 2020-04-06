## Transform

abstract type TransformedObject <: Object end

export AffineTransform
struct AffineTransform <:TransformedObject
    object :: Object
    A :: Array{Float64,2}
    b :: Array{Float64,1}
    function AffineTransform(object,A,b)
        if object == Empty()
            Empty()
        elseif (rank(A) ≤ 1)
            Empty()
        else
            new(object,A,b)
        end
    end
    function AffineTransform(object::Object, A::RealMatrix; fixedpoint=[0.0,0.0,0.0])
        b = fixedpoint-A*fixedpoint
        return AffineTransform(object,A,b)
    end
end

export ParallelTranslation
struct ParallelTranslation <:TransformedObject
    object :: Object
    b :: Array{Float64,1}
    function ParallelTranslation(object,b)
        if object == Empty()
            Empty()
        else
            new(object,b)
        end
    end
end

export Scaling
struct Scaling <:TransformedObject
    object :: Object
    k :: Float64
    function Scaling(object, k; fixedpoint=[0.0,0.0,0.0])
        if object == Empty()
            Empty()
        elseif k == 0
            Empty()
        elseif norm(fixedpoint) ≠ 0
            new(ParallelTranslation(object,fixedpoint*(1-k)/k),k)
        else
            new(object,k)
        end
    end
    function Scaling(object::Object,k)
        b=fixedpoint-k*fixedpoint
        return ParallelTranslation(Scaling(object,k),b)
    end
end

function translate2pov(affinetransform :: AffineTransform)
    return "object{"*translate2pov(affinetransform.object)*" matrix "*translate2pov(vcat(reshape(affinetransform.A,9),affinetransform.b))*"}"
end
function translate2pov(paralleltranslation :: ParallelTranslation)
    return "object{"*translate2pov(paralleltranslation.object)*" translate "*translate2pov(paralleltranslation.b)*"}"
end
function translate2pov(scaling :: Scaling)
    return "object{"*translate2pov(scaling.object)*" scale "*translate2pov(scaling.k)*"}"
end

function Rotate(object::Object, v::RealVector, θ::Real, fixedpoint::RealVector=[0.0,0.0,0.0])
    return AffineTransform(object,rotatematrix(v,θ),fixedpoint=fixedpoint)
end

export Mirror
function Mirror(object::Object, p1::RealVector, p2::RealVector, p3::RealVector)
    e₃ = NormalVector(p1,p2,p3)
    e₁ = OrthogonalVector(e₃)
    e₂ = cross(e₃,e₁)
    A = hcat(e₁,e₂,e₃)*transpose(hcat(e₁,e₂,-e₃))
    return AffineTransform(object,A,fixedpoint=(p1+p2+p3)/3)
end
