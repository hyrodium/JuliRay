## Transform

abstract type TransformedObject <: Object end

export AffineTransform
struct AffineTransform <: TransformedObject
    object :: Object
    matrix :: Array{Float64,2}
    vector :: Array{Float64,1}
    function AffineTransform(object,matrix,vector)
        if object == Empty()
            Empty()
        elseif (rank(matrix) ≤ 1)
            Empty()
        else
            new(object,matrix,vector)
        end
    end
    function AffineTransform(object::Object, matrix::RealMatrix; fixedpoint=[0.0,0.0,0.0])
        vector = fixedpoint - matrix*fixedpoint
        return AffineTransform(object,matrix,vector)
    end
end

export ParallelTranslation
struct ParallelTranslation <: TransformedObject
    object :: Object
    vector :: Array{Float64,1}
    function ParallelTranslation(object,vector)
        if object == Empty()
            Empty()
        else
            new(object,vector)
        end
    end
end

export Scaling
struct Scaling <: TransformedObject
    object :: Object
    scalar :: Float64
    function Scaling(object, scalar; fixedpoint=[0.0,0.0,0.0])
        if object == Empty()
            Empty()
        elseif scalar == 0
            Empty()
        elseif norm(fixedpoint) ≠ 0
            new(ParallelTranslation(object,fixedpoint*(1-scalar)/scalar),scalar)
        else
            new(object,scalar)
        end
    end
end

function povray_script(affinetransform :: AffineTransform; preindent = 0)
    script = "object{\n" * "  "^(preindent+1)
    script *= povray_script(affinetransform.object, preindent=preindent+1) * "\n" * "  "^(preindent+1)
    script *= "matrix " * povray_script(vcat(reshape(affinetransform.matrix,9),affinetransform.vector)) * "\n" * "  "^preindent
    script *= "}"
    return script
end

function povray_script(paralleltranslation :: ParallelTranslation; preindent = 0)
    script = "object{\n" * "  "^(preindent+1)
    script *= povray_script(paralleltranslation.object, preindent=preindent+1) * "\n" * "  "^(preindent+1)
    script *= "translate "*povray_script(paralleltranslation.vector) * "\n" * "  "^preindent
    script *= "}"
    return script
end

function povray_script(scaling :: Scaling; preindent = 0)
    script = "object{\n" * "  "^(preindent+1)
    script *= povray_script(scaling.object, preindent=preindent+1) * "\n" * "  "^(preindent+1)
    script *= "scale "*povray_script(scaling.scalar) * "\n" * "  "^preindent
    script *= "}"
    return script
end

export Rotate
function Rotate(object::Object, v::RealVector, θ::Real; fixedpoint::RealVector=[0.0,0.0,0.0])
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
