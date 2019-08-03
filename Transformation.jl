## Transform

abstract type TransformedObject <: Object end


struct AffineTransform <:TransformedObject
    object :: Object
    A :: Array{Float64,2}
    b :: Array{Float64,1}
    AffineTransform(object,A,b) = (object ≠ Empty && rank(A) ≥ 2) ? new(object,A,b) : Empty
end

struct ParallelTranslation <:TransformedObject
    object :: Object
    b :: Array{Float64,1}
end

function translate2pov(affinetransform :: AffineTransform)
    return "object{"*translate2pov(affinetransform.object)*" matrix"*translate2pov(vcat(reshape(affinetransform.A,9),affinetransform.b))*"}"
end

function translate2pov(affinetransform :: ParallelTranslation)
    return "object{"*translate2pov(affinetransform.object)*" matrix"*translate2pov(vcat(reshape(affinetransform.A,9),affinetransform.b))*"}"
end


