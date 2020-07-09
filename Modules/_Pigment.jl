using Colors

abstract type ColoredObject <: Object end
abstract type Pigment <: JR end

export FT
struct FT <: Pigment
    filter::Float64
    transmit::Float64

    function FT(filter,transmit)
        if filter<0 || transmit<0 || filter+transmit>1
            error("だめです")
        else
            new(filter,transmit)
        end
    end
end

export rgbColor
struct rgbColor <: ColoredObject
    object :: Object
    color :: Color

    function rgbColor(object,color)
        if object == Empty()
            Empty()
        else
            new(object,color)
        end
    end
end

export rgbftColor
struct rgbftColor <:ColoredObject
    object :: Object
    color :: Color
    transparence :: FT

    function rgbftColor(object,color,ft)
        if object == Empty()
            Empty()
        else
            new(object,color,ft)
        end
    end
end

export rgbftColor2
struct rgbftColor2 <:ColoredObject
    object :: Object
    color :: Color
    transparence :: FT

    function rgbftColor2(object,color,ft)
        if object == Empty()
            Empty()
        else
            new(object,color,ft)
        end
    end
end



# Color
function povray_script(color :: Color)
    r = string(Float64(color.r))
    g = string(Float64(color.g))
    b = string(Float64(color.b))
    return "rgb<"*r*","*g*","*b*">"
end

function povray_script(rgbcolor :: rgbColor; preindent = 0)
    script = "object{\n" * "  "^(preindent+1)
    script *= povray_script(rgbcolor.object, preindent=preindent+1) * "\n" * "  "^(preindent+1)
    script *= "pigment{" * povray_script(rgbcolor.color) * "}" * "\n" * "  "^preindent
    script *= "}"
    return script
end

function povray_script(rgbftcolor :: rgbftColor; preindent = 0)
    r = string(Float64(rgbftcolor.color.r))
    g = string(Float64(rgbftcolor.color.g))
    b = string(Float64(rgbftcolor.color.b))
    f = string(Float64(rgbftcolor.transparence.filter))
    t = string(Float64(rgbftcolor.transparence.transmit))

    script = "object{\n" * "  "^(preindent+1)
    script *= povray_script(rgbftcolor.object, preindent=preindent+1) * "\n" * "  "^(preindent+1)
    script *= "pigment{rgbft<"*r*","*g*","*b*","*f*","*t*">}" * "\n" * "  "^preindent
    script *= "}"
    return script
end

function povray_script(rgbftcolor :: rgbftColor2; preindent = 0)
    r = string(Float64(rgbftcolor.color.r))
    g = string(Float64(rgbftcolor.color.g))
    b = string(Float64(rgbftcolor.color.b))
    f = string(Float64(rgbftcolor.transparence.filter))
    t = string(Float64(rgbftcolor.transparence.transmit))

    script = "object{\n" * "  "^(preindent+1)
    script *= povray_script(rgbftcolor.object, preindent=preindent+1) * "\n" * "  "^(preindent+1)
    script *= "texture{pigment{rgbft<"*r*","*g*","*b*","*f*","*t*">}}" * "\n" * "  "^preindent
    script *= "texture{pigment{rgbft<"*r*","*g*","*b*","*f*","*t*">}}" * "\n" * "  "^preindent
    script *= "}"
    return script
end



export Transparent
function Transparent(object::PrimitiveObject,ft::FT)
    return object
end
function Transparent(csg::csgObject,ft::FT)
    csgtype = typeof(csg)
    return csgtype(map(obj->Transparent(obj,ft),csg.objects))
end

function Transparent(rgbcolor::rgbColor,ft::FT)
    return rgbftColor(rgbcolor.object,rgbcolor.color,ft)
end
function Transparent(rgbftcolor::rgbftColor,ft::FT)
    f0 = rgbftcolor.transparence.filter
    t0 = rgbftcolor.transparence.transmit
    f1 = ft.filter
    t1 = ft.transmit
    f2 = 1-(1-f0)*(1-f1)
    t2 = 1-(1-t0)*(1-t1)
    return rgbftColor(rgbftcolor.object,rgbftcolor.color,FT(f2,t2))
end
function Transparent(affinetransform::AffineTransform,ft::FT)
    return AffineTransform(Transparent(affinetransform.object,ft),affinetransform.matrix, affinetransform.vector)
end
function Transparent(paralleltranslation::ParallelTranslation,ft::FT)
    return ParallelTranslation(Transparent(paralleltranslation.object,ft),paralleltranslation.vector)
end
function Transparent(scaling::Scaling,ft::FT)
    return Scaling(Transparent(scaling.object,ft),scaling.scalar)
end
#TODO: Transparent for rgbftColors
