using Colors

abstract type ColoredObject <: Object end
abstract type Pigment <: JR end

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

# Color
function translate2pov(color :: Color)
    r = string(Float64(color.r))
    g = string(Float64(color.g))
    b = string(Float64(color.b))
    return "rgb<"*r*","*g*","*b*">"
end
function translate2pov(rgbcolor :: rgbColor)
    return "object{"*translate2pov(rgbcolor.object)*" pigment{"*translate2pov(rgbcolor.color)*"}}"
end
function translate2pov(rgbftcolor :: rgbftColor)
    r = string(Float64(rgbftcolor.color.r))
    g = string(Float64(rgbftcolor.color.g))
    b = string(Float64(rgbftcolor.color.b))
    f = string(Float64(rgbftcolor.transparence.filter))
    t = string(Float64(rgbftcolor.transparence.transmit))
    return "object{"*translate2pov(rgbftcolor.object)*" pigment{rgbft<"*r*","*g*","*b*","*f*","*t*">}}"
end

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
    return AffineTransform(Transparent(affinetransform.object,ft),affinetransform.A,affinetransform.b)
end
function Transparent(affinetransform::ParallelTranslation,ft::FT)
    return ParallelTranslation(Transparent(affinetransform.object,ft),affinetransform.b)
end
