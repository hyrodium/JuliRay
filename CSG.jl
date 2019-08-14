## csgObject

abstract type csgObject <: Object end

struct csgUnion <: csgObject
    objects :: Array{Object,1}
    csgUnion(objects) =
    if(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 0)
        Empty
    elseif(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 1)
        objects[1]
    else
        new(DeleteDuplicates(deleteat!(objects, objects.== Empty)))
    end
end
struct csgIntersection <: csgObject
    objects :: Array{Object,1}
    csgIntersection(objects) =
    if(Empty ∈ objects)
        Empty
    elseif(length(DeleteDuplicates(objects)) == 0)
        Empty
    elseif(length(DeleteDuplicates(objects)) == 1)
        objects[1]
    else
        new(DeleteDuplicates(objects))
    end
end
struct csgMerge <: csgObject
    objects :: Array{Object,1}
    csgMerge(objects) =
    if(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 0)
        Empty
    elseif(length(DeleteDuplicates(deleteat!(objects, objects.== Empty))) == 1)
        objects[1]
    else
        new(DeleteDuplicates(deleteat!(objects, objects.== Empty)))
    end
end
struct csgDifference <: csgObject
    objects :: Array{Object,1}
    csgDifference(objects) =
    if(length(objects) ≠ 2)
        error("Too many objects.")
    elseif(objects[1] == Empty)
        Empty
    elseif(objects[2] == Empty)
        objects[1]
    else
        new(objects)
    end
end
struct csgClip <: csgObject
    objects :: Array{Object,1}
    csgDifference(objects) =
    if(length(objects) ≠ 2)
        error("Too many objects.")
    elseif(objects[1] == Empty)
        Empty
    elseif(objects[2] == Empty)
        objects[1]
    else
        new(objects)
    end
end

function csgUnion(objects::Object...)
    return csgUnion([objects...])
end
function csgIntersection(objects::Object...)
    return csgIntersection([objects...])
end
function csgMerge(objects::Object...)
    return csgMerge([objects...])
end
function csgDifference(object1::Object,object2::Object)
    return csgDifference([object1,object2])
end
function csgClip(object1::Object,object2::Object)
    return csgClip([object1,object2])
end

function translate2pov(csg :: csgUnion)
    return "union{"* *(translate2pov.(csg.objects)...)*"}"
end
function translate2pov(csg :: csgIntersection)
    return "intersection{"* *(translate2pov.(csg.objects)...)*"}"
end
function translate2pov(csg :: csgMerge)
    return "merge{"* *(translate2pov.(csg.objects)...)*"}"
end
function translate2pov(csg :: csgDifference)
    return "difference{"* *(translate2pov.(csg.objects)...)*"}"
end
function translate2pov(csg :: csgClip)
    return "object{"*translate2pov(csg.objects[1])*"clipped_by{"*translate2pov(csg.objects[2])*"}}"
end
