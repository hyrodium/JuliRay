## csgObject

abstract type csgObject <: Object end

export csgUnion
struct csgUnion <: csgObject
    objects :: Array{Object,1}

    function csgUnion(objects)
        if length(DeleteDuplicates(deleteat!(objects, isempty.(objects)))) == 0
            Empty()
        elseif length(DeleteDuplicates(deleteat!(objects, isempty.(objects)))) == 1
            objects[1]
        else
            new(DeleteDuplicates(deleteat!(objects, isempty.(objects))))
        end
    end
    function csgUnion(objects::Object...)
        return csgUnion([objects...])
    end
    function csgUnion(itr::Base.Generator)
        return csgUnion(collect(itr))
    end
end

export csgIntersection
struct csgIntersection <: csgObject
    objects :: Array{Object,1}

    function csgIntersection(objects)
        if Empty() ∈ objects
            Empty()
        elseif length(DeleteDuplicates(objects)) == 0
            Empty()
        elseif length(DeleteDuplicates(objects)) == 1
            objects[1]
        else
            new(DeleteDuplicates(objects))
        end
    end
    function csgIntersection(objects::Object...)
        return csgIntersection([objects...])
    end
    function csgIntersection(itr::Base.Generator)
        return csgIntersection(collect(itr))
    end
end

export csgMerge
struct csgMerge <: csgObject
    objects :: Array{Object,1}
    function csgMerge(objects)
        if length(DeleteDuplicates(deleteat!(objects, isempty.(objects)))) == 0
            Empty()
        elseif length(DeleteDuplicates(deleteat!(objects, isempty.(objects)))) == 1
            objects[1]
        else
            new(DeleteDuplicates(deleteat!(objects, isempty.(objects))))
        end
    end
    function csgMerge(objects::Object...)
        return csgMerge([objects...])
    end
    function csgMerge(itr::Base.Generator)
        return csgMerge(collect(itr))
    end
end

export csgDifference
struct csgDifference <: csgObject
    objects :: Array{Object,1}
    function csgDifference(objects)
        if length(objects) ≠ 2
            error("Too many objects.")
        elseif objects[1] == Empty()
            Empty()
        elseif objects[2] == Empty()
            objects[1]
        else
            new(objects)
        end
    end
    function csgDifference(object1::Object,object2::Object)
        return csgDifference([object1,object2])
    end
end

export csgClip
struct csgClip <: csgObject
    objects :: Array{Object,1}
    function csgClip(objects)
        if length(objects) ≠ 2
            error("Too many objects.")
        elseif objects[1] == Empty()
            Empty()
        elseif objects[2] == Empty()
            Empty()
        else
            new(objects)
        end
    end
    function csgClip(object1::Object,object2::Object)
        return csgClip([object1,object2])
    end
end

export csgBound
struct csgBound <: csgObject
    objects :: Array{Object,1}
    function csgBound(objects)
        if length(objects) ≠ 2
            error("Too many objects.")
        elseif objects[1] == Empty()
            Empty()
        elseif objects[2] == Empty()
            Empty()
        else
            new(objects)
        end
    end
    function csgBound(object1::Object,object2::Object)
        return csgBound([object1,object2])
    end
end

function povray_script(csg :: csgUnion; preindent = 0)
    script = "union{\n" * "  "^(preindent+1)
    for object in csg.objects
        script *= povray_script(object; preindent=preindent+1) * "\n" * "  "^(preindent+1)
    end
    script = script[1:end-2] * "}"
    return script
end

function povray_script(csg :: csgIntersection; preindent = 0)
    script = "intersection{\n" * "  "^(preindent+1)
    for object in csg.objects
        script *= povray_script(object; preindent=preindent+1) * "\n" * "  "^(preindent+1)
    end
    script = script[1:end-2] * "}"
    return script
end

function povray_script(csg :: csgMerge; preindent = 0)
    script = "merge{\n" * "  "^(preindent+1)
    for object in csg.objects
        script *= povray_script(object; preindent=preindent+1) * "\n" * "  "^(preindent+1)
    end
    script = script[1:end-2] * "}"
    return script
end

function povray_script(csg :: csgDifference; preindent = 0)
    script = "difference{\n" * "  "^(preindent+1)
    for object in csg.objects
        script *= povray_script(object; preindent=preindent+1) * "\n" * "  "^(preindent+1)
    end
    script = script[1:end-2] * "}"
    return script
end

function povray_script(csg :: csgClip; preindent = 0)
    return "object{"*povray_script(csg.objects[1])*"clipped_by{"*povray_script(csg.objects[2])*"}}"
end
function povray_script(csg :: csgBound; preindent = 0)
    return "object{"*povray_script(csg.objects[1])*"bounded_by{"*povray_script(csg.objects[2])*"}}"
end
