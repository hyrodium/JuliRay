using IntervalSets

using Random
using ForwardDiff
using Statistics

export ParametricCurve
function ParametricCurve(parametrization::Function, D::Interval, r::Real; mesh::Int=50)
    n = mesh
    a, b = endpoints(D)
    𝒑 = parametrization
    ts = range(a,b,length=n+1)
    pts = csgUnion(Sphere(𝒑(t),r) for t in ts)
    segts = csgUnion(Cylinder(𝒑(ts[i]),𝒑(ts[i+1]),r) for i in 1:n)
    return csgUnion(pts,segts)
end

export ParametricSurface
struct ParametricSurface <: PrimitiveObject
    parametrization :: Function
    D1 :: Interval{:closed,:closed,Float64}
    D2 :: Interval{:closed,:closed,Float64}
    mesh :: Tuple{Int,Int}
    smooth :: Bool
    function ParametricSurface(parametrization::Function, D1::Interval, D2::Interval; mesh::Tuple{Int,Int}=(10,10), smooth=true)
        new(parametrization, ClosedInterval(D1), ClosedInterval(D2), mesh, smooth)
    end
end

export ColoredParametricSurface
struct ColoredParametricSurface <: PrimitiveObject
    parametrization :: Function
    D1 :: Interval{:closed,:closed,Float64}
    D2 :: Interval{:closed,:closed,Float64}
    mesh :: Tuple{Int,Int}
    smooth :: Bool
    rgb :: Function
    ft :: Function
    smooth_texture :: Bool
    function ColoredParametricSurface(parametrization::Function, D1::Interval, D2::Interval; mesh::Tuple{Int,Int}=(10,10), smooth=true, rgb=(t->RGB(1,0,0)), ft=(t->FT(0,0)), smooth_texture=false)
        new(parametrization, ClosedInterval(D1), ClosedInterval(D2), mesh, smooth, rgb, ft, smooth_texture)
    end
end

function povray_script(parametricsurface::ParametricSurface; preindent=0)
    𝒑 = parametricsurface.parametrization
    D1, D2 = parametricsurface.D1, parametricsurface.D2
    mesh = parametricsurface.mesh
    smooth = parametricsurface.smooth
    return mesh2(𝒑, D1, D2; mesh=mesh, smooth=smooth, preindent=preindent)
end

function povray_script(parametricsurface::ColoredParametricSurface; preindent=0)
    𝒑 = parametricsurface.parametrization
    D1, D2 = parametricsurface.D1, parametricsurface.D2
    mesh = parametricsurface.mesh
    smooth = parametricsurface.smooth
    rgb = parametricsurface.rgb
    ft = parametricsurface.ft
    smooth_texture = parametricsurface.smooth_texture
    return mesh2_with_texture(𝒑, D1, D2; mesh=mesh, smooth=smooth, rgb=rgb, ft=ft, smooth_texture=smooth_texture, preindent=preindent)
end

function mesh2(𝒑, D1, D2; mesh=(10,10), smooth=true, preindent=0)
    n1, n2 = mesh
    D1₋, D1₊ = endpoints(D1)
    D2₋, D2₊ = endpoints(D2)
    Δ1, Δ2 = Δ = (D1₊-D1₋)/n1, (D2₊-D2₋)/n2

    ts = [[t1, t2] for t1 in range(D1₋, D1₊, length = n1 +1), t2 in range(D2₋, D2₊, length = n2 +1)]
    tc = [mean([ts[i1,i2], ts[i1+1,i2], ts[i1,i2+1], ts[i1+1,i2+1]]) for i1 in 1:n1, i2 in 1:n2]

    𝒑′(t)=ForwardDiff.jacobian(𝒑,t)
    𝒆(t) = normalize(cross(𝒑′(t)[1:3,1],𝒑′(t)[1:3,2]))

    𝒑s = 𝒑.(ts)
    𝒑c = 𝒑.(tc)

    𝒆s = 𝒆.(ts)
    𝒆c = 𝒆.(tc)

    Ns(i1, i2) = i1 + (n1+1) * (i2-1)
    Nc(i1, i2) = i1 + n1 * (i2-1)

    F1 = [[Ns(i1,i2)-1, Ns(i1+1,i2)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F2 = [[Ns(i1,i2)-1, Ns(i1,i2+1)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F3 = [[Ns(i1+1,i2+1)-1, Ns(i1+1,i2)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F4 = [[Ns(i1+1,i2+1)-1, Ns(i1,i2+1)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]

    np = (n1+1)*(n2+1) + n1*n2
    nf = 4*n1*n2

    script = ""
    script *= "mesh2{\n" * "  "^(preindent)
    script *= "  vertex_vectors{\n" * "  "^(preindent)
    script *= "    " * povray_script(np) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([𝒑s...]) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([𝒑c...]) * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)

    if smooth
        script *= "  normal_vectors{\n" * "  "^(preindent)
        script *= "    " * povray_script(np) * ", \n" * "  "^(preindent)
        script *= "    " * povray_script([𝒆s...]) * ", \n" * "  "^(preindent)
        script *= "    " * povray_script([𝒆c...]) * "\n" * "  "^(preindent)
        script *= "  }\n" * "  "^(preindent)
    end

    script *= "  face_indices{\n" * "  "^(preindent)
    script *= "    " * povray_script(nf) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([F1...]) * "\n" * "  "^(preindent)
    script *= "    " * povray_script([F2...]) * "\n" * "  "^(preindent)
    script *= "    " * povray_script([F3...]) * "\n" * "  "^(preindent)
    script *= "    " * povray_script([F4...]) * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)
    script *= "}"

    return script
end


function mesh2_with_texture(𝒑, D1, D2; mesh=(10,10), smooth=true, smooth_texture=true, rgb=(t->RGB(1,0,0)), ft=(t->FT(0,0)), preindent=0)
    n1, n2 = mesh
    D1₋, D1₊ = endpoints(D1)
    D2₋, D2₊ = endpoints(D2)
    Δ1, Δ2 = Δ = (D1₊-D1₋)/n1, (D2₊-D2₋)/n2

    ts = [[t1, t2] for t1 in range(D1₋, D1₊, length = n1 +1), t2 in range(D2₋, D2₊, length = n2 +1)]
    tc = [mean([ts[i1,i2], ts[i1+1,i2], ts[i1,i2+1], ts[i1+1,i2+1]]) for i1 in 1:n1, i2 in 1:n2]

    𝒑s = 𝒑.(ts)
    𝒑c = 𝒑.(tc)

    if smooth
        𝒑′(t)=ForwardDiff.jacobian(𝒑,t)
        𝒆(t) = normalize(cross(𝒑′(t)[1:3,1],𝒑′(t)[1:3,2]))

        𝒆s = 𝒆.(ts)
        𝒆c = 𝒆.(tc)
    end

    Cs = rgb.(ts)
    Cc = rgb.(tc)
    Ts = ft.(ts)
    Tc = ft.(tc)

    Ns(i1, i2) = i1 + (n1+1) * (i2-1)
    Nc(i1, i2) = i1 + n1 * (i2-1)

    F1 = [[Ns(i1,i2)-1, Ns(i1+1,i2)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F2 = [[Ns(i1,i2)-1, Ns(i1,i2+1)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F3 = [[Ns(i1+1,i2+1)-1, Ns(i1+1,i2)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F4 = [[Ns(i1+1,i2+1)-1, Ns(i1,i2+1)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]

    np = (n1+1)*(n2+1) + n1*n2
    nf = 4*n1*n2

    script = ""
    script *= "mesh2{\n" * "  "^(preindent)
    script *= "  vertex_vectors{\n" * "  "^(preindent)
    script *= "    " * povray_script(np) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([𝒑s...]) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([𝒑c...]) * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)

    if smooth
        script *= "  normal_vectors{\n" * "  "^(preindent)
        script *= "    " * povray_script(np) * ", \n" * "  "^(preindent)
        script *= "    " * povray_script([𝒆s...]) * ", \n" * "  "^(preindent)
        script *= "    " * povray_script([𝒆c...]) * "\n" * "  "^(preindent)
        script *= "  }\n" * "  "^(preindent)
    end

    script *= "  texture_list{\n" * "  "^(preindent)
    script *= "    " * povray_script(np) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([Cs...],[Ts...]) * ", \n" * "  "^(preindent)
    script *= "    " * povray_script([Cc...],[Tc...]) * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)

    F1_ = [povray_script(f)*","*string(f[end]) for f in F1]
    F2_ = [povray_script(f)*","*string(f[end]) for f in F2]
    F3_ = [povray_script(f)*","*string(f[end]) for f in F3]
    F4_ = [povray_script(f)*","*string(f[end]) for f in F4]

    script *= "  face_indices{\n" * "  "^(preindent)
    script *= "    " * povray_script(nf) * ", \n" * "  "^(preindent)
    script *= "    " * join(F1_,",") * "\n" * "  "^(preindent)
    script *= "    " * join(F2_,",") * "\n" * "  "^(preindent)
    script *= "    " * join(F3_,",") * "\n" * "  "^(preindent)
    script *= "    " * join(F4_,",") * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)
    script *= "}"

    return script
end
