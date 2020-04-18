using IntervalSets

using Random
using ForwardDiff
using Statistics


export ParametricSurface
struct ParametricSurface <: PrimitiveObject
    parametrization :: Function
    D1 :: Interval{:closed,:closed,Float64}
    D2 :: Interval{:closed,:closed,Float64}
    mesh :: Tuple{Int,Int}
    function ParametricSurface(parametrization::Function, D1::Interval, D2::Interval; mesh::Tuple{Int,Int}=(10,10))
        new(parametrization, ClosedInterval(D1), ClosedInterval(D2), mesh)
    end
end

function povray_script(parametricsurface::ParametricSurface)
    name = randstring()
    name = "mesh2.inc"
    path = BASE_DIR * name
    𝒑 = parametricsurface.parametrization
    D1, D2 = parametricsurface.D1, parametricsurface.D2
    mesh = parametricsurface.mesh
    return mesh2(𝒑, D1, D2; mesh=mesh, smooth = true)
end

function mesh2(𝒑, D1, D2; mesh=(10,10), smooth=true, preindnet=0)
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
    # [ts...][Ns(8,11)] == ts[8,11]
    # [tc...][Nc(8,11)] == tc[8,11]

    F1 = [[Ns(i1,i2)-1, Ns(i1+1,i2)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F2 = [[Ns(i1,i2)-1, Ns(i1,i2+1)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F3 = [[Ns(i1+1,i2+1)-1, Ns(i1+1,i2)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]
    F4 = [[Ns(i1+1,i2+1)-1, Ns(i1,i2+1)-1, (n1+1)*(n2+1)+Nc(i1,i2)-1] for i1 in 1:n1, i2 in 1:n2]

    np = (n1+1)*(n2+1) + n1*n2
    nf = 4*n1*n2

    Ns = [1,2,3]

    mesh2 = ""
    mesh2 *= "mesh2{\n"
    mesh2 *= "  vertex_vectors{\n"
    mesh2 *= "    " * povray_script(np) * ", \n"
    mesh2 *= "    " * povray_script([𝒑s...]) * ", \n"
    mesh2 *= "    " * povray_script([𝒑c...]) * "\n"
    mesh2 *= "  }\n"
    mesh2 *= "  normal_vectors{\n"
    mesh2 *= "    " * povray_script(np) * ", \n"
    mesh2 *= "    " * povray_script([𝒆s...]) * ", \n"
    mesh2 *= "    " * povray_script([𝒆c...]) * "\n"
    mesh2 *= "  }\n"
    mesh2 *= "  face_indices{\n"
    mesh2 *= "    " * povray_script(nf) * ", \n"
    mesh2 *= "    " * povray_script([F1...]) * "\n"
    mesh2 *= "    " * povray_script([F2...]) * "\n"
    mesh2 *= "    " * povray_script([F3...]) * "\n"
    mesh2 *= "    " * povray_script([F4...]) * "\n"
    mesh2 *= "  }\n"
    mesh2 *= "}"

    return mesh2
end
