push!(LOAD_PATH, "Modules")
using Colors
using ColorVectorSpace
using IntervalSets
using Revise
using JuliRay
using BasicBSpline
Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

function Smooth(x::Real)
    f(x) = ifelse(x>0,exp(-1/x),0)
    return f(2x/√3)/(f(2x/√3)+f(2/√3-2x/√3))
end
function Smooth(a::Real, b::Real, x::Real)
    return Smooth((x-a)/(b-a))
end
function Smooth(dict::Dict{A,B} where A <: Real where B <: Any, t::Real)
    T=sort(collect(keys(dict)))
    L=length(T)
    V=[dict[T[i]] for i in 1:L]

    if t ≤ T[1]
        return V[1]
    end
    for i in 2:L
        if t ≤ T[i]
            return V[i-1]+(V[i]-V[i-1])*Smooth(T[i-1],T[i],t)
        end
    end
    return V[end]
end

_p = 3
_k = Knots(1:13)
_P = BSplineSpace(_p,_k)
_a = [[[rand()-1/2,rand()-1/2,rand()-1/2] for _ in 1:dim(_P)-_p] for i in 1:4, j in 1:4]
_a = [vcat(_a[i,j],_a[i,j][1:_p]) for i in 1:4, j in 1:4]
_M = [BSplineManifold([_P], _a[i,j]) for i in 1:4, j in 1:4]

N=20
for i in 0:N-1
    t = _k[1+_p] + (_k[end-_p] - _k[1+_p]) * i/N
    k = 4*Knots([0,1])
    p = 3
    P = BSplineSpace(p,k)
    a0 = [[i-2.5,j-2.5,0]+mapping(_M[i,j],[t]) for i in 1:4, j in 1:4]
    a2 = [[0,0,Float64((i∈2:3) & (j∈2:3))] for i in 1:4, j in 1:4]
    a = a0+0*a2
    M = BSplineManifold([P,P],a)

    𝒑(u) = mapping(M,u)
    D1 = 0..1
    D2 = 0..1
    mesh = (40,40)
    nn = 5 # checker mesh
    rgb(u) = Float64(xor(sin(2π*nn*u[1])<0,sin(2π*nn*u[2])>0))*RGB(1,0.6,0.2) + Float64(xor(sin(2π*nn*u[1])>0,sin(2π*nn*u[2])>0))*RGB(1,0.2,0.2)
    ft(u) = FT(0.1,0.1)
    surface = ColoredParametricSurface(𝒑, D1, D2, mesh=mesh, rgb=rgb, ft=ft, smooth = false)

    ctrlpts = rgbColor(csgUnion([Sphere(p,0.05) for p in [a...]]),RGB(0.02,0.02,0.02))
    ctrlsgts1 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[2,1],0.02),Cylinder(a[2,1],a[3,1],0.02),Cylinder(a[3,1],a[4,1],0.02),
        Cylinder(a[1,2],a[2,2],0.02),Cylinder(a[2,2],a[3,2],0.02),Cylinder(a[3,2],a[4,2],0.02),
        Cylinder(a[1,3],a[2,3],0.02),Cylinder(a[2,3],a[3,3],0.02),Cylinder(a[3,3],a[4,3],0.02),
        Cylinder(a[1,4],a[2,4],0.02),Cylinder(a[2,4],a[3,4],0.02),Cylinder(a[3,4],a[4,4],0.02),
        ]), RGB(0.1,0.1,0.1))
    ctrlsgts2 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[1,2],0.02),Cylinder(a[1,2],a[1,3],0.02),Cylinder(a[1,3],a[1,4],0.02),
        Cylinder(a[2,1],a[2,2],0.02),Cylinder(a[2,2],a[2,3],0.02),Cylinder(a[2,3],a[2,4],0.02),
        Cylinder(a[3,1],a[3,2],0.02),Cylinder(a[3,2],a[3,3],0.02),Cylinder(a[3,3],a[3,4],0.02),
        Cylinder(a[4,1],a[4,2],0.02),Cylinder(a[4,2],a[4,3],0.02),Cylinder(a[4,3],a[4,4],0.02),
        ]), RGB(0.1,0.1,0.1))
    ctrlsgts = csgUnion(ctrlsgts1, ctrlsgts2)

    lng = 20°
    lat = 30°
    lookat = [0,0,0.25]*0
    zoom = 0.25
    render(csgUnion(surface,ctrlpts,ctrlsgts),camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat, width=500,height=300),name="DancingBezierSurface",index=i+1)
end

T = 0.2

N=40
for i in 0:N
    t = i/N
    tt = _k[1+_p] + (_k[end-_p] - _k[1+_p]) * T
    k = 4*Knots([0,1])
    p = 3
    P = BSplineSpace(p,k)
    a0 = [[i-2.5,j-2.5,0]+0.6*mapping(_M[i,j],[tt]) for i in 1:4, j in 1:4]
    a2 = [[0,0,Float64((i∈2:3) & (j∈2:3))] for i in 1:4, j in 1:4]
    a = a0+1.5*a2
    M = BSplineManifold([P,P],a)

    𝒑(u) = mapping(M,u)
    D1 = 0..1
    D2 = 0..1
    mesh = (40,40)
    nn = 5 # checker mesh
    rgb(u) = Float64(xor(sin(2π*nn*u[1])<0,sin(2π*nn*u[2])>0))*RGB(0.3,0.3,0.3) + Float64(xor(sin(2π*nn*u[1])>0,sin(2π*nn*u[2])>0))*RGB(0.5,0.5,0.5)
    ft(u) = FT(0.1,0.1)
    surface = ColoredParametricSurface(𝒑, D1, D2, mesh=mesh, rgb=rgb, ft=ft, smooth = false)

    ctrlpts = rgbColor(csgUnion([Sphere(p,0.05) for p in [a...]]),RGB(0,0.01,0.02))
    ctrlsgts1 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[2,1],0.02),Cylinder(a[2,1],a[3,1],0.02),Cylinder(a[3,1],a[4,1],0.02),
        Cylinder(a[1,2],a[2,2],0.02),Cylinder(a[2,2],a[3,2],0.02),Cylinder(a[3,2],a[4,2],0.02),
        Cylinder(a[1,3],a[2,3],0.02),Cylinder(a[2,3],a[3,3],0.02),Cylinder(a[3,3],a[4,3],0.02),
        Cylinder(a[1,4],a[2,4],0.02),Cylinder(a[2,4],a[3,4],0.02),Cylinder(a[3,4],a[4,4],0.02),
        ]), RGB(0,0.1,0.2))
    ctrlsgts2 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[1,2],0.02),Cylinder(a[1,2],a[1,3],0.02),Cylinder(a[1,3],a[1,4],0.02),
        Cylinder(a[2,1],a[2,2],0.02),Cylinder(a[2,2],a[2,3],0.02),Cylinder(a[2,3],a[2,4],0.02),
        Cylinder(a[3,1],a[3,2],0.02),Cylinder(a[3,2],a[3,3],0.02),Cylinder(a[3,3],a[3,4],0.02),
        Cylinder(a[4,1],a[4,2],0.02),Cylinder(a[4,2],a[4,3],0.02),Cylinder(a[4,3],a[4,4],0.02),
        ]), RGB(0,0.1,0.2))
    ctrlsgts = csgUnion(ctrlsgts1)

    _k_ = 4*Knots(0,1)
    _p_ = 3
    _P_ = BSplineSpace(_p_, _k_)
    _M1_ = [BSplineManifold([_P_], [a[1,i],a[2,i],a[3,i],a[4,i]]) for i in 1:4]
    c1 = rgbColor(csgUnion(ParametricCurve(t->mapping(_M1_[i],[t]),0..1,0.02) for i in 1:4), RGB(0,0.8,1))

    la1 = [mapping(_M1_[i],[t]) for i in 1:4]
    lM1 = BSplineManifold([_P_], la1)
    ca1 = rgbColor(csgUnion([Sphere(p,0.05) for p in [la1...]]),RGB(0.02,0.01,0))
    ctr1 = rgbColor(csgUnion([
        Cylinder(la1[1],la1[2],0.02),Cylinder(la1[2],la1[3],0.02),Cylinder(la1[3],la1[4],0.02),
        ]), RGB(0.2,0.1,0))
    lc1 = rgbColor(ParametricCurve(t->mapping(lM1,[t]),0..1,0.02), RGB(1,0.5,0.0))


    lng = 20°
    lat = 30°
    lookat = [0,0,0.1]
    zoom = 0.25
    render(csgUnion(surface,ctrlpts,ctrlsgts, c1, lc1, ca1, ctr1),camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat, width=600,height=400),name="BezierSurfaceLocus1",index=i+1)
end


for i in 0:N
    t = i/N
    tt = _k[1+_p] + (_k[end-_p] - _k[1+_p]) * T
    k = 4*Knots([0,1])
    p = 3
    P = BSplineSpace(p,k)
    a0 = [[i-2.5,j-2.5,0]+0.6*mapping(_M[i,j],[tt]) for i in 1:4, j in 1:4]
    a2 = [[0,0,Float64((i∈2:3) & (j∈2:3))] for i in 1:4, j in 1:4]
    a = a0+1.5*a2
    M = BSplineManifold([P,P],a)

    𝒑(u) = mapping(M,u)
    D1 = 0..1
    D2 = 0..1
    mesh = (40,40)
    nn = 5 # checker mesh
    rgb(u) = Float64(xor(sin(2π*nn*u[1])<0,sin(2π*nn*u[2])>0))*RGB(0.3,0.3,0.3) + Float64(xor(sin(2π*nn*u[1])>0,sin(2π*nn*u[2])>0))*RGB(0.5,0.5,0.5)
    ft(u) = FT(0.1,0.1)
    surface = ColoredParametricSurface(𝒑, D1, D2, mesh=mesh, rgb=rgb, ft=ft, smooth = false)

    ctrlpts = rgbColor(csgUnion([Sphere(p,0.05) for p in [a...]]),RGB(0,0.01,0.02))
    ctrlsgts1 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[2,1],0.02),Cylinder(a[2,1],a[3,1],0.02),Cylinder(a[3,1],a[4,1],0.02),
        Cylinder(a[1,2],a[2,2],0.02),Cylinder(a[2,2],a[3,2],0.02),Cylinder(a[3,2],a[4,2],0.02),
        Cylinder(a[1,3],a[2,3],0.02),Cylinder(a[2,3],a[3,3],0.02),Cylinder(a[3,3],a[4,3],0.02),
        Cylinder(a[1,4],a[2,4],0.02),Cylinder(a[2,4],a[3,4],0.02),Cylinder(a[3,4],a[4,4],0.02),
        ]), RGB(0,0.1,0.2))
    ctrlsgts2 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[1,2],0.02),Cylinder(a[1,2],a[1,3],0.02),Cylinder(a[1,3],a[1,4],0.02),
        Cylinder(a[2,1],a[2,2],0.02),Cylinder(a[2,2],a[2,3],0.02),Cylinder(a[2,3],a[2,4],0.02),
        Cylinder(a[3,1],a[3,2],0.02),Cylinder(a[3,2],a[3,3],0.02),Cylinder(a[3,3],a[3,4],0.02),
        Cylinder(a[4,1],a[4,2],0.02),Cylinder(a[4,2],a[4,3],0.02),Cylinder(a[4,3],a[4,4],0.02),
        ]), RGB(0,0.1,0.2))
    ctrlsgts = csgUnion(ctrlsgts2)

    _k_ = 4*Knots(0,1)
    _p_ = 3
    _P_ = BSplineSpace(_p_, _k_)
    _M1_ = [BSplineManifold([_P_], [a[i,1],a[i,2],a[i,3],a[i,4]]) for i in 1:4]
    c1 = rgbColor(csgUnion(ParametricCurve(t->mapping(_M1_[i],[t]),0..1,0.02) for i in 1:4), RGB(0,0.8,1))

    la1 = [mapping(_M1_[i],[t]) for i in 1:4]
    lM1 = BSplineManifold([_P_], la1)
    ca1 = rgbColor(csgUnion([Sphere(p,0.05) for p in [la1...]]),RGB(0.02,0.01,0))
    ctr1 = rgbColor(csgUnion([
        Cylinder(la1[1],la1[2],0.02),Cylinder(la1[2],la1[3],0.02),Cylinder(la1[3],la1[4],0.02),
        ]), RGB(0.2,0.1,0))
    lc1 = rgbColor(ParametricCurve(t->mapping(lM1,[t]),0..1,0.02), RGB(1,0.5,0.0))


    lng = 20°
    lat = 30°
    lookat = [0,0,0.1]
    zoom = 0.25
    render(csgUnion(surface,ctrlpts,ctrlsgts, c1, lc1, ca1, ctr1),camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat, width=600,height=400),name="BezierSurfaceLocus2",index=i+1)
end



for i in 0:0
    tt = _k[1+_p] + (_k[end-_p] - _k[1+_p]) * T
    k = 4*Knots([0,1])
    p = 3
    P = BSplineSpace(p,k)
    a0 = [[i-2.5,j-2.5,0]+0.6*mapping(_M[i,j],[tt]) for i in 1:4, j in 1:4]
    a2 = [[0,0,Float64((i∈2:3) & (j∈2:3))] for i in 1:4, j in 1:4]
    a = a0+1.5*a2
    M = BSplineManifold([P,P],a)

    𝒑(u) = mapping(M,u)
    D1 = 0..1
    D2 = 0..1
    mesh = (40,40)
    nn = 5 # checker mesh
    rgb(u) = Float64(xor(sin(2π*nn*u[1])<0,sin(2π*nn*u[2])>0))*RGB(0.3,0.3,0.3) + Float64(xor(sin(2π*nn*u[1])>0,sin(2π*nn*u[2])>0))*RGB(0.5,0.5,0.5)
    ft(u) = FT(0.1,0.1)
    surface = ColoredParametricSurface(𝒑, D1, D2, mesh=mesh, rgb=rgb, ft=ft, smooth = false)

    ctrlpts = rgbColor(csgUnion([Sphere(p,0.05) for p in [a...]]),RGB(0,0.01,0.02))
    ctrlsgts1 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[2,1],0.02),Cylinder(a[2,1],a[3,1],0.02),Cylinder(a[3,1],a[4,1],0.02),
        Cylinder(a[1,2],a[2,2],0.02),Cylinder(a[2,2],a[3,2],0.02),Cylinder(a[3,2],a[4,2],0.02),
        Cylinder(a[1,3],a[2,3],0.02),Cylinder(a[2,3],a[3,3],0.02),Cylinder(a[3,3],a[4,3],0.02),
        Cylinder(a[1,4],a[2,4],0.02),Cylinder(a[2,4],a[3,4],0.02),Cylinder(a[3,4],a[4,4],0.02),
        ]), RGB(0,0.1,0.2))
    ctrlsgts2 = rgbColor(csgUnion([
        Cylinder(a[1,1],a[1,2],0.02),Cylinder(a[1,2],a[1,3],0.02),Cylinder(a[1,3],a[1,4],0.02),
        Cylinder(a[2,1],a[2,2],0.02),Cylinder(a[2,2],a[2,3],0.02),Cylinder(a[2,3],a[2,4],0.02),
        Cylinder(a[3,1],a[3,2],0.02),Cylinder(a[3,2],a[3,3],0.02),Cylinder(a[3,3],a[3,4],0.02),
        Cylinder(a[4,1],a[4,2],0.02),Cylinder(a[4,2],a[4,3],0.02),Cylinder(a[4,3],a[4,4],0.02),
        ]), RGB(0,0.1,0.2))
    ctrlsgts = csgUnion(ctrlsgts1, ctrlsgts2)

    lng = 20°
    lat = 30°
    lookat = [0,0,0.1]
    zoom = 0.25
    render(csgUnion(surface,ctrlpts,ctrlsgts),camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat, width=600,height=400),name="BezierSurfaceLocus")
end
