push!(LOAD_PATH, "Modules")
using Colors
using IntervalSets
using Revise
using JuliRay
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


N=80
for i in 1:N
    t = π*i/N
    catenoid(u) = [cos(u[2])*cosh(u[1]),sin(u[2])*cosh(u[1]),u[1]]
    helicoid(u) = [cos(u[2])*sinh(u[1]),sin(u[2])*sinh(u[1]),u[2]]

    # 𝒑(u) = cos(t)*catenoid(u) + sin(t)*helicoid
    𝒑(u) = cos(t)*catenoid(u)+sin(t)*helicoid(u)
    D1 = -π/2..π/2
    D2 = -π..π
    mesh = (40,40)
    obj = rgbColor(ParametricSurface(𝒑, D1, D2, mesh=mesh), RGB(0.1,1,0.1))
    lng = 30°+90°*i/60
    lat = 20°
    lookat = [0,0,1]
    zoom = 0.1
    render(obj,camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat),name="CatenoidHelicoid",index=i)
end
