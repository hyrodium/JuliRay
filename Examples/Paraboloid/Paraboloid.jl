push!(LOAD_PATH, "Modules")
using Colors
using IntervalSets
using Revise
using JuliRay
Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

N=60
for i in 1
    𝒑(u) = [u...,u[1]^2+u[2]^2]
    D1 = -1..1
    D2 = -1..1
    mesh = (m,m)
    obj = rgbColor(ParametricSurface(𝒑, D1, D2, mesh=mesh, smooth=false), RGB(0.1,1,0.1))
    lng = 30°
    lat = 20°
    lookat = [0,0,0]
    zoom = 0.3
    render(obj,camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat),name="Paraboloid",index=i)
end
