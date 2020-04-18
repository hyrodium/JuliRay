push!(LOAD_PATH, "Modules")
using Colors
using Revise
using IntervalSets
using JuliRay
Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

begin
    𝒑(u) = [u...,u[1]^2+u[2]^2]
    D1 = -1..1
    D2 = -1..1
    mesh = (50,50)
    obj = rgbColor(ParametricSurface(𝒑, D1, D2, mesh=mesh), RGB(0.8,1,0.5))
    lng = 30°
    lat = 20°
    lookat = [0,0,1]
    zoom = 0.3
    render(obj,camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat),name="Paraboloid",index=0)
end
