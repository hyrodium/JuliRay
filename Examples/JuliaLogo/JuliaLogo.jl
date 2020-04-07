push!(LOAD_PATH, "Modules")
using Colors
using JuliRay

a=rgbColor(Sphere([2,-1,-1]/3,0.6),RGB(1,0,0))
b=rgbColor(Sphere([-1,2,-1]/3,0.6),RGB(0.5,0,0.5))
c=rgbColor(Sphere([-1,-1,2]/3,0.6),RGB(0,1,0))
abc=csgMerge(a,b,c)
render(Transparent(abc,FT(0.1,0.2)),name="julialogo",index=0,camera=LngLatCamera(lng=π/5,lat=π/5,pers=0.5,zoom=0.3,height=1000,width=1000))
