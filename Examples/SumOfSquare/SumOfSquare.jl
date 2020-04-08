push!(LOAD_PATH, "Modules")
using Colors
using Revise
using JuliRay
Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

function f(x)
    ifelse(x>0,exp(-1/x),0)
end
function Smooth(x::Real)
    f(2x/√3)/(f(2x/√3)+f(2/√3-2x/√3))
end
function Smooth(a::Real, b::Real, x::Real)
    Smooth((x-a)/(b-a))
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

n=5
colors=[
    RGB(0.8,0.05,0.15), # Red
    RGB(0.8,0.3,0), # Orange
    RGB(0.8,0.7,0), # Yellow
    RGB(0,0.8,0.3), # Mint Green
    RGB(0,0.6,0.8)  # Cyan Blue
    ]

# S(i)=rgbColor(Box([0,0,n-i].+0.03,[i,i,n-i+1].-0.03), colors[i])
S(i)=rgbColor(RoundedBox([0,0,n-i].+0.001,[i,i,n-i+1].-0.001, 0.1), colors[i])

Ss1(t)=csgUnion(S(i) for i in 1:n)
Ss2(t)=ParallelTranslation(Rotate(Ss1(t),[1,1,-1], 2*π/3, fixedpoint=[1,1,n]), [0,Smooth(Dict(0=>2n, 1=>0),t),Smooth(Dict(0=>1.5n, 1=>-1),t)])
Ss3(t)=ParallelTranslation(Rotate(Ss1(t),[1,1,-1], -2*π/3, fixedpoint=[1,1,n]), [-1,0,Smooth(Dict(1=>2n,2=>0),t)])
SSS(t)=ParallelTranslation(
    Rotate(
        csgUnion(Ss1(t), Ss2(t), Ss3(t)),[0,1,0], Smooth(Dict(2=>0°, 4=>-90°),t), fixedpoint=[n,n,n]/2
    ),  Smooth(Dict(2=>[0,0,0], 4=>[4,-4,0], 5=>[-3,-3,0]),t)
)

SSS′(t)=Rotate(
    ParallelTranslation(
        Rotate(
            ParallelTranslation(
                Scaling(SSS(4),-1), [0,0,n]
            ), [0,0,1], -90°
        ), [-1,1,0]*Smooth(Dict(2=>3n, 4=>0),t)+[-1,0,0]
    ), [0,0,1], Smooth(Dict(4=>0, 5=>90°),t)
)

N = 20
for i ∈ 1:5N
    t=i/N
    lng=Smooth(Dict(0=>60°, 1=>30°, 2=>-20°, 4=>-150°, 5=>-120°),t)
    lat=20°
    lookat=[n,n,n]/2+Smooth(Dict(4=>[0,0,0], 5=>[-5.5,0.5-2,0]),t)
    zoom=Smooth(Dict(2=>0.1, 4=>0.055, 5=>0.07),t)
    render(csgUnion(SSS(t),SSS′(t)),camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat),name="SumOfSquare",index=i)
end

for i ∈ 5N:7N
    t=i/N
    s1=ParallelTranslation(Rotate(Ss1(0), [0,0,1], 180°),[-3,3,0])
    s2=ParallelTranslation(Rotate(s1,[1,1,1],120°, fixedpoint=[-4,2,5]),[0,0,-1])
    s3=ParallelTranslation(Rotate(s1,[1,1,1],-120°, fixedpoint=[-4,2,5]),[-1,0,-1])
    s4=Scaling(s1,-1,fixedpoint=[-3.5, 0, 2.5])
    s5=Scaling(s2,-1,fixedpoint=[-3.5, 0, 2.5])
    s6=Scaling(s3,-1,fixedpoint=[-3.5, 0, 2.5])

    t1=s1
    t2=Transparent(s2,FT(0,Smooth(Dict(6.0=>0, 6.5=>1),t)))
    t3=Transparent(s3,FT(0,Smooth(Dict(5.0=>0, 5.5=>1),t)))
    t4=Transparent(s4,FT(0,Smooth(Dict(5.5=>0, 6.0=>1),t)))
    t5=Transparent(s5,FT(0,Smooth(Dict(6.0=>0, 6.5=>1),t)))
    t6=Transparent(s6,FT(0,Smooth(Dict(6.5=>0, 7.0=>1),t)))
    obj=csgUnion(t1,t2,t3,t4,t5,t6)
    lng=Smooth(Dict(0=>60°, 1=>30°, 2=>-20°, 4=>-150°, 5=>-120°, 7=>-120°),t)
    lat=20°
    lookat=[n,n,n]/2+Smooth(Dict(4=>[0,0,0], 5=>[-5.5,0.5-2,0], 7=>[-8,-2,0]),t)
    zoom=Smooth(Dict(2=>0.1, 4=>0.055, 5=>0.07, 7=>0.1),t)
    render(obj,camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat),name="SumOfSquare",index=i)
end
