push!(LOAD_PATH, "Modules")
using Colors
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
    RGB(0.8,0,0.2), # Red
    RGB(0.8,0.3,0), # Orange
    RGB(0.8,0.7,0), # Yellow
    RGB(0,0.8,0.3), # Mint Green
    RGB(0,0.6,0.8)  # Cyan Blue
    ]
S(i)=rgbColor(Box([0,0,n-i+1],[i,i,n-i]), colors[i])

Ss1(t)=csgUnion(S(i) for i in 1:n)
Ss2(t)=ParallelTranslation(Rotate(Ss1(t),[1,1,-1], 2*π/3, fixedpoint=[1,1,n]), [0,Smooth(Dict(0=>2n, 1=>0),t),Smooth(Dict(0=>1.5n, 1=>-1),t)])
Ss3(t)=ParallelTranslation(Rotate(Ss1(t),[1,1,-1], -2*π/3, fixedpoint=[1,1,n]), [-1,0,Smooth(Dict(1=>2n,2=>0),t)])
SSS(t)=ParallelTranslation(
    Rotate(
        csgUnion(Ss1(t), Ss2(t), Ss3(t)),[0,1,0], Smooth(Dict(2=>0°, 3=>-90°),t), fixedpoint=[n,n,n]/2
    ),  Smooth(Dict(2=>[0,0,0], 3=>[4,-4,0], 4=>[-3,-3,0]),t)
)

SSS′(t)=Rotate(
    ParallelTranslation(
        Rotate(
            ParallelTranslation(
                Scaling(SSS(3),-1), [0,0,n]
            ), [0,0,1], -90°
        ), [-1,1,0]*Smooth(Dict(2=>3n, 3=>0),t)+[-1,0,0]
    ), [0,0,1], Smooth(Dict(3=>0, 4=>90°),t)
)

N = 20
for i ∈ 0:4N
    t=i/N
    lng=Smooth(Dict(0=>60°, 1=>30°, 2=>-20°, 3=>-150°, 4=>-120°),t)
    lat=20°
    lookat=[n,n,n]/2+Smooth(Dict(3=>[0,0,0], 4=>[-5.5,0.5-2,0]),t)
    zoom=Smooth(Dict(2=>0.1, 3=>0.055, 4=>0.07),t)
    render(csgUnion(SSS(t),SSS′(t)),camera=LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=zoom,lookat=lookat),name="SumOfSquare",index=i)
end
