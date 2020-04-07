push!(LOAD_PATH, "Modules")
using LinearAlgebra
using Colors
using JuliRay

Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

function Plane2Sphere(p::Array{T,1}) where T<:Real
    if length(p)≠2
        error("No Point on ℝ²")
    end
    return [2p[1]/(1+p[1]^2+p[2]^2),2p[2]/(1+p[1]^2+p[2]^2),(-1+p[1]^2+p[2]^2)/(1+p[1]^2+p[2]^2)]
end

function Plane2Sphere(θ::Real,p::Array{T,1}) where T<:Real
    if length(p)≠2
        error("No Point on ℝ²")
    end
    c=[0,0,tan(θ)]
    r=1/cos(θ)
    P=p*cot(θ/2+π/4)
    return Plane2Sphere(P)*r+c
end

function POLYGON(vertices;r1=1.0,r2=0.05,r3=0.025)
    n=length(vertices)
    V=rgbColor(csgUnion((p->Sphere(p,r2)).(vertices)),RGB(0.1,0.1,0.1))
    E=rgbColor(csgUnion([Cylinder(vertices[i],vertices[mod(i,n)+1],r3) for i in 1:n]),RGB(0.1,0.1,0.1))
    F=rgbftColor(Polygon(vertices),RGB(0.5,1,1),FT(0.1,0.1))
    return csgUnion(V,E,F)
end

function Mirror(q,p1,p2,p3)
    e₃=JuliRay.NormalVector(p1,p2,p3)
    e₁=JuliRay.OrthogonalVector(e₃)
    e₂=cross(e₃,e₁)
    A=hcat(e₁,e₂,e₃)*transpose(hcat(e₁,e₂,-e₃))
    fixedpoint=(p1+p2+p3)/3
    b=fixedpoint-A*fixedpoint
    return A*q+b
end

function Put(F,i,O)
    p1=O
    n=length(F)
    p2,p3=F[mod(i-1,n)+1],F[mod(i,n)+1]
    return (p->Mirror(p,p1,p2,p3)).(F)
end

function f(x)
    ifelse(x>0,exp(-1/x),0)
end
function g(x)
    f(2x/√3)/(f(2x/√3)+f(2/√3-2x/√3))
end
function Smooth(a,b,x)
    g((x-a)/(b-a))
end

## Cube

h=1/√3 # 単位球面に内接させた際の内接球半径
n=4 # n角形
m=6 # m面体
M=80 # アニメーションの刻み数

r1=√(1-h^2) # n角形の外接円半径
for i ∈ 1:M
    t=3(i-1)/M
    θ=(1-Smooth(0.5,1.4,t)+Smooth(1.6,2.5,t))*π/2 # パラメータθは0 → π/2 → 0の順で滑らかに変化
    R=1/cos(θ)
    O=[0,0,tan(θ)]
    N=O+[0,0,R]
    H=O-[0,0,√(R^2-r1^2)]
    F=[[[r1*cos(2π*i/n),r1*sin(2π*i/n),0.0]+H for i in 1:n] for i ∈ 1:m] # 最初のn角形の頂点
    NET(F)=Scaling(Rotate(csgUnion([POLYGON(F[i],r2=0.02,r3=0.02) for i ∈ 1:m]),[0,0,1],45°),1.5/√3)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[1],2,O)
    F[4]=Put(F[1],3,O)
    F[5]=Put(F[1],4,O)
    F[6]=Put(F[4],1,O)
    object01=NET(F)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[2],2,O)
    F[5]=Put(F[2],4,O)
    F[6]=Put(F[3],1,O)
    object02=NET(F)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[1],2,O)
    F[5]=Put(F[2],4,O)
    F[6]=Put(F[3],1,O)
    object03=NET(F)

    F[2]=Put(F[1],4,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[3],2,O)
    F[5]=Put(F[2],1,O)
    F[6]=Put(F[4],1,O)
    object04=NET(F)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[3],2,O)
    F[5]=Put(F[1],4,O)
    F[6]=Put(F[3],1,O)
    object05=NET(F)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[3],1,O)
    F[5]=Put(F[2],2,O)
    F[6]=Put(F[4],4,O)
    object06=NET(F)

    F[2]=Put(F[1],4,O)
    F[3]=Put(F[1],2,O)
    F[4]=Put(F[1],3,O)
    F[5]=Put(F[2],1,O)
    F[6]=Put(F[4],1,O)
    object07=NET(F)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[3],1,O)
    F[5]=Put(F[2],4,O)
    F[6]=Put(F[3],2,O)
    object08=NET(F)

    F[2]=Put(F[1],4,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[3],1,O)
    F[5]=Put(F[2],1,O)
    F[6]=Put(F[3],2,O)
    object09=NET(F)

    F[2]=Put(F[1],4,O)
    F[3]=Put(F[1],3,O)
    F[4]=Put(F[3],1,O)
    F[5]=Put(F[2],1,O)
    F[6]=Put(F[4],2,O)
    object10=NET(F)

    F[2]=Put(F[1],1,O)
    F[3]=Put(F[2],3,O)
    F[4]=Put(F[1],2,O)
    F[5]=Put(F[4],3,O)
    F[6]=Put(F[5],1,O)
    object11=NET(F)

    object=csgUnion(
        ParallelTranslation(object01,[-8,4,0]),
        ParallelTranslation(object02,[-3,4,0]),
        ParallelTranslation(object03,[2,4,0]),
        ParallelTranslation(object04,[7,4,0]),
        ParallelTranslation(object05,[-8,0,0]),
        ParallelTranslation(object06,[-3,0,0]),
        ParallelTranslation(object07,[2,0,0]),
        ParallelTranslation(object08,[7,0,0]),
        ParallelTranslation(object09,[-8,-4,0]),
        ParallelTranslation(object10,[-3,-4,0]),
        ParallelTranslation(object11,[5,-3.5,0])
    )
    lng = -90°-15°*(Smooth(0.1,0.5,t)-Smooth(2.5,2.9,t))
    lat = 90°-50°*(Smooth(0.1,0.5,t)-Smooth(2.5,2.9,t))
    camera = LngLatCamera(lng=lng,lat=lat,pers=0.2,zoom=0.065,width=640*2,height=360*2)
    render(object, camera=camera, name="NetOfCubes", index=i)
end
