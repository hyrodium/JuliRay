push!(LOAD_PATH, "Modules")
using Colors
using JuliRay

φ=MathConstants.φ

Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

function POLYGON(vertices;r1=1.0,r2=0.05,r3=0.025)
    n=length(vertices)
    V=rgbColor(csgUnion((p->Sphere(p,r2)).(vertices)),RGB(0.1,0.1,0.1))
    E=rgbColor(csgUnion([Cylinder(vertices[i],vertices[mod(i,n)+1],r3) for i in 1:n]),RGB(0.1,0.1,0.1))
    F=rgbftColor(Polygon(vertices),RGB(0.5,1,1),FT(0.1,0.1))
    return csgUnion(V,E,F)
end

function Mirror(q,p1,p2,p3)
    e₃=NormalVector(p1,p2,p3)
    e₁=OrthogonalVector(e₃)
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

## Tetrahedron

h=1/3 # 単位球面に内接させた際の内接球半径
n=3 # n角形
m=4 # m面体
M=60 # アニメーションの刻み数

r1=√(1-h^2) # n角形の外接円半径
for i ∈ 1:M
    θ=(1-Smooth(0,1,2(i-1)/M)+Smooth(1,2,2(i-1)/M))*π/2 # パラメータθは0→1→0の順で滑らかに変化
    R=1/cos(θ)
    O=[0,0,tan(θ)]
    N=O+[0,0,R]
    H=O-[0,0,√(R^2-r1^2)]

    F=[[[r1*cos(2π*i/n),r1*sin(2π*i/n),0.0]+H for i in 1:n] for i ∈ 1:m] # 最初のn角形の頂点
    F[2]=Put(F[1],1,O) # 以下はn角形を付け加える操作
    F[3]=Put(F[2],2,O)
    F[4]=Put(F[2],3,O)

    object=csgUnion([POLYGON(F[i],r2=0.03,r3=0.01) for i ∈ 1:m]) # 頂点座標をもとにオブジェクトを生成
    render(object,camera=LngLatCamera(lng=30°,lat=30°,pers=0.2,zoom=0.2,width=640,height=360),name="Np"*string(m),index=i)
end
