push!(LOAD_PATH, "Modules")
using Colors
using JuliRay

φ=MathConstants.φ

Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

function POLYGON2(vertices,O,R;r1=1.0,r2=0.05,r3=0.025)
    n=length(vertices)
    C=O+R*normalize(+(vertices...)/n-O)
    V=rgbColor(csgUnion((p->Sphere(p,r2)).(vertices)),RGB(0.1,0.1,0.1))
    E=rgbColor(csgUnion([Arc(vertices[i],O+R*normalize((vertices[i]+vertices[mod(i,n)+1])/2-O),vertices[mod(i,n)+1],r3) for i in 1:n]),RGB(0.1,0.1,0.1))
    if(R<10^2)
        F=csgClip(
            rgbftColor(
                Sphere(O,R),RGB(1,0.5,1),FT(0.1,0.1)
            ),
            csgIntersection(
                [
                    (v=JuliRay.NormalVector(O,vertices[i],vertices[mod(i,n)+1]);s=sign(dot(v,C-O));Cylinder(O,O+s*R*v,R))
                for i ∈ 1:n]
            )
        )
    else
        F=rgbftColor(Polygon(vertices),RGB(1,0.5,1),FT(0.1,0.1))
    end
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

    object=csgUnion([POLYGON2(F[i],O,R,r2=0.03,r3=0.01) for i ∈ 1:m]) # 頂点座標をもとにオブジェクトを生成
    render(csgUnion(object),camera=LngLatCamera(lng=30°,lat=30°,pers=0.2,zoom=0.2,width=640,height=360),name="Ns"*string(m),index=i)
end

## Cube
h=1/√3 # 単位球面に内接させた際の内接球半径
n=4 # n角形
m=6 # m面体
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
    F[3]=Put(F[2],3,O)
    F[4]=Put(F[1],2,O)
    F[5]=Put(F[4],3,O)
    F[6]=Put(F[5],1,O)

    object=csgUnion([POLYGON2(F[i],O,R,r2=0.03,r3=0.01) for i ∈ 1:m]) # 頂点座標をもとにオブジェクトを生成
    render(object,camera=LngLatCamera(lng=30°,lat=30°,pers=0.2,zoom=0.2,width=640,height=360),name="Ns"*string(m),index=i)
end

## Octahedron
h=1/√3 # 単位球面に内接させた際の内接球半径
n=3 # n角形
m=8 # m面体
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
    F[3]=Put(F[2],3,O)
    F[4]=Put(F[3],2,O)
    F[5]=Put(F[4],3,O)
    F[6]=Put(F[1],2,O)
    F[7]=Put(F[6],3,O)
    F[8]=Put(F[7],2,O)

    object=csgUnion([POLYGON2(F[i],O,R,r2=0.03,r3=0.01) for i ∈ 1:m]) # 頂点座標をもとにオブジェクトを生成
    render(object,camera=LngLatCamera(lng=30°,lat=30°,pers=0.2,zoom=0.2,width=640,height=360),name="Ns"*string(m),index=i)
end

## Dodecahedron
φ=2cos(π/5)
ξ=2sin(π/5)

h=φ/(√3*ξ) # 単位球面に内接させた際の内接球半径
n=5 # n角形
m=12 # m面体
M=60 # アニメーションの刻み数

r1=√(1-h^2) # n角形の外接円半径
for i ∈ 1:M
    θ=(1-Smooth(0,1,2(i-1)/M)+Smooth(1,2,2(i-1)/M))*π/2 # パラメータθは0→1→0の順で滑らかに変化
    R=1/cos(θ)
    O=[0,0,tan(θ)]
    N=O+[0,0,R]
    H=O-[0,0,√(R^2-r1^2)]

    F=[[[r1*cos(2π*i/n),r1*sin(2π*i/n),0.0]+H for i in 1:n] for i ∈ 1:m] # 最初のn角形の頂点
    F[2]=Put(F[1],4,O) # 以下はn角形を付け加える操作
    F[3]=Put(F[2],1,O)
    F[4]=Put(F[2],2,O)
    F[5]=Put(F[2],3,O)
    F[6]=Put(F[2],5,O)
    F[7]=Put(F[1],2,O)
    F[8]=Put(F[7],5,O)
    F[9]=Put(F[8],1,O)
    F[10]=Put(F[8],2,O)
    F[11]=Put(F[8],3,O)
    F[12]=Put(F[8],4,O)

    object=csgUnion([POLYGON2(F[i],O,R,r2=0.03,r3=0.01) for i ∈ 1:m]) # 頂点座標をもとにオブジェクトを生成
    render(object,camera=LngLatCamera(lng=30°,lat=30°,pers=0.2,zoom=0.2,width=640,height=360),name="Ns"*string(m),index=i)
end

## Icosahedron
φ=2cos(π/5)
ξ=2sin(π/5)

h=φ/(√3*ξ) # 単位球面に内接させた際の内接球半径
n=3 # n角形
m=20 # m面体
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
    F[3]=Put(F[2],3,O)
    F[4]=Put(F[3],2,O)
    F[5]=Put(F[4],1,O)
    F[6]=Put(F[5],3,O)
    F[7]=Put(F[1],2,O)
    F[8]=Put(F[7],3,O)
    F[9]=Put(F[8],1,O)
    F[10]=Put(F[9],2,O)
    F[11]=Put(F[1],3,O)
    F[12]=Put(F[2],2,O)
    F[13]=Put(F[3],1,O)
    F[14]=Put(F[4],3,O)
    F[15]=Put(F[5],2,O)
    F[16]=Put(F[6],1,O)
    F[17]=Put(F[7],1,O)
    F[18]=Put(F[8],2,O)
    F[19]=Put(F[9],3,O)
    F[20]=Put(F[10],1,O)

    object=csgUnion([POLYGON2(F[i],O,R,r2=0.03,r3=0.01) for i ∈ 1:m]) # 頂点座標をもとにオブジェクトを生成
    render(object,camera=LngLatCamera(lng=30°,lat=30°,pers=0.2,zoom=0.2,width=640,height=360),name="Ns"*string(m),index=i)
end
