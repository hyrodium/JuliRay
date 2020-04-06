include("../JuliRay.jl")
Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))
φ=MathConstants.φ
include("../PlatonicSolids.jl")

##%
ε=10^(-1.2)

function f(x)
    ifelse(x>0,exp(-1/x),0)
end
function g(x)
    f(2x/√3)/(f(2x/√3)+f(2/√3-2x/√3))
end
function Smooth(a,b,x)
    g((x-a)/(b-a))
end

function ℝ³⭢S³(p::Array{T,1}) where T<:Real
    if length(p)≠3
        error("No Point on ℝ³")
    end
    return [2p[1]/(1+p[1]^2+p[2]^2+p[3]^2),2p[2]/(1+p[1]^2+p[2]^2+p[3]^2),2p[3]/(1+p[1]^2+p[2]^2+p[3]^2),(-1+p[1]^2+p[2]^2+p[3]^2)/(1+p[1]^2+p[2]^2+p[3]^2)]
end
function ℝ³⭢S³(p::Array{T,1},θ::Real) where T<:Real
    if length(p)≠3
        error("No Point on ℝ³")
    elseif θ<ε
        return [p[1:3]...,0.0]
    else
        c=[0,0,0,cot(θ)]
        r=1/sin(θ)
        P=p*cot((π/2-θ)/2+π/4)
        return ℝ³⭢S³(P)*r+c
    end
end
function S³⭢ℝ³(q::Array{T,1}) where T<:Real
    if length(q)≠4
        error("No Point on S³")
    elseif !(norm(q) ≈ 1.0)
        print(norm(q))
        error("No Point on unit S³, radius: $(norm(q))")
    end
    return [q[1]/(1-q[4]),q[2]/(1-q[4]),q[3]/(1-q[4])]
end
function S³⭢ℝ³(q::Array{T,1},θ::Real) where T<:Real
    if length(q)≠4
        error("No Point on S³")
    elseif θ<ε
        return q[1:3]
    else
        c=[0,0,0,cot(θ)]
        r=1/sin(θ)
        return S³⭢ℝ³((q-c)/r)*tan((π/2-θ)/2+π/4)
    end
end
function ℝ⁴⭢S³(p::Array{T,1},θ::Real) where T<:Real
    if θ<ε
        return [p[1:3]...,0.0]
    else
        R=1/sin(θ)
        O=[0,0,0,cot(θ)]
        return O+R*normalize(p-O)
    end
end

function NormalVector(p₁::RealVector,p₂::RealVector,p₃::RealVector,p₄::RealVector)
    A=hcat(p₁-p₄,p₂-p₄,p₃-p₄)
    return normalize([(-1)^i*det(A[deleteat!(collect(1:4),i),:]) for i ∈ 1:4])
end

function Mirror(q::RealVector,p₁::RealVector,p₂::RealVector,p₃::RealVector)
    𝒏=NormalVector(p₁,p₂,p₃)
    return q-2*dot(𝒏,q-p₁)*𝒏
end
function Mirror(q::RealVector,p₁::RealVector,p₂::RealVector,p₃::RealVector,p₄::RealVector)
    𝒏=NormalVector(p₁,p₂,p₃,p₄)
    return q-2*dot(𝒏,q-p₁)*𝒏
end
function Mirror(q::RealVector,p₁::RealVector,p₂::RealVector,p₃::RealVector,p₄::RealVector, θ)
    𝒏=NormalVector(p₁,p₂,p₃,p₄)
    return ℝ⁴⭢S³(q-2*dot(𝒏,q-p₁)*𝒏,θ)
end

function PickFace(cell::CELL, v::RealVector,POINTS³)
    if norm(v) ≈ 0
        error("vector v must be non-zero")
    end
    faces=vertices.(cell)
    fpts=(i->POINTS³[i]).(faces)
    return cell[findmax([dot(v,+(pts...)) for pts ∈ fpts])[2]]
end

function NewCell(cell::CELL,v::RealVector,θ,POINTS³, POINTS⁴)
    n=length(POINTS³)
    c=copy(cell)
    face=PickFace(cell,v,POINTS³)
    IND_face=vertices(face)
    IND_cell=Int[]
    for f ∈ c for e ∈ f for v ∈ e push!(IND_cell,v) end end end
    IND_cell=union(IND_cell)
    IND_cell2=[
        if i ∈ IND_face
            i
        else
            push!(POINTS³,Mirror(POINTS³[i],POINTS³[IND_face[1]],POINTS³[IND_face[2]],POINTS³[IND_face[3]]));
            if θ<ε
                push!(POINTS⁴,[POINTS³[end]...,0]);
            else
                O=[0,0,0,cot(θ)]
                push!(POINTS⁴,Mirror(POINTS⁴[i],O,POINTS⁴[IND_face[1]],POINTS⁴[IND_face[2]],POINTS⁴[IND_face[3]],θ));
            end
            n=n+1
        end
        for i ∈ IND_cell]
    return [[[IND_cell2[findfirst(w->w==v,IND_cell)] for v ∈ e] for e ∈ f] for f ∈ c]
end

function SphericalSphere(v,r::Real,θ::Real) where T<:RealVector
    V=S³⭢ℝ³(v,θ)
    return Sphere(V,r)
end

function SphericalCylinder(v₁,v₂,r::Real,θ::Real) where T<:RealVector
    w₁=ℝ⁴⭢S³((v₁+v₂)/2,θ)
    V₁=S³⭢ℝ³(v₁,θ)
    V₂=S³⭢ℝ³(v₂,θ)
    W₁=S³⭢ℝ³(w₁,θ)
    return Arc(V₁,W₁,V₂,r)
end

function SphericalPolygon(v::Array{T,1},θ::Real) where T<:RealVector
    n=length(v)
    u=ℝ⁴⭢S³(+(v...)/n,θ)
    v₁=v[1]
    v₂=v[2]
    v₃=v[3]
    V=(q->S³⭢ℝ³(q,θ)).(v)
    U=S³⭢ℝ³(u,θ)
    w=[ℝ⁴⭢S³((v[i]+v[mod(i,length(v))+1])/2,θ) for i ∈ 1:length(v)]
    W=(q->S³⭢ℝ³(q,θ)).(w)
    if θ<ε
        return Polygon(V)
    elseif rank(hcat(V...),atol=1.0e-4)==2
    # elseif rank(hcat(V...),atol=1.0e-12)==2
        m=4
        vw=copy(v)
        for _ ∈ 1:m
            l=length(vw)
            vw=[ℝ⁴⭢S³((vw[(i+1)÷2]+vw[mod(i÷2,l)+1])/2,θ) for i ∈ 1:2l]
        end
        VW=(q->S³⭢ℝ³(q,θ)).(vw)
        return csgUnion(Polygon(VW))
    else
        O=Circumcenter(U,V[1],V[2],V[3])
        sphere=Sphere(O,norm(U-O))
        N=NormalVector(V[1],W[1],V[2]);
        direction=sign(dot(N,O-U))

        cylinders=csgIntersection([
                (V₁=V[i];
                V₂=V[mod(i,n)+1];
                V₃=V[mod(i+1,n)+1];
                C=Circumcenter(V₁,W[i],V₂);
                N=NormalVector(V₁,W[i],V₂);
                direction=sign(dot(N,U-W[i]));
                cylinder=Cylinder(C,C+2*direction*norm(U-O)*N,norm(U-O)))
                for i ∈ 1:n
                ])
        return csgClip(sphere,cylinders)
    end
end

function Cells2Object(cells::Array{CELL,1},θ,POINTS⁴;rᵥ=0.05,rₑ=0.025,color=RGB(1,1,1))
    cs=copy(cells)
    fs=DeleteDuplicates(vcat(cs...))
    es=DeleteDuplicates(vcat(fs...))
    vs=DeleteDuplicates(vcat(es...))
    V=rgbColor(csgUnion([SphericalSphere(POINTS⁴[v],rᵥ,θ) for v ∈ vs]),RGB(0.1,0.1,0.1))
    E=rgbColor(csgUnion([SphericalCylinder(POINTS⁴[e[1]],POINTS⁴[e[2]],rₑ,θ) for e ∈ es]),RGB(0.2,0.2,0.2))
    # F=rgbftColor(csgUnion([SphericalPolygon([POINTS⁴[i] for i ∈ vertices(f)],θ) for f ∈ fs]),color,FT(0,0))
    F=rgbftColor(csgUnion([SphericalPolygon([POINTS⁴[i] for i ∈ vertices(f)],θ) for f ∈ fs]),color,FT(0.4,0.2))
    return csgUnion(V,E,F)
end

##%

function C₈(θ)
    h=1/2
    r=√(1-h^2)

    R=1/sin(θ)
    O=[0,0,0,cot(θ)]
    N=O+[0,0,0,R]
    H=O-[0,0,0,√(R^2-r^2)]

    POINTS, c₁=F₆
    POINTS=r*copy(POINTS)
    POINTS³=copy(POINTS)
    POINTS⁴=[[𝒑...,0]+H for 𝒑 ∈ POINTS]
    c₂=NewCell(c₁,[-1,0,0],θ,POINTS³, POINTS⁴)
    c₃=NewCell(c₂,[-1,0,0],θ,POINTS³, POINTS⁴)
    c₄=NewCell(c₁,[1,0,0],θ,POINTS³, POINTS⁴)
    c₅=NewCell(c₁,[0,0,-1],θ,POINTS³, POINTS⁴)
    c₆=NewCell(c₅,[0,-1,0],θ,POINTS³, POINTS⁴)
    c₇=NewCell(c₆,[0,-1,0],θ,POINTS³, POINTS⁴)
    c₈=NewCell(c₅,[0,1,0],θ,POINTS³, POINTS⁴)

    cells=[c₁,c₂,c₃,c₄,c₅,c₆,c₇,c₈]
    return cells, POINTS³, POINTS⁴
end


##%

function C₁₆(θ)
    h=1/2
    r=√(1-h^2)

    R=1/sin(θ)
    O=[0,0,0,cot(θ)]
    N=O+[0,0,0,R]
    H=O-[0,0,0,√(R^2-r^2)]

    POINTS, c₁=F₄
    POINTS=r*copy(POINTS)
    POINTS³=copy(POINTS)
    POINTS⁴=[[𝒑...,0]+H for 𝒑 ∈ POINTS]
    V=[1,-2,0]
    c₂=NewCell(c₁,V,θ,POINTS³, POINTS⁴)
    c₃=NewCell(c₂,V,θ,POINTS³, POINTS⁴)
    c₄=NewCell(c₃,V,θ,POINTS³, POINTS⁴)
    c₅=NewCell(c₁,-V,θ,POINTS³, POINTS⁴)
    c₆=NewCell(c₅,-V,θ,POINTS³, POINTS⁴)
    c₇=NewCell(c₆,-V,θ,POINTS³, POINTS⁴)
    c₈=NewCell(c₇,-V,θ,POINTS³, POINTS⁴)
    c₉=NewCell(c₁,[-2,-1,0],θ,POINTS³, POINTS⁴)
    W=[-1,-0.5,1]
    c₁₀=NewCell(c₉,W,θ,POINTS³, POINTS⁴)
    c₁₁=NewCell(c₁₀,W,θ,POINTS³, POINTS⁴)
    c₁₂=NewCell(c₁₁,W,θ,POINTS³, POINTS⁴)
    c₁₃=NewCell(c₉,-W,θ,POINTS³, POINTS⁴)
    c₁₄=NewCell(c₁₃,-W,θ,POINTS³, POINTS⁴)
    c₁₅=NewCell(c₁₄,-W,θ,POINTS³, POINTS⁴)
    c₁₆=NewCell(c₁₅,-W,θ,POINTS³, POINTS⁴)

    cells=[c₁,c₂,c₃,c₄,c₅,c₆,c₇,c₈,c₉,c₁₀,c₁₁,c₁₂,c₁₃,c₁₄,c₁₅,c₁₆]
    return cells, POINTS³, POINTS⁴
end

##%

function C₂₄(θ)
    h=1/√2
    r=√(1-h^2)

    R=1/sin(θ)
    O=[0,0,0,cot(θ)]
    N=O+[0,0,0,R]
    H=O-[0,0,0,√(R^2-r^2)]

    POINTS, c₁=F₈
    POINTS=r*copy(POINTS)
    POINTS³=copy(POINTS)
    POINTS⁴=[[𝒑...,0]+H for 𝒑 ∈ POINTS]
    c₂=NewCell(c₁,[1,1,1],θ,POINTS³, POINTS⁴)
    c₃=NewCell(c₂,[1,1,1],θ,POINTS³, POINTS⁴)
    c₄=NewCell(c₃,[1,1,1],θ,POINTS³, POINTS⁴)
    c₅=NewCell(c₁,-[1,1,1],θ,POINTS³, POINTS⁴)
    c₆=NewCell(c₅,-[1,1,1],θ,POINTS³, POINTS⁴)
    c₇=NewCell(c₁,[1,-1,-1],θ,POINTS³, POINTS⁴)
    c₈=NewCell(c₁,[-1,1,-1],θ,POINTS³, POINTS⁴)
    c₉=NewCell(c₁,[-1,-1,1],θ,POINTS³, POINTS⁴)
    c₁₀=NewCell(c₂,[-1,0,0],θ,POINTS³, POINTS⁴)
    c₁₁=NewCell(c₂,[0,-1,0],θ,POINTS³, POINTS⁴)
    c₁₂=NewCell(c₂,[0,0,-1],θ,POINTS³, POINTS⁴)
    c₁₃=NewCell(c₃,[1,-1,-1],θ,POINTS³, POINTS⁴)
    c₁₄=NewCell(c₃,[-1,1,-1],θ,POINTS³, POINTS⁴)
    c₁₅=NewCell(c₃,[-1,-1,1],θ,POINTS³, POINTS⁴)
    c₁₆=NewCell(c₄,[-1,0,0],θ,POINTS³, POINTS⁴)
    c₁₇=NewCell(c₄,[0,-1,0],θ,POINTS³, POINTS⁴)
    c₁₈=NewCell(c₄,[0,0,-1],θ,POINTS³, POINTS⁴)
    c₁₉=NewCell(c₆,[1,-1,-1],θ,POINTS³, POINTS⁴)
    c₂₀=NewCell(c₆,[-1,1,-1],θ,POINTS³, POINTS⁴)
    c₂₁=NewCell(c₆,[-1,-1,1],θ,POINTS³, POINTS⁴)
    c₂₂=NewCell(c₅,[-1,0,0],θ,POINTS³, POINTS⁴)
    c₂₃=NewCell(c₅,[0,-1,0],θ,POINTS³, POINTS⁴)
    c₂₄=NewCell(c₅,[0,0,-1],θ,POINTS³, POINTS⁴)


    cells=[c₁,c₂,c₃,c₄,c₅,c₆,c₇,c₈,c₉,c₁₀,c₁₁,c₁₂,c₁₃,c₁₄,c₁₅,c₁₆,c₁₇,c₁₈,c₁₉,c₂₀,c₂₁,c₂₂,c₂₃,c₂₄]
    return cells, POINTS³, POINTS⁴
end

##%

function C₁₂₀(θ)
    h=norm([φ^3,0,0,0])/norm([φ^3,1,1,1])
    r=√(1-h^2)

    R=1/sin(θ)
    O=[0,0,0,cot(θ)]
    N=O+[0,0,0,R]
    H=O-[0,0,0,√(R^2-r^2)]

    POINTS, c₁=F₁₂
    POINTS=r*copy(POINTS)
    POINTS³=copy(POINTS)
    POINTS⁴=[[𝒑...,0]+H for 𝒑 ∈ POINTS]
    Va,_=F₂₀
    Vb=[v-2*Va[1]*dot(Va[1],v) for v ∈ Va]
    c₂=NewCell(c₁,Va[1],θ,POINTS³, POINTS⁴)
    c₃=NewCell(c₂,Va[1],θ,POINTS³, POINTS⁴)
    c₄=NewCell(c₃,Va[1],θ,POINTS³, POINTS⁴)
    c₅=NewCell(c₄,Va[1],θ,POINTS³, POINTS⁴)
    c₆=NewCell(c₅,Va[1],θ,POINTS³, POINTS⁴)
    c₇=NewCell(c₁,-Va[1],θ,POINTS³, POINTS⁴)
    c₈=NewCell(c₇,-Va[1],θ,POINTS³, POINTS⁴)
    c₉=NewCell(c₈,-Va[1],θ,POINTS³, POINTS⁴)
    c₁₀=NewCell(c₉,-Va[1],θ,POINTS³, POINTS⁴)

    cells=[c₁,c₂,c₃,c₄,c₅,c₆,c₇,c₈,c₉,c₁₀]

    push!(cells,NewCell(cells[1],Va[3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[1],Va[7],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[1],Va[8],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[1],Va[10],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[1],Va[12],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[2],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[2],Vb[5],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[2],Vb[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[2],Vb[9],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[2],Vb[11],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[3],Va[3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[3],Va[7],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[3],Va[8],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[3],Va[10],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[3],Va[12],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[4],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[4],Vb[5],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[4],Vb[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[4],Vb[9],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[4],Vb[11],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[5],Va[3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[5],Va[7],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[5],Va[8],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[5],Va[10],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[5],Va[12],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[6],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[6],Vb[5],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[6],Vb[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[6],Vb[9],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[6],Vb[11],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[7],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[7],Vb[5],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[7],Vb[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[7],Vb[9],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[7],Vb[11],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[8],Va[3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[8],Va[7],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[8],Va[8],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[8],Va[10],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[8],Va[12],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[9],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[9],Vb[5],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[9],Vb[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[9],Vb[9],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[9],Vb[11],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[10],Va[3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[10],Va[7],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[10],Va[8],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[10],Va[10],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[10],Va[12],θ,POINTS³, POINTS⁴))


    push!(cells,NewCell(cells[18],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],Va[5],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],-Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],-Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],-Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],-Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end-4],Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],Va[11],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[end],Va[11],θ,POINTS³, POINTS⁴))


    push!(cells,NewCell(cells[17+45],Va[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[17+45],Va[4],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[17+45],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[17+45],Vb[12],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[18+45],Va[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[18+45],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[18+45],Va[1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[18+45],Va[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[18+45],Va[8],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[19+45],[-1,-2,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[19+45],[1,-1,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[19+45],Vb[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[19+45],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[19+45],[-1,-1,-1],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[20+45],Va[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[20+45],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[20+45],Va[1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[20+45],Va[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[20+45],Va[8],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[21+45],[-1,-2,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[21+45],[1,-1,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[21+45],Vb[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[21+45],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[21+45],[-1,-1,-1],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[22+45],Va[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[22+45],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[22+45],Va[1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[22+45],Va[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[22+45],Va[8],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[23+45],[-1,-2,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[23+45],[1,-1,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[23+45],[2,-2,-1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[23+45],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[23+45],[-1,-1,-1],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[24+45],Va[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[24+45],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[24+45],Va[1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[24+45],Va[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[24+45],Va[8],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[25+45],[-1,-2,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[25+45],[1,-1,3],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[25+45],[2,-2,-1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[25+45],Vb[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[25+45],[-1,-1,-1],θ,POINTS³, POINTS⁴))

    push!(cells,NewCell(cells[26+45],Va[2],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[26+45],Va[12],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[26+45],Va[1],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[26+45],Va[6],θ,POINTS³, POINTS⁴))
    push!(cells,NewCell(cells[26+45],Va[8],θ,POINTS³, POINTS⁴))



    return cells, POINTS³, POINTS⁴
end

##%
M=120
for i ∈ 0:2M
    θ=π/2*(Smooth(0,1,i/M)-Smooth(1,2,i/M))
    cells, POINTS³, POINTS⁴=C₈(θ)
    render(Cells2Object(cells,θ,POINTS⁴,color=RGB(0.2,1,1),rₑ=0.012,rᵥ=0.03),camera=LngLatCamera(lng=180°+360°*i/M,lat=25°,pers=0.2,zoom=0.15,width=1200,height=900),name="8-Cell",index=i+1)
end

M=120
for i ∈ 0:2M
    θ=π/2*(Smooth(0,1,i/M)-Smooth(1,2,i/M))
    cells, POINTS³, POINTS⁴=C₁₆(θ)
    render(Cells2Object(cells,θ,POINTS⁴,color=RGB(0.2,1,1)),camera=LngLatCamera(lng=180°+360°*i/M,lat=25°,pers=0.2,zoom=0.12,width=600,height=450),name="16-Cell",index=i+1)
end

M=120
for i ∈ 0:2M
    θ=π/2*(Smooth(0,1,i/M)-Smooth(1,2,i/M))
    cells, POINTS³, POINTS⁴=C₂₄(θ)
    render(Cells2Object(cells,θ,POINTS⁴,color=RGB(0.2,1,1)),camera=LngLatCamera(lng=180°+360°*i/M,lat=25°,pers=0.2,zoom=0.12,width=600,height=450),name="24-Cell",index=i+1)
end

M=120
i=24
θ=π/2*(Smooth(0,1,i/M)-Smooth(1,2,i/M))
cells, POINTS³, POINTS⁴=C₁₂₀(θ)
render(Cells2Object(cells,θ,POINTS⁴,color=RGB(0.2,1,1)),camera=LngLatCamera(lng=180°+360°*i/M,lat=25°,pers=0.2,zoom=0.12,width=600,height=450),name="120-CellB",index=i+1)

M=120
for i ∈ 0:2M
    θ=π/2*(Smooth(0,1,i/M)-Smooth(1,2,i/M))
    cells, POINTS³, POINTS⁴=C₁₂₀(θ)
    try
        render(Cells2Object(cells,θ,POINTS⁴,color=RGB(0.2,1,1)),camera=LngLatCamera(lng=180°+360°*i/M,lat=25°,pers=0.2,zoom=0.07,width=500,height=500),name="120-CellD",index=i+1)
    catch
        println(i)
    end
end
