VERTEX=Int
EDGE=Array{VERTEX,1}
FACE=Array{EDGE,1}
CELL=Array{FACE,1}

function vertices(edge::EDGE)
    return edge
end
function vertices(face::FACE)
    f=copy(face)
    v=union(vcat(f...))
    w=[v[1]]
    for j ∈ 1:(length(v)-1)
        i=findfirst(e->w[end] ∈ e, f)
        push!(w,filter(v->v≠w[end],f[i])[1])
        deleteat!(f,i)
    end
    return w
end
function vertices(cell::CELL)
    return unique(vcat(vcat(cell...)...))
end
function edges(face::FACE)
    return face
end
function edges(cell::CELL)
    return unique(sort.(vcat(cell...)))
end
function faces(cell::CELL)
    return cell
end

import Base.∈
function ∈(vertex::VERTEX, face::FACE)
    return vertex ∈ vertices(face)
end
function ∈(vertex::VERTEX, cell::CELL)
    return vertex ∈ vertices(cell)
end
function ∈(edge::EDGE, cell::CELL)
    return edge ∈ edges(cell)
end

function DualRegularPlolyhedron(POINTS,c)
    fs=copy(c)
    es=edges(fs)
    vs=vertices(fs)
    POINTS_=normalize.([+(POINTS[vertices(f)]...) for f ∈ fs])
    vs_=collect(1:length(faces(fs)))
    es_=[findall(f->e ∈ f, fs) for e ∈ es]
    fs_=[es_[findall(e->issubset(e,findall(f->v ∈ f, c)),es_)] for v ∈ vs]
    return POINTS_, fs_
end

# F₄
𝒑₁=[1,1,1]
𝒑₂=[1,-1,-1]
𝒑₃=[-1,1,-1]
𝒑₄=[-1,-1,1]
POINTS=[𝒑₁,𝒑₂,𝒑₃,𝒑₄]
v₁=1
v₂=2
v₃=3
v₄=4
e₁=[v₁,v₂]
e₂=[v₁,v₃]
e₃=[v₁,v₄]
e₄=[v₂,v₃]
e₅=[v₂,v₄]
e₆=[v₃,v₄]

f₁=[e₄,e₅,e₆]
f₂=[e₂,e₃,e₆]
f₃=[e₁,e₃,e₅]
f₄=[e₁,e₂,e₄]
c₁=[f₁,f₂,f₃,f₄]
F₄=normalize.(POINTS), c₁

# F₆
𝒑₁=[1,1,1]
𝒑₂=[1,1,-1]
𝒑₃=[1,-1,1]
𝒑₄=[1,-1,-1]
𝒑₅=[-1,1,1]
𝒑₆=[-1,1,-1]
𝒑₇=[-1,-1,1]
𝒑₈=[-1,-1,-1]
POINTS=[𝒑₁,𝒑₂,𝒑₃,𝒑₄,𝒑₅,𝒑₆,𝒑₇,𝒑₈]
v₁=1
v₂=2
v₃=3
v₄=4
v₅=5
v₆=6
v₇=7
v₈=8
e₁=[v₁,v₂]
e₂=[v₁,v₃]
e₃=[v₁,v₅]
e₄=[v₂,v₄]
e₅=[v₂,v₆]
e₆=[v₃,v₄]
e₇=[v₃,v₇]
e₈=[v₄,v₈]
e₉=[v₅,v₆]
e₁₀=[v₅,v₇]
e₁₁=[v₆,v₈]
e₁₂=[v₇,v₈]
f₁=[e₁,e₂,e₄,e₆]
f₂=[e₁,e₃,e₅,e₉]
f₃=[e₂,e₃,e₇,e₁₀]
f₄=[e₄,e₅,e₈,e₁₁]
f₅=[e₆,e₇,e₈,e₁₂]
f₆=[e₉,e₁₀,e₁₁,e₁₂]
c₁=[f₁,f₂,f₃,f₄,f₅,f₆]
F₆=normalize.(POINTS), c₁

F₈=DualRegularPlolyhedron(F₆...)

# F₁₂
𝒑₁=[1,1,1]
𝒑₂=[1,1,-1]
𝒑₃=[1,-1,1]
𝒑₄=[1,-1,-1]
𝒑₅=[-1,1,1]
𝒑₆=[-1,1,-1]
𝒑₇=[-1,-1,1]
𝒑₈=[-1,-1,-1]
𝒑₉=[0,1/φ,φ]
𝒑₁₀=[0,1/φ,-φ]
𝒑₁₁=[0,-1/φ,φ]
𝒑₁₂=[0,-1/φ,-φ]
𝒑₁₃=[φ,0,1/φ]
𝒑₁₄=[φ,0,-1/φ]
𝒑₁₅=[-φ,0,1/φ]
𝒑₁₆=[-φ,0,-1/φ]
𝒑₁₇=[1/φ,φ,0]
𝒑₁₈=[1/φ,-φ,0]
𝒑₁₉=[-1/φ,φ,0]
𝒑₂₀=[-1/φ,-φ,0]
POINTS=[𝒑₁,𝒑₂,𝒑₃,𝒑₄,𝒑₅,𝒑₆,𝒑₇,𝒑₈,𝒑₉,𝒑₁₀,𝒑₁₁,𝒑₁₂,𝒑₁₃,𝒑₁₄,𝒑₁₅,𝒑₁₆,𝒑₁₇,𝒑₁₈,𝒑₁₉,𝒑₂₀]
v₁=1
v₂=2
v₃=3
v₄=4
v₅=5
v₆=6
v₇=7
v₈=8
v₉=9
v₁₀=10
v₁₁=11
v₁₂=12
v₁₃=13
v₁₄=14
v₁₅=15
v₁₆=16
v₁₇=17
v₁₈=18
v₁₉=19
v₂₀=20
e₁=[v₁,v₉]
e₂=[v₅,v₉]
e₃=[v₃,v₁₁]
e₄=[v₇,v₁₁]
e₅=[v₉,v₁₁]
e₆=[v₂,v₁₀]
e₇=[v₆,v₁₀]
e₈=[v₄,v₁₂]
e₉=[v₈,v₁₂]
e₁₀=[v₁₀,v₁₂]
e₁₁=[v₁,v₁₃]
e₁₂=[v₃,v₁₃]
e₁₃=[v₂,v₁₄]
e₁₄=[v₄,v₁₄]
e₁₅=[v₁₃,v₁₄]
e₁₆=[v₅,v₁₅]
e₁₇=[v₇,v₁₅]
e₁₈=[v₆,v₁₆]
e₁₉=[v₈,v₁₆]
e₂₀=[v₁₅,v₁₆]
e₂₁=[v₁,v₁₇]
e₂₂=[v₂,v₁₇]
e₂₃=[v₅,v₁₉]
e₂₄=[v₆,v₁₉]
e₂₅=[v₁₇,v₁₉]
e₂₆=[v₃,v₁₈]
e₂₇=[v₄,v₁₈]
e₂₈=[v₇,v₂₀]
e₂₉=[v₈,v₂₀]
e₃₀=[v₁₈,v₂₀]
f₁=[e₁,e₃,e₅,e₁₁,e₁₂]
f₂=[e₂,e₄,e₅,e₁₆,e₁₇]
f₃=[e₆,e₈,e₁₀,e₁₃,e₁₄]
f₄=[e₇,e₉,e₁₀,e₁₈,e₁₉]
f₅=[e₁₁,e₁₃,e₁₅,e₂₁,e₂₂]
f₆=[e₁₂,e₁₄,e₁₅,e₂₆,e₂₇]
f₇=[e₁₆,e₁₈,e₂₀,e₂₃,e₂₄]
f₈=[e₁₇,e₁₉,e₂₀,e₂₈,e₂₉]
f₉=[e₁,e₂,e₂₁,e₂₃,e₂₅]
f₁₀=[e₆,e₇,e₂₂,e₂₄,e₂₅]
f₁₁=[e₃,e₄,e₂₆,e₂₈,e₃₀]
f₁₂=[e₈,e₉,e₂₇,e₂₉,e₃₀]
c₁=[f₁,f₂,f₃,f₄,f₅,f₆,f₇,f₈,f₉,f₁₀,f₁₁,f₁₂]
F₁₂=normalize.(POINTS), c₁

F₂₀=DualRegularPlolyhedron(F₁₂...)
