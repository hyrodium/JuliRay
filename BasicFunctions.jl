## Basic functions

function DeleteDuplicates(v::Array{T,1})::Array{T,1} where T <: Any
    w=Array{eltype(v),1}()
    if (length(v)>0)
        for e ∈ v
            if e ∉ w
                push!(w,e)
            end
        end
    end
    return w
end

function OrthogonalVector(v::Array{T,1}) where T <: Real
    if(norm(v)==0)
        error("zero vector")
    else
        v₁=normalize(v)
        i=argmin(abs.(v))
        v₂=zeros(length(v))
        v₂[i]=1.0
        v₂=normalize(v₂-dot(v₁,v₂)*v₁)
        return v₂
    end
end

function isorthogonal(M::AbstractMatrix{T}) where T<:Real
    size(M)
    return M*M'≈oneunit(M)
end

function rotatematrix(v::AbstractVector,θ::Real)
    n=normalize(v)
    N=[0 -n[3] n[2];n[3] 0 -n[1];-n[2] n[1] 0]
    return (oneunit(N)+N^2+(-cos(θ)*N^2+sin(θ)*N))
end

function Circumcenter(point₁,point₂,point₃)
    l₁²=(norm(point₂-point₃))^2
    l₂²=(norm(point₃-point₁))^2
    l₃²=(norm(point₁-point₂))^2
    c₁=l₁²*(l₂²+l₃²-l₁²)
    c₂=l₂²*(l₃²+l₁²-l₂²)
    c₃=l₃²*(l₁²+l₂²-l₃²)
    center=(c₁*point₁+c₂*point₂+c₃*point₃)/(c₁+c₂+c₃)
    return center
end

function NormalVector(v₁::Array{Float64,1},v₂::Array{Float64,1},v₃::Array{Float64,1})
    return normalize(cross(v₁,v₂)+cross(v₂,v₃)+cross(v₃,v₁))
end

