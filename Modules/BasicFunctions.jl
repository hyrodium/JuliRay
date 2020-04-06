## Basic functions

# RealVector=Union{Array{T,1} where T <: Real}
RealVector = Vector{<:Real}
# Base.@irrational ° 0.0174532925199432957692369076848861271344 (big(pi)/big(180))

function DeleteDuplicates(v::Array{T,1})::Array{T,1} where T <: Any
    w = Array{eltype(v),1}()
    if length(v) > 0
        for e ∈ v
            if e ∉ w
                push!(w,e)
            end
        end
    end
    return w
end

function OrthogonalVector(v::RealVector) where T <: Real
    if norm(v)==0
        error("zero vector")
    else
        v₁ = normalize(v)
        i = argmin(abs.(v))
        v₂ = zeros(length(v))
        v₂[i] = 1.0
        v₂ = normalize(v₂-dot(v₁,v₂)*v₁)
        return v₂
    end
end

function isorthogonal(M::AbstractMatrix{T}) where T<:Real
    size(M)
    return M*M'≈oneunit(M)
end

function rotatematrix(v::RealVector,θ::Real)
    n = normalize(v)
    N = [0 -n[3] n[2];n[3] 0 -n[1];-n[2] n[1] 0]
    return oneunit(N)+N^2+(-cos(θ)*N^2+sin(θ)*N)
end

function Circumcenter(p₁::RealVector,p₂)
    c₁ = 1
    c₂ = 1
    center = (c₁*p₁+c₂*p₂)/(c₁+c₂)
    return center
end

function Circumcenter(p₁::RealVector,p₂::RealVector,p₃::RealVector)
    l₁² = (norm(p₂-p₃))^2
    l₂² = (norm(p₃-p₁))^2
    l₃² = (norm(p₁-p₂))^2
    c₁ = l₁²*(l₂²+l₃²-l₁²)
    c₂ = l₂²*(l₃²+l₁²-l₂²)
    c₃ = l₃²*(l₁²+l₂²-l₃²)
    center = (c₁*p₁+c₂*p₂+c₃*p₃)/(c₁+c₂+c₃)
    return center
end

function Circumcenter(p₁::RealVector,p₂::RealVector,p₃::RealVector,p₄::RealVector)
    n = 3
    𝒑 = hcat(p₁,p₂,p₃,p₄)
    p² = [dot(𝒑[1:n,i],𝒑[1:n,i]) for i in 1:(n+1)]
    P = vcat(𝒑,reshape([1. for i in 1:(n+1)],1,n+1))
    Q = inv(P)
    𝒄 = (Q'*p²)[1:n]/2
    r = sqrt(dot(𝒄,𝒄)+dot(p²,Q[:,end]))
    return 𝒄
end

function NormalVector(p₁::RealVector,p₂::RealVector,p₃::RealVector)
    return normalize(cross(p₁,p₂)+cross(p₂,p₃)+cross(p₃,p₁))
end

function Angles(p₁::RealVector,p₂::RealVector,p₃::RealVector)
    l₁ = norm(p₂-p₃)
    l₂ = norm(p₃-p₁)
    l₃ = norm(p₁-p₂)
    cosα₁ = dot(p₂-p₁,p₃-p₁)/(l₂*l₃)
    cosα₂ = dot(p₃-p₂,p₁-p₂)/(l₃*l₁)
    cosα₃ = dot(p₁-p₃,p₂-p₃)/(l₁*l₂)
    α₁ = acos(clamp(cosα₁,-1,1))
    α₂ = acos(clamp(cosα₂,-1,1))
    α₃ = acos(clamp(cosα₃,-1,1))
    return α₁,α₂,α₃
end
