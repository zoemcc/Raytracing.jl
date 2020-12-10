abstract type AbstractMaterial end

@inline function scatter(mat::AbstractMaterial, ray_in::Ray{T}, hitrecord::HitRecord{T}) where {T<:Real} 
     error("Scatter for $(typeof(mat)) has not been defined.")
end

struct Lambertian{T<:Real} <: AbstractMaterial
    albedo::RGB{T}
end

struct NormalShading <: AbstractMaterial
end

@inline function scatter(::NormalShading, ray_in::Ray{T}, hitrecord::HitRecord{T}) where {T<:Real}
    color = RGB((((normal(hitrecord) .* (front_face(hitrecord) ? 1 : -1)) .+ 1) ./ 2)...)
    ray_out = Ray(at(ray_in, t(hitrecord)), direction(ray_in)) # direction could be the normal?
    absorbed = true
    (color=color, ray_out=ray_out, absorbed=absorbed)
end

