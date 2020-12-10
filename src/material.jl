abstract type AbstractMaterial end

@inline function scatter(mat::AbstractMaterial, ray_in::Ray{T}, hitrecord::HitRecord{T}, rng::MersenneTwister) where {T<:Real} 
     error("Scatter for $(typeof(mat)) has not been defined.")
end

struct Lambertian{T<:Real} <: AbstractMaterial
    albedo::Vec{3, T}
end

albedo(lambertian::Lambertian{T}) where {T<:Real} = lambertian.albedo

@inline function scatter(lambertian::Lambertian{T}, ray_in::Ray{T}, hitrecord::HitRecord{T}, rng::MersenneTwister) where {T<:Real}
    color = albedo(lambertian)
    facing_normal = normal(hitrecord)
    new_dir = normalize(facing_normal + random_unit_vector3(rng, T))
    ray_out = Ray(point(hitrecord), new_dir) # direction could be the normal?
    absorbed = false
    (color=color, ray_out=ray_out, absorbed=absorbed)
end

struct NormalShading <: AbstractMaterial
end

@inline function scatter(::NormalShading, ray_in::Ray{T}, hitrecord::HitRecord{T}, rng::MersenneTwister) where {T<:Real}
    color = (((normal(hitrecord) .* (front_face(hitrecord) ? 1 : -1)) .+ 1) ./ 2)
    ray_out = Ray(point(hitrecord), direction(ray_in)) # direction could be the normal?
    absorbed = true
    (color=color, ray_out=ray_out, absorbed=absorbed)
end

