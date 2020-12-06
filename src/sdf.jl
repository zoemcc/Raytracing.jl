abstract type AbstractSignedDistanceField{T<:Real} end

(sdf::AbstractSignedDistanceField{T})(::Point{3, T}) where {T<:Real} = error("Signed Distance Field for $(typeof(sdf)) has not been defined.")
normal(sdf::AbstractSignedDistanceField{T}, ::Point{3, T}) where {T<:Real} = error("Normal of Signed Distance Field for $(typeof(sdf)) has not been defined.")

struct SphereSignedDistanceField{T<:Real} <: AbstractSignedDistanceField{T}
    radius::T
end

radius(sdf::SphereSignedDistanceField) = sdf.radius
(sdf::SphereSignedDistanceField{T})(point::Point{3, T}) where {T<:Real} = norm(point) - radius(sdf)

# Currently always pointing out 
normal(::SphereSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = Vec3{T}(normalize(point))
normalzygote(sdf::SphereSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(gradient(sdf, point)[1]))
