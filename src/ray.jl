struct Ray{T<:Real} 
    origin::Point{3, T}
    direction::Vec{3, T}
end

@inline origin(ray::Ray) = ray.origin
@inline direction(ray::Ray) = ray.direction
@inline at(ray::Ray{T}, t::T) where {T<:Real} = origin(ray) .+ t .* direction(ray)

import Base.≈
≈(ray1::Ray, ray2::Ray; kwargs...) = ≈(origin(ray1), origin(ray2); kwargs...) && ≈(direction(ray1), direction(ray2); kwargs...)


struct HitRecord{T<:Real}
    point::Point{3, T}
    normal::Vec{3, T}
    t::T
    front_face::Bool
end

point(hitrecord::HitRecord) = hitrecord.point
normal(hitrecord::HitRecord) = hitrecord.normal
t(hitrecord::HitRecord) = hitrecord.t
front_face(hitrecord::HitRecord) = hitrecord.front_face

function face_normal(ray::Ray{T}, outward_normal::Vec{3, T}) where {T<:Real}
    front_face = dot(direction(ray), outward_normal) < 0
    normal = front_face ? outward_normal : -outward_normal
    (front_face, normal)
end

function HitRecord(ray::Ray{T}, outward_normal::Vec{3, T}, t::T) where {T<:Real}
    (front_face, normal) = face_normal(ray, outward_normal)
    HitRecord(at(ray, t), normal, t, front_face)
end






