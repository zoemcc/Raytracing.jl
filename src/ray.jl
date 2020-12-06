struct Ray{T<:Real} 
    origin::Point{3, T}
    direction::Vec{3, T}
end

origin(ray::Ray) = ray.origin
direction(ray::Ray) = ray.direction

import Base.≈
≈(ray1::Ray, ray2::Ray; kwargs...) = ≈(origin(ray1), origin(ray2); kwargs...) && ≈(direction(ray1), direction(ray2); kwargs...)

