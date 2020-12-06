struct Camera{T<:Real}
    origin::Point3{T}
    rotation::RotMatrix3{T}
    horizontal::Vec3{T}
    vertical::Vec3{T}
    lower_left_corner::Vec3{T}
    aspect_ratio::T
end

origin(cam::Camera) = cam.origin
rotation(cam::Camera) = cam.rotation
horizontal(cam::Camera) = cam.horizontal
vertical(cam::Camera) = cam.vertical
lower_left_corner(cam::Camera) = cam.lower_left_corner
aspect_ratio(cam::Camera) = cam.aspect_ratio

function generate_ray(camera::Camera{T}, s::T, t::T)::Ray{T} where {T<:Real}
    Ray(origin(camera), normalize(lower_left_corner(camera) + s*horizontal(camera) + t*vertical(camera) - origin(camera)))
end


function make_camera(lookfrom::Point3{T}, lookat::Point3{T}, vup::Vec3{T}, vfov::T, aspect_ratio::T)::Camera{T} where {T<:Real}
    theta = deg2rad(vfov)
    h = tan(theta/2)
    viewport_height = 2 * h
    viewport_width = aspect_ratio * viewport_height

    w = normalize(Vec3{T}(lookfrom - lookat))
    u = normalize(cross(vup, w))
    v = cross(w, u)

    origin = lookfrom
    rotation = RotMatrix3{T}([w;u;v])
    horizontal = viewport_width * u
    vertical = viewport_height * v
    lower_left_corner = Vec3{T}(origin - horizontal/2 - vertical/2 - w)

    camera = Camera{T}(origin, rotation, horizontal, vertical, lower_left_corner, aspect_ratio)
end
