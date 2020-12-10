
function raytrace_image(scene::Scene, cam::Camera{T}, image_height::Integer, samples_per_pixel::Integer, max_bounces_per_ray::Integer, max_steps_per_bounce::Integer, distance_tolerance::T) where {N<:Integer, T<:Real}
    aspect_ratio_cam = aspect_ratio(cam)
    image_width = Int64(floor(aspect_ratio_cam * image_height))
    out_image = zeros(RGB{T}, image_height, image_width)
    for i in 1:image_height, j in 1:image_width
        s = T((j - 1)/(image_width - 1))
        t = T((i - 1)/(image_height - 1))
        ray = generate_ray(cam, s, t)
        ray_color = raymarch(scene, ray, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance)
        out_image[i, j] = ray_color
    end
    out_image
end

function raymarch(scene::Scene, ray::Ray{T}, max_bounces_left::Integer, max_steps_per_bounce::Integer, distance_tolerance::T) where {N<:Integer, T<:Real}
    if max_bounces_left > 0
        t = T(0)
        for i in 1:max_steps_per_bounce
            current_point = origin(ray) + t*direction(ray)
            current_min_dist, current_min_index = min_dist_index(scene, current_point)
            if current_min_dist < t * distance_tolerance
                min_shape = shapes(scene)[current_min_index]
                outward_normal = normal(sdf(min_shape), current_point)
                hitrecord = HitRecord(ray, outward_normal, t)
                color, ray_out, absorbed = scatter(mat(min_shape), ray, hitrecord)
                if absorbed
                    return color
                else
                    return color * raymarch(scene, ray_out, max_bounces_left - 1, max_steps_per_bounce, distance_tolerance)
                end
                #return RGB(((normal(sdf(shapes(scene)[current_min_index]), current_point) .+ 1) ./ 2)...)
                #return albedo(shapes(scene)[current_min_index])
            else
                # Assumes the ray direction is normalized
                t += current_min_dist
            end
        end
    end
    # max bounces <= 0 or ran out of march steps
    zero(RGB{T})
end

saveimage(savename::AbstractString, image) = save("./images/$(savename)", map(clamp01nan, image))


