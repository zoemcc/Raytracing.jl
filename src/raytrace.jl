
function raytrace_image(scene::Scene, cam::Camera{T}, image_height::Integer, samples_per_pixel::Integer, max_bounces_per_ray::Integer, max_steps_per_bounce::Integer, distance_tolerance::T, rng::MersenneTwister) where {N<:Integer, T<:Real}
    aspect_ratio_cam = aspect_ratio(cam)
    image_width = Int64(floor(aspect_ratio_cam * image_height))
    out_image = zeros(RGB{T}, image_height, image_width)
    for i in 1:image_height, j in 1:image_width
        pixel_color = zero(Vec{3, T})
        for rands in 1:samples_per_pixel
            s = T((j - 1 + rand(rng))/(image_width - 1))
            t = T((i - 1 + rand(rng))/(image_height - 1))
            ray = generate_ray(cam, s, t)
            ray_color = raymarch(scene, ray, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance, rng::MersenneTwister)
            pixel_color += ray_color ./ samples_per_pixel 
        end
        out_image[i, j] = RGB{T}(pixel_color...)
    end
    out_image
end

function raytrace_image_cuda(scene::Scene, cam::Camera{T}, image_height::Integer, samples_per_pixel::Integer, max_bounces_per_ray::Integer, max_steps_per_bounce::Integer, distance_tolerance::T, rng::MersenneTwister) where {N<:Integer, T<:Real}
    aspect_ratio_cam = aspect_ratio(cam)
    image_width = Int64(floor(aspect_ratio_cam * image_height))
    out_image = zeros(RGB{T}, image_height, image_width)
    for i in 1:image_height, j in 1:image_width
        pixel_color = zero(Vec{3, T})
        for rands in 1:samples_per_pixel
            s = T((j - 1 + rand(rng))/(image_width - 1))
            t = T((i - 1 + rand(rng))/(image_height - 1))
            ray = generate_ray(cam, s, t)
            ray_color = raymarch(scene, ray, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance, rng::MersenneTwister)
            pixel_color += ray_color ./ samples_per_pixel 
        end
        out_image[i, j] = RGB{T}(pixel_color...)
    end
    out_image
end



#=

CUDA.@sync begin
        #@cuda threads=numthreads blocks=numblocks mandelbrotandregiongpu!(cudaimage,
            #centerx, xstart, xrangeextent, centery, ystart, yrangeextent,
            #numiters, height, width, cosrot, sinrot, a, b, c, d)
        @cuda threads=numthreads blocks=numblocks mandelbrotandregiongpurational!(cudaimage,
            centerx, xstart, xrangeextent, centery, ystart, yrangeextent,
            numiters, height, width, cosrot, sinrot, conformal)
    end

    CUDA.copyto!(outimage[], cudaimage

    function mandelbrotandregiongpurational!(escape_color, centerx::T, xstart::T, xrangeextent::T, centery::T, ystart::T, yrangeextent::T, 
        numiters, height, width, cosrot::T, sinrot::T, rational::ComplexRational{T, N}) where {T <: Real, N}
    indexx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    indexy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stridex = blockDim().x * gridDim().x
    stridey = blockDim().y * gridDim().y
    #@cuprintln("thread $indexx, $indexy, stride $stridex, $stridey, weight $width $height")
    for i in indexx:stridex:width, j in indexy:stridey:height
        prex_ij = (i - 1) / (width - 1) * xrangeextent + xstart
        prey_ij = (j - 1) / (height - 1) * yrangeextent + ystart
        z_ij = Complex{T}(cosrot * prex_ij - sinrot * prey_ij + centerx, sinrot * prex_ij + cosrot * prey_ij + centery)
        z_ij_rational = transform(rational, z_ij)
        escape_ij = mandelbrot(z_ij_rational, numiters) / numiters
        escape_ij_f32 = Float32(escape_ij)
        @inbounds escape_color[i, j] = ColorTypes.RGBA{Float32}(escape_ij_f32, escape_ij_f32, escape_ij_f32)
    end
    return nothing
end
=#


function raymarch(scene::Scene, ray::Ray{T}, max_bounces_left::Integer, max_steps_per_bounce::Integer, distance_tolerance::T, rng::MersenneTwister) where {T<:Real}
    if max_bounces_left > 0
        t = T(0)
        for i in 1:max_steps_per_bounce
            current_point = origin(ray) + t*direction(ray)
            current_min_dist, current_min_index = min_dist_index(scene, current_point)
            if current_min_dist < t * distance_tolerance
                min_shape = shapes(scene)[current_min_index]
                outward_normal = normal(sdf(min_shape), current_point)
                hitrecord = HitRecord(ray, outward_normal, t)
                color, ray_out, absorbed = scatter(mat(min_shape), ray, hitrecord, rng)
                if absorbed
                    return color
                else
                    return color .* raymarch(scene, ray_out, max_bounces_left - 1, max_steps_per_bounce, distance_tolerance, rng)
                end
            else
                # Assumes the ray direction is normalized
                t += current_min_dist
            end
        end
    end
    # max bounces <= 0 or ran out of march steps
    zero(Vec{3, T})
end

saveimage(savename::AbstractString, image) = save("./images/$(savename)", map(clamp01nan, image))


