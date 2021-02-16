
abstract type AbstractRayTracer end

struct CudaRayTracer{T<:Real, N<:Integer, C<:NamedTuple} <: AbstractRayTracer
    aspect_ratio::T
    image_height::N
    image_width::N
    samples_per_pixel::N
    max_steps_per_bounce::N
    distance_tolerance::T
    compute_buffer_size::N
    _array_cache::C
end

aspect_ratio(raytracer::CudaRayTracer) = raytracer.aspect_ratio
image_height(raytracer::CudaRayTracer) = raytracer.image_height
image_width(raytracer::CudaRayTracer) = raytracer.image_width
samples_per_pixel(raytracer::CudaRayTracer) = raytracer.samples_per_pixel
max_steps_per_bounce(raytracer::CudaRayTracer) = raytracer.max_steps_per_bounce
distance_tolerance(raytracer::CudaRayTracer) = raytracer.distance_tolerance
compute_buffer_size(raytracer::CudaRayTracer) = raytracer.compute_buffer_size
_array_cache(raytracer::CudaRayTracer) = raytracer._array_cache

function CudaRayTracer(aspect_ratio::T, image_height::N, samples_per_pixel::N, max_steps_per_bounce::N, distance_tolerance::T, compute_buffer_size::N) where {T<:Real, N<:Integer}
    image_width = N(floor(aspect_ratio * image_height))

    cube_width_axis = (0.0f0:1f0/(image_width - 1):1.0f0)
    cube_height_axis = (0.0f0:1f0/(image_height - 1):1.0f0)
    cube_coords = cu(cat([y for y in cube_height_axis, x in cube_width_axis], [x for y in cube_height_axis, x in cube_width_axis], dims=3))

    _array_cache = (image_accumulated = cu(zeros(Float32, image_height, image_width, 3)),
                    image_output = cu(zeros(Float32, image_height, image_width, 3)),
                    compute_buffer = cu(zeros(Float32, compute_buffer_size, 3)),
                    cube_coords = cube_coords)
    CudaRayTracer(aspect_ratio, image_height, image_width, samples_per_pixel, max_steps_per_bounce, distance_tolerance, compute_buffer_size, _array_cache)


end

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

function raymarch_iterative(scene::Scene, ray::Ray{T}, max_bounces::Integer, max_steps_per_bounce::Integer, distance_tolerance::T, rng::MersenneTwister) where {T<:Real}
    if max_bounces <= 0
        return zero(Vec{3, T})
    else
        accumulated_color = MVector(T(1.0), T(1.0), T(1.0))
        max_bounces_left = max_bounces
        while max_bounces_left > 0
            t = T(0)
            for i in 1:max_steps_per_bounce
                current_point = origin(ray) + t*direction(ray)
                current_min_dist, current_min_index = min_dist_index(scene, current_point)
                if current_min_dist < t * distance_tolerance
                    min_shape = shapes(scene)[current_min_index]
                    outward_normal = normal(sdf(min_shape), current_point)
                    hitrecord = HitRecord(ray, outward_normal, t)
                    color, ray_out, absorbed = scatter(mat(min_shape), ray, hitrecord, rng)
                    accumulated_color .*= color
                    if absorbed
                        return Vec{3, T}(accumulated_color) 
                    else
                        max_bounces_left -= 1
                        break
                    end
                else
                    # Assumes the ray direction is normalized
                    t += current_min_dist
                end
            end
        end
        return Vec{3, T}(accumulated_color)
    end
end

@inline function raymarch_iterative_simplecore(ray::Ray{T}, scene::Scene, max_steps_per_bounce::Integer, distance_tolerance::T) where {T<:Real}
    t = T(0)
    current_point = origin(ray)
    for i in 1:max_steps_per_bounce
        current_point = (t .* direction(ray)) .+ origin(ray)
        current_min_dist, current_min_index = min_dist_index(scene, current_point)
        if current_min_dist < t * distance_tolerance
            return (current_point, current_min_index)
        else
            # Assumes the ray direction is normalized
            t += current_min_dist
        end
    end
    return (current_point, 0)
end

function get_color_from_hit((point, index), ray::Ray{T}, scene::Scene) where {T<:Real}
    min_shape = shapes(scene)[index]
    outward_normal = normal(sdf(min_shape), point)
    hitrecord = HitRecord(ray, outward_normal, t)
    color, ray_out, absorbed = scatter(mat(min_shape), ray, hitrecord, rng)
    accumulated_color .*= color
    if absorbed
        return Vec{3, T}(accumulated_color) 
    else
        max_bounces_left -= 1
    end
end




@kernel function raytrace_from_pixel_s_t!(output_image::AbstractMatrix, pixelxs::AbstractVector, pixelys::AbstractVector, scene::Scene, camera::Camera{T}, max_steps_per_bounce::Integer, distance_tolerance::T)  where {T<:Real}
    i, j = @index(Global, NTuple)
    #ray[I] = Raytracing.Ray(orig, )
    @inbounds coord_y, coord_x = pixelys[i], pixelxs[j]
    #output_image[i, j] = (coord_y, coord_x)
    orig = Raytracing.origin(camera)
    llc = Raytracing.lower_left_corner(camera)
    horiz = Raytracing.horizontal(camera)
    vert = Raytracing.vertical(camera)



    raydir = normalize((coord_x .* horiz) .+ (coord_y .* vert) + (llc .- orig))
    ray = Raytracing.Ray{T}(orig, raydir)
    #ray = Raytracing.generate_ray(camera, coord_y, coord_x)
    #Raytracing.Ray(orig, normalize(llc .+ coord_x .* horiz .+ coord_y .* vert .- orig))
    #orig, orig)
    t = T(0)
    current_point = copy(orig)
    marched_point = copy(orig)
    shape_index = 0
    num_iters = 0
    for i in 1:max_steps_per_bounce
        current_point = (t .* raydir) .+ orig
        current_min_dist, current_min_index = min_dist_index(scene, current_point)
        #(marched_point, shape_index) = raymarch_iterative_simplecore(ray, scene, 50, 1e-3f0)
        if current_min_dist < t .* distance_tolerance
            (marched_point, shape_index, num_iters) = (current_point, current_min_index, i) # for shading for 
            #(marched_point, shape_index) = (current_point, current_min_index) # actual
            break
        else
            # Assumes the ray direction is normalized
            t += current_min_dist
        end
    end

    normal_vec = zeros(Point3{T})
    for k in 1:length(scene.shapes)
        normal_k = normalforwarddiff(scene.shapes[k].sdf, current_point)
        if shape_index == k
            normal_vec = normal_k
        end
    end
    rescaled_normal = (normal_vec .+ 1) ./ 2
    @inbounds output_image[i, j] = RGB{T}(rescaled_normal[1], rescaled_normal[2], rescaled_normal[3])
    #output_image[i, j] = gray_to_RGB(T(num_iters) / T(max_steps_per_bounce))
    #color = (((normal(hitrecord) .* (front_face(hitrecord) ? 1 : -1)) .+ 1) ./ 2)
    nothing
end



@inline function gray_to_RGB(gray::T) where {T<:Real}
    RGB(gray, gray, gray)
end



function raytrace_image_cuda(scene::Scene, cam::Camera{T}, image_height::Integer, samples_per_pixel::Integer, max_bounces_per_ray::Integer, max_steps_per_bounce::Integer, distance_tolerance::T, rng::MersenneTwister) where {T<:Real}
    aspect_ratio_cam = aspect_ratio(cam)
    image_width = Int64(floor(aspect_ratio_cam * image_height))
    out_image = zeros(RGB{T}, image_height, image_width)
    pixelxs = collect(Float32.((0.0:1/(image_width - 1):1.0)))
    pixelys = collect(Float32.((0.0:1/(image_height - 1):1.0)))
    combined_pixels = collect([(0.0f0, 0.0f0) for y in pixelys, x in pixelxs])
    pixelxs_cu = CuArray(pixelxs)
    pixelys_cu = CuArray(pixelys)
    combined_pixels_cu = CuArray(combined_pixels)
    shaded = collect([Raytracing.gray_to_RGB(T(0)) for y in pixelys, x in pixelxs])
    shaded_cu = CuArray(shaded)


    kernel! = Raytracing.raytrace_from_pixel_s_t!(CPU(), 1)
    @time begin
        event = kernel!(combined_pixels, pixelxs, pixelys, cam, ndrange=(image_height, image_width))
        wait(event)

    end



    pixels = [(y, x) for y in pixelys, x in pixelxs]
    pixels = CuArray([(y, x) for y in pixelys, x in pixelxs])

    orig = Raytracing.origin(cam)
    llc = Raytracing.lower_left_corner(cam)
    horiz = Raytracing.horizontal(cam)
    vert = Raytracing.vertical(cam)

    #@inline raymap((y, x)) = Raytracing.generate_ray(cam, x, y)
    raymap((y, x)) = Raytracing.Ray(orig, normalize(llc .+ x .* horiz .+ y .* vert .- orig))
    rays = map(raymap, pixels)

    @inline hitmap(ray) = Raytracing.raymarch_iterative_simplecore(ray, scene, max_steps_per_bounce, distance_tolerance)

    @time hits = hitmap.(rays)
    hit_indices = map(z->z[2], hits)


end

function raytrace(scene::Scene, cam::Camera{T}, raytracer::CudaRayTracer{T, N, C}) where {T<:Real, N<:Integer, C<:NamedTuple}
    aspect_ratio_tracer = aspect_ratio(raytracer)
    image_height = image_height(raytracer)
    image_width = image_width(raytracer)
    out_image = zeros(RGB{T}, image_height, image_width)


    #### buffers needed: output_image::RGB, cube_coords::(2), 

    orig = Raytracing.origin(cam)
    llc = Raytracing.lower_left_corner(cam)
    horiz = Raytracing.horizontal(cam)
    vert = Raytracing.vertical(cam)

    #@inline raymap((y, x)) = Raytracing.generate_ray(cam, x, y)
    raymap((y, x)) = Raytracing.Ray(orig, normalize(llc .+ x .* horiz .+ y .* vert .- orig))
    rays = map(raymap, pixels)

    @inline hitmap(ray) = Raytracing.raymarch_iterative_simplecore(ray, scene, max_steps_per_bounce, distance_tolerance)

    @time hits = hitmap.(rays)
    hit_indices = map(z->z[2], hits)


end




saveimage(savename::AbstractString, image) = save("./images/$(savename)", map(clamp01nan, image))


