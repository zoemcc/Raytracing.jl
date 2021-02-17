begin
    using Raytracing

    using GeometryBasics
    using Random
    using CoordinateTransformations
    using ColorTypes
    using ColorVectorSpace
    using StaticArrays
    using Images
    using ImageIO
    using FileIO
    using ImageMagick
    using CUDA
    using KernelAbstractions
    using BenchmarkTools
    using GLMakie
    using Observables
    using AbstractPlotting
end



#function main()
begin
    T = Float32
    origin = Point3{T}(0,0,1.5)
    lookat = Point3{T}(0,0,0)
    vup = Point3{T}(0,1,0)
    aspectratio = T(16/9)
    vfov = T(40)
    camera = Raytracing.make_camera(origin, lookat, vup, vfov, aspectratio)

    spheresdf = Raytracing.SphereSignedDistanceField{T}(0.5)

    #mandelsdf = Raytracing.MandelBulbSignedDistanceField(50)
    #boundingmandelspheresdf = Raytracing.SphereSignedDistanceField{T}(1.25)
    #overshootfactor = T(1.5)
    #boundingmandelsdf = Raytracing.BoundingVolumeSignedDistanceField(overshootfactor, boundingmandelspheresdf, mandelsdf)

    point = Point3{T}(0, 0, 1)
    #mandelsdf(point)
    Raytracing.normalforwarddiff(spheresdf, point)

    normalmat = Raytracing.NormalShading()

    id = IdentityTransformation()
    #mandelalbedo = Vec{3, T}(0.8, 0.25, 0.25)
    spherealbedo = Point{3, T}(0.8, 0.25, 0.25)

    #lambertianmat = Raytracing.Lambertian(mandelalbedo)
    
    #mandelshape = Raytracing.Shape(id, mandelsdf, lambertianmat)
    sphereshape = Raytracing.Shape(id, spheresdf, normalmat)
    #mandelshape = Raytracing.Shape(id, mandelsdf, normalmat)


    skyboxspheresdf = Raytracing.SphereSignedDistanceField{T}(100)
    skyboxsphereshape = Raytracing.Shape(id, skyboxspheresdf, normalmat)

    sceneshapes = (skyboxsphereshape, sphereshape)
    scene = Raytracing.Scene(sceneshapes)

    image_height = 1080
    aspect_ratio = aspectratio
    image_width = Int(floor(aspectratio * image_height))
    samples_per_pixel = 3
    max_bounces_per_ray = 1
    max_steps_per_bounce = 400
    distance_tolerance = T(1e-2)
    compute_buffer_size = image_height * image_width
    rng = MersenneTwister(1234)
    rat = 2.0 / 0.70
    len = 1.2
    camera3 = Raytracing.make_camera(Point3{T}(0.0, rat * len, 1/rat * len), lookat, vup, vfov, aspectratio)
end
begin
    @time image = Raytracing.raytrace_image(scene, camera3, image_height, samples_per_pixel, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance, rng)

    savename = "lofi_sphere_normalshading.png"
    lofibig = Raytracing.supersample_crisp(image, 100)
    #Raytracing.saveimage("$(savename)", lofibig)


    cudaRayTracer = Raytracing.CudaRayTracer(aspect_ratio, image_height, samples_per_pixel, max_steps_per_bounce, distance_tolerance, compute_buffer_size)

    cudaRayTracer._array_cache.cube_coords 



    j = 5
    i = 5
    s = T((j - 1)/(image_width - 1))
    t = T((i - 1)/(image_height - 1))
    ray = Raytracing.generate_ray(camera, s, t)
    @time Raytracing.raymarch_iterative_simplecore(ray, scene, max_steps_per_bounce, distance_tolerance)

    Raytracing.get_color_from_hit()

    Raytracing.raytrace_image_cuda(scene, camera3, image_height, samples_per_pixel, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance, rng)


end




for i in 1:3
    for j in 4:6
        @show (i, j)
        if j - 2 == i
            break
        end
    end
end


begin
    out_image = zeros(RGB{T}, image_height, image_width)
    pixelxs = collect(Float32.((0.0:1/(image_width - 1):1.0)))
    pixelys = collect(Float32.((0.0:1/(image_height - 1):1.0)))
    combined_pixels = collect([(0.0f0, 0.0f0) for y in pixelys, x in pixelxs])
    pixelxs_cu = CuArray(pixelxs)
    pixelys_cu = CuArray(pixelys)
    combined_pixels_cu = CuArray(combined_pixels)

    shaded = collect([Raytracing.gray_to_RGB(T(0)) for y in pixelys, x in pixelxs])
    shaded_cu = CuArray(shaded)

    pixelxs

end

@btime begin
    kernel! = Raytracing.raytrace_from_pixel_s_t!(CPU(), 1)
    event = kernel!(shaded, pixelxs, pixelys, scene, camera, max_steps_per_bounce, distance_tolerance, ndrange=(image_height, image_width))
    wait(event)

end

shaded_cu_cpu = Array(shaded_cu)
GLMakie.image(shaded_cu_cpu)

@btime begin
    kernel_cu! = Raytracing.raytrace_from_pixel_s_t!(CUDADevice(), 256)
    event = kernel_cu!(shaded_cu, pixelxs_cu, pixelys_cu, scene, camera, max_steps_per_bounce, distance_tolerance, ndrange=(image_height, image_width))
    wait(event)
    CUDA.copyto!(shaded_cu_cpu, shaded_cu)

end


begin
    theme = AbstractPlotting.Attributes(show_axis=false, raw=false, scale_plot=true, 
                        padding=Point3(0.0f0, 0.0f0, 0.0f0), 
                        align= (:middle, :middle),
                        #resolution=(640, 360), #limits=FRect(0, 0, 640, 360),
                        panbutton=nothing, update_limits=true)
    AbstractPlotting.set_theme!(theme)
    shaded_cu_cpu = Array(shaded_cu)
    makie_scene = AbstractPlotting.image(collect(shaded_cu_cpu'); show_axis=false)
    makie_image = makie_scene.plot[:image]
    kernel_cu! = Raytracing.raytrace_from_pixel_s_t!(CUDADevice(), 256)
    display(makie_scene)
end
begin

    for i in 1:100
        @time begin 
            origin_move = Point3{T}(0,0,3.5 + 0.5 * sin((i / 10)))
            lookat_move = Point3{T}(0 + sin(i/10),0 + cos(i/10),0)
            vup_move = Point3{T}(0,1,0 + sin(i/100))
            aspectratio = T(16/9)
            vfov = T(40)
            camera_move = Raytracing.make_camera(origin_move, lookat_move, vup_move, vfov, aspectratio)

            event = kernel_cu!(shaded_cu, pixelxs_cu, pixelys_cu, scene, camera_move, max_steps_per_bounce, distance_tolerance, ndrange=(image_height, image_width))
            wait(event)
            CUDA.copyto!(shaded_cu_cpu, shaded_cu)
            makie_image[] .= shaded_cu_cpu'
            Observables.notify!(makie_image)
        end
    end
end
