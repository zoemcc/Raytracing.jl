using Raytracing
using GeometryBasics
using CoordinateTransformations
using ColorTypes
using ColorVectorSpace
using StaticArrays
using Images
using ImageIO
using FileIO
using ImageMagick
using Test

@testset "Raytracing.jl" begin
    T = Float64
    origin = Point3{T}(0,0,1.5)
    lookat = Point3{T}(0,0,0)
    vup = Vec3{T}(0,1,0)
    aspectratio = T(16/9)
    vfov = T(40)
    camera = Raytracing.make_camera(origin, lookat, vup, vfov, aspectratio)
    genray = Raytracing.generate_ray(camera, 0.5, 0.5)
    raydir = Vec3{T}(0,0,-1)
    @test Raytracing.origin(genray) ≈ origin
    @test Raytracing.direction(genray) ≈ raydir
    @test genray ≈ Raytracing.Ray(origin, raydir)
    @test ≈(genray, Raytracing.Ray(origin, raydir); atol=1e-2)

    origin2 = Point3{T}(0,0.5,1)
    camera2 = Raytracing.make_camera(origin2, lookat, vup, vfov, aspectratio)
    #@show dump(camera2)
    #@show Raytracing.generate_ray(camera2, 0.5, 0.5)
    #@show Raytracing.generate_ray(camera2, 0.0, 0.0)

    spheresdf = Raytracing.SphereSignedDistanceField{T}(0.5)
    @test spheresdf(origin) ≈ 1.0
    @test spheresdf(lookat) ≈ -0.5
    @test spheresdf(Point3{T}(0,0,0.5)) ≈ 0.0

    #@show Raytracing.normal(spheresdf, origin)
    #@show Raytracing.normalzygote(spheresdf, origin)
    @test Raytracing.normal(spheresdf, origin) ≈ Raytracing.normalzygote(spheresdf, origin)

    normalmat = Raytracing.NormalShading()

    id = IdentityTransformation()
    spherealbedo = RGB{T}(0.1, 0.8, 0.8)
    
    sphereshape = Raytracing.Shape(id, spheresdf, normalmat)


    skyboxspheresdf = Raytracing.SphereSignedDistanceField{T}(100)
    skyboxsphereshape = Raytracing.Shape(id, skyboxspheresdf, normalmat)

    sceneshapes = [skyboxsphereshape, sphereshape]
    #sceneshapes = @SVector [sphereshape]
    scene = Raytracing.Scene(sceneshapes)

    image_height = 360
    samples_per_pixel = 1
    max_bounces_per_ray = 1
    max_steps_per_bounce = 200
    distance_tolerance = 1e-4
    image = Raytracing.raytrace_image(scene, camera, image_height, samples_per_pixel, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance)

    tmpdir = mktempdir()
    saved_image_path = "$(tmpdir)/lofi_sphere_skybox.png"
    save(saved_image_path, image)
    image_reloaded = load(saved_image_path)

    reference_image = load("lofi_sphere_skybox.png")

    # Regression test
    @test image_reloaded ≈ reference_image


end
