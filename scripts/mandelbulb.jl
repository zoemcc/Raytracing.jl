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



function main()
    T = Float64
    origin = Point3{T}(0,0,1.5)
    lookat = Point3{T}(0,0,0)
    vup = Vec3{T}(0,1,0)
    aspectratio = T(16/9)
    vfov = T(40)
    camera = Raytracing.make_camera(origin, lookat, vup, vfov, aspectratio)

    spheresdf = Raytracing.SphereSignedDistanceField{T}(0.5)

    mandelsdf = Raytracing.MandelBulbSignedDistanceField{T}(3)

    point = Point3{T}(0, 0, 1)
    mandelsdf(point)
    Raytracing.normalforwarddiff(mandelsdf, point)

    normalmat = Raytracing.NormalShading()

    id = IdentityTransformation()
    mandelalbedo = Vec{3, T}(0.8, 0.25, 0.25)

    #lambertianmat = Raytracing.Lambertian(mandelalbedo)
    
    #mandelshape = Raytracing.Shape(id, mandelsdf, lambertianmat)
    mandelshape = Raytracing.Shape(id, mandelsdf, normalmat)


    skyboxspheresdf = Raytracing.SphereSignedDistanceField{T}(100)
    skyboxsphereshape = Raytracing.Shape(id, skyboxspheresdf, normalmat)

    sceneshapes = [skyboxsphereshape, mandelshape]
    scene = Raytracing.Scene(sceneshapes)

    image_height = 360
    samples_per_pixel = 3
    max_bounces_per_ray = 1
    max_steps_per_bounce = 400
    distance_tolerance = 1e-4
    rng = MersenneTwister(1234)
    rat = 2.0 / 0.70
    len = 1.2
    camera3 = Raytracing.make_camera(Point3{T}(0.0, rat * len, 1/rat * len), lookat, vup, vfov, aspectratio)
    @time image = Raytracing.raytrace_image(scene, camera3, image_height, samples_per_pixel, max_bounces_per_ray, max_steps_per_bounce, distance_tolerance, rng)

    savename = "me_mandelbulb_3_normalshading.png"
    Raytracing.saveimage("$(savename)", image)


end
