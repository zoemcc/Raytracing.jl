abstract type AbstractSignedDistanceField{T<:Real} end

(sdf::AbstractSignedDistanceField{T})(::Point{3, T}) where {T<:Real} = error("Signed Distance Field for $(typeof(sdf)) has not been defined.")
normal(sdf::AbstractSignedDistanceField{T}, ::Point{3, T}) where {T<:Real} = error("Normal of Signed Distance Field for $(typeof(sdf)) has not been defined.")
normalzygote(sdf::AbstractSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(Zygote.gradient(sdf, point)[1]))
normalfinitediff(sdf::AbstractSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(FiniteDifferences.grad(central_fdm(5, 1), sdf, point)[1]))
normalforwarddiff(sdf::AbstractSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(ForwardDiff.gradient(sdf, point)))

struct SphereSignedDistanceField{T<:Real} <: AbstractSignedDistanceField{T}
    radius::T
end

radius(sdf::SphereSignedDistanceField) = sdf.radius
(sdf::SphereSignedDistanceField{T1})(point::Point{3, T2}) where {T1<:Real, T2<:Real} = 
    norm(point) - T2(radius(sdf))

# Currently always pointing out 
normal(::SphereSignedDistanceField{T1}, point::Point{3, T2}) where {T1<:Real, T2<:Real} = 
    Vec3{T2}(normalize(point))


struct MandelBulbSignedDistanceField{T<:Real} <: AbstractSignedDistanceField{T}
    numfractaliter::Int64
end

numfractaliter(sdf::MandelBulbSignedDistanceField) = sdf.numfractaliter
function (sdf::MandelBulbSignedDistanceField{T1})(point::Point{3, T2}) where {T1<:Real, T2<:Real} 
    px, py, pz = point
    dz = one(T2)
    w = copy(point)
    m = dot(w, w)
    for i in 1:numfractaliter(sdf)
        w2 = w .* w
        w4 = w2 .* w2
        x, y, z = w
        x2, y2, z2 = w2
        x4, y4, z4 = w4
        m2 = m .* m
        m4 = m2 .* m2

        dz = 8*sqrt(m4*m2*m)*dz + 1

        k3 = x2 + z2
        k2 = 1 / sqrt(k3*k3*k3*k3*k3*k3*k3)
        k1 = x4 + y4 + z4 - 6y2*z2 - 6x2*y2 + 2z2*x2
        k4 = x2 - y2 + z2

        newwx = px + 64x*y*z*k4*k1*k2 * (x2 - z2) * (x4 - 6x2*z2 + z4)
        newwy = py - 16y2*k3*k4*k4 + k1*k1
        newwz = pz - 8y*k4*k1*k2 * (x4*x4 - 28x4*x2*z2 + 70x4*z4 - 28x2*z2*z4 + z4*z4)

        w = Point{3, T2}(newwx, newwy, newwz)
        m = dot(w, w)
        if m > 256
            break
        end
    end

    T2(0.25) * log(m) * sqrt(m) / dz
end

normal(sdf::MandelBulbSignedDistanceField{T1}, point::Point{3, T2}) where {T1<:Real, T2<:Real} = 
    normalforwarddiff(sdf, point)




