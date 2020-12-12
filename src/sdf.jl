abstract type AbstractSignedDistanceField{T<:Real} end

(sdf::AbstractSignedDistanceField{T})(::Point{3, T}) where {T<:Real} = error("Signed Distance Field for $(typeof(sdf)) has not been defined.")
normal(sdf::AbstractSignedDistanceField{T}, ::Point{3, T}) where {T<:Real} = error("Normal of Signed Distance Field for $(typeof(sdf)) has not been defined.")

struct SphereSignedDistanceField{T<:Real} <: AbstractSignedDistanceField{T}
    radius::T
end

radius(sdf::SphereSignedDistanceField) = sdf.radius
(sdf::SphereSignedDistanceField{T})(point::Point{3, T}) where {T<:Real} = norm(point) - radius(sdf)

# Currently always pointing out 
normal(::SphereSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = Vec3{T}(normalize(point))
normalzygote(sdf::SphereSignedDistanceField{T}, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(gradient(sdf, point)[1]))


    

#=  old mandelbulb for reference

function mandelbulbdistfunc(numfractaliter::Integer) 
    function mandelbulbdist(p::Point{3, T}) where {T <: Real}
        px, py, pz = p[1], p[2], p[3]
        dz = 1.0
        w = Array{T, 1}(undef, 3)
        w .= 0

        stop_iter = 0

        w .+= p
        #@show typeof(w)
        m = normsq(w)
        #trap = vcat(abs.(w), m)

        for i = 1:numfractaliter

            w2 = w .* w
            w4 = w2 .* w2

            x, y, z = w[1], w[2], w[3]
            x2, y2, z2 = w2[1], w2[2], w2[3]
            x4, y4, z4 = w4[1], w4[2], w4[3]

            m2 = m * m
            m4 = m2 * m2

            dz *= 8.0 * sqrt(m4 * m2 * m)
            dz += 1.0

            k3 = x2 + z2
            k2 = inversesqrt(k3 * k3 * k3 * k3 * k3 * k3 * k3)
            k1 = x4 + y4 + z4 - 
                6.0 * y2 * z2 - 
                6.0 * x2 * y2 + 
                2.0 * z2 * x2
            k4 = x2 - y2 + z2

            neww = [
                px + 
                64.0 * x * y * z * (x2 - z2) * k4 * 
                (x4 - 6.0 * x2 * z2 + z4) * k1 * k2,
                py +
                -16.0 * y2 * k3 * k4 * k4 + k1 * k1,
                pz +
                -8.0 * y * k4 * (x4 * x4 - 
                28.0 * x4 * x2 * z2 +
                70.0 * x4 * z4 -
                28.0 * x2 * z2 * z4 +
                z4 * z4) * k1 * k2
                ]
            w[:] = neww[:]

            #trap[:, :, :] = min.(trap, vcat(abs.(w), m)) .* (1.0 .- stop_iter) .+ stop_iter .* trap
            m = normsq(w)
            if m > 256.0
                return -0.0005
            end
            #m[:, :, :] = min.(normsq(w), 256.0)
            #m[map(isnan, m)] .= 256.0
            #dz[map(isnan, dz)] .= 1.0

            #@show i
            #@show w
            #@show m
            #@show stop_iter
            #@show dz
            stop_iter = i

        end

        dist = 0.25 * log(m) * sqrt(m) / dz
        return dist
    end
    return mandelbulbdist
end

=#