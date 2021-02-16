struct Scene{T<:Tuple}
    shapes::T
end

shapes(scene::Scene) = scene.shapes

@inline function min_dist_index(scene::Scene, point::AbstractArray{T}) where {T<:Real} 
    curmin = T(Inf)
    index = 0
    for (i, shape) in enumerate(shapes(scene))
        dist = abs(shape(point))
        if dist < curmin
            index = i
            curmin = dist
        end
    end
    (curmin, index)
end



@inline function normalforwarddiff(scene::Scene, point::AbstractArray{T}) where {T<:Real}
    (curmin, index) = min_dist_index(scene, point)
    normalforwarddiff(scene.shapes[index].sdf, point)
end