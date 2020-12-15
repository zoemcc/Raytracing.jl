struct Scene{T<:Tuple}
    shapes::T
end

shapes(scene::Scene) = scene.shapes

min_dist_index(scene::Scene, point::Point{3, T}) where {T<:Real} = findmin([abs(shape(point)) for shape in shapes(scene)])
