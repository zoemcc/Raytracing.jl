struct Shape{T<:Real, SDF<:AbstractSignedDistanceField{T}, TRANSFORM<:Transformation}
    # Store the inverse pose since that's what we use to transform input points
    invpose::TRANSFORM 
    sdf::SDF
    albedo::RGB{T}
end    

invpose(shape::Shape) = shape.invpose
sdf(shape::Shape) = shape.sdf
albedo(shape::Shape) = shape.albedo
# Calculate the signed distance within the sdf's coordinate frame by using the inverse pose
(shape::Shape{T, SDF, TRANSFORM})(point::Point{3, T}) where {T<:Real, SDF<:AbstractSignedDistanceField{T}, TRANSFORM<:Transformation} = sdf(shape)(invpose(shape)(point))


