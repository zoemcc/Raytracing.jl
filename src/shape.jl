struct Shape{Transform<:Transformation, SDF<:AbstractSignedDistanceField, Mat<:AbstractMaterial}
    # Store the inverse pose since that's what we use to transform input points
    invpose::Transform 
    sdf::SDF
    mat::Mat
end    

@inline invpose(shape::Shape) = shape.invpose
@inline sdf(shape::Shape) = shape.sdf
@inline mat(shape::Shape) = shape.mat
# Calculate the signed distance within the sdf's coordinate frame by using the inverse pose
@inline function (shape::Shape{Transform, SDF, Mat})(point::Point{3, T}) where 
        {T<:Real, Transform<:Transformation, SDF<:AbstractSignedDistanceField, Mat<:AbstractMaterial} 
    sdf(shape)(invpose(shape)(point))
end


