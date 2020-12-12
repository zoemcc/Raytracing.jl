
function random_unit_vector3(rng::MersenneTwister, T::DataType)
    normalize(Vec{3, T}(randn(rng, T, 3)))
end

function supersample_crisp(image::AbstractArray{C, 2}, upsize) where {C<:Color}
    bigimage = zeros(C, (size(image) .* upsize)...)
    height, width = size(image)
    for i in 1:height, j in 1:width
        bigimage[(i-1) * upsize + 1:i * upsize, (j-1) * upsize + 1:j * upsize] .= image[i, j]
    end
    bigimage
end

