
function random_unit_vector3(rng::MersenneTwister, T::DataType)
    normalize(Vec{3, T}(randn(rng, T, 3)))
end

