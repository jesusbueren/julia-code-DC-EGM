struct ModelParameter
    nkk::Int
    A_min::Float64
    A_max::Float64
    A_grid::Vector{Float64}
    T::Int
    R::Float64
    y::Float64
    δ::Float64
    ρ::Float64
    β::Float64
    ϵ_c::Float64
end

function ModelParameter()
    nkk = 500
    A_min = 0.001
    A_max = 400.0
    A_grid = collect(range(A_min, stop=A_max, length=nkk))
    T = 21
    R = 1.0
    y = 20.0
    δ = 1.0
    ρ = 1.0
    β = 0.98
    ϵ_c = 1e-5

    return ModelParameter(nkk, A_min, A_max, A_grid, T, R, y, δ, ρ, β, ϵ_c)
end
