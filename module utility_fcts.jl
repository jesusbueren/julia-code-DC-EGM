module utility_fcts
export u, u_prime,inv_u_prime,u_prime_prime,u_prime_prime_prime
    include("module model_parameter.jl")
    using .model_parameter
    function u(c::Float64,dₜ::Int64)
        if ρ==1.0
            if dₜ==1
                return log(max(c,ϵ_c))-δ
            else
                return log(max(c,ϵ_c))
            end
        else
            if dₜ==1
                return (max(c,ϵ_c)^(1.0-ρ)-1)/(1.0-ρ)-δ
            else
                return (max(c,ϵ_c)^(1.0-ρ)-1)^(1.0-ρ)/(1.0-ρ)
            end
        end
    end
    function u_prime(c::Float64)
            return c^(-ρ)
    end
    function inv_u_prime(c::Float64)
        return c^(-1.0/ρ)
    end
end