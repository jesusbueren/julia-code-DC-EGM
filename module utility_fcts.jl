
    function u(c::Float64,dₜ::Int64, params::ModelParameter)
        ϵ_c=params.ϵ_c
        δ=params.δ
        ρ=params.ρ    
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
    function u_prime(c::Float64,params::ModelParameter)
        ρ=params.ρ
            return c^(-ρ)
    end
    function inv_u_prime(c::Float64,params::ModelParameter)
        ρ=params.ρ
        return c^(-1.0/ρ)
    end