
using Dierckx, Plots, Debugger

include("module model_parameter.jl")
params=ModelParameter()

include("module utility_fcts.jl")
include("module EGM_algorithm.jl")


function policy_fct(params::ModelParameter)
    # Unpack parameters
    T=params.T
    nkk=params.nkk
    A_grid=params.A_grid

    cʳ=zeros(T,nkk) #consumption policy function for retired indv (fct of age, assets)
    Vʳ=zeros(T,nkk) #value function for retired indv (fct of age, assets)
    cʷ=zeros(T,nkk,2) #consumption policy function for workers (fct of age, assets, ret decision)
    Vʷ=zeros(T,nkk,2) #value policy function for workers (fct of age, assets, ret decision)
    M̄_d=zeros(T,nkk,2)
    M̄_r=zeros(nkk)
    v_r=zeros(nkk)

    #Policy and value fct of the retiree
    c_fct_r=Array{Any}(undef,T)
    v_fct_r=Array{Any}(undef,T)

    #Policy and value fct of the worker
    c_fct_d=Array{Any}(undef,T,2)
    v_fct_d=Array{Any}(undef,T,2)

    for t=T:-1:1 #backward induction
        println(t)
        for sₜ ∈ [0,1] #retired or not
            if sₜ==0  #retired  
                if t==T 
                    cʳ[t,:]=A_grid #consume everything at t=T
                    local d_v=Array{Int64}(undef, nkk)
                    d_v=0 #not work at t=T
                    Vʳ[t,:]=u.(A_grid,d_v,Ref(params))
                    c_fct_r[t]=Spline1D(A_grid,cʳ[t,:],k=1,bc="extrapolate") 
                    v_fct_r[t]=Spline1D(A_grid,Vʳ[t,:],k=1,bc="extrapolate")
                else
                    global M̄_r
                    global v_r
                    (c_fct_r[t],v_fct_r[t],M̄_r,v_r)=EGM_retired(c_fct_r[t+1],v_fct_r[t+1],params)
                end
            else #not retired
                for dₜ ∈ [0,1] #can decide to work or not
                    if t==T
                        cʷ[t,:,dₜ+1]=A_grid #consume everything at t=T
                        M̄_d[t,:,dₜ+1]=A_grid #consume everything at t=T
                        local d_v=Array{Int64}(undef, nkk)
                        if dₜ==0 
                            d_v.=0
                        else
                            d_v.=1
                        end
                        Vʷ[t,:,dₜ+1]=u.(A_grid,d_v,Ref(params)) 
                        c_fct_d[t,dₜ+1]=Spline1D(A_grid,cʷ[t,:,dₜ+1],k=1,bc="extrapolate")
                        v_fct_d[t,dₜ+1]=Spline1D(A_grid,Vʷ[t,:,dₜ+1],k=1,bc="extrapolate")
                    else
                        if dₜ==0
                            c_fct_d[t,dₜ+1]=c_fct_r[t]
                            v_fct_d[t,dₜ+1]=v_fct_r[t]
                            Vʷ[t,:,dₜ+1]=v_fct_d[t,dₜ+1].(A_grid)
                            cʷ[t,:,dₜ+1]=c_fct_d[t,dₜ+1].(A_grid)
                        else
                            (cʷ[t,:,dₜ+1],Vʷ[t,:,dₜ+1],M̄_d[t,:,dₜ+1])=EGM_worker(c_fct_d[t+1,:],v_fct_d[t+1,:],params)
                            #(c_fct_d[t,dₜ+1],v_fct_d[t,dₜ+1])=upper_envelope(cʷ[t,:,dₜ+1],Vʷ[t,:,dₜ+1],M̄_d[t,:,dₜ+1])
                            (c_fct_d[t,dₜ+1],v_fct_d[t,dₜ+1])=upper_envelope(cʷ[t,:,dₜ+1],Vʷ[t,:,dₜ+1],M̄_d[t,:,dₜ+1],params)
                        end
                    end
                end 
            end
        end
    end
    return c_fct_r,v_fct_r,c_fct_d,v_fct_d
end

c_fct_r,v_fct_r,c_fct_d,v_fct_d=policy_fct(params)

t=15
dₜ=1
scatter(params.A_grid,c_fct_d[t,dₜ+1].(params.A_grid),label="v_d",color=:black,markersize=2)
t=1
scatter!(params.A_grid,c_fct_d[t,dₜ+1].(params.A_grid),label="v_d",color=:red,markersize=2)










