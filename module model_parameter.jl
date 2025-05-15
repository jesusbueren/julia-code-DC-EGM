module model_parameter
export A_grid,T,nkk,R,y,δ,ρ,β,ϵ_c
    const nkk=500 #number of points in the grid of savings
    const A_min=0.001 #minimum savings
    const A_max=400.0 # maximum savings
    const A_grid=collect(Float64,range(A_min,stop=A_max,length=nkk)) #define exogenous grid of savings
    const T=21 # time horizon
    const R=1.0 # interest rate
    const y=20.0 # income
    const δ=1.0 # disutility of work 
    const ρ=1.0 # risk-aversion
    const β=0.98 # discount factor
    const ϵ_c=1e-5 #mininum consumption to avoid numerical problems
end 