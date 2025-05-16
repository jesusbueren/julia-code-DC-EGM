
    
    function EGM_retired(pol_c:: Spline1D,V:: Spline1D,params::ModelParameter)
        # Unpack parameters
        nkk=params.nkk
        A_grid = params.A_grid
        R= params.R
        β= params.β

        c_j=zeros(nkk)
        v_j=zeros(nkk)
        m_grid=zeros(nkk)
        c_j=inv_u_prime.(β*R*u_prime.(pol_c(A_grid.*R),Ref(params)),Ref(params))
        v_j=u.(c_j,zeros(Int64,nkk),Ref(params))+β*V(A_grid.*R)
        m_grid=c_j.+A_grid
        pol_c_new=Spline1D(m_grid,c_j,k=1,bc="extrapolate")
        V_new=Spline1D(m_grid,v_j,k=1,bc="extrapolate")
        return pol_c_new,V_new,m_grid,v_j
    end
    
    function EGM_worker(pol_c_d::Vector{Any},V_d::Vector{Any},params::ModelParameter)
        # Unpack parameters
        nkk=params.nkk
        A_grid = params.A_grid
        R= params.R
        β= params.β
        y= params.y

        c_j=zeros(nkk)
        v_j=zeros(nkk)
        m_grid=zeros(nkk)
        for (j_l,a) in enumerate(A_grid)
            if V_d[2](R*a+y)>V_d[1](R*a+y) 
                c_j[j_l]=inv_u_prime(β*R*u_prime(pol_c_d[2](a*R+y),params),params)
                v_j[j_l]=u(c_j[j_l],1,params)+β*R*(V_d[2](a*R+y))
            else
                c_j[j_l]=inv_u_prime(β*R*u_prime(pol_c_d[1](a*R+y),params),params)
                v_j[j_l]=u(c_j[j_l],1,params)+β*(V_d[1](a*R+y))
            end
            m_grid[j_l]=c_j[j_l]+a
        end
        return c_j,v_j,m_grid
    end



    function upper_envelope(c::Vector{Float64}, v::Vector{Float64}, m_grid::Vector{Float64},params::ModelParameter)
        # Unpack parameters
        nkk=params.nkk
        ind = ones(Int, nkk)  # 1 means keep, 0 means discard
        ϵ_c=params.ϵ_c


        # Arrays to store extra M_star points and their interpolated values
        m_star_list = Float64[]
        c_star_list = Float64[]
        v_star_list = Float64[]

        limit_left=0
        limit_right=0
        eps = 1e-6

        # Helper: perform first-order Taylor approximation at `target` using closest to `bound`
        function taylor_approx(bound, target, direction)
            step = direction == :right ? 1 : -1
            closest = bound
            
            while 1 <= closest + step <= nkk &&
                abs(m_grid[target] - m_grid[closest]) >
                abs(m_grid[target] - m_grid[closest + step])
                closest += step
            end

            α1 = u_prime(c[closest],params)
            α2 = v[closest] - α1 * m_grid[closest]
            return α1 * m_grid[target] + α2
        end

        function trim!(i, direction::Symbol)
            bound_left = i
            bound_right = i + 1

            # Set the initial bound and update step based on direction
            bound = direction == :right ? bound_left : bound_right
            step  = direction == :right ? -1 : 1

            while 1 <= bound <= nkk
                # Arguments are flipped depending on direction
                v_approx = direction == :right ?
                    taylor_approx(bound_right, bound, :right) :
                    taylor_approx(bound_left, bound, :left)

                # Value comparison and deletion
                if v_approx > v[bound]
                    ind[bound] = 0
                    bound += step
                else
                    if direction == :right
                        limit_left = bound
                    else
                        limit_right = bound
                    end
                    break
                end
            end
        end

        # Linear extrapolation helper
        function linear_extrap(x1, y1, x2, y2, x_star)
            slope = (y2 - y1) / (x2 - x1)
            return y1 + slope * (x_star - x1)
        end


        # Main loop
        for i in 1:nkk-1
            if m_grid[i+1] < m_grid[i]

                #Trim right
                trim!(i,:right)
                #Trim left
                trim!(i,:left)

                if limit_left>0
                    #Compute intersection
                    α1_1=u_prime(c[limit_left],params)
                    α1_2=u_prime(c[limit_right],params)
                    α2_1=v[limit_left] - α1_1 * m_grid[limit_left]
                    α2_2=v[limit_right] - α1_2 * m_grid[limit_right]
                    M_star=(α2_2-α2_1)/(α1_1-α1_2)
                    v_star = α1_1 * M_star + α2_1

                    # Left extrapolation for c
                    if limit_left >= 2
                        c_L = linear_extrap(m_grid[limit_left-1], c[limit_left-1],
                                            m_grid[limit_left], c[limit_left],
                                            M_star - eps)
                    else
                        c_L = c[limit_left]  # fallback if no left point
                    end

                    # Right extrapolation for c
                    if limit_right + 1 <= nkk
                        c_R = linear_extrap(m_grid[limit_right], c[limit_right],
                                            m_grid[limit_right+1], c[limit_right+1],
                                            M_star + eps)
                    else
                        c_R = c[limit_right]  # fallback if no right point
                    end

                    # Store both sides
                    push!(m_star_list, M_star - eps)
                    push!(c_star_list, c_L)
                    push!(v_star_list, v_star)

                    push!(m_star_list, M_star + eps)
                    push!(c_star_list, c_R)
                    push!(v_star_list, v_star)
                end
                #println(M_star)
            end
        end

        # Rebuild splines
        kept_indices = ind .== 1
        m_vals = m_grid[kept_indices]
        c_vals = c[kept_indices]
        v_vals = v[kept_indices]

        # Combine original points with M_star jump points
        m_all = vcat(m_vals, m_star_list)
        c_all = vcat(c_vals, c_star_list)
        v_all = vcat(v_vals, v_star_list)

        # Sort by m
        perm = sortperm(m_all)
        m_all_sorted = m_all[perm]
        c_all_sorted = c_all[perm]
        v_all_sorted = v_all[perm]

        # Final splines
        pol_c_new = Spline1D(vcat(0.0, m_all_sorted), vcat(0.0, c_all_sorted), k=1, bc="extrapolate")
        V_new     = Spline1D(vcat(0.0,m_all_sorted), vcat(u(ϵ_c,1,params),v_all_sorted), k=1, bc="extrapolate")

        return pol_c_new, V_new
    end

