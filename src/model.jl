using LinearAlgebra
using Graphs
using JuMP
# using GLPK
# using CPLEX
using HiGHS

# include("data.jl")

# auxiliary stuff
SPINLOCK = Threads.SpinLock()

# custom types

CapacityConstraintRef = JuMP.Containers.DenseAxisArray#{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 1, Tuple{Vector{Tuple{Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}}}

# include("toy.jl")


"""fills used_arcs with arcs used by new_lambda_vertices and returns the cost"""
function get_used_arcs_and_cost(new_lambda_vertices::Vector{Int64}, used_arcs::Vector{Tuple{Int64, Int64}}, arc_to_cost_k::Dict{Tuple{Int64, Int64}, Float64})::Float64

    # get set of used arcs and cost
    # used_arcs = Tuple{Int64, Int64}[]
    lambda_cost::Float64 = 0.0
    path_length::Int64 = length(new_lambda_vertices)

    for i in 1:path_length-1
        j = i+1
        a = (new_lambda_vertices[i], new_lambda_vertices[j])
        lambda_cost += arc_to_cost_k[a] 

        push!(used_arcs, a)
    end
    return lambda_cost
end

# capacity_con::JuMP.Containers.DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 1, Tuple{Vector{Tuple{Int64, Int64}}}, Tuple{JuMP.Containers._AxisLookup{Dict{Tuple{Int64, Int64}, Int64}}}}

"""adds the lambda to the rmlp and ref array"""
function add_lambda_to_rmlp(rmlp::GenericModel{Float64}, k::Int64, lambda_cost::Float64, new_lambda_vertices::Vector{Int64}, used_arcs::Vector{Tuple{Int64, Int64}}, S::Array{Array{Int64}}, T::Array{Array{Int64}}, lambdas_ref::Vector{Vector{VariableRef}}, capacity_con, last_used_col::Vector{Vector{Int64}})::VariableRef
        
    # create variable in master
    new_lambda_ref::VariableRef = @variable(rmlp, base_name="λ$(k)_$(length(lambdas_ref[k])+1)", lower_bound=0.0)

    # set coefficient in capacity constraints
    for a in used_arcs
        set_normalized_coefficient(capacity_con[a], new_lambda_ref, 1.0) # set new column coeff. for each capacity constraint
    end

    # set coefficient in demand constraints
    if new_lambda_vertices[1] in S[k]
        set_normalized_coefficient(constraint_by_name(rmlp, "source_k$(k)_v$(new_lambda_vertices[1])"), new_lambda_ref, 1.0)
    end
    if new_lambda_vertices[end] in T[k]
        set_normalized_coefficient(constraint_by_name(rmlp, "sink_k$(k)_v$(new_lambda_vertices[end])"), new_lambda_ref, 1.0) 
    end

    set_objective_coefficient(rmlp, new_lambda_ref, lambda_cost) # set new column cost on rmlp

    # save variable ref
    push!(lambdas_ref[k], new_lambda_ref)

    # register tracker for new var
    push!(last_used_col[k], 0)

    return new_lambda_ref
end

function update_lambdas_in_use(last_used_col::Vector{Vector{Int64}}, l_val::Vector{Vector{Float64}}, lambdas_ref::Vector{Vector{VariableRef}}, max_counter::Int64)::Int64
    
    # is lambda not being used?
    is_zero::Vector{BitVector} = [lvs .< 1e-4 for lvs in l_val]
    
    # increment not used's counter 
    last_used_col .+= is_zero

    # get indices of used to reset counter
    in_use_list = findall.(x -> !(x), is_zero)

    # reset counter of used
    for (k, in_use_list_k) in enumerate(in_use_list)
        last_used_col[k][in_use_list_k] .= 0
    end

    # lambdas that will not remain in the rmlp
    newly_cleaned::Int64 = 0

    to_clean_list = findall.(x -> x > max_counter, last_used_col)
    for (k, list_k) in enumerate(to_clean_list)
        for i in list_k
            fix(lambdas_ref[k][i], 0.0, force=true)
            newly_cleaned += 1
        end
    end

    return newly_cleaned
end

function mnfp_lp(instance::MnfpData; optimizer::Union{DataType, Nothing}=nothing)::Tuple{GenericModel{Float64},Dict{Any, Any},Float64}

    if isnothing(optimizer)
        optimizer = HiGHS.Optimizer
    end

    # lp::Model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_THREADS" => Threads.nthreads()))
    lp::Model = Model(optimizer)
    set_silent(lp)

    K = instance.K
    vertex_amount = instance.vertex_amount
    A = instance.A
    u = instance.u
    d = instance.d
    c = instance.c
    

    # x variables
    x = Vector{VariableRef}[VariableRef[] for k in K]
    for k in K
        for a in instance.A
            new_x::VariableRef = @variable(lp, base_name="x_k$(k)_a$(a)", lower_bound=0.0)
            push!(x[k], new_x)
        end
    end

    for v in 1:vertex_amount
        for k in K
            @constraint(lp, 
                sum([x[k][n] for (n, a) in enumerate(A) if a[2] == v]) - 
                sum([x[k][n] for (n, a) in enumerate(A) if a[1] == v]) == d[k, v],
                base_name="flow_con_k$(k)_v$(v)"
            )
        end
    end

    for (n, a) in enumerate(A)
        @constraint(lp, 
            sum([x[k][n] for k in K]) <= u[n],
            base_name="capacity_con_a$(a)"
        )
    end

    @objective(lp, Min, sum([c[k,n]*x[k][n] for k in K, (n, a) in enumerate(A)]))

    optimize!(lp)

    x_val = [value.(x[k]) for k in K]

    solution = Dict()
    for k in K
        flow = Dict{Tuple{Int64, Int64}, Float64}()
        for (n, a) in enumerate(A)
            val = x_val[k][n]
            if val > 1e-6
                flow[a] = val
            end
        end
        solution[k] = flow
    end

    for k in K
        println("\ncommodity $(k):")
        println("x_sol: ", solution[k], "\n")
    end


    return lp, solution, objective_value(lp)
end


# MnfpData(vertex_amount::Int64, A::Vector{Tuple{Int64, Int64}}, u::Matrix{Float64}, c::Matrix{Float64}, d::Matrix{Float64})
# function mnfp_cga(
#     vertex_amount::Int64, A::Vector{Tuple{Int64, Int64}}, u::Matrix{Float64}, c::Matrix{Float64}, d::Matrix{Float64}; 
#     quiet::Bool=false,
#     verbose::Bool=false,
#     full_output::Bool=false,
#     max_time::Int64=1200,
#     optimizer::Union{DataType,
#     Nothing}=nothing)::Union{
#         Tuple{Model, Vector{Dict{Tuple{Int64,Int64}, Float64}}, Float64, Bool},
#         Tuple{Model, Vector{Dict{Tuple{Int64,Int64}, Float64}}, Float64, Bool, Bool, Bool, Bool, Vector{Float64}}
#     }

#     return mnfp_cga()
# end

function mnfp_cga(
    instance::MnfpData; 
    quiet::Bool=false,
    verbose::Bool=false,
    max_time::Int64=1200,
    it_print_interval::Int64=1,
    max_unused_iterations::Int64=10,
    optimizer::Union{DataType, Nothing}=nothing)::@NamedTuple{
        rmlp::Model,
        solution::Vector{Dict{Tuple{Int64,Int64}, Float64}},
        rmlp_obj::Float64,
        
        opt_failed::Bool,
        art_vars_in_sol::Bool,
        globally_optimal::Bool,
        hit_time_limit::Bool,
        
        pi_bar::Vector{Float64},
        v_bar::Vector{Vector{Float64}},
        n_bar::Vector{Vector{Float64}},

        elapsed_time::Float64,
        heur_time::Float64,
        rmlp_time::Float64,
        subp_time::Float64,
        
        lambdas_ref::Vector{Vector{VariableRef}},
        capacity_con::CapacityConstraintRef,
        sources_con::Vector{Vector{ConstraintRef}},
        sinks_con::Vector{Vector{ConstraintRef}},

        cols_gen::Int64,
        available_cols::Int64
    }

    # M::Float64 = 1e6
    M::Float64 = 1e9
    max_iter::Int64=1e9

    ################################################ read instance

    K::UnitRange{Int64} = instance.K
    V::UnitRange{Int64} = instance.V
    vertex_amount::Int64 = instance.vertex_amount
    s::Int64 = instance.s
    t::Int64 = instance.t
    V_st::Vector{Int64} = instance.V_st
    A::Vector{Tuple{Int64, Int64}} = instance.A
    arc_amount::Int64 = instance.arc_amount
    u::Matrix{Float64} = instance.u
    c::Matrix{Float64} = instance.c
    arc_to_cost::Vector{Dict{Tuple{Int64, Int64}, Float64}} = instance.arc_to_cost
    d::Matrix{Float64} = instance.d
    S::Vector{Array{Int64}} = instance.S
    T::Vector{Array{Int64}} = instance.T
    A_s::Vector{Vector{Tuple{Int64, Int64}}} = instance.A_s
    A_t::Vector{Vector{Tuple{Int64, Int64}}} = instance.A_t
    base_graph::SimpleDiGraph = instance.base_graph
    sub_graphs::Vector{SimpleDiGraph{Int64}} = instance.sub_graphs

    ################################################### make graphs auxiliary structures and take base metrics

    sub_cost_mod::Vector{Matrix{Float64}} = Matrix{Float64}[zeros(t,t) for k in K]

    done::Vector{Bool} = Bool[false for k in K]

    # auxiliary data structure
    lambdas_ref::Vector{Vector{VariableRef}} = Array{VariableRef}[VariableRef[] for k in K]
    lambdas_used_arcs::Vector{Vector{Vector{Tuple{Int64,Int64}}}} = Array{Array{Tuple{Int64,Int64}}}[Array{Tuple{Int64,Int64}}[] for k in K]
    lambdas_vertices::Vector{Vector{Vector{Int64}}} = Array{Array{Int64}}[Vector{Vector{Int64}}[] for k in K]
    lambdas_cost::Vector{Vector{Float64}} = [Float64[] for k in K]

    heur_time::Float64 = time()  
    ################################################### Make inicial columns

    # infinity::Vector{Float64} = Float64[Inf for k in K]

    ########################################### set infinity to diameter(max_g)

    # # build max g
    # max_g = SimpleGraph(vertex_amount)
    # for (i, j) in A
    #     add_edge!(max_g, i, j)
    # end

    # max_g_dists = zeros(vertex_amount, vertex_amount)
    # for (n, (i, j)) in enumerate(A)
    #     max_g_dists[i,j] = maximum(c[:,n])
    # end
    # for (n, (j, i)) in enumerate(A)
    #     m = max(max_g_dists[i,j], max_g_dists[j,i])
    #     max_g_dists[i,j] = m
    #     max_g_dists[j,i] = m
    # end

    # # effectively infinity for base path calculations
    # max_g_diameter = diameter(max_g, max_g_dists)
    # for k in K
    #     infinity[k] = max_g_diameter
    # end

    ########################################### set infinity to diameter(base_graph)

    # for k in K
    #     base_costs::Matrix{Float64} = zeros(t,t)

    #     # original arcs cost
    #     for (n, (i, j)) in enumerate(A)
    #         # add_edge!(g, i, j)
    #         base_costs[i,j] = c[k,n]
    #     end

    #     infinity[k] = diameter(base_graph, base_costs)
    # end

    ###########################################

    # # !quiet && println("infinity: $(infinity)")

    biased_K::Vector{Int64} = sort(K, by = x -> (length(S[x])^4)*length(T[x]))
    
    dealers_amount::Int64 = sum([length(S[k])+length(T[k]) for k in K])

    max_initial::Int64 = dealers_amount
    # max_initial::Int64 = 99999999999
    # max_initial::Int64 = -1 
    current_initial::Int64 = 1

    Threads.@threads for k in biased_K 

        if current_initial >= max_initial
            break
        end

        # take measures
        for (n, (i, j)) in enumerate(A)
            sub_cost_mod[k][i,j] = c[k,n]
        end
        
        ################### any s -> every t
        # st_path_state = dijkstra_shortest_paths(sub_graphs[k], s, sub_cost_mod[k])
        # # st_path_state = dijkstra_shortest_paths(sub_graphs[k], s, sub_cost_mod[k]; maxdist=infinity[k])

        # st_paths = enumerate_paths(st_path_state, T[k])
        # for (n, v) in enumerate(T[k])
        #     new_lambda_vertices = st_paths[n][2:end] # cut beginning (source). The end is not t
        #     # !quiet && println("$(k), $(n): $(st_paths[n])")

        #     used_arcs::Array{Tuple{Int64, Int64}} = Tuple{Int64, Int64}[]
        #     lambda_cost::Float64 = get_used_arcs_and_cost(new_lambda_vertices, used_arcs, arc_to_cost[k])

        #     push!(lambdas_used_arcs[k], used_arcs)
        #     push!(lambdas_cost[k], lambda_cost)
        #     push!(lambdas_vertices[k], new_lambda_vertices)
        # end

        ################### every s -> every t
        # for (m, x) in enumerate(S[k])
        for m in 1:length(S[k])
            if current_initial >= max_initial
                break
            end
            
            x::Int64 = S[k][m]

            # st_path_state = dijkstra_shortest_paths(sub_graphs[k], x, sub_cost_mod[k])
            st_path_state = bellman_ford_shortest_paths(sub_graphs[k], x, sub_cost_mod[k])


            # st_path_state = dijkstra_shortest_paths(sub_graphs[k], s, sub_cost_mod[k]; maxdist=infinity[k])

            st_paths::Vector{Vector{Int64}} = enumerate_paths(st_path_state, T[k])
            for (n, v) in enumerate(T[k])

                new_lambda_vertices::Vector{Int64} = st_paths[n]
                # !quiet && println("$(k), ($(x), $(v)): $(st_paths[n])")

                verbose && println("adding path between $(x) and $(v) for commodity $(k): ", new_lambda_vertices)
                if new_lambda_vertices == []
                    println("warning: no path between $(x) and $(v) for commodity $(k)")
                    continue
                end

                used_arcs::Array{Tuple{Int64, Int64}} = Tuple{Int64, Int64}[]
                lambda_cost::Float64 = get_used_arcs_and_cost(new_lambda_vertices, used_arcs, arc_to_cost[k])

                Threads.lock(SPINLOCK) do 
                    push!(lambdas_used_arcs[k], used_arcs)
                    push!(lambdas_cost[k], lambda_cost)
                    push!(lambdas_vertices[k], new_lambda_vertices)
    
                    current_initial += 1
                end
                if current_initial >= max_initial
                    break
                end

            end            
        end

        # !quiet && println("initial arcs for $(k): $(lambdas_used_arcs[k])")
    end

    heur_time = time() - heur_time  

    #####################################################

    ## master
    if isnothing(optimizer)
        optimizer = HiGHS.Optimizer
    end
    # rmlp::Model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPX_PARAM_THREADS" => Threads.nthreads()))
    rmlp::Model = Model(optimizer)
    set_silent(rmlp)

    capacity_con::CapacityConstraintRef = @constraint(rmlp, capacity[a in A], 0.0 <= 0.0)
    for (n,a) in enumerate(A)
        set_normalized_rhs(capacity_con[a], u[n]) # constante do lado direito
    end

    artificial_vars::Vector{Vector{VariableRef}} = Vector{VariableRef}[VariableRef[] for k in K]
    sources::Vector{Vector{ConstraintRef}} = Vector{ConstraintRef}[ConstraintRef[] for k in K]
    sinks::Vector{Vector{ConstraintRef}} = Vector{ConstraintRef}[ConstraintRef[] for k in K]

    for k in K
        
        # source constraints and artificial vars
        for (n, v) in enumerate(S[k])
            new_cons::ConstraintRef = @constraint(rmlp, 0 == 0, base_name="source_k$(k)_v$(v)")
            new_var::VariableRef = @variable(rmlp, base_name="art_k$(k)_v$(v)", lower_bound=0.0) 
            
            set_normalized_coefficient(new_cons, new_var, 1)
            set_normalized_rhs(new_cons, abs(d[k, v]))
            
            push!(sources[k], new_cons)
            push!(artificial_vars[k], new_var)
        end

        # sink constraints and artificial vars
        for (n, v) in enumerate(T[k])
            new_cons::ConstraintRef = @constraint(rmlp, 0 == 0, base_name="sink_k$(k)_v$(v)")
            new_var::VariableRef = @variable(rmlp, base_name="art_k$(k)_v$(v)", lower_bound=0.0) 
            
            set_normalized_coefficient(new_cons, new_var, 1)
            set_normalized_rhs(new_cons, d[k, v])
            
            push!(sinks[k], new_cons)
            push!(artificial_vars[k], new_var)
        end
    end

    @objective(rmlp, Min, sum([sum(M*artificial_vars[k]) for k in K]))

    # tracker of how many iterations since a lambda was used
    last_used_col::Vector{Vector{Int64}} = Vector{Int64}[Int64[] for k in K]
    
    # add initial paths
    cols_gen::Int64 = 0
    for k in K
        for (used_arcs, lambda_cost, new_lambda_vertices) in zip(lambdas_used_arcs[k], lambdas_cost[k], lambdas_vertices[k])
            
            verbose && println("adding initial lambda to rmlp for commodity $(k) with cost $(lambda_cost) and vertices $(new_lambda_vertices)")
            add_lambda_to_rmlp(rmlp, k, lambda_cost, new_lambda_vertices, used_arcs, S, T, lambdas_ref, capacity_con, last_used_col)
            cols_gen += 1
        end
    end
    initial_paths_amount::Int64 = sum([length(lambdas_vertices[k]) for k in K])
    !quiet && println("initial paths amount: $(initial_paths_amount)")

    !quiet && verbose && println("rmlp:\n", rmlp)

    !quiet && println("starting CGA")

    rmlp_obj::Float64 = Inf
    solution::Vector{Dict{Tuple{Int64,Int64}, Float64}} = Vector{Dict{Tuple{Int64,Int64}, Float64}}[]

    # CGA:
    # max_iter = 30

    mod_const::Int64 = it_print_interval
    if verbose
        mod_const = 1
    end

    globally_optimal::Bool = false
    hit_time_limit::Bool = false

    elapsed_time::Float64 = 0.0
    
    rmlp_time::Float64 = 0.0
    subp_time::Float64 = 0.0
    
    subp_time_aux::Float64 = 0.0
    
    cleaned_lambdas::Int64 = 0
    
    l_val::Vector{Vector{Float64}} = Vector{Float64}[Float64[0.0 for _ in lambdas_ref[k]] for k in K]
    pi_bar::Vector{Float64} = Float64[0.0 for _ in A]
    v_bar::Vector{Vector{Float64}} = Vector{Float64}[Float64[0.0 for _ in sinks[k]] for k in K]
    n_bar::Vector{Vector{Float64}} = Vector{Float64}[Float64[0.0 for _ in sources[k]] for k in K]
    art_vals::Vector{Vector{Float64}} = Vector{Float64}[Float64[0.0 for _ in k] for k in artificial_vars]
    art_sum::Float64 = 0.0

    elapsed_time_start::Float64 = time()

    for it_num in 1:max_iter
        
        !quiet && verbose && println("######################################## ITERATION $(it_num)")
        
        # optimize master
        optimize!(rmlp)
        optimal = termination_status(rmlp) in [MOI.OPTIMAL, MOI.ALMOST_INFEASIBLE, MOI.ALMOST_DUAL_INFEASIBLE]
        rmlp_time += solve_time(rmlp)
        rmlp_obj = objective_value(rmlp)

        elapsed_time = time() - elapsed_time_start
        if elapsed_time > max_time
            !quiet && println("reached time limit")
            hit_time_limit = true
            break
        end

        # check if it is feasible
        if !(optimal)
            !quiet && println("rmlp infeasible")
            break
        end

        ############################################# extract model data
        aux::Float64 = time()

        # get solution lambdas
        l_val = Vector{Float64}[value.(lambdas_ref[k]) for k in K]

        # get artificial values
        art_vals = Vector{Float64}[value.(k) for k in artificial_vars]

        # get duals
        pi_bar = dual.(capacity_con).data
        v_bar = Vector{Float64}[dual.(sinks[k]) for k in K]
        n_bar = Vector{Float64}[dual.(sources[k]) for k in K]

        #############################################


        # clean lambdas
        # println(last_used_col)
        # println(lambdas_ref)
        cleaned_lambdas = update_lambdas_in_use(last_used_col, l_val, lambdas_ref, max_unused_iterations)
        
        # auxiliary variables for subproblem solving
        s_red_cost::Vector{Vector{Float64}} = .- n_bar 
        t_red_cost::Vector{Vector{Float64}} = .- v_bar
        arc_red_cost::Matrix{Float64} = c .- transpose(pi_bar) 
        
        # function get_lambda_reduced_cost(lambdas_ref, s_red_cost, t_red_cost, arc_red_cost, lambdas_used_arcs)
        
        # end
        
        if !quiet && verbose
            println("\n\n")        
            println(rmlp)
            println(artificial_vars)
            println(art_vals)
            for k in K
                println("lambdas $(k): ", value.(lambdas_ref[k]))
                println("lvertices $(k): ", lambdas_vertices[k])
            end
            println("pi_bar: $(pi_bar)")
            println("v_bar: $(v_bar)")
            println("n_bar: $(n_bar)")
        end

        if !quiet && mod(it_num, mod_const) == 0
            println("aux time: ", time()-aux)
            println(
                "it_num: ", it_num,
                " | elapsed_time: ", elapsed_time,
                " | rmlp_time: ", rmlp_time,
                " | subp_time: ", subp_time,
                " | rmlp_obj: ", rmlp_obj,
                " | cols_gen: ", cols_gen,
                " | available_cols: ", cols_gen-cleaned_lambdas,
            )
        end

        subp_time_aux = time()

        # for (k, g) in enumerate(sub_graphs)
        Threads.@threads for k in K
            g::SimpleDiGraph = sub_graphs[k] 

            # c_k::Vector{Float64} = c[k,:]
            # A_cost::Vector{Float64} = c_k - pi_bar

            # A_s_k_cost::Vector{Float64} = -n_bar[k]
            # A_t_k_cost::Vector{Float64} = -v_bar[k]

            # original arcs reduced cost
            for (n, (i, j)) in enumerate(A)
                # sub_cost_mod[k][i,j] = A_cost[n]
                sub_cost_mod[k][i,j] = arc_red_cost[k,n]
            end

            # source arcs reduced cost
            for (n, (i, j)) in enumerate(A_s[k])
                # sub_cost_mod[k][i,j] = A_s_k_cost[n]
                sub_cost_mod[k][i,j] = s_red_cost[k][n]
            end

            # sink arcs reduced cost
            for (n, (i, j)) in enumerate(A_t[k])
                # sub_cost_mod[k][i,j] = A_t_k_cost[n]
                sub_cost_mod[k][i,j] = t_red_cost[k][n]
            end

            # st_path_state = dijkstra_shortest_paths(g, s, sub_cost_mod[k])
            st_path_state = bellman_ford_shortest_paths(g, s, sub_cost_mod[k])

            done[k] = st_path_state.dists[t] > 0 - 1e-4
        
            if !(done[k])

                # vertices being visited by the new lambda
                new_lambda_vertices::Vector{Int64} = enumerate_paths(st_path_state, t)[2:end-1] # cut beginning and end (source and sink)
                push!(lambdas_vertices[k], new_lambda_vertices)

                !quiet && verbose && println("new lambda$(k): λ$(k)_$(length(lambdas_ref[k])+1) = ", new_lambda_vertices, " subobj: $(st_path_state.dists[t])")
                            
                used_arcs::Array{Tuple{Int64, Int64}} = Tuple{Int64, Int64}[]
                lambda_cost::Float64 = get_used_arcs_and_cost(new_lambda_vertices, used_arcs, arc_to_cost[k])
                
                push!(lambdas_used_arcs[k], used_arcs)
                push!(lambdas_cost[k], lambda_cost)

                Threads.lock(SPINLOCK) do 
                    add_lambda_to_rmlp(rmlp, k, lambda_cost, new_lambda_vertices, used_arcs, S, T, lambdas_ref, capacity_con, last_used_col)
                    cols_gen += 1
                end
            
            else
                !quiet && verbose && println("g $(k) did not generate a new column")
            end
        
        end

        subp_time += time() - subp_time_aux

        # se todos os subproblemas acabaram, entao terminou
        if all(done)
            !quiet && println("done")
            elapsed_time = time() - elapsed_time_start
            globally_optimal = true
            break
        end
    end

    art_sum= sum([sum(vals) for vals in art_vals])
    
    # get solution
    for k in K
        x_sol::Dict{Tuple{Int64,Int64}, Float64} = Dict{Tuple{Int64,Int64}, Float64}()
        
        # get solution in x variables
        for (val, used_arcs) in zip(l_val[k], lambdas_used_arcs[k])
            
            for a in used_arcs
                if val > .00001
                    if haskey(x_sol, a)
                        x_sol[a] += val
                    else
                        x_sol[a] = val
                    end
                end
                
            end
        end
        
        !quiet && verbose && println("\ncommodity $(k):")
        !quiet && verbose && println("lambda_sol:", l_val[k])
        !quiet && verbose && println("x_sol: ", x_sol, "\n")
        push!(solution, x_sol)
    end

    !quiet && println("artificial vars sum: ", art_sum)

    art_vars_in_sol = art_sum > 1e-12
    opt_failed = art_vars_in_sol || !(globally_optimal)


    return (
        rmlp=rmlp,
        solution=solution,
        rmlp_obj=rmlp_obj,

        opt_failed=opt_failed,
        art_vars_in_sol=art_vars_in_sol,
        globally_optimal=globally_optimal,
        hit_time_limit=hit_time_limit,
        
        pi_bar=pi_bar,
        v_bar=v_bar,
        n_bar=n_bar,

        elapsed_time=elapsed_time,
        heur_time=heur_time,
        rmlp_time=rmlp_time,
        subp_time=subp_time,

        lambdas_ref=lambdas_ref,
        capacity_con=capacity_con,
        sources_con=sources,
        sinks_con=sinks,

        cols_gen=cols_gen,
        available_cols=cols_gen-cleaned_lambdas
    )
end 

function get_solution_cost(instance::MnfpData, solution::Vector{Dict{Tuple{Int64,Int64}, Float64}})::Float64
    total_cost::Float64 = 0.0

    K::UnitRange{Int64} = instance.K
    c::Matrix{Float64} = instance.c

    for k in K
        arc_to_cost_k::Dict{Tuple{Int64, Int64}, Float64} = Dict{Tuple{Int64, Int64}, Float64}()
        for (n, a) in enumerate(instance.A)
            arc_to_cost_k[a] = c[k,n]
        end

        x_sol = solution[k]
        for (a, val) in x_sol
            total_cost += arc_to_cost_k[a] * val
        end
    end

    return total_cost
end

# first_time = @elapsed mnfp_cga(toy_instance, quiet=false, verbose=true)
# first_time = @elapsed mnfp_cga(toy_instance, quiet=true, it_print_interval=1)
# passed_time = @elapsed mnfp_cga(toy_instance)
# println("time: $(first_time)")
# println("time: $(passed_time)")
