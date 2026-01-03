using Random
using Graphs
using GraphPlot
using Compose
using Cairo
using Fontconfig
using Colors
using Serialization

# LongCoord = Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}
LongCoord = Tuple{Int64, Int64}


struct MnfpData
    K::UnitRange{Int64}
    k_amount::Int64
    V::UnitRange{Int64}
    vertex_amount::Int64
    s::Int64
    t::Int64
    V_st::Vector{Int64}
    A::Vector{Tuple{Int64, Int64}}
    arc_amount::Int64
    u::Matrix{Float64}
    c::Matrix{Float64}
    arc_to_cost::Vector{Dict{Tuple{Int64, Int64}, Float64}}
    d::Matrix{Float64}
    S::Vector{Array{Int64}}
    T::Vector{Array{Int64}}
    A_s::Vector{Vector{Tuple{Int64, Int64}}}
    A_t::Vector{Vector{Tuple{Int64, Int64}}}
    base_graph::SimpleDiGraph
    sub_graphs::Vector{SimpleDiGraph{Int64}}
    locs_x::Union{Vector{Int64}, Nothing}
    locs_y::Union{Vector{Int64}, Nothing}
    bhratio::Union{Float64, Nothing}
    guaranteed_planar::Union{Bool, Nothing}
    cols::Union{Int64, Nothing}
    rows::Union{Int64, Nothing}

    function MnfpData(vertex_amount::Int64, A::Vector{Tuple{Int64, Int64}}, u::Matrix{Float64}, c::Matrix{Float64}, d::Matrix{Float64}; locs_x::Union{Vector{Int64}, Nothing}=nothing, locs_y::Union{Vector{Int64}, Nothing}=nothing, bhratio::Union{Float64, Nothing}=nothing, guaranteed_planar::Union{Bool, Nothing}=nothing, cols::Union{Int64, Nothing}=nothing, rows::Union{Int64, Nothing}=nothing)

        V::UnitRange{Int64} = 1:vertex_amount
        K::UnitRange{Int64} = 1:length(d[:, 1])

        k_amount::Int64 = collect(K)[end]

        s::Int64 = max(V...)+1
        t::Int64 = max(V...)+2
        V_st::Vector{Int64} = vcat(s, V, t)

        arc_amount::Int64 = length(A)

        arc_to_cost::Vector{Dict{Tuple{Int64, Int64}, Float64}} = Dict{Tuple{Int64, Int64}, Int64}[Dict() for k in K]
        for k in K
            for (n, a) in enumerate(A)
                arc_to_cost[k][a] = c[k,n]
            end
        end

        S::Vector{Array{Int64}} = Vector{Int64}[Int64[] for k in K]
        T::Vector{Array{Int64}} = Vector{Int64}[Int64[] for k in K]

        # # demanda < 0: source
        # S = [["A", "B"], ["A"]]
        # # demanda > 0: sink
        # T = [["D"], ["C", "D"]]
        for k in K
            for v in V
                if (d[k,v] < 0)
                    push!(S[k], v)
                elseif (d[k,v] > 0)
                    push!(T[k], v)
                end
            end
        end

        # G_k = {V_st, vcat(A, A_st[k])}
        A_s::Vector{Vector{Tuple{Int64, Int64}}} = [[] for k in K]
        A_t::Vector{Vector{Tuple{Int64, Int64}}} = [[] for k in K]
        for k in K
            for v in S[k]
                push!(A_s[k], (s,v))
            end

            for v in T[k]
                push!(A_t[k], (v,t))
            end
        end

        # build base graph
        base_graph::SimpleDiGraph = SimpleDiGraph(vertex_amount)
        for (i, j) in A
            add_edge!(base_graph, i, j)
        end

        # make subgraphs with source and sink
        sub_graphs::Vector{SimpleDiGraph{Int64}} = SimpleDiGraph{Int64}[]
        for k in K
            
            # copy base graph and add sources and sinks
            g::SimpleDiGraph = deepcopy(base_graph)
            add_vertices!(g, 2)
            
            for v in S[k]
                add_edge!(g, s, v)
            end
            
            for v in T[k]
                add_edge!(g, v, t)
            end
            push!(sub_graphs, g)
        end

        new(K, k_amount, V, vertex_amount, s, t, V_st, A, arc_amount, u, c, arc_to_cost, d, S, T, A_s, A_t, base_graph, sub_graphs, locs_x, locs_y, bhratio, guaranteed_planar, cols, rows) 
    end
end

"""
    print_mnfp_instance(mnfp::MnfpData)

Print a summary of the MNFP instance.
"""
function print_mnfp_instance(mnfp::MnfpData; max_n_print::Int64=5)

    
    println("Planar MNFP Instance:")
    println("  Vertices: $(mnfp.vertex_amount)")
    println("  Arcs: $(mnfp.arc_amount)")
    println("  Commodities: $(mnfp.k_amount)")

    println("\nSources and Sinks (first $(max_n_print)):")
    for k in 1:min(max_n_print, mnfp.k_amount)
        println("  Commodity $(k):")
        println("    Sources: $(mnfp.S[k])")
        println("    Sinks: $(mnfp.T[k])")
    end

    println("\nArc capacities (first $(max_n_print)):")
    for i in 1:min(max_n_print, mnfp.arc_amount)
        println("  Arc $(mnfp.A[i]): capacity = $(mnfp.u[i])")
    end
    println("\nCommodity costs (first $(max_n_print) arcs):")
    for i in 1:min(max_n_print, mnfp.arc_amount)
        println("  Arc $(mnfp.A[i]): cost = $(mnfp.c[:, i])")
    end
    println("\nDemands (non-zero only):")
    for k in 1:min(max_n_print,mnfp.k_amount)
        for i in 1:min(max_n_print, mnfp.vertex_amount)
            if abs(mnfp.d[k, i]) > 1e-12
                println("  Vertex $i: demand in $(k) = $(mnfp.d[k, i])")
            end
        end
    end
end

function get_edge_labels(u, A; us::Union{Nothing, Matrix}=nothing, c::Union{Nothing, Matrix}=nothing, specific_k::Union{Int64, Nothing}=nothing, show_cost::Bool=false)

    cap::Vector = Int64.(u[1, :])

    got_int::Bool = true

    used_cap::Vector = []
    if us !== nothing
        try
            used_cap = Int64.(us[1, :])
        catch e
            got_int = false
            used_cap = round.(100*deepcopy(us[1, :]))./100
        end
    end

    if c !== nothing
        if specific_k !== nothing
            acost = Int64.(c[specific_k, :])
            scost = ["$(i)" for i in acost]
        else
            scost = ["[$(Int64(minimum(c[:,i]))), $(Int64(maximum(c[:,i])))]" for i in 1:length(A)]
        end
    else
        show_cost=false
    end

    edge_label_dict = Dict()
    for (n, (i,j)) in enumerate(A)

        x = n

        # is the forward arc?
        f = i < j

        # does (j,i) exist?
        inv_exists = (j,i) in A

        if inv_exists
            y = findfirst(x -> x == (j,i), A)

            if !(f)
                x, y = y, x
            end
            label = "$(x) | $(y)"
            label = "$(cap[x]) | $(cap[y])"
            if us !== nothing
                label = "$(used_cap[x])/$(cap[x]) | $(used_cap[y])/$(cap[y])"
            end
            
            if show_cost
                label *= "\n$(scost[x]) | $(scost[y])"
            end

        else
            label = "$(x)"
            label = "$(cap[x])"
            if us !== nothing
                label = "$(used_cap[x])/$(cap[x])"
            end
            
            if show_cost
                label *= "\n$(scost[x])"
            end
        end

        edge_label_dict[(i,j)] = label
    end
    # return String[edge_label_dict[a] for a in A]
    # for a in A
    #     if 4 in a
    #         println("Arc $(a): Label = $(edge_label_dict[a])")
    #     end
    # end

    return edge_label_dict
end

function get_used_capacity(mnfp::MnfpData, solution; specific_k::Union{Int64, Nothing}=nothing)

    us=zeros(1, length(mnfp.A))

    arc_pos = Dict()
    for (n,a) in enumerate(mnfp.A)
        arc_pos[a] = n
    end

    if specific_k != nothing
        for (a, val) in solution[specific_k]
            us[1, arc_pos[a]] += val 
        end
    else
        for k in mnfp.K
            for (a, val) in solution[k]
                us[1, arc_pos[a]] += val 
            end
        end
    end

    return us
end

"""saves a graph plot of the MNFP instance to a file"""
function print_mnfp_graph(
    mnfp::MnfpData; 
    file_path::Union{String, Nothing}="mnfp.png", 
    solution::Union{Vector{Dict{Tuple{Int64, Int64}, Float64}}, Nothing}=nothing,
    multi::Bool=false,
    multi_folder::String="mnfp_multi",
    specific_k::Union{Nothing, Int}=nothing,
    show_cost::Bool=true,
    locs_x::Union{Vector, Nothing}=nothing,
    locs_y::Union{Vector, Nothing}=nothing,
    )

    if isnothing(locs_x) || isnothing(locs_x)

        if isnothing(mnfp.locs_x) || isnothing(mnfp.locs_y)
            locs_x, locs_y = spring_layout(mnfp.base_graph)
        else
            locs_x = Float64.(mnfp.locs_x)
            locs_y = Float64.(mnfp.locs_y)
        end

    end


    if multi
        # make folder
        if !isdir(multi_folder)
            mkdir(multi_folder)
        end
    
        for k in mnfp.K
            print_mnfp_graph(
                mnfp,
                file_path="$(multi_folder)/k$(k).png",
                solution=solution,
                specific_k=k,
                show_cost=show_cost,
                locs_x=locs_x,
                locs_y=locs_y,
            )
        end
    end

    arc_pos = Dict()
    for (n,a) in enumerate(mnfp.A)
        arc_pos[a] = n
    end

    if solution == nothing
        us=nothing
        total_us=nothing
        nodelabel=1:mnfp.vertex_amount
        nodefillc=colorant"lightgrey"

        us = nothing
        total_us = nothing

    else
        us=zeros(1, length(mnfp.A))
        total_us=zeros(1, length(mnfp.A))

        for k in mnfp.K
            for (a, val) in solution[k]
                total_us[1, arc_pos[a]] += val 
            end
        end

        if specific_k != nothing
            for (a, val) in solution[specific_k]
                us[1, arc_pos[a]] += val 
            end
        else
            us = total_us
        end

    end

    if specific_k != nothing
        total_flow = round.(Int64, [sum(mnfp.d[specific_k, i]) for i in 1:mnfp.vertex_amount])
    else    
        total_flow = round.(Int64, [sum(mnfp.d[:, i]) for i in 1:mnfp.vertex_amount])
    end
    
    nodelabel = []
    for i in 1:mnfp.vertex_amount
        if total_flow[i] != 0
            push!(nodelabel, "$(i) ($(total_flow[i]))")
        else
            push!(nodelabel, "$(i)") 
        end
    end
    
    nodefillc = []
    for i in total_flow
        if i > 0
            # push!(nodefillc, colorant"magenta") 
            # push!(nodefillc, colorant"purple") 
            push!(nodefillc, colorant"orange") 
        elseif i < 0
            push!(nodefillc, colorant"cyan") 
        else
            push!(nodefillc, colorant"lightgray")
        end
    end

    if specific_k != nothing
        show_cost = show_cost
    else
        show_cost = false
    end

    edge_label_dict = get_edge_labels(mnfp.u, mnfp.A; us=us, c=mnfp.c, specific_k=specific_k, show_cost=show_cost)
    edgelabel = String[edge_label_dict[a] for a in mnfp.A]

    edgelabelc = [colorant"lightblue" for i in 1:mnfp.arc_amount]
    
    if specific_k !== nothing && solution !== nothing
        # used_perc = total_us[1,:] ./ mnfp.u[1,:]
        # used_perc = us[1,:] ./ mnfp.u[1,:]
        
        for (i, a) in enumerate(mnfp.A)
            
            if haskey(arc_pos, (a[2], a[1]))
                val = max(us[1,i], us[1,arc_pos[(a[2], a[1])]])
            else
                val = us[1,i]
            end

            if val > 0
                edgelabelc[i] = colorant"chartreuse"
            end
        end
    else
        if total_us !== nothing
            used_perc = total_us[1,:] ./ mnfp.u[1,:]
        else
            used_perc = 0 .* mnfp.u[1,:]
        end

        for (i, a) in enumerate(mnfp.A) 
            if haskey(arc_pos, (a[2], a[1]))
                val = max(used_perc[i], used_perc[arc_pos[(a[2], a[1])]])
            else
                val = used_perc[i]
            end

            if val < .8
                edgelabelc[i] = colorant"white"
            elseif val >= 1 - 1e-12
                edgelabelc[i] = colorant"brown1"
            else
                edgelabelc[i] = colorant"firebrick"
            end
        end
    end
    
    if isnothing(mnfp.cols) || isnothing(mnfp.rows)
        cols = ceil(sqrt(mnfp.vertex_amount))
        rows = cols
    else
        cols = mnfp.cols
        rows = mnfp.rows     
    end


    if !(isnothing(file_path))
        # draw(PNG(file_path, 5*cols*cm, 5*rows*cm), gplot(mnfp.base_graph, Float64.(locs_x), Float64.(locs_y), nodelabel=1:n, edgelinewidth=2,  nodelabelsize=40, nodesize=10, arrowlengthfrac=.02, edgelabel=edgelabel, edgelabelc=colorant"lightblue", edgelabelsize=6, edgestrokec=colorant"black", nodefillc=colorant"lightgrey"))
        # draw(PDF(file_path, 5*cols*cm, 5*rows*cm), gplot(mnfp.base_graph, Float64.(locs_x), Float64.(locs_y), nodelabel=1:n, edgelinewidth=2,  nodelabelsize=40, nodesize=10, arrowlengthfrac=.02, edgelabel=edgelabel, edgelabelc=colorant"lightblue", edgelabelsize=6, edgestrokec=colorant"black", nodefillc=colorant"lightgrey"))
        
        if locs_x !== nothing
            draw(PNG(file_path, 7*cols*cm, 7*rows*cm), gplot(mnfp.base_graph, locs_x, locs_y, nodelabel=nodelabel, edgelinewidth=2,  nodelabelsize=40, nodesize=12, arrowlengthfrac=.02, edgelabel=edgelabel, edgelabelc=edgelabelc, edgelabelsize=6, edgestrokec=colorant"gray20", nodefillc=nodefillc))
        else
            draw(PNG(file_path, 7*cols*cm, 7*rows*cm), gplot(mnfp.base_graph, locs_x, locs_y, nodelabel=nodelabel, edgelinewidth=2,  nodelabelsize=40, nodesize=12, edgelabel=edgelabel, edgelabelc=edgelabelc, edgelabelsize=6, edgestrokec=colorant"gray20", nodefillc=nodefillc))
        end
    end
end

"""returns how many times each long edge crosses another long edge"""
function get_cross_count(hlong_edges_grid, vlong_edges_grid)
    hcross_count = Dict{Tuple{Int64, LongCoord}, Int64}()
    vcross_count = Dict{Tuple{Int64, LongCoord}, Int64}()
    
    #initialize the dicts with zero
    for (i, edge_list) in enumerate(hlong_edges_grid)
        for e in edge_list
            hcross_count[(i, e)] = 0
        end
    end
    
    for (i, edge_list) in enumerate(vlong_edges_grid)
        for e in edge_list
            vcross_count[(i, e)] = 0
        end
    end
    
    # doing horizontal
    for (i, edge_list_h) in enumerate(hlong_edges_grid)
        for e_h in edge_list_h
            
            for k in e_h[1]+1:e_h[2]-1
                for e_v in vlong_edges_grid[k]
                    if e_v[1] <= i <= e_v[2]
                        hcross_count[(i, e_h)] += 1
                        vcross_count[(k, e_v)] += 1
                    end
                end
            end
        end
    end

    return hcross_count, vcross_count
end

function get_dist_on_grid(v1, v2, v_grid_loc)
    (i1, j1) = v_grid_loc[v1]
    (i2, j2) = v_grid_loc[v2]

    return sqrt((i1 - i2)^2 + (j1 - j2)^2)
end

function get_closest_arc(component1, component2, v_grid_loc)
    min_dist = Inf
    best_arc::Tuple{Int64, Int64} = (0,0)

    for v1 in component1
        for v2 in component2
            dist = get_dist_on_grid(v1, v2, v_grid_loc)
            if dist < min_dist
                min_dist = dist
                best_arc = (v1, v2)
            end
        end
    end

    return best_arc
end

function random_split(value::Int, n::Int, allow_zero::Bool=false)
    
    if allow_zero
        # Generate n-1 random split points between 0 and value
        split_points = sort(rand(0:value, n-1))
        
        # Add boundaries
        boundaries = [0; split_points; value]
        
        # Calculate differences between consecutive boundaries
        splits = diff(boundaries)
        
        return splits

    else
        if n > value
            error("Cannot split $value into $n positive parts")
        end
        
        # Start with 1 for each part
        result = ones(Int, n)
        remaining = value - n
        
        # Distribute remaining value randomly
        for _ in 1:remaining
            result[rand(1:n)] += 1
        end
        
        return result
    end
end

function get_distance_on_grid(v1, v2, v_grid_loc)
    (i1, j1) = v_grid_loc[v1]
    (i2, j2) = v_grid_loc[v2]

    return sqrt((i1 - i2)^2 + (j1 - j2)^2)
end

"""
    generate_planar_grid_mnfp()

Generate a planar MNFP instance based on a grid graph. Returns a named Tuple with output data.
"""
function generate_planar_grid_mnfp(; 
    vertex_amount::Int64 = 20,
    base_mult::Float64 = 3.0,
    bhratio::Float64 = 1.0,
    K::Int64 = 5,
    capacity_range::Tuple{Float64,Float64}=(10.0, 50.0),
    supply_range::Tuple{Float64,Float64}=(1.0, 20.0),
    cost_range::Tuple{Float64,Float64}=(1.0, 2.0),
    demand_density::Float64=0.1,
    random_capacity_slack::Float64 = 0.2, # how much extra capacity to add randomly
    base_capacity_slack::Float64 = 1.0, # how much capacity to add overall
    seed::Union{Int64,Nothing}=nothing,
    raw_data::Bool=false,
    )

    #################################################################################################### data setup
    guaranteed_planar::Bool = true

    if !isnothing(seed)
        Random.seed!(seed)
    end

    h::Int64 = 1
    b::Int64 = ceil(h/bhratio)  

    while (b*h < base_mult*vertex_amount)
        h += 1
        b = ceil(h/bhratio)
    end

    cols::Int64 = b
    rows::Int64 = h
    println("making graph with $(vertex_amount) vertices in area $(b*h)")

    A::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[]

    # extra A sets
    A_basic::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[]

    vlong_edges_grid::Vector{Vector{LongCoord}} = Vector{LongCoord}[ LongCoord[] for i in 1:b]
    hlong_edges_grid::Vector{Vector{LongCoord}} = Vector{LongCoord}[ LongCoord[] for i in 1:h]

    locs_x = Int64[]
    locs_y = Int64[]

    kept = randperm(b*h)[1:vertex_amount]
    sort!(kept)


    #################################################################################################### make base arc set


    current_v = 1

    # make grid positions
    grid_pos::Matrix{Int64} = zeros(rows,cols)
    v_grid_loc::Dict{Int64, Tuple{Int64,Int64}} = Dict()
    for j in 1:cols
        for i in 1:rows

            base_v = i + (j-1)*rows 
            v = current_v

            if base_v in kept
                grid_pos[i, j] = v
                v_grid_loc[v] = (i,j)
        
                push!(locs_x, j)
                push!(locs_y, i)

                current_v += 1
            end
        end
    end

    # Create grid edges (planar by construction)
    # Vertex (i,j) has index: i + (j-1)*rows
    for j in 1:cols
        for i in 1:rows

            base_v = i + (j-1)*rows 

            a1 = grid_pos[i, j]

            if a1 == 0
                continue
            end

            # Right neighbor
            if j < cols
                for k in j+1:cols
                    a2 = grid_pos[i, k]

                    if a2 != 0
                        push!(A, (a1, a2))
                        push!(A, (a2, a1))  # reverse arc

                        if k == j+1
                            push!(A_basic, (a1, a2))
                            push!(A_basic, (a2, a1))
                        else
                            # push!(hlong_edges_grid[i], (min(a1, a2), max(a1,a2)))
                            # push!(hlong_edges_grid[i], (base_v, i+(k-1)*rows))
                            push!(hlong_edges_grid[i], (j, k))
                        end

                        break
                    end
                end
            end

            # Down neighbor
            if i < rows
                for k in i+1:rows
                    a2 = grid_pos[k, j]

                    if a2 != 0
                        push!(A, (a1, a2))
                        push!(A, (a2, a1))  # reverse arc

                        if k == i+1
                            push!(A_basic, (a1, a2))
                            push!(A_basic, (a2, a1))
                        else
                            # push!(vlong_edges_grid[j], (min(a1, a2), max(a1,a2)))
                            # push!(vlong_edges_grid[j], (base_v, k+(j-1)*rows))
                            # push!(vlong_edges_grid[j], ((i,j), (k,j)))
                            push!(vlong_edges_grid[j], (i, k))
                        end

                        break
                    end
                end
            end

            
            # diagonal (down) neighbor
            if i < rows && j < cols && j > 1
                # toss a coin to choose direction ()
                if rand() > .5
                    x=1
                else
                    x=-1
                end

                # test if new arc might cross another, if it does, dont add it
                down = grid_pos[i+1, j]
                side = grid_pos[i, j+x]
                if down != 0 && side != 0
                    if (down,side) in A || (side, down) in A
                        continue
                    end
                end

                a3 = grid_pos[i+1, j+x]
                if a3 != 0
                    
                    # toss a coin to make a directed arc
                    if rand() > .5
                        a = (a1, a3)
                    else
                        a = (a3, a1)
                    end
                    push!(A, a)
                    push!(A_basic, a)
                end
                
            end
        end
    end

    ############################################################################################### clean unwanted long A

    arcs_to_remove = Tuple{Int64, Int64}[]

    hcross_count, vcross_count = get_cross_count(hlong_edges_grid, vlong_edges_grid)

    max_cross = max(maximum(values(hcross_count)), maximum(values(vcross_count)))

    done = false
    for iter in 1:1e9 # while max_cross > 0

        h_to_remove = Vector{LongCoord}[LongCoord[] for _ in 1:length(hlong_edges_grid)]
        for ((k, l), vertex_amount) in hcross_count
            if vertex_amount == max_cross
                push!(h_to_remove[k], l)
            end
        end

        v_to_remove = Vector{LongCoord}[LongCoord[] for _ in 1:length(vlong_edges_grid)]
        for ((k, l), vertex_amount) in vcross_count
            if vertex_amount == max_cross
                push!(v_to_remove[k], l)
            end
        end

        # add correspondent A to arcs_to_remove
        for (k, edge_list) in enumerate(h_to_remove)
            for e in edge_list
                i = k

                a1 = (grid_pos[i,e[1]], grid_pos[i,e[2]])
                a2 = (grid_pos[i,e[2]], grid_pos[i,e[1]])
                push!(arcs_to_remove, a1)
                push!(arcs_to_remove, a2)

            end
            filter!(x -> !(x in edge_list), hlong_edges_grid[k])
        end
        for (k, edge_list) in enumerate(v_to_remove)
            for e in edge_list
                j = k

                a1 = (grid_pos[e[1],j], grid_pos[e[2],j])
                a2 = (grid_pos[e[2],j], grid_pos[e[1],j])
                push!(arcs_to_remove, a1)
                push!(arcs_to_remove, a2)

            end
            filter!(x -> !(x in edge_list), vlong_edges_grid[k])
        end

        hcross_count, vcross_count = get_cross_count(hlong_edges_grid, vlong_edges_grid)

        if isempty(hcross_count)
            max_cross_h = 0
        else
            max_cross_h = maximum(values(hcross_count))
        end
        if isempty(vcross_count)
            max_cross_v = 0
        else
            max_cross_v = maximum(values(vcross_count))
        end

        max_cross = max(max_cross_h, max_cross_v)

        if max_cross == 0
            done = true
            break
        end
    end
    
    if !(done)
        println("ARC REMOVAL FAILED")
        guaranteed_planar=false
    end

    A = setdiff(A, arcs_to_remove)
    arcs_pruned = deepcopy(A)

    ############################################################################################### add arcs to orphaned vertex

    v_in = Vector{Tuple{Int64, Int64}}[Tuple{Int64, Int64}[] for i in 1:vertex_amount]
    v_out = Vector{Tuple{Int64, Int64}}[Tuple{Int64, Int64}[] for i in 1:vertex_amount]

    for a in A
        v1, v2 = a
        push!(v_in[v2], a)
        push!(v_out[v1], a)
    end

    multi = Tuple{Int64,Int64}[(x,y) for x in [1,-1,0] for y in [1,0,-1]]
    for v in 1:vertex_amount
        # needs_in = length(v_in[v]) == 0
        # needs_out = length(v_out[v]) == 0
        
        needs_in = length(v_in[v]) <= 1
        needs_out = length(v_out[v]) <= 1
        
        if needs_in || needs_out
        
            println("$(v) is alone")
            
            i, j = v_grid_loc[v]

            done = false
            # search for a friend
            for i_mod in 1:h
                for j_mod in 1:b
                    for (i_mult, j_mult) in multi
                        if i_mult == 0 && j_mult == 0
                            continue
                        end

                        i_k = i_mult*i_mod + i
                        j_k = j_mult*j_mod + j

                        if i_k > h || i_k < 1 || j_k > b || j_k < 1 || i_mod + j_mod == 0
                            continue
                        end
                        
                        v_k = grid_pos[i_k, j_k]
                        if v_k != 0
                            
                            # new friend
                            
                            if needs_out && !((v, v_k) in v_out[v])
                                println("adding arc $((v,v_k))")
                                push!(A, (v, v_k))
                                push!(v_in[v_k], (v, v_k))
                                push!(v_out[v], (v, v_k))
                            end
                            if needs_in && !((v_k, v) in v_in[v])
                                println("adding arc $((v_k, v))")
                                push!(A, (v_k, v))
                                push!(v_in[v], (v_k, v))
                                push!(v_out[v_k], (v_k, v))
                            end

                            done = true

                            if i_mod*j_mod > 1
                                guaranteed_planar=false
                            end
                            break
                        end
                    end

                    if done
                        break
                    end
                end

                if done
                    break
                end
            end
            
        end
    end

    ############################################################################################### ensure connectivity 1
    
    for v in 1:vertex_amount
        if isempty(v_in[v]) || isempty(v_out[v])
            a = get_closest_arc([v], [x for x in 1:vertex_amount if x != v], v_grid_loc)
        end

        if isempty(v_in[v])
            push!(A, (a[2], a[1]))
            println("adding edge $((a[2], a[1])) to ensure connectivity for $(v)")
        end
        if isempty(v_out[v])
            push!(A, a)
            println("adding edge $(a) to ensure connectivity for $(v)") 
        end
        
    end
    
    ############################################################################################### ensure connectivity 2

    done=false
    # build base graph
    base_graph::SimpleDiGraph = SimpleDiGraph(vertex_amount)
    for (i, j) in A
        add_edge!(base_graph, i, j)
    end

    for iter in 1:1e5
        components = connected_components(base_graph)
        
        done = length(components) == 1
        if done
            break
        end
        
        guaranteed_planar=false

        println("graph has $(length(components)) components, connecting...")
        
        component1, component2 = components[1:2]
        a = get_closest_arc(component1, component2, v_grid_loc)
        println("adding edge $(a)")
        push!(A, a)
        push!(A, (a[2], a[1]))

        add_edge!(base_graph, a[1], a[2])
        add_edge!(base_graph, a[2], a[1])

    end

    if !(done)
        println("FAILED TO ENSURE CONNECTIVITY")
    end

    ############################################################################################## extra arcs


    # extra_arc_amount = 0

    # arcs_not_in_graph = []
    # for i in 1:vertex_amount
    #     for j in 1:vertex_amount
    #         if i == j
    #             continue
    #         end
    #         if !((i, j) in A) 
    #             push!(arcs_not_in_graph, (i,j))
    #         end
    #     end
    # end

    # new_arcs = randperm(arcs_not_in_graph)[1:extra_arc_amount]

    # for a in new_arcs
    #     push!(A, a)
    # end

    
    ############################################################################################### generate capacities, costs and demands
    
    sort!(A)
    m::Int64 = length(A)

    # Generate demands per commodity
    d::Matrix{Float64} = zeros(K, vertex_amount)

    total_vertex_demand::Vector{Float64} = zeros(vertex_amount)
    total_vertex_supply::Vector{Float64} = zeros(vertex_amount)

    single_s_or_t::Bool = true

    for k in 1:K

        # either source or sinks
        num_of_dealers::Int64 = max(2, round(Int64, demand_density * vertex_amount))

        num_sources::Int64 = 1
        num_sinks::Int64 = 1

        # randomly choose sources and sinks
        if single_s_or_t
            
            if rand() < 0.5
                num_sources += num_of_dealers-2
            else
                num_sinks += num_of_dealers-2
            end
        
        else
            for i in 1:num_of_dealers-2
                if rand() < 0.5
                    num_sources += 1
                else
                    num_sinks += 1
                end
            end
        end

        # num_sinks, num_sources = num_sources, num_sinks

        sources = randperm(vertex_amount)[1:num_sources]
        remaining = setdiff(1:vertex_amount, sources)
        # sinks = remaining[randperm(length(remaining))[1:min(num_sinks, length(remaining))]]
        sinks = remaining[randperm(length(remaining))[1:num_sinks]]
        
        # Total supply
        # total_supply = rand(supply_range[1]:0.5:supply_range[2])
        total_supply::Float64 = rand(supply_range[1]:supply_range[2])*max(num_sinks, num_sources)

        # split total supply among sources and sinks
        supply_split = random_split(round(Int64, total_supply), num_sources)
        demand_split = random_split(round(Int64, total_supply), num_sinks)
        
        for (i, s) in enumerate(sources)
            d[k, s] = -supply_split[i]
            total_vertex_supply[s] += supply_split[i]
        end

        for (i, t) in enumerate(sinks)
            d[k, t] = demand_split[i]
            total_vertex_demand[t] += demand_split[i]
        end
    end

    # get minimum capacity needed per arc
    min_capacity_per_arc::Dict{Tuple{Int64, Int64}, Float64} = Dict{Tuple{Int64, Int64}, Float64}()
    for a in A
        min_capacity_per_arc[a] = 0.0
    end

    for v in 1:vertex_amount
        if isempty(v_in[v])
            println("vertex $(v) has no incoming arcs")
            continue
        end
        if isempty(v_out[v])
            println("vertex $(v) has no outgoing arcs")
        end 
    end

    # for each vertex, split its demand among incoming arcs, and its supply among outgoing arcs
    for v in 1:vertex_amount
        if total_vertex_demand[v] > 0 
            incoming_arcs = v_in[v]
            num_incoming = length(incoming_arcs)

            if num_incoming == 0
                continue
            end

            # demand_split = random_split(round(Int64, total_vertex_demand[v]), num_incoming, true)
            demand_split = ceil(total_vertex_demand[v]/num_incoming)

            for (i, a) in enumerate(incoming_arcs)
                min_capacity_per_arc[a] += demand_split
            end
        end
        if total_vertex_supply[v] > 0
            outgoing_arcs = v_out[v]
            num_outgoing = length(outgoing_arcs)

            if num_outgoing == 0
                continue
            end

            # supply_split = random_split(round(Int64, total_vertex_supply[v]), num_outgoing, true)
            supply_split = ceil(total_vertex_supply[v]/num_outgoing)

            for (i, a) in enumerate(outgoing_arcs)
                min_capacity_per_arc[a] += supply_split
            end
        end
    end

    demand_sum::Float64 = sum(total_vertex_demand)

    # Generate arc costs per commodity
    c::Matrix{Float64} = rand(cost_range[1]:0.1:cost_range[2], K, m)

    base_cost::Vector{Float64} = Float64[]
    for a in A
        push!(base_cost, 10*get_distance_on_grid(a[1], a[2], v_grid_loc))
    end

    for k in 1:K
        c[k, :] = ceil.(base_cost.*c[k,:])
    end
    
    # Generate arc capacities
    u::Matrix{Float64} = zeros(1, m)
   
    # how much capacity to divide equally among arcs
    base_capacity::Float64 = ceil(base_capacity_slack*demand_sum / m)
    
    for (i, a) in enumerate(A)
        u[1, i] = min_capacity_per_arc[a] + base_capacity
    end

    for i in rand(1:m, round(Int64, random_capacity_slack*demand_sum))
        u[1, i] += 1.0
    end

    if raw_data
        return (
            vertex_amount=vertex_amount,
            bhratio=bhratio,
            K=K,
            capacity_range=capacity_range,
            supply_range=supply_range,
            cost_range=cost_range,
            demand_density=demand_density,
            base_mult=base_mult,
            h=h,
            b=b,
            cols=cols,
            rows=rows,
            A=A,
            A_basic=A_basic,
            vlong_edges_grid=vlong_edges_grid,
            hlong_edges_grid=hlong_edges_grid,
            locs_x=locs_x,
            locs_y=locs_y,
            v_in=v_in,
            v_out=v_out,
            kept=kept,
            grid_pos=grid_pos,
            v_grid_loc=v_grid_loc,
            m=m,
            u=u,
            c=c,
            random_capacity_slack=random_capacity_slack,
            base_capacity_slack=base_capacity_slack,
            d=d,
            arcs_pruned=arcs_pruned
        )
    else
        MnfpData(
            vertex_amount,
            A,
            u,
            c,
            d;
            locs_x=locs_x,
            locs_y=locs_y,
            bhratio=bhratio,
            guaranteed_planar=guaranteed_planar,
            cols=cols,
            rows=rows,
        )
    end
end

function save_instance(mnfp_instance::MnfpData, mod::String, instance_base_path::String, instance_base_img_path::String, instance_counter::Int64; multi_k::Bool=true)
    suffix::String = "$(instance_counter)$(mod)"
    
    save_mnfp(mnfp_instance, "$(instance_base_path)_$(suffix).jls")

    if !isdir("$(instance_base_img_path)_$(suffix)")
        mkdir("$(instance_base_img_path)_$(suffix)")
    end

    print_mnfp_graph(mnfp_instance, file_path="$(instance_base_img_path)_$(suffix)/all_sum.png", multi_folder="$(instance_base_img_path)_$(suffix)", multi=multi_k)

end

function harden_instance(mnfp_instance::MnfpData, instance_base_path, instance_base_img_path; max_iter::Int64=200, max_time::Int64=5, multi_k::Bool=True, instance_counter::Int64=1)


    ############################################# ITERATING

    best_u = deepcopy(mnfp_instance.u)
    best_score = -Inf

    frozen = Int64[]
    hard_found = false

    println("initial u: ", mnfp_instance.u)

    iter = 0
    while !(hard_found) && iter < max_iter
        iter += 1

        ##################### solve the instance
        # save_mnfp(mnfp_instance, "current_instance.jls")

        passed_time = @elapsed rmlp, solution, rmlp_obj, opt_failed, art_vars_in_sol, globally_optimal, hit_time_limit, pi_bar = mnfp_cga(mnfp_instance, quiet=false, verbose=false, max_time=max_time)

        # println("solution: ", solution)

        us = get_used_capacity(mnfp_instance, solution)
        diff = mnfp_instance.u - us

        println("\n",iter, " last diffsum: ", sum(diff), " time: ", passed_time, " seconds")

        ###################### check usefulness
        if opt_failed

            # instance simply needed more time 
            if hit_time_limit && !(art_vars_in_sol) && termination_status(rmlp) == MOI.OPTIMAL

                # save instance
                save_instance(mnfp_instance, "h", instance_base_path, instance_base_img_path, instance_counter, multi_k=multi_k)
                

                # easy it to generate medium and easy
                for (i, a) in enumerate(mnfp_instance.A)
                    mnfp_instance.u[1, i] = ceil(1.1*mnfp_instance.u[1, i])
                end
                save_instance(mnfp_instance, "m", instance_base_path, instance_base_img_path, instance_counter, multi_k=multi_k)

                
                for (i, a) in enumerate(mnfp_instance.A)
                    mnfp_instance.u[1, i] = ceil(1.1*mnfp_instance.u[1, i])
                end
                save_instance(mnfp_instance, "e", instance_base_path, instance_base_img_path, instance_counter, multi_k=multi_k)

                hard_found = true            
                
                continue
                
            else # instance is infeasible, or we are not sure if it is feasible (treat as infeasible)

                println("treating infeasibility")
                
                for (i, a) in enumerate(mnfp_instance.A)
                    mnfp_instance.u[1, i] += 1.0

                    if pi_bar[i] < -1e10 # too tight, causing infeasibility
                        println("increased u_$(a) to avoid infeasibility")
                        
                        # if !(i in frozen)
                        #     push!(frozen, i)
                        # end

                        # mnfp_instance.u[1, i] = mnfp_instance.u[1, i]+rand(1:10)
                        mnfp_instance.u[1, i] = 2*(mnfp_instance.u[1, i]+rand(1:10))
                    end
                end

            end
                                                    
        else # we found a feasible instance that can be solved within max_time 

            if length(frozen) == mnfp_instance.arc_amount
                println("froze all arcs")
                break
            end
            
            # neither -inf nor 0
            useful = filter(x -> x < -1e-10 && x > -1e10, pi_bar) # filtering x < -1e10 is redundant, but no harm
            score = 0.0
            if useful != []
                # score = -sum(useful)*std(useful)
                score = length(useful)
            end
            
            if score > best_score 
                println("new best score: $(score) / $(mnfp_instance.arc_amount), new u: $(mnfp_instance.u)")
                best_u = deepcopy(mnfp_instance.u)
                best_score = score

                save_instance(mnfp_instance, "b", instance_base_path, instance_base_img_path, instance_counter, multi_k=multi_k)

                # save_mnfp(mnfp_instance, "best_instance.jls")
                # save_solution(solution, "best_solution.jls")

                # if passed_time > max_time
                #     println("instance long enough")
                #     break
                # end

            else # return to the best found instance
                println("rolling back to previous best instance")
                for (i, a) in enumerate(mnfp_instance.A)
                    mnfp_instance.u[1, i] = best_u[1, i] 
                end
            end

            println("adjusting looseness")
            # for (i, a) in enumerate(mnfp_instance.A)
            #     if us[1, i] > 0 

            #         if rand() < 0.5
            #             println("tightened ", a)
        
            #             if rand() < 0.2
            #                 mnfp_instance.u[1, i] = max(ceil(.5*us[1, i]),1)
            #             else
            #                 mnfp_instance.u[1, i] = max(us[1, i]-1,1)
            #             end
            #         end
        
            #     end
            # end

            # sort arc flows from largest to minimum
            neg_flow = -us[1, :]
            largest_flow = sort(1:mnfp_instance.arc_amount, by=x -> neg_flow[x])
            total_flow = sum(max.(mnfp_instance.d, 0))
            
            # dislocate cut_perc of arcs flow from largest flow to smallest until flow_target is reached
            cut_perc = 0.5
            flow_target = 0.05 * total_flow
            
            cur_flow = 0.0
            for i in largest_flow
                if i in frozen
                    continue
                else
                    push!(frozen, i)
                end

                flow = us[1, i]

                new_u_i = max( ceil( flow*(1-cut_perc) ), 1 )
                mnfp_instance.u[1, i] = new_u_i


                cur_flow += flow - new_u_i
                if cur_flow > flow_target
                    break
                end
            end
        end
        
        if max_iter == iter
            println("maximum iter reached")
        end
    end
    println("done searching")

    # passed_time = @elapsed rmlp, solution, rmlp_obj, opt_failed, art_vars_in_sol, globally_optimal, hit_time_limit, pi_bar = mnfp_cga(mnfp_instance, quiet=false, verbose=false, full_output=true)

    # println("\n\ntime: $(passed_time)")
    # # println("solution: $(solution)")
    # println("rmlp_obj: $(rmlp_obj)")
    # println("opt_failed: $(opt_failed)")
    # println("pi_bar: $(pi_bar)")

    #############################################

    for (i, a) in enumerate(mnfp_instance.A)
        mnfp_instance.u[1, i] = best_u[1, i] 
    end
    # return rmlp
end



# Serialize an instance struct to a file
function save_mnfp(obj, filename::String)
    open(filename, "w") do io
        serialize(io, obj)
    end
    println("Saved to $filename")
end

# Serialize a solution struct to a file
function save_solution(obj, filename::String)
    open(filename, "w") do io
        serialize(io, obj)
    end
    println("Saved to $filename")
end

# read instance from txt
function read_mnfp(filename::String)
    # Initialize data structures using basic types
    num_nodes = 0
    num_arcs = 0
    density = 0.0
    num_commodities = 0
    
    # Using Dicts and Arrays instead of custom structures
    sources = Dict{Int, Dict{Int, Float64}}()  # commodity => node => supply
    sinks = Dict{Int, Dict{Int, Float64}}()    # commodity => node => demand
    costs = Dict{Int, Dict{Tuple{Int,Int}, Float64}}()  # commodity => (i,j) => cost
    capacities = Dict{Tuple{Int,Int}, Float64}()  # (i,j) => capacity
    
    open(filename, "r") do f
        section = ""
        
        for line in eachline(f)
            line = strip(line)
            isempty(line) && continue
            
            # Check for section headers
            if occursin("SOURCES_SECTION", line)
                section = "SOURCES"
                continue
            elseif occursin("SINKS_SECTION", line)
                section = "SINKS"
                continue
            elseif occursin("COSTS_SECTION", line)
                section = "COSTS"
                continue
            elseif occursin("CAPACITY_SECTION", line)
                section = "CAPACITY"
                continue
            end
            
            # Parse header information
            if startswith(line, "num_nodes:")
                num_nodes = parse(Int, split(line, ":")[2])
            elseif startswith(line, "num_arcs:")
                num_arcs = parse(Int, split(line, ":")[2])
            elseif startswith(line, "density:")
                density = parse(Float64, split(line, ":")[2])
            elseif startswith(line, "num_commodities:")
                num_commodities = parse(Int, split(line, ":")[2])
            # Parse section data
            elseif section == "SOURCES"
                parts = split(line)
                commodity = parse(Int, parts[1])
                node = parse(Int, parts[2])
                supply = parse(Float64, parts[3])
                if !haskey(sources, commodity)
                    sources[commodity] = Dict{Int, Float64}()
                end
                sources[commodity][node] = supply
            elseif section == "SINKS"
                parts = split(line)
                commodity = parse(Int, parts[1])
                node = parse(Int, parts[2])
                demand = parse(Float64, parts[3])
                if !haskey(sinks, commodity)
                    sinks[commodity] = Dict{Int, Float64}()
                end
                sinks[commodity][node] = demand
            elseif section == "COSTS"
                parts = split(line)
                commodity = parse(Int, parts[1])
                from_node = parse(Int, parts[2])
                to_node = parse(Int, parts[3])
                cost = parse(Float64, parts[4])
                if !haskey(costs, commodity)
                    costs[commodity] = Dict{Tuple{Int,Int}, Float64}()
                end
                costs[commodity][(from_node, to_node)] = cost
            elseif section == "CAPACITY"
                parts = split(line)
                from_node = parse(Int, parts[1])
                to_node = parse(Int, parts[2])
                capacity = parse(Float64, parts[3])
                capacities[(from_node, to_node)] = capacity
            end
        end
    end

    # Return MnfpData

    vertex_amount::Int64 = num_nodes

    A::Vector{Tuple{Int64, Int64}} = [a for (a, u) in capacities]
    sort!(A)
    
    u::Matrix{Float64} = zeros(1, length(A))
    for (i, a) in enumerate(A)
        u[1, i] = capacities[a]
    end

    c::Matrix{Float64} = zeros(num_commodities, num_arcs)
    for k in 1:num_commodities
        for (i, a) in enumerate(A)
            c[k, i] = costs[k][a]
        end
    end

    d::Matrix{Float64} = zeros(num_commodities, num_nodes)
    for dealer in [sinks, sources]
        for (k, k_dealer) in dealer
            for (v, d_v) in k_dealer
                d[k, v] = d_v
            end
        end
    end

    return MnfpData(vertex_amount, A, u, c, d)
    
    # return (
    #     num_nodes = num_nodes,
    #     num_arcs = num_arcs,
    #     density = density,
    #     num_commodities = num_commodities,
    #     sources = sources,
    #     sinks = sinks,
    #     costs = costs,
    #     capacities = capacities
    # )
end

# Load a struct from a file
function load_mnfp(filename::String)

    if filename[end-3:end] == ".txt"
        return read_mnfp(filename)
    else

        obj = open(filename, "r") do io
            deserialize(io)
        end
        println("Loaded from $filename")
        return obj

    end
end

function write_mnfp_txt(filename::String, mnfp_instance::MnfpData)
    
    vertex_amount::Int = mnfp_instance.vertex_amount
    A::Vector{Tuple{Int,Int}} = mnfp_instance.A
    u::Matrix{Float64} = mnfp_instance.u
    c::Matrix{Float64} = mnfp_instance.c
    d::Matrix{Float64} = mnfp_instance.d
    
    num_commodities = size(c, 1)
    num_arcs = length(A)
    density = 0.0
    
    open(filename, "w") do f
        # Write header
        println(f, "num_nodes: $vertex_amount")
        println(f, "num_arcs: $num_arcs")
        println(f, "density: $density")
        println(f, "num_commodities: $num_commodities")
        
        # Write SOURCES_SECTION (negative values in d matrix)
        println(f, "SOURCES_SECTION")
        for k in 1:num_commodities
            for node in 1:vertex_amount
                if d[k, node] < 0
                    println(f, "$k $node $(d[k, node])")
                end
            end
        end
        
        # Write SINKS_SECTION (positive values in d matrix)
        println(f, "SINKS_SECTION")
        for k in 1:num_commodities
            for node in 1:vertex_amount
                if d[k, node] > 0
                    println(f, "$k $node $(d[k, node])")
                end
            end
        end
        
        # Write COSTS_SECTION
        println(f, "COSTS_SECTION")
        for k in 1:num_commodities
            for (arc_idx, (i, j)) in enumerate(A)
                println(f, "$k $i $j $(c[k, arc_idx])")
            end
        end
        
        # Write CAPACITY_SECTION
        println(f, "CAPACITY_SECTION")
        for (arc_idx, (i, j)) in enumerate(A)
            println(f, "$i $j $(u[arc_idx])")
        end
    end
    
    println("MNFP instance written to $filename")
end