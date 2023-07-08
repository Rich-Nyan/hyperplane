using CSV, DataFrames
using Plots 

function reader(str::DataFrame)
    matrix = Float64[]
    for row in eachcol(str)
        for value in row
            valuer = eval(Meta.parse(value))
            push!(matrix, valuer)
        end
    end
    matrix = reshape(matrix, nrow(str), nrow(str))
    return matrix
end

function hyperplane(states::AbstractString,params::Tuple{AbstractString, AbstractString, AbstractString, AbstractString})
    # Input
    mat = readlines(states)
    num = rsplit(mat[1], ",")
    mati = zeros(Float64,floor(Int, length(mat)/4), 4 , length(num))
    for i in 1:length(mat)
        index = ceil(Int,i/4)
        input = map(x -> parse(Float64, x), rsplit(mat[i],","))
        for j in 1:length(num)
            mati[index, i - 4*(index-1), j] = input[j]
        end
    end

    adj_mat =  Matrix(CSV.read(params[1], DataFrame, header = false, delim = ' ', types = Bool))
    αs = reader(CSV.read(params[2], DataFrame, header = false, delim = ' ', types = String))
    ωs = reader(CSV.read(params[3], DataFrame, header = false, delim = ' ', types = String))
    ρs = reader(CSV.read(params[4], DataFrame, header = false, delim = ' ', types = String))
   
    # Constants
    players = floor(Int,length(mat)/4)
    colors = palette(:default)[1:players] # Players can be max at 16
    # Boundaries
    x_domain = extrema(mati[:,1,:]) .+ (-0.01, 0.01)
    y_domain = extrema(mati[:,2,:]) .+ (-0.01, 0.01)
    domain  = [minimum([x_domain[1],y_domain[1]]),maximum([x_domain[2],y_domain[2]])]


    numofconstraints = sum(adj_mat, init = 0)
    indexvalues = findall(x -> x == true, adj_mat)

    # Plot
    gify = @animate for i=1:1:length(num)
        # Setup
        plot(
        linewidth = 4,
        label = false,
        xlabel = 'x',
        ylabel = 'y',
        title = string(players)*" Player Hyperplane",
        aspectratio = 1,
        ylimits = domain,
        xlimits = domain,
        legend = false
        )

        # Trajectory
        for j in 1:players
            plot!(
                [mati[j,1,k] for k in 1:1:i],
                [mati[j,2,k] for k in 1:1:i],
                j,
                linewidth = 1,
                linecolor = colors[j],
                label = false
            )
        end

        # Points
        for j in 1:players
            scatter!(
                [mati[j,1,i]],
                [mati[j,2,i]],
                markersize = 5,
                color = colors[j],
            )
            annotate!(mati[j,1,i], mati[j,2,i], text(string(j), 8, :black))
        end

        # Hyperplane
        for j in 1:numofconstraints
            constrain = indexvalues[j][1]
            center = indexvalues[j][2]
            line = (x, y) -> (tan(pi/2 .+ αs[constrain,center] .+ ωs[constrain,center] .* x) .* (y .- mati[center,1,x] .- ρs[constrain,center] .* cos(αs[constrain,center] .+ ωs[constrain,center] .* x)) .+ mati[center,2,x] .+ ρs[constrain,center] .* sin(αs[constrain,center] .+ ωs[constrain,center] .* x))
            plot!(
                [k for k in domain],
                [line(i,k) for k in domain],
                players + j,
                linewidth = 1,
                seriesalpha = 0.1,
                linecolor = colors[constrain],
                label = false
            )
        end

        # Circle
        for j in 1:numofconstraints
            constrain = indexvalues[j][1]
            center = indexvalues[j][2]
            c = x -> (mati[center,1,i] .+ ρs[constrain, center] .* cos.(x), mati[center,2,i] .+ ρs[constrain, center] .* sin.(x)) 
            
            plot!(
                [c(t) for t = range(0, stop=2pi, length=100)],
                linewidth = 1,
                linecolor = colors[constrain],
                seriesalpha = 1,
                # fill = true,
                # fillalpha = 0.1,
                label = false
            )
        end

        # Arrow
        for j in 1:numofconstraints
            constrain = indexvalues[j][1]
            center = indexvalues[j][2]
            u = [ρs[constrain, center] * cos(αs[constrain, center] + ωs[constrain, center] * i)]
            v = [ρs[constrain, center] * sin(αs[constrain, center] + ωs[constrain, center] * i)]
            quiver!([mati[center,1,i]], [mati[center,2,i]], quiver =(u, v), arrow = true, arrowhead_length = 0.05, arrowhead_width = 0.03, linecolor = colors[constrain], linewidth = 2)
        end
    end
    gif(gify, "hyperplane_" * string(players) * "p.gif", fps = 15)
end

hyperplane("input/KKT_trajectory_state.csv",("params/2players/adj_mat_2p.txt","params/2players/alpha_2p.txt","params/2players/omega_2p.txt","params/2players/rho_2p.txt"))
hyperplane("input/KKT_trajectory_state_3p.csv",("params/3players/adj_mat_3p.txt","params/3players/alpha_3p.txt","params/3players/omega_3p.txt","params/3players/rho_3p.txt"))
