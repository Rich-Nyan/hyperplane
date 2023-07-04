using Plots

# Read from Input
# file = open("input/KKT_trajectory_state.csv")
file = open("input/KKT_trajectory_state_3p.csv")
mat = readlines(file)
num = rsplit(mat[1], ",")
mati = zeros(Float64,floor(Int, length(mat)/4), 4 , length(num))
for i in 1:length(mat)
    index = ceil(Int,i/4)
    input = map(x -> parse(Float64, x), rsplit(mat[i],","))
    for j in 1:length(num)
        mati[index, i - 4*(index-1), j] = input[j]
    end
end

# Constants
players = floor(Int,length(mat)/4)
colors = palette(:default)[1:players] # Players can be max at 16
# Boundaries
x_domain = extrema(mati[:,1,:]) .+ (-0.01, 0.01)
y_domain = extrema(mati[:,2,:]) .+ (-0.01, 0.01)
domain  = [minimum([x_domain[1],y_domain[1]]),maximum([x_domain[2],y_domain[2]])]

# KoZ matrices
# 2 Players
# ωs = [0.0 0.03;
#       0.0 0.0]

# ρs = [0.0 0.25;
#       0.0 0.0]

# αs = [0.0 3/4*pi;
#       0.0 0.0]

# adj_mat = [false true;
#            false false]

# 3 Players
ωs = [0.0 0.03 0.03;
      0.0 0.0 -0.03;
      0.0 0.0  0.0]

ρs = [0.0 0.25 0.25;
      0.0 0.0  0.1
      0.0 0.0  0.0]

αs = [0.0 3/4*pi pi;
      0.0 0.0 5/4*pi;
      0.0 0.0 0.0]

adj_mat = [false true true;
           false false true;
           false false false]

numofconstraints = sum(adj_mat)
indexvalues = findall(x -> x == true, adj_mat)

# Plot
gify = @animate for i=1:1:length(num)
    # Setup
    p = plot(
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
# gif(gify, "hyperplane_2p.gif", fps = 15)
gif(gify, "hyperplane_3p.gif", fps = 15)
