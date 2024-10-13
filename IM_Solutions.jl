# Packages for plotting, formatting, storing data, basic statistics, and random number generation
using Plots, PGFPlotsX, LaTeXStrings, DataFrames, Statistics, StableRNGs

include("IM_Graphs.jl")

# Plots backend
pgfplotsx()

# Random seed
rng = StableRNG(12345)

# Default plotting parameters
ticksize = 10
legendsize = 10
labelsize = 14
legendpos = :topright
	
# Edge object
struct Edge
    dst::Tuple{Int, Int, Int}
    weight::Float64
end

# Object for exact enumeration including distribution, observables
struct Ising_exact
	vertices::Array{Tuple{Int, Int, Int}}
	β::Float64
	J::Dict{Tuple{Int, Int, Int}, Set{Edge}}
	h::Dict{Tuple{Int, Int, Int}, Float64}
	P::Vector{Float64}
	σ::Vector{Dict{Tuple{Int, Int, Int}, Int64}}
	e::Vector{Dict{Tuple{Int, Int, Int}, Float64}}
	F::Float64
end 

# Function for exact enumeration
function Ising_exact(vertices::Array{Tuple{Int, Int, Int}}, β::Float64, J::Dict{Tuple{Int, Int, Int}, Set{Edge}}, h::Dict{Tuple{Int, Int, Int}, Float64})
	N = length(vertices)
	σ_N = 2^N
	σ = Vector{Dict{Tuple{Int, Int, Int}, Int64}}()
	for i in 0:σ_N-1
		binary_string = string(i, base=2)
		padded_binary_string = lpad(binary_string, N, '0')
		σ_i = parse.(Int, collect(padded_binary_string))
		σ_i[σ_i .== 0] .= -1
		σ_i_dict = Dict(vertex => σ_i[j] for (j, vertex) in enumerate(vertices))
		push!(σ, σ_i_dict)
	end
	e = [Dict(vertex => -0.5 * σ_i[vertex] * sum([σ_i[edge.dst] * edge.weight for edge in J[vertex]]) + h[vertex] * σ_i[vertex] for vertex in vertices) for σ_i in σ]			 
	E = [sum(values(e_i)) for e_i in e]
	W = exp.(-β*E)
	Z = sum(W)
	P = W/Z
	F = -(1/β)*log(Z)
	Ising_exact(vertices, β, J, h, P, σ, e, F) 
end

# Object for MCMC sample observations
mutable struct Ising_MC
    vertices::Array{Tuple{Int, Int, Int}}
    β::Float64
    J::Dict{Tuple{Int, Int, Int}, Set{Edge}}
    h::Dict{Tuple{Int, Int, Int}, Float64}
    σ::Dict{Tuple{Int, Int, Int}, Int}
	E::Float64
	E_min::Float64
end 

# Object for collecting MCMC samples
struct Ising_MC_observables
    β::Float64 
	s::Vector{Float64}
    E::Vector{Float64}
end

# Function for generating initial MCMC configuration
function Ising_MC(vertices::Array{Tuple{Int, Int, Int}}, β::Float64, J::Dict{Tuple{Int, Int, Int}, Set{Edge}}, h::Dict{Tuple{Int, Int, Int}, Float64})
    σ = Dict(i => rand(rng, [-1 1]) for i in vertices)
	e = sum([-0.5 * σ[i] * sum([σ[edge.dst] * edge.weight for edge in J[i]]) + h[i] * σ[i] for i in vertices])
	E = sum(e)
	σ_min = Dict(i => 1 for i in vertices)
	e_min = sum([-0.5 * σ_min[i] * sum([σ_min[edge.dst] * edge.weight for edge in J[i]]) + h[i] * σ_min[i] for i in vertices])
	E_min = e_min
    Ising_MC(vertices, β, J, h, σ, E, E_min)
end

# MCMC step functions
function vertex_energy(ising::Ising_MC, vertex::Tuple{Int, Int, Int})
    e = -0.5 * ising.σ[vertex] * sum([ising.σ[edge.dst] * edge.weight for edge 
	in ising.J[vertex]]) + ising.h[vertex] * ising.σ[vertex]
    return e
end 
   
   	function flip_energy(ising::Ising_MC, vertex::Tuple{Int, Int, Int})
       	E = 2*vertex_energy(ising, vertex) - ising.h[vertex] * ising.σ[vertex]
       	dE = -2E
       	return dE
   	end 
   
   	function metropolis_step!(ising::Ising_MC)
       	rand_vertex = rand(rng, ising.vertices)
       	dE = flip_energy(ising, rand_vertex)
       	if dE < 0 || rand(rng) < exp(-ising.β * dE)
        	ising.σ[rand_vertex] *= -1
			ising.E += dE
		end
   	end

	# commented out: nonzero field adaptation for Swendsen-Wang
  	function identify_clusters(ising::Ising_MC)
       	clusters = []
       	visited = Dict(i => false for i in ising.vertices)
       	#ghosts = []

       	for vertex in ising.vertices
           	if !visited[vertex]
               	cluster = []
               	queue = [vertex]
               	visited[vertex] = true
               	#ghost = false

               	while !isempty(queue)
                   	current_vertex = popfirst!(queue)
                   	push!(cluster, current_vertex)

                   	for neighbour in ising.J[current_vertex]
						neighbour_vertex = neighbour.dst
						if !visited[neighbour_vertex]
							if ising.σ[current_vertex] == ising.σ[neighbour_vertex] && rand(rng) < 1 - exp(-2 * ising.β * neighbour.weight)
                           		push!(queue, neighbour_vertex)
                           		visited[neighbour_vertex] = true
							end
                       	end
                   	end

                   #=	if ising.σ[current_vertex] == -sign(ising.h[current_vertex]) && 
					!ghost
                       	if rand(rng) < 1 - exp(-2 * ising.β * 
						abs(ising.h[current_vertex]))
                           	ghost = true
                       	end
                   	end=#
               	end

               	push!(clusters, cluster)
               #	push!(ghosts, ghost)
           	end
       	end
           
       	#return [cluster for (cluster, ghost) in zip(clusters, ghosts) if !ghost]
		return clusters
   	end
   
   	function swendsen_wang!(ising::Ising_MC)
   	clusters = identify_clusters(ising)   
       	for cluster in clusters
           	if rand(rng) < 0.5
               	for vertex in cluster
					dE = flip_energy(ising, vertex)
                   	ising.σ[vertex] *= -1
					ising.E += dE
               	end #could improve efficency by calculating changes due to each bond, neglecting within-cluster changes
           	end
       	end
   	end

   	function hybrid_step!(ising::Ising_MC, p::Float64)
		n_v = length(ising.vertices)
		if rand(rng) < p/(n_v+p)
       		swendsen_wang!(ising)
		else
       		metropolis_step!(ising)
       	end
   	end

# Function for recording MCMC observations
function MCMC_observables!(ising::Ising_MC, N::Int, p::Float64; 
	burn::Int = 5000)
	mc_steps = 0
	good_steps = N - burn
	s = Vector{Float64}(undef, good_steps)
    E = Vector{Float64}(undef, good_steps)
	for i in 1:burn
		mc_steps += 1
		if p >= 0
			if p == 1 && i == 1
				swendsen_wang!(ising)
			else
				hybrid_step!(ising, p)
			end
		else
			swendsen_wang!(ising)
		end
	end
	for i in 1:good_steps
		mc_steps += 1
		if p >= 0
			hybrid_step!(ising, p)
		else
			swendsen_wang!(ising)
		end
		s[i] = sum(values((ising.σ)))
        E[i] = ising.E
	end
	return Ising_MC_observables(ising.β, s, E)	
end

# Functions for computing observables
function thermodynamics_global_exa(ising::Ising_exact)
	U = sum([sum(values(ising.e[j]))*ising.P[j] for j in 1:length(ising.P)])
	U_sq = sum([(sum(values(ising.e[j])))^2*ising.P[j] for j in 1:length(ising.P)])
	Cᵥ = k_B*ising.β^2*(U_sq-U^2)
	return U, Cᵥ
end

function magnetics_global_exa(ising::Ising_exact)
	M = sum([sum(values(ising.σ[j]))*ising.P[j] for j in 1:length(ising.P)])
	M_abs = sum([abs(sum(values(ising.σ[j]))*ising.P[j]) for j in 1:length(ising.P)])
	return M, M_abs
end

function thermodynamics_global_mc(obs::Ising_MC_observables)
	E_steps = obs.E
    E = mean(E_steps)
    E_sq = mean(E_steps.^2)
    Cᵥ = k_B*obs.β^2*(E_sq-E^2)
    return E, Cᵥ
end

function magnetics_global_mc(obs::Ising_MC_observables)
    M_steps = obs.s
    M = mean(M_steps)
    M_abs = mean(abs.(M_steps))
    return M, M_abs
end

# Functions for running multiple simulations and outputting plots
function phasediagrams(vertices::Array{Tuple{Int, Int, Int}}, β_vec, J_dict::Dict{Tuple{Int, Int, Int}, Set{Edge}}, h::Dict{Tuple{Int, Int, Int}, Float64}; N = 10000, N_multi = false, p=1.0, burn = 5000, burn_multi = false)
	n_v = length(vertices)
	n_e = sum(J_dict_to_adj_matrix(vertices, J_dict))/2
	if n_v <= 14
		ising_exa_arr = Ising_exact.(Ref(vertices), β_vec, Ref(J_dict), Ref(h)) 
		magnetics_data = magnetics_global_exa.(ising_exa_arr)
		thermodynamics_data = thermodynamics_global_exa.(ising_exa_arr)
		
		c_vec = getindex.(thermodynamics_data, 2)/(k_B*n_e)
		c_plt = plot(β_vec, c_vec, xlabel = L"K", ylabel = L"c", label = false, tickfont = ticksize, guidefont = labelsize, legendfontsize = legendsize, legend = legendpos, grid = false, framestyle=:box, color = :black)

		m_vec = getindex.(magnetics_data, 2)/n_v
		m_plt = plot(β_vec, m_vec, xlabel = L"K", ylabel = L"\chi", label = false, tickfont = ticksize, guidefont = labelsize, legendfontsize = legendsize, legend = legendpos, grid = false, framestyle=:box, color = :black)

		u_vec = getindex.(thermodynamics_data, 1)/n_e
		u_plt = plot(β_vec, u_vec, xlabel = L"K", ylabel = L"u", label = false, tickfont = ticksize, guidefont = labelsize, legendfontsize = legendsize, legend = legendpos, grid = false, framestyle=:box, color = :black)
		
		plt_data = [β_vec, c_vec, m_vec, u_vec]
		return c_plt, m_plt, u_plt, plt_data
	else	
		if burn_multi
			burn = 200*n_v
		end
		if N_multi
			N = 2000*n_v
		end
		ising_mc_vec = Ising_MC.(Ref(vertices), β_vec, Ref(J_dict), Ref(h)) 	
		obs_mc_vec = MCMC_observables!.(ising_mc_vec, Ref(N), Ref(p), burn = burn)
		thermodynamics_data = thermodynamics_global_mc.(obs_mc_vec)
		
		c_vec = getindex.(thermodynamics_data, 2)/(k_B*n_e)
		c_plt = plot(β_vec, c_vec, xlabel = L"K", ylabel = L"c", label = false, tickfont = ticksize, guidefont = labelsize, legendfontsize = legendsize, legend = legendpos, grid = false, framestyle=:box, color = :black)

		m_vec =  getindex.(magnetics_global_mc.(obs_mc_vec), 2)/n_v
		m_plt = plot(β_vec, m_vec, xlabel = L"K", ylabel = L"\chi", label = false, tickfont = ticksize, guidefont = labelsize, legendfontsize = legendsize, legend = legendpos, grid = false, framestyle=:box, color = :black)

		u_vec = getindex.(thermodynamics_data, 1)/n_e
		u_plt = plot(β_vec, u_vec, xlabel = L"K", ylabel = L"u", label = false, tickfont = ticksize, guidefont = labelsize, legendfontsize = legendsize, legend = legendpos, grid = false, framestyle=:box, color = :black)
		plt_data = [β_vec, c_vec, m_vec, u_vec]
		return c_plt, m_plt, u_plt, plt_data
	end
end

function isingsequence(verticesseq::Vector{Array{Tuple{Int, Int, Int}, 3}}, graphseq::Vector{Matrix{Float64}}, β_vec, h_0::Float64)
	runs = length(graphseq)
	J_dict_vec = [adj_matrix_to_J_dict(verticesseq[i], graphseq[i]) for i in 1:runs]
	sequence_df = DataFrame(heatcapacityplots = [], absmplots = [], energyplots = [], plots_data = [])

   	for j in 1:runs
    	J_dict = J_dict_vec[j]
		vertices = verticesseq[j]
		h = h_dict(vertices, h_0)
		
		c_plot, m_plot, u_plot, plot_data = phasediagrams(vertices, β_vec, J_dict, h; N_multi = true, burn_multi = true)		
		sequence_row = (heatcapacityplots = c_plot, absmplots = m_plot, energyplots = u_plot, plots_data = plot_data)
		push!(sequence_df, sequence_row)
	end
	return sequence_df
end

function isingsequencevisualiser(sequencedata::DataFrame)
	overlays_data = sequencedata.plots_data
	c_overlays_plot = plot(xlabel = L"K", ylabel = L"c", guidefont = labelsize, legendfontsize = legendsize, tickfont = ticksize, legend = legendpos, grid = false, framestyle=:box)
	m_overlays_plot = plot(xlabel = L"K", ylabel = "abs(m)", guidefont = labelsize, legendfontsize = legendsize, tickfont = ticksize, legend = legendpos, grid = false, framestyle=:box)
	u_overlays_plot = plot(xlabel = L"K", ylabel = L"u", guidefont = labelsize, legendfontsize = legendsize, tickfont = ticksize, legend = legendpos, grid = false, framestyle=:box)
	overlays = length(overlays_data)
	overlay_cmap = cgrad(:blues, overlays)
	
	for i in 1:overlays
		colour = overlay_cmap[i]
		β_data, c_data, m_data, u_data = overlays_data[i]
		plot!(c_overlays_plot, β_data, c_data, label = "Graph $i", color = colour)
		plot!(m_overlays_plot, β_data, m_data, label = "Graph $i", color = colour)
		plot!(u_overlays_plot, β_data, u_data, label = "Graph $i", color = colour)
	end

	return c_overlays_plot, m_overlays_plot, u_overlays_plot
end	

# Example output for some complete graphs
complete_n_vec = 4:4:20
complete_dims_vec = [(1, n) for n in complete_n_vec]
complete_vertices_vec = vertex_indices.(complete_dims_vec)
complete_adj_matrices = adjacency_matrix_complete.(complete_n_vec)
complete_data = isingsequence(complete_vertices_vec, complete_adj_matrices, 0:0.01:8, 0.0)
complete_output = isingsequencevisualiser(complete_data)
