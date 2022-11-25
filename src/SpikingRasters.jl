module SpikingRasters
using Colors
using Random

mutable struct SpikeTrainCollection{R}
  N::Int
  row_indexes::Vector{Int}
  t_start::R
  t_end::R
  trains::Vector{Vector{R}}
end

function SpikeTrainCollection(trains::Vector{Vector{R}};
    t_start::R=zero(R), t_end::R=zero(R)) where {R} 
  N = length(trains)
  row_indexes = collect(1:N)
  t_start = zero(R)
  if iszero(t_end) 
    t_end = maximum(last.(trains))
  end
  return  SpikeTrainCollection(N, row_indexes, t_start, t_end, trains)
end

function SpikeTrainCollection(binary_trains::BitMatrix,times::AbstractVector{R}) where R
  t_start = first(times)
  t_end = last(times)
  spiketimes = Vector{R}[]
  for bin_row in eachrow(binary_trains)
    push!(spiketimes, times[findall(bin_row)])
  end
  return   SpikeTrainCollection(spiketimes;t_start=t_start,t_end=t_end)
end

"""
    count_spikes(Y::Vector{R},dt::R,Tend::R;Tstart::R=0.0) where R

# Arguments
  + `Y::Vector{<:Real}` : vector of spike times
  + `dt::Real` : time bin size
  + `Tend::Real` : end time for the raster
# Optional argument
  + `Tstart::Real=0.0` : start time for the raster

# Returns   
  + `binned_spikes::Vector{<:Integer}` : `binned_spikes[k]` is the number of spikes that occur 
      in the timebin `k`  (i.e. between `Tstart + (k-1)*dt` and `Tstart + k*dt`)
"""
function count_spikes(Y::Vector{R},dt::R,t_end::R;t_start::R=0.0) where R
  times = range(t_start,t_end;step=dt)  
  ret = fill(0,length(times)-1)
  for y in Y
    if t_start < y <= last(times)
      k = searchsortedfirst(times,y)-1
      ret[k] += 1
    end
  end
  return ret
end


"""
    detect_spikes(Y::Vector{R},dt::R,Tend::R;Tstart::R=0.0) where R

# Arguments
  + `Y::Vector{<:Real}` : vector of spike times
  + `dt::Real` : time bin size
  + `Tend::Real` : end time for the raster
# Optional argument
  + `Tstart::Real=0.0` : start time for the raster

# Returns   
  + `detected_spikes::Vector{<:Integer}` : `detected_spike[k]` is true if 
    at least one spike has been found.
"""
function detect_spikes(Y::Vector{R},dt::R,t_end::R;t_start::R=0.0) where R
  times = range(t_start,t_end;step=dt)  
  nt = length(times)-1
  ret = falses(nt)
  times = range(t_start,t_end;step=dt)  
  for y in Y
    if t_start < y <= last(times)
      k = searchsortedfirst(times,y)-1
      ret[k] = true
    end
  end
  return ret
end


"""
  generate_spike_raster(trains::SpikeTrainCollection,
      dt::Real;
      spike_size::Integer = 5,
      spike_separator::Integer = 1,
      background_color::Color=RGB(1.,1.,1.),
      spike_colors::Union{C,Vector{C}}=RGB(0.,0.0,0.0),
      max_size::Real=1E4) where C<:Color

Draws a matrix that contains the raster plot of the spike train.

# Arguments
  + `Trains` :  The spike trains. 
  + `dt` : time interval representing one pixel horizontally  

# Optional arguments
  + `spike_size::Integer` : heigh of spike (in pixels)
  + `spike_separator::Integer` : space between spikes, and vertical padding
  + `background_color::Color` : self-explanatory
  + `spike_colors::Union{Color,Vector{Color}}` : if a single color, color of all spikes, if vector of colors, 
     color for each neuron (length should be same as number of neurons)
  + `max_size::Integer` : throws an error if image is larger than this number (in pixels)

# Returns
  + `raster_matrix::Matrix{Color}` you can save it as a png file
"""
function generate_spike_raster(trains::SpikeTrainCollection{R},dt::R;
    spike_size::Integer = 5,
    spike_separator::Integer = 1,
    background_color::Color=RGB(1.,1.,1.),
    spike_colors::Union{C,Vector{C}}=RGB(0.,0.0,0.0),
    max_size::R=5E4) where {R,C<:Color}
  @assert trains.t_end > trains.t_start "t_end should be greater than t_start. But found: t_end=$(trains.t_end), t_start=$(trains.t_start)"
  nneus = trains.N
  if typeof(spike_colors) <: Color
    spike_colors = repeat([spike_colors,];outer=nneus)
  else
    @assert length(spike_colors) == nneus "error in setting colors"
  end
  binned_binary = detect_spikes.(trains.trains,dt,trains.t_end;t_start=trains.t_start)
  ntimes = length(binned_binary[1])
  ret = fill(background_color,
    (nneus*spike_size + # spike sizes
      spike_separator*nneus + # spike separators (incl. top/bottom padding) 
      spike_separator),ntimes)
  @assert all(size(ret) .< max_size ) "The image is too big! Please change time limits"  
  for (neu,binv,col) in zip(trains.row_indexes,binned_binary,spike_colors)
    spk_idx = findall(binv)
    _idx_pre = (neu-1)*(spike_size+spike_separator)+spike_separator
    y_color = _idx_pre+1:_idx_pre+spike_size
    ret[y_color,spk_idx] .= col
  end
  return ret
end


function make_poisson_train(rate::R,t_end::R ; t_start=zero(R)) where R
  n_estimate = round(Integer,1.2*(rate+2.0)*t_end+10.0) # preallocate heuristically
  ret = Vector{R}(undef,n_estimate)
  t_curr = zero(R)
  k_curr = 1
  while t_curr <= t_end
    Δt = -log(rand())/rate
    t_curr += Δt
    ret[k_curr] = t_curr
    k_curr += 1
  end
  keepat!(ret,1:(k_curr-2))
  return ret .+ t_start
end

function make_poisson_trains(rates::Vector{R},t_end::R ; t_start=zero(R)) where R
  return make_poisson_train.(rates,t_end;t_start=t_start)
end





# # do draw spikes on existing raster plot
# function add_to_spike_raster!(raster::Matrix,
#     trains::Vector{Vector{Float64}},neu_idxs::Vector{Int64},
#     dt::Real,Tend::Real, spike_color::C;
#       Tstart::Real=0.0,
#       spike_size::Integer = 5,
#       spike_separator::Integer = 1,
#       ) where C<:Color
#   @assert length(neu_idxs) == length(trains)    
#   binned_binary  = map(trains) do train
#     .! iszero.(bin_spikes(train,dt,Tend;Tstart=Tstart))
#   end
#   ntimes = length(binned_binary[1])
#   @assert size(raster,2) == ntimes "Mismatch in time binning!"
#   for (neu,binv) in zip(neu_idxs,binned_binary)
#     spk_idx = findall(binv)
#     _idx_pre = (neu-1)*(spike_size+spike_separator)+spike_separator
#     y_color = _idx_pre+1:_idx_pre+spike_size
#     raster[y_color,spk_idx] .= spike_color
#   end
#   return nothing
# end




end # of module
