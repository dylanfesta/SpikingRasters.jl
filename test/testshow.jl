using SpikingRasters
global const SR = SpikingRasters


t_end = 50.0
trains_test = SR.make_poisson_trains([4., 50., 3. , 4.],t_end)


collection_test = SR.SpikeTrainCollection(trains_test;t_end=t_end)

collection_test.t_start = 40.0


bao = SR.generate_spike_raster(collection_test, 5E-2;spike_size=5, spike_separator= 1)

t_end = 10.0
times = range(0.0, t_end; step=0.5)
bin_test = rand(5,length(times)) .< 0.03

collection_test = SR.SpikeTrainCollection(bin_test,times)

bao = SR.generate_spike_raster(collection_test, 0.5;spike_size=5, spike_separator= 1)