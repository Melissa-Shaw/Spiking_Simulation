function [spikeMat,spike_raster,time_array] = create_poisson_spiking(FR,sim_length)
    dt = 1/1000; % ms bins in seconds
    nBins = floor(sim_length/dt); % round down to nearest integer
    spikeMat = poissrnd(FR,1,nBins); % nBins length vector of poisson distributed firing rates around lambda = FR
    spike_raster = spikeMat < FR*dt; % spike train 
    time_array = 0:dt:sim_length-dt; 
end