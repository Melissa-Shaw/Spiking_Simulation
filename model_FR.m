function [spikeMat,M_FR,FR_accuracy] = model_FR(sim_length,actual_M_FR,num_samples,interval,removal)
    % create poisson distribution of firing rates
    spikeMat = poissrnd(actual_M_FR,1,sim_length);

    % randomly remove percentage of data
    if removal > 0
        num_removals = removal*numel(spikeMat);
        for r = 1:num_removals
            idx = randi([1 sim_length],1);
            spikeMat(idx) = 0;
        end
    end
    
    % find sample mean FR at randomised start points
    M_FR = []; count = 1;
    while count <= num_samples
        start = randi([1 sim_length],1);
        if start+interval <= sim_length
            M_FR = [M_FR mean(spikeMat(start:(start+interval)))];
            count = count + 1;
        end
    end

    % calc ratio of sample FR to actual FR
    FR_accuracy = M_FR./actual_M_FR;
end