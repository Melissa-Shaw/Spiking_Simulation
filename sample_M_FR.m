
function [M_FR] = sample_M_FR(FR,num_samples,interval)
    
    sim_length = size(FR,2);

    % find sample mean FR at randomised start points
    M_FR = []; count = 1;
    while count <= num_samples
        start = randi([1 sim_length],1);
        if start+interval <= sim_length
            M_FR = [M_FR mean(FR(start:(start+interval)))];
            count = count + 1;
        end
    end
    
end