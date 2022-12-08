function [FR] = find_thinned_FR(FR,removal)

    num_units = size(FR,1);
    cond_length = size(FR,2);

    % find thinned mean FR
    num_removals = removal*cond_length;
    for n = 1:num_units
        idx_record = [];
        r = 1;
        while r <= num_removals
            idx = randi([1 cond_length],1);
            new_idx = true;
            for i = 1:numel(idx_record)
                if idx == idx_record(i)
                    new_idx = false;
                end
            end
            if new_idx == true
                FR(n,idx) = 0;
                idx_record = [idx_record idx];
                r = r + 1;
            end
        end
    end
    
    %M_thin_FR = mean(FR,2);
    
    % correct for left shift
    %M_thin_FR = M_thin_FR./(1-removal);
    
end