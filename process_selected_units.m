function [GROUP] = process_selected_units(pre_FR,post_FR,units,removal)
    % extract FR
    GROUP.FR{1} = pre_FR(units,:);
    GROUP.FR{2} = post_FR(units,:);

    % find M FR
    GROUP.M_FR{1} = mean(GROUP.FR{1},2);
    GROUP.M_FR{2} = mean(GROUP.FR{2},2);

    % find modulation index
    GROUP.mod_idx = log10(GROUP.M_FR{2}./GROUP.M_FR{1});
    GROUP.mod_idx(GROUP.mod_idx == -Inf) = NaN;

    % find thinned FR for pre and post
    for cond = 1:2
        [GROUP.thin_FR{cond}] = find_thinned_FR(GROUP.FR{cond},removal);
        GROUP.M_thin_FR{cond} = mean(GROUP.thin_FR{cond},2);
        GROUP.M_thin_FR{cond} = GROUP.M_thin_FR{cond}./(1-removal); % correct for left shift
    end

    % find thinned mod index
    GROUP.thin_mod_idx = log10(GROUP.M_thin_FR{2}./GROUP.M_thin_FR{1});
    GROUP.thin_mod_idx(GROUP.thin_mod_idx == -Inf) = NaN;
    
    % log transform M_FR
    %for cond = 1:2
    %    GROUP.M_FR{cond} = log10(GROUP.M_FR{cond});
    %    GROUP.M_thin_FR{cond} = log10(GROUP.M_thin_FR{cond}); 
    %end
     
    
end