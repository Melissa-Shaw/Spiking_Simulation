addpath('X:\cortical_dynamics\User\ms1121\Code\General');
run('makedb_TCB2_MS');
clear Batch2V1 Batch3PFC AnaesV1 AnaesPFC AwakeV1

% set parameters
gap = 60*1000; % 60 seconds in ms
section = 8*60*1000; % 8 min in ms
bin_edges = [-1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 1.75 2.25]; % log FR bins
removal = 0.9; % 90% thinning

% find exp for each drug condition
tcb_exp = []; con_exp = [];
for exp = [Batch1PFC Batch2PFC]
    if db(exp).dose > 0 && ~strcmp(db(exp).animal,'M191023_A_BD')
        tcb_exp = [tcb_exp exp];
    elseif db(exp).dose == 0 && ~strcmp(db(exp).animal,'M191023_A_BD')
        con_exp = [con_exp exp];
    end
end

%% Extract FR for all selected exp
pre_FR = []; post_FR = [];
pre_FR_sample2 = []; post_FR_sample2 = [];
for exp = tcb_exp 
        % load spikestruct
        disp(['Exp: ' num2str(exp) ' Loading spikestruct...']);
        [spikestruct] = load_spikestruct('X:',db,exp);
        
        % check length of conditions
        for c = 1:2
            if size(spikestruct.condspikevector{db(exp).cond(c)},2) <= (2*gap + 2*section)
                disp(['NOTE: Exp: ' num2str(exp) ' Cond: ' num2str(c) ' is too short for set gap and section length.']);
            end
        end
        
        % find sample 1 FR for pre
        FR = [];
        end_cond = size(spikestruct.condspikevector{db(exp).cond(1)},2); disp(num2str(end_cond));
        pre_raster = spikestruct.condspikevector{db(exp).cond(1)}(:,end_cond-gap-section+1:end_cond-gap);
        num_units = size(pre_raster,1);
        for n = 1:num_units
            FR(n,:) = calc_running_sum(pre_raster(n,:),1000); % get FR in sp/s for each neuron
        end
        pre_FR = [pre_FR; FR];
        
        % find sample 1 FR for post
        FR = [];
        end_cond = size(spikestruct.condspikevector{db(exp).cond(2)},2); disp(num2str(end_cond));
        post_raster = spikestruct.condspikevector{db(exp).cond(2)}(:,end_cond-gap-section+1:end_cond-gap);
        for n = 1:num_units
            FR(n,:) = calc_running_sum(post_raster(n,:),1000); % get FR in sp/s for each neuron
        end
        post_FR = [post_FR; FR];
        
        % find sample 2 FR for pre
        FR = [];
        pre_raster = spikestruct.condspikevector{db(exp).cond(1)}(:,gap+1:gap+section);
        for n = 1:num_units
            FR(n,:) = calc_running_sum(pre_raster(n,:),1000); % get FR in sp/s for each neuron
        end
        pre_FR_sample2 = [pre_FR_sample2; FR];
        
        % find sample 2 FR for post
        FR = [];
        post_raster = spikestruct.condspikevector{db(exp).cond(2)}(:,gap+1:gap+section);
        for n = 1:num_units
            FR(n,:) = calc_running_sum(post_raster(n,:),1000); % get FR in sp/s for each neuron
        end
        post_FR_sample2 = [post_FR_sample2; FR];
        
        clear end_cond pre_raster post_raster FR num_units spikestruct n c
end
clear exp

% find up and down regulated units
pre_M_FR_sample2 = mean(pre_FR_sample2,2);
post_M_FR_sample2 = mean(post_FR_sample2,2);
up_units = pre_M_FR_sample2 <= post_M_FR_sample2;
down_units = pre_M_FR_sample2 > post_M_FR_sample2;

%% Run for up regulated units

% find thinned FR for pre and post
pre_M_FR = mean(pre_FR(up_units,:),2);
[pre_thin_FR,pre_M_thin_FR] = find_thinned_FR(pre_FR(up_units,:),removal);
post_M_FR = mean(post_FR(up_units,:),2);
[post_thin_FR,post_M_thin_FR] = find_thinned_FR(post_FR(up_units,:),removal);

% find modulation index
mod_idx = log10(post_M_FR./pre_M_FR);
mod_idx(mod_idx == -Inf) = NaN;
pre_M_FR = log10(pre_M_FR); post_M_FR = log10(post_M_FR);

% find thinned mod index
thin_mod_idx = log10(post_M_thin_FR./pre_M_thin_FR);
thin_mod_idx(thin_mod_idx == -Inf) = NaN;
pre_M_thin_FR = log10(pre_M_thin_FR);
post_M_thin_FR = log10(post_M_thin_FR);

% set up figure
figure; T = tiledlayout(2,2);
title(T,'Up-regulated Units');
set(gcf,'color','w'); clear T;

% plot FR vs mod_idx for actual data
nexttile;
scatter(pre_M_FR,mod_idx,[],[0.8500 0.3250 0.0980],'.'); 
xlabel('Baseline log FR'); ylabel('Mod Idx'); title('Actual');
for e = 1:numel(bin_edges)
   xline(bin_edges(e),'--');
end
xlim([-2 2]); ylim([-0.5 2.5]);

% plot swarmchart of binned FR for actual data
nexttile;
[mod_idx_bin] = bin_dataY_by_dataX(bin_edges,pre_M_FR,mod_idx);
mod_std = NaN(1,numel(mod_idx_bin)); mod_M = NaN(1,numel(mod_idx_bin));
for b = 1:numel(mod_idx_bin)
    mod_std(b) = nanstd(mod_idx_bin{b});
    mod_M(b) = nanmean(mod_idx_bin{b});
end
plot_swarmchart(bin_edges,mod_idx_bin,[0.8500 0.3250 0.0980]);
hold on
errorbar(0.5*(bin_edges(1:end-1) + bin_edges(2:end)),mod_M,mod_std,'k','LineStyle','none','LineWidth',1); 
hold off;
xlabel('log FR'); ylabel('Mod Idx'); title('Actual');
ylim([-0.5 2.5]);

% plot FR vs mod_idx for thinned data
nexttile;
scatter(pre_M_thin_FR,thin_mod_idx,'b.'); 
xlabel('Baseline log FR'); ylabel('Mod Idx'); title('Thinned');
for e = 1:numel(bin_edges)
   xline(bin_edges(e),'--');
end
xlim([-2 2]); ylim([-0.5 2.5]);

% plot swarmchart of binned FR for thinned data
nexttile
[thin_mod_idx_bin] = bin_dataY_by_dataX(bin_edges,pre_M_thin_FR,thin_mod_idx);
thin_mod_std = NaN(1,numel(thin_mod_idx_bin)); thin_mod_M = NaN(1,numel(thin_mod_idx_bin));
for b = 1:numel(thin_mod_idx_bin)
    thin_mod_std(b) = nanstd(thin_mod_idx_bin{b});
    thin_mod_M(b) = nanmean(thin_mod_idx_bin{b});
end
plot_swarmchart(bin_edges,thin_mod_idx_bin,'b');
hold on
errorbar(0.5*(bin_edges(1:end-1) + bin_edges(2:end)),thin_mod_M,thin_mod_std,'k','LineStyle','none','LineWidth',1); 
hold off;
xlabel('log FR'); ylabel('Mod Idx'); title('Thinned');
ylim([-0.5 2.5]);

clear e

%% Run for down regulated units

% find thinned FR for pre and post
pre_M_FR = mean(pre_FR(down_units,:),2);
[pre_thin_FR,pre_M_thin_FR] = find_thinned_FR(pre_FR(down_units,:),removal);
post_M_FR = mean(post_FR(down_units,:),2);
[post_thin_FR,post_M_thin_FR] = find_thinned_FR(post_FR(down_units,:),removal);

% find modulation index
mod_idx = log10(post_M_FR./pre_M_FR);
mod_idx(mod_idx == -Inf) = NaN;
pre_M_FR = log10(pre_M_FR); post_M_FR = log10(post_M_FR);

% find thinned mod index
thin_mod_idx = log10(post_M_thin_FR./pre_M_thin_FR);
thin_mod_idx(thin_mod_idx == -Inf) = NaN;
pre_M_thin_FR = log10(pre_M_thin_FR);
post_M_thin_FR = log10(post_M_thin_FR);

% set up figure
figure; T = tiledlayout(2,2);
title(T,'Down-regulated Units');
set(gcf,'color','w'); clear T;

% plot FR vs mod_idx for actual data
nexttile;
scatter(pre_M_FR,mod_idx,[],[0.4940 0.1840 0.5560],'.'); 
xlabel('Baseline log FR'); ylabel('Mod Idx'); title('Actual');
for e = 1:numel(bin_edges)
   xline(bin_edges(e),'--');
end
xlim([-2 2]); ylim([-2.5 0.5]);

% plot swarmchart of binned FR for actual data
nexttile;
[mod_idx_bin] = bin_dataY_by_dataX(bin_edges,pre_M_FR,mod_idx);
mod_std = NaN(1,numel(mod_idx_bin)); mod_M = NaN(1,numel(mod_idx_bin));
for b = 1:numel(mod_idx_bin)
    mod_std(b) = nanstd(mod_idx_bin{b});
    mod_M(b) = nanmean(mod_idx_bin{b});
end
plot_swarmchart(bin_edges,mod_idx_bin,[0.4940 0.1840 0.5560]);
hold on
errorbar(0.5*(bin_edges(1:end-1) + bin_edges(2:end)),mod_M,mod_std,'k','LineStyle','none','LineWidth',1); 
hold off;
xlabel('log FR'); ylabel('Mod Idx'); title('Actual');
ylim([-2.5 0.5]);

% plot FR vs mod_idx for thinned data
nexttile;
scatter(pre_M_thin_FR,thin_mod_idx,'b.'); 
xlabel('Baseline log FR'); ylabel('Mod Idx'); title('Thinned');
for e = 1:numel(bin_edges)
   xline(bin_edges(e),'--');
end
xlim([-2 2]); ylim([-2.5 0.5]);

% plot swarmchart of binned FR for thinned data
nexttile
[thin_mod_idx_bin] = bin_dataY_by_dataX(bin_edges,pre_M_thin_FR,thin_mod_idx);
thin_mod_std = NaN(1,numel(thin_mod_idx_bin)); thin_mod_M = NaN(1,numel(thin_mod_idx_bin));
for b = 1:numel(thin_mod_idx_bin)
    thin_mod_std(b) = nanstd(thin_mod_idx_bin{b});
    thin_mod_M(b) = nanmean(thin_mod_idx_bin{b});
end
plot_swarmchart(bin_edges,thin_mod_idx_bin,'b');
hold on
errorbar(0.5*(bin_edges(1:end-1) + bin_edges(2:end)),thin_mod_M,thin_mod_std,'k','LineStyle','none','LineWidth',1); 
hold off;
xlabel('log FR'); ylabel('Mod Idx'); title('Thinned');
ylim([-2.5 0.5]);

clear e

%% LOCAL FUNCTIONS

function [binned_data] = bin_dataY_by_dataX(bin_edges,x,y)
    [~,~,bin] = histcounts(x,bin_edges);
    num_bins = numel(bin_edges)-1;
    binned_data = cell(1,num_bins);
    for b = 1:num_bins % for each bin
        binned_data{b} = y(bin==b);
    end
end

function plot_swarmchart(bin_edges,bin_mod_idx,marker_colour)
    mid_bin = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));
    x = []; y = [];
    for b = 1:(numel(bin_edges)-1)
        x = [x ; repelem(mid_bin(b),numel(bin_mod_idx{b}))'];
        y = [y ; bin_mod_idx{b}];
    end
    swarmchart(x,y,10,marker_colour,'filled');
    for e = 1:numel(bin_edges)
        xline(bin_edges(e),'--');
    end
end

