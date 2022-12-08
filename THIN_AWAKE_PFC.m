addpath('X:\cortical_dynamics\Shared\Code\matlib\stats');
addpath('X:\cortical_dynamics\User\ms1121\Code\General');
run('makedb_TCB2_MS');

% set parameters
gap = 1*60*1000; % 1 min in ms
section = 5*60*1000; %8*60*1000; % 8 min in ms
bin_edges = [-1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 1.75 2.25]; % log FR bins
removal = 0.9; % 90% thinning

% find exp for each drug condition
tcb_exp = []; con_exp = [];
for exp = Batch3PFC%[Batch1PFC Batch2PFC]
    if db(exp).dose > 0 && ~strcmp(db(exp).animal,'M191023_A_BD')
        tcb_exp = [tcb_exp exp];
    elseif db(exp).dose == 0 && ~strcmp(db(exp).animal,'M191023_A_BD')
        con_exp = [con_exp exp];
    end
end
ket_control = false;

%% Extract FR for all selected exp
pre_FR = []; post_FR = [];
pre_FR_sample2 = []; post_FR_sample2 = [];
for exp = tcb_exp 
        % load spikestruct
        disp(['Exp: ' num2str(exp) ' Loading spikestruct...']);
        [spikestruct] = load_spikestruct('X:',db,exp);
        
        % ketamine control needs saline condition
        if ket_control == true
            db(exp).cond(2) = 2; % TEMP override for saline condition
        end
        
        % check length of conditions
        for c = 1:2
            if size(spikestruct.condspikevector{db(exp).cond(c)},2) < (2*gap + 2*section)
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
clear pre_M_FR_sample2 pre_FR_sample2 post_M_FR_sample2 post_FR_sample2

%% find mean FR
pre_M_FR = mean(pre_FR,2);
pre_M_FR(pre_M_FR==0) = NaN;
post_M_FR = mean(post_FR,2);
post_M_FR(post_M_FR==0) = NaN;
mod_idx = log10(post_M_FR./pre_M_FR);

%% 
figure
T = tiledlayout(3,4);
set(gcf,'color','w');

nexttile(1)
plot_scatter_changeFR_all(pre_M_FR,post_M_FR,'r.');
xlabel('Pre FR'); ylabel('Post FR');
nexttile(2)
plot_cumulative_dist(pre_M_FR,'r','--');
hold on; plot_cumulative_dist(post_M_FR,'r','-');
hold off;
nexttile(3)
plot_gamma_dist(pre_M_FR,'--','r');
hold on; plot_gamma_dist(post_M_FR,'-','r');
hold off;


% Run for up regulated units
[UP_REG] = process_selected_units(pre_FR,post_FR,up_units,removal);
marker_colour = [0.8500 0.3250 0.0980]; y_limits = [-0.5 2.5];

nexttile(5)
plot_mod_idx(log10(UP_REG.M_FR{1}),UP_REG.mod_idx,bin_edges,marker_colour,y_limits); title('Actual');
nexttile(6)
plot_binned_mod_idx(log10(UP_REG.M_FR{1}),UP_REG.mod_idx,bin_edges,marker_colour,y_limits); title('Actual');
nexttile(9)
plot_mod_idx(log10(UP_REG.M_thin_FR{1}),UP_REG.thin_mod_idx,bin_edges,marker_colour,y_limits); title('Thinned');
nexttile(10)
plot_binned_mod_idx(log10(UP_REG.M_thin_FR{1}),UP_REG.thin_mod_idx,bin_edges,marker_colour,y_limits); title('Thinned');


% Run for down regulated units
[DOWN_REG] = process_selected_units(pre_FR,post_FR,down_units,removal);
marker_colour = [0.4940 0.1840 0.5560]; y_limits = [-2.5 0.5];

nexttile(7)
plot_mod_idx(log10(DOWN_REG.M_FR{1}),DOWN_REG.mod_idx,bin_edges,marker_colour,y_limits); title('Actual');
nexttile(8)
plot_binned_mod_idx(log10(DOWN_REG.M_FR{1}),DOWN_REG.mod_idx,bin_edges,marker_colour,y_limits); title('Actual');
nexttile(11)
plot_mod_idx(log10(DOWN_REG.M_thin_FR{1}),DOWN_REG.thin_mod_idx,bin_edges,marker_colour,y_limits); title('Thinned');
nexttile(12)
plot_binned_mod_idx(log10(DOWN_REG.M_thin_FR{1}),DOWN_REG.thin_mod_idx,bin_edges,marker_colour,y_limits); title('Thinned');

% Plot distribution of up and down units
nexttile(4)
plot_gamma_dist(UP_REG.M_FR{1},'--',[0.8500 0.3250 0.0980]);
hold on
plot_gamma_dist(UP_REG.M_FR{2},'-',[0.8500 0.3250 0.0980]);
plot_gamma_dist(DOWN_REG.M_FR{1},'--',[0.4940 0.1840 0.5560]);
plot_gamma_dist(DOWN_REG.M_FR{2},'-',[0.4940 0.1840 0.5560]);
title([num2str((sum(up_units)/numel(up_units))*100) ' % up units, ' num2str((sum(down_units)/numel(down_units))*100) ' % down units']);


%% Local Functions
function plot_scatter_changeFR_all(neuronFR_cond1,neuronFR_cond2,marker_style)
loglog(neuronFR_cond1,neuronFR_cond2,marker_style);
hline = refline(1,0); hline.Color = [0.3 0.3 0.3]; hline.LineStyle = '--';
xlim([10^-2 10^2]); ylim([10^-2 10^2]);
box off; axis square;
end

function plot_cumulative_dist(FR,marker_colour,marker_style)
f = cdfplot(FR);
f.Color = marker_colour; f.LineStyle = marker_style;
xlabel('Firing Rate'); ylabel('Cumulative Probability'); title('');
set(gca,'Xscale','log'); xlim([10^-2 10^2]);
grid off; box off; axis square;
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.0;
end
end

function plot_gamma_dist(FR,marker_style,marker_colour)
e = [-3:0.02:2];
FR_params = fitdist(FR, 'Gamma');
plot(e, log(10)*logammaPDF(FR_params.b, FR_params.a,e*log(10)),marker_style,'Color',marker_colour,'LineWidth', 1); % logammaPDF gives natural log, need log(10)*logammaPDF) to convert to log10
xlabel('Log Firing Rate'); ylabel('Probability');
xlim([-2 2])
box off; axis square;
end

function plot_mod_idx(pre_M_FR,mod_idx,bin_edges,marker_colour,y_limits)
    scatter(pre_M_FR,mod_idx,[],marker_colour,'.'); 
    xlabel('Baseline log FR'); ylabel('Mod Idx');
    for e = 1:numel(bin_edges)
       xline(bin_edges(e),'--');
    end
    xlim([-2 2]); ylim(y_limits);
    axis square; box off;
end

% plot swarmchart of binned FR for actual data
function plot_binned_mod_idx(pre_M_FR,mod_idx,bin_edges,marker_colour,y_limits)
    [mod_idx_bin] = bin_dataY_by_dataX(bin_edges,pre_M_FR,mod_idx);
    mod_std = NaN(1,numel(mod_idx_bin)); mod_M = NaN(1,numel(mod_idx_bin));
    for b = 1:numel(mod_idx_bin)
        mod_std(b) = std(mod_idx_bin{b},'omitnan');
        mod_M(b) = mean(mod_idx_bin{b},'omitnan');
    end
    plot_swarmchart(bin_edges,mod_idx_bin,marker_colour);
    hold on
    errorbar(0.5*(bin_edges(1:end-1) + bin_edges(2:end)),mod_M,mod_std,'k','LineStyle','none','LineWidth',1); 
    hold off;
    xlabel('log FR'); ylabel('Mod Idx');
    ylim(y_limits);
    axis square; box off;
end

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

