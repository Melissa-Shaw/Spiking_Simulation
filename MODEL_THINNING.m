%% Set parameters
sim_length = 10*60*60; % length of simulation in seconds
num_samples = 200; % number of means to be taken
sample_length = 8*60; % sample length
mod_factor = 5; % set multiplicative modulation factor
removal = 0.5; % 90% data removal
actual_M_FR = [0.1 1 10 100]; % array of set mean FRs

%% Model and plot sample FR to actual FR accuracy
figure
tiledlayout(1,4)
set(gcf,'color','w');

axes = [];
for interval = [5 10 20 60]
    ax1 = nexttile;
    for FR = actual_M_FR
        [~,~,FR_accuracy] = model_FR(sim_length,FR,num_samples,interval*60,0);

        scatter(repmat(FR,[1 numel(FR_accuracy)]),FR_accuracy,'k.');
        hold on;
    end
    yline(1,'--'); xlabel('Sample'); ylabel('Sample FR to Actual FR Ratio'); title(['Sample length: ' num2str(interval) ' min']);
    box off; set(gca,'XScale','log'); xlim([0.01 100]);
    axes = [axes ax1];
end
linkaxes(axes,'y');
clear axes ax1

%% Model and plot modulated firing rates
figure
set(gcf,'color','w');
tiledlayout(1,2);

ax1 = nexttile;
i = 1;
for FR = actual_M_FR
    [pre_spikeMat(i,:),pre_FR(i,:),~] = model_FR(sim_length,FR,num_samples,sample_length,0);
    [post_spikeMat(i,:),post_FR(i,:),~] = model_FR(sim_length,FR*mod_factor,num_samples,sample_length,0); % increased M_FR by factor of 5
    mod_ind(i,:) = post_FR(i,:)./pre_FR(i,:);
    scatter(pre_FR(i,:),mod_ind(i,:),'k.');
    hold on;
    i = i+1;
end
yline(mod_factor,'--');
hold off; xlabel('Mean FR'); ylabel('Modulation Index'); title('No Thinning');
set(gca,'XScale','log'); xlim([10^-2 10^3]); ylim([0 10]);

ax2 = nexttile;
[thin_pre_FR] = find_thinned_FR(pre_spikeMat,removal);
[thin_post_FR] = find_thinned_FR(post_spikeMat,removal);


for i = 1:numel(actual_M_FR)
    %[thin_pre_FR(i,:)] = thin_FR(pre_spikeMat(i,:),num_samples,interval,removal);
    [M_thin_pre_FR(i,:)] = sample_M_FR(thin_pre_FR(i,:),num_samples,sample_length);
    M_thin_pre_FR = M_thin_pre_FR./(1-removal); % correct for left shift
    
    [M_thin_post_FR(i,:)] = sample_M_FR(thin_post_FR(i,:),num_samples,sample_length);
    M_thin_post_FR = M_thin_post_FR./(1-removal); 
    
    thin_mod_ind(i,:) = M_thin_post_FR(i,:)./M_thin_pre_FR(i,:);
    scatter(M_thin_pre_FR(i,:),thin_mod_ind(i,:),'k.');
    hold on;
end
yline(mod_factor,'--');
hold off; xlabel('Mean FR'); ylabel('Modulation Index'); title('Thinning');
set(gca,'XScale','log'); xlim([10^-2 10^3]); ylim([0 10]);
%linkaxes([ax1 ax2],'y');
clear ax1 ax2 i 









