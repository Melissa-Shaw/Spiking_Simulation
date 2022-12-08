% Function to model and plot accuracy of measured FR from poisson modelled
% FRs over different time intervals.
    % Inputs:
        % sim_length --> length in seconds of model (min 10 hours) 
                        %(e.g. sim_length = 10*60*60 for 10 hours)
        % actual_M_FR --> array of values for actual mean FR for the model (lambda) 
                        %(e.g. actual_M_FR = [0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 50 100];)
    % Outputs:
        % spikeMat --> firing rates per second (column) for poisson neuron with corresponding actual_M_FR (row)
        % measured_M_FR --> matrix of measured mean firing rates (row) for different time intervals (column)
        % interval --> array of time intervals used for comparison
        
function [spikeMat,measured_M_FR,interval] = model_time_intervals(sim_length,actual_M_FR)

% create poisson distribution of firing rates
i = 1;
for FR = actual_M_FR
    spikeMat(i,:) = poissrnd(FR,1,sim_length); % 1 second bins of sim_length vector of poisson distributed firing rates around lambda (FR)
    i = i+1;
end

% plot measured mean FRs compared to actual mean FRs for different time intervals
figure; tiledlayout('flow');
set(gcf,'color','w'); set(gcf, 'Position', get(0, 'Screensize'));
interval = [1*60 5*60 20*60 60*60 5*60*60 10*60*60]; % 1min, 5min, 20min, 1hr, 5hr, 10hr intervals for measuring mean FR
interval_names = {'1min' '5min' '20min' '1hr' '5hr' '10hr'};
symbol = {'*','square','o','+','^','x'};
i = 1;
measured_M_FR = [];
for int = interval
    loglog(actual_M_FR,mean(spikeMat(:,1:int),2),symbol{i},'LineWidth',1,'MarkerSize',10);
    hold on
    i = i+1;
    measured_M_FR = [measured_M_FR mean(spikeMat(:,1:int),2)];
end
legend(interval_names,'location','southeast'); legend box off;
hline = refline(1,0); hline.Color = [0.3 0.3 0.3]; hline.LineStyle = '--';
xlabel('Actual Mean FR'); ylabel('Measured Mean FR');
axis square; box off;

end