% EXAMPLE CODE FOR PUBLICATION

clear
close all
clc

addpath('fit_functions');
addpath('fit_functions/third_party')
rng(10) %100 10

% PERSONALIZE THE SIMULATION
csi = 0.2; % noise scale factor
N_trials = 20; % number of noise realizations

% TIME VECTOR
dt = [  ones(12,1)*10; ones(2,1)*30;  ones(3,1)*60; ...
    ones(2,1)*120; ones(5,1)*300; ones(3,1)*600];
t = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
time = mean(t,2);

% RADIOACTIVE DECAY CONSTANT
dk = log(2)/(109.8);

% INPUT FUNCTION SIMULATION
inputFunParam = [0.535  ...                       % delay
    158.58819 7.61367 5.28450 ...    % A1 A2 A3
    -2.065 -0.0114 -0.048];          % lambda1 lambda2 lambda2
[Cp, ppCp,ifParams] = FengIF_simulation(time, inputFunParam);

% KINETIC PARAMETERS
%      K1     k2      k3     k4      fv
K = [.1104  .1910    .1024  .0094  .0446   % grey atter
    .0876  .1621    .0879  .0097  .0412   % whole brain
    .0622  .1248    .0700  .0097  .0270   % White matter
    .0640  .1272    .0738  .0081  .0415]; % tumor

% TISSUE TACs SIMULATION
N_tac = size(K,1);
TAC = zeros(length(time),N_tac);

for n_tac =1:N_tac
    k = K(n_tac,:);
    % the function simulation compartmental model is a function of the
    % auxiliary parameter set, as discussed in the paper, so we need to
    % convert the kinetic parameters K into the auxiliry set [alpha1,
    % beta1, alpha2, beta2]
    d    = abs(sqrt((k(2)+k(3)+k(4)).^2 - 4*k(2)*k(4)));
    p(2) = (k(2) + k(3) + k(4) + d) ./ 2;   %  beta1
    p(4) = (k(2) + k(3) + k(4) - d) ./ 2;   %  beta2
    p(1) = (k(1)*(p(2) - k(3) - k(4)))./ d; %  alpha1
    p(3) = (k(1)*(-p(4) + k(3) + k(4)))./ d; % alpha2
    p(5) = k(5); % fv = k5
    
    [TAC(:,n_tac), J] = TwoTissueModel_simulation(p, t, inputFunParam, Cp);
end

% ADDING NOISE
Noisy_TAC = zeros(length(time),N_tac,N_trials);
Std = csi * sqrt( TAC./repmat(dt/60,[1,N_tac]) );  %  Time in min
for n_trial = 1:N_trials
    Noisy_TAC(:,:,n_trial) = TAC + Std.*randn(size(Std));
end
Noisy_TAC(Noisy_TAC<=0) = eps;

% VISUALIZE SIMULATION
figure
subplot(2,4,[2,3]),plot(time/60,Cp,'-','linewidth', 1.5)
xlim([0 60]),
xlabel('time [sec]'),
ylabel('emission activity [KBq/cc]')
title('Simulated Arterial Input Function'),
subplot(2,4,[5,6]),plot(time/60,TAC,'-','linewidth', 1.5)
xlim([0 60]), ylim([0 13.5])
xlabel('time [sec]'),
ylabel('emission activity [KBq/cc]')
legend({'Gray matter','Brain','White matter', 'Tumor'},'Location','best')
title('Simulated tissue TACs')
subplot(2,4,[7,8]),plot(time/60,squeeze(Noisy_TAC(:,:,1)),'-','linewidth', 1.5)
xlim([0 60]), ylim([0 13.5])
xlabel('time [sec]'),
ylabel('emission activity [KBq/cc]')
legend({'Gray matter','Brain','White matter', 'Tumor'},'Location','best')
title('Noisy tissue TACs')

% FITTING
time_NN   = zeros(N_tac,N_trials);
fit_NN    = zeros(size(time,1),N_tac,N_trials);
params_NN = zeros(N_tac,6,N_trials);
time_AN   = time_NN;
fit_AN    = fit_NN;
params_AN = params_NN;
time_AA   = time_NN;
fit_AA    = fit_NN;
params_AA = params_NN;

for n_tac = 1:N_tac % for each curve,
    for n_trial = 1:N_trials % fit each noise realization
        
        tac = squeeze(Noisy_TAC(:,n_tac,n_trial));
        
        % NM-NC implementation
        tic
        [fit, params] = fit_NMNC(tac,time,Cp);
        time_NN(n_tac,n_trial) = toc;
        fit_NN(:,n_tac,n_trial) = fit;
        params_NN(n_tac,:,n_trial) = params;
        
        % AM-NC implementation
        tic
        [fit, params] = fit_AMNC(tac,time,Cp);
        time_AN(n_tac,n_trial) = toc;
        fit_AN(:,n_tac,n_trial) = fit;
        params_AN(n_tac,:,n_trial) = params;
        
        % AM-AC implementation
        tic
        [fit, params] = fit_AMAC(tac,time,inputFunParam,Cp);
        time_AA(n_tac,n_trial) = toc;
        fit_AA(:,n_tac,n_trial) = fit;
        params_AA(n_tac,:,n_trial) = params;
    end
end

tit = {'GREY MATTER','WHOLE BRAIN','WHITE MATTER', 'TUMOR'};
figure,
for c = 1:size(Noisy_TAC,2)
    subplot(2,2,c)
    hold on
    plot(time/60,TAC(:,c),'k')
    plot(time/60,Noisy_TAC(:,c,1),'*-k')
    errorbar(time/60,mean(fit_NN(:,c,:),3),std(fit_NN(:,c,:),0,3),'b')
    errorbar(time/60,mean(fit_AN(:,c,:),3),std(fit_AN(:,c,:),0,3),'r')
    errorbar(time/60,mean(fit_AA(:,c,:),3),std(fit_AA(:,c,:),0,3),'g')
    ylim([0,15])
    xlabel('time [sec]','FontSize',12),
    ylabel('emission activity [KBq/cc]','FontSize',12)
    legend({'real','noisy','NM-NC','AM-NC','AM-AC'},'Location','best')
    title(tit(c))
end
set(gcf,'Position',[0 0 900 700])

fprintf('\n\n============================================================================\n')
fprintf('COMPUTATIONAL TIME TO FIT ONE CURVE [s]\n')
fprintf('============================================================================\n')
fprintf('                GM               Br                WM              Tumor        \n')
fprintf('  NM-NC:   %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)      \n',...
    [mean(time_NN,2) , std(time_NN,0,2)]')
fprintf('  AM-NC:   %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)      \n',...
    [mean(time_AN,2) , std(time_AN,0,2)]')
fprintf('  AM-AC:   %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)      \n',...
    [mean(time_AA,2) , std(time_AA,0,2)]')
fprintf('============================================================================\n')


for c = 1:size(Noisy_TAC,2)
    fprintf('\n\n===============================================================================================================\n')
    fprintf('PARAMETERS FOR %s REGION\n',tit{c})
    fprintf('===============================================================================================================\n')
    fprintf('                K1               k2               k3                k4              fv              Ki \n')
    fprintf('  NM-NC:   %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)      \n',...
        [mean(params_NN(c,:,:),3,'omitnan') ; std(params_NN(c,:,:),0,3,'omitnan')])
    fprintf('  AM-NC:   %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)      \n',...
        [mean(params_AN(c,:,:),3,'omitnan') ; std(params_AN(c,:,:),0,3,'omitnan')])
    fprintf('  AM-AC:   %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)    %.4f(%.3f)      \n',...
        [mean(params_AA(c,:,:),3,'omitnan') ; std(params_AA(c,:,:),0,3,'omitnan')])
    fprintf('===============================================================================================================\n')
    
end








