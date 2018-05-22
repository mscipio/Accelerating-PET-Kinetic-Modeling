%|=================================================================================
%|   PART-ANALYTIC (MODEL) PART-NUMERIC (CONVOLUTION) IMPLEMENTATION 
%|   OF FITTING OF A TWO-TISSUE COMPARTMENT MODEL
%|
%|   INPUTS:
%|       tac:       time series to fit [time_points x 1]
%|       scanTime:  vector of start and end times of each frame [time_points x 2]
%|       IF:        arterial inputn function [time_points x 1]
%|
%|   OUTPUTS:
%|       fit:       fitted time curve [time_points x 1]
%|       params:    parameter vector (K1,k2,k3,k4,fv,Ki) [6 x 1]
%|       JAC:       jacobian [time_points x 6]
%|       resnorm:   sum of squared residuals [1 x 1]
%|       residual:  residual vector [time_points x 1]
%|       output:    structure with summary info about the fitting process
%|
%|  Last revision:
%|  22 May 2018
%|  Michele Scipioni, Univeristy of Pisa
%|
%|=================================================================================


function [fit, par_out, JAC, resnorm, residual, output] = fit_AMNC(tac,scanTime,IF) 

warning off;

if max(scanTime(:))>180
    scanTime = scanTime./60; % time has to be in minutes
end

% create options
options = optimoptions('lsqcurvefit','Display', 'none');
% options = optimoptions('lsqcurvefit');

options.TolFun=1e-12;  %1e-6;
options.TolX=1e-12;          % 1e-6;
options.MaxFunEvals=1e6;
options.MaxIter=1e6;
% options.PlotFcn = @optimplotresnorm;
lb    = [0.   0.   0.   0.  0.]; % all parameters should be at least positive
ub    = [1.   1.8   1.   1.  1.]; 
par0  = [0.1  0.1  0.01  0.01  0.01 ];
fixed = [0    0    0     0     0];

func = @(k,t) numeric_conv(k, scanTime, IF);
[params,JAC,resnorm,residual,~,output] = fit_nl(func, par0, scanTime, tac, fixed, lb, ub, options);

n  = params(1)*params(2) + params(3)*params(4);
d  = params(1) + params(3);
par_out(5) = params(5);
par_out(1) = d;
par_out(2) = n / d;
par_out(3) = (params(1) * params(3) * (params(2) - params(4))^2 ) / (d * n);
par_out(4) = (params(2) * params(4) * (params(1) + params(3))) / n;

Ki = par_out(1)*par_out(3)/(par_out(2)+par_out(3)); %k1*k3/k2+k3
par_out = [par_out Ki]; % k1 k2 k3 k4 vB Ki

fit = numeric_conv(params, scanTime, IF);

end

%% MODELLO COMPARTIMENTALE COMPLETO
function Ct = numeric_conv(k, scanTime, Cp)

dk = log(2)/109.8;  % radioactive decay constant for F-18

% use fine time sampling and implement the convolution
% using discrete mathematics.
if max(scanTime(:))>180
    scanTime = scanTime./60; % time has to be in minutes
end
time = mean(scanTime,2); %[scanTime(1,1); scanTime(:,2)];
dt = .01;  % time step, in minutes    
t = (0:dt:time(end))';
% t = (scanTime(1,1):dt:scanTime(end,2))';
IF = max(0,interp1(time,Cp,t,'linear',0) );

% Convolution with input function
sol = conv(modello_exp(k, t), IF)*dt;
tsol = (0:dt:2*max(t))';

% Taking into account natural radioactive decay
sol = sol .* exp(-dk * tsol);

% Taking into account blood fraction in tissue
sol = max(0,interp1(tsol,sol,time,'linear',0) );
Ct = (1-k(5))*sol + k(5)* Cp;
Ct(Ct<=0) = eps;

end

%% SOLUZIONE ANALITICA SISTEMA ODE
function Ct = modello_exp(k, time)


% d      = sqrt((k(2)+k(3)+k(4)).^2 - 4*k(2)*k(4));
% beta1  = ((k(2)+k(3)+k(4))+d)/2;
% beta2  = ((k(2)+k(3)+k(4))-d)/2;
% alpha1 = (k(3)+k(4)-beta1)*(k(1)/(beta2-beta1));
% alpha2 = (beta2-k(3)-k(4))*(k(1)/(beta2-beta1));

% Ct = (time >= 0) .* ( alpha1*exp(-beta1*time)+alpha2*exp(-beta2*time) );
Ct = (time >= 0) .* ( k(1)*exp(-k(2)*time)+k(3)*exp(-k(4)*time) );

end

