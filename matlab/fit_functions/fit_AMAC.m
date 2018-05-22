%|=================================================================================
%|   FULLY-ANALYTIC (model + convolution) IMPLEMENTATION 
%|   OF FITTING OF A TWO-TISSUE COMPARTMENT MODEL
%|
%|   INPUTS:
%|       tac:       time series to fit [time_points x 1]
%|       scanTime:  vector of start and end times of each frame [time_points x 2]
%|       IFparams:  input function model parameters 
%|                  (delay, A1, A2, A3, lambda1, lambda2, lambda3) [1 x 7]
%|       IFfit:     fitted arterial inputn function [time_points x 1]
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


function [fit, par_out, JAC, resnorm, residual, output] = fit_AMAC(tac,scanTime,IFparams, IFfit) 

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
options.SpecifyObjectiveGradient    =	false;
% options.PlotFcn = @optimplotresnorm;
lb    = [0.   0.   0.   0.  0.]; % all parameters should be at least positive
ub    = [1.   1.8   1.   1.  1.]; 
par0  = [0.1  0.1  0.01  0.01  0.01 ];
fixed = [0    0    0     0     0];

func = @(k,t) analytic_model(k, scanTime, IFparams, IFfit);
%  tic
[params,JAC,resnorm,residual,~,output] = fit_nl(func, par0, scanTime, tac, fixed, lb, ub, options);
%  toc
n  = params(1)*params(2) + params(3)*params(4);
d  = params(1) + params(3);
par_out(5) = params(5);
par_out(1) = d;
par_out(2) = n / d;
par_out(3) = (params(1) * params(3) * (params(2) - params(4))^2 ) / (d * n);
par_out(4) = (params(2) * params(4) * (params(1) + params(3))) / n;

Ki = par_out(1)*par_out(3)/(par_out(2)+par_out(3)); %k1*k3/k2+k3
par_out = [par_out Ki]; % k1 k2 k3 k4 vB Ki

fit = analytic_model(params, scanTime,  IFparams, IFfit);

end

%% MODELLO COMPARTIMENTALE COMPLETO
function [Ct, J] = analytic_model(p, scanTime, IFparams, IFfit)

if nargout > 1
    jacBool = 1;
else
    jacBool = 0;
end

if max(scanTime(:))>180
    scanTime = scanTime./60; % time has to be in minutes
end
time = mean(scanTime,2);%[scanTime(1,1); scanTime(:,2)];

delay = IFparams(1);
idx = find(time > delay);
t=time(idx)-delay;

A = IFparams(2:4); % [A1 A2 A3]
l = -IFparams(5:7); % [lambda1 lambda2 lambda3] has to be > 0

irf = zeros(size(time,1),1);

J = zeros(size(time,1),5);

for i = 1:2:4 % 2 compartimenti
    sumTerm = zeros(size(t));
    gradSumTerm_Bi = zeros(size(t));
	gradSumTerm_Li = zeros(size(t));
    Bi = p(i);
    Li = p(i+1);
    delta0 = 1.0/(Li-l(1));
    Ahat = [-A(2)-A(3)-(A(1)/(Li -l(1)))   A(2)  A(3)];
    Abar = [-A(2)-A(3) , A(2) , A(3)];
    
    for j = 1:3 % IF modeled as Feng
        Ahat_j = Ahat(j);
        Abar_j = Abar(j);
        lj = l(j);
        delta = 1./(Li-lj);
        
        sumTerm = sumTerm + Ahat_j*delta*(exp(-lj*t)-exp(-Li*t));
        
        if jacBool
            gradSumTerm_Bi = gradSumTerm_Bi + ...
                                Ahat_j * delta * (exp(-lj*t)-exp(-Li*t));
            gradSumTerm_Li = gradSumTerm_Li + ...
                                Abar_j * (delta^2) .* (exp(-Li*t)-exp(-lj*t)) +...
                                Abar_j * t * delta .* exp(-Li*t);
        end
    end
    sumTerm = sumTerm * Bi + ((A(1)*Bi*t) / (Li-l(1))) .* exp(-l(1)*t);
    irf(idx) = irf(idx) + sumTerm;
    
    if jacBool
        gradSumTerm_Bi = gradSumTerm_Bi + (A(1) * t * delta0) .* exp(-l(1)*t);
        gradSumTerm_Li = gradSumTerm_Li * Bi + ...
                        (A(1) * t * Bi * (delta0^2)) .* (exp(-Li*t)-exp(-l(1)*t)) + ...
                        (2 * A(1) * Bi * (delta0^3)) .* (exp(-Li*t)-exp(-l(1)*t));
        
        J(idx,i) = gradSumTerm_Bi./2;
        J(idx,i+1) = gradSumTerm_Li./2;
    end
    
end

dk = log(2)/109.8;  % radioactive decay constant for F-18
irf(idx) = irf(idx) .* exp(-dk * t);
Ct = (1-p(5))*irf + p(5)* IFfit; % vB = k5

Ct(Ct<=0) = eps;

if jacBool
    J(:,5) = IFfit - irf;
end

end
