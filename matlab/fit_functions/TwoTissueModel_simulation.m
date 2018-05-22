%|=================================================================================
%|   TWO-TISSUE COMPARTMENT MODEL SIMULATION
%|
%|   INPUTS:
%|       p:         model parameters in auxiliary domain 
%|                  (alpha1, beta1, aòpha2, beta2, fv) [5 x 1]
%|       scanTime:  vector of start and end times of each frame [time_points x 2]
%|       IFparams:  input function model parameters 
%|                  (delay, A1, A2, A3, lambda1, lambda2, lambda3) [1 x 7]
%|       IFfit:     fitted arterial inputn function [time_points x 1]
%|
%|   OUTPUTS:
%|       Ct:        simulated time curve [time_points x 1]
%|       J:         jacobian [time_points x 6]
%|
%|  Last revision:
%|  22 May 2018
%|  Michele Scipioni, Univeristy of Pisa
%|
%|=================================================================================


function [Ct, J] = TwoTissueModel_simulation(p, scanTime, IFparams, IFfit)

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