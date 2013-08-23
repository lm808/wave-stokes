function [k] = fDispersionV5(d,T,H,order)
% Dispersion equation, with variable order (w.r.t. Hk/2) from 1st to 5th
% (return flow is incoprated, current effects still to be included)
% -----------------------
% d - mean water depth
% T - wave period
% H - trough-to-crest wave height
% order - order of the stokes theory to be applied
% -----------------------
% Li Ma, July 2013

    switch nargin
        case {0,1}
            error('Insufficient input arguments.')
        case 2
            order = 1;
            H = NaN; % does not need H for 1st order calculations
        case 3
            warning('Calling 5th order dispersion equation by default.')
            order = 5;
        case 4
        otherwise
            error('Too many input arguments.')
    end

    omega = 2*pi/T;
    k = fzero(@(k) calcLHS(k,d,omega,H,order),omega^2/9.81);
end

function [lhs] = calcLHS(k,d,omega,H,order)
    g = 9.81;
    kd = k*d;
    epsilon = H/2*k;
    S = 1/cosh(2*kd);
    
    ur = -((H/2)^2)*k*pi*(1/tanh(kd))/((2*pi/omega)*kd);
    
    C0 = sqrt(tanh(kd));
    C2 = ((sqrt(tanh(kd)))*(2+7*S^2))/(4*(1-S)^2);
    C4 = ((sqrt(tanh(kd)))*(4+32*S - 116*S^2 - 400*S^3- 71*S^4 +146*S^5))/(32*(1-S)^5);
    
    D2 = -0.5*sqrt(1/tanh(kd));
    D4 = (sqrt(1/tanh(kd))*(2 + 4*S + S^2 + 2*S^3))/(8*(1-S)^3);
    
    switch order
        case 1
            lhs = (omega/k) * sqrt(k/g) - C0;
        case 2
            lhs = (omega/k-ur) * sqrt(k/g) - C0;
        case {3,4}
            lhs = (omega/k-ur) * sqrt(k/g) - C0 - (epsilon^2)*(C2+D2/kd);
        case 5
            lhs = (omega/k-ur) * sqrt(k/g) - C0 - (epsilon^2)*(C2+D2/kd) - (epsilon^4)*(C4+D4/kd);
    end
end

