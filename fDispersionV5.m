function [k] = fDispersionV5(d,T,H,order,varargin)
% Dispersion equation, with variable order (w.r.t. Hk/2) from 1st to 5th
% Fenton, 1985
% ------------------------------------------------------------------------
% d - mean water depth
% T - wave period
% H - trough-to-crest wave height
% order - order of the stokes theory to be applied
% varargin:
%   'ReturnFlow' : ['on' | 'off']
%   'DTerms' : ['on' | 'off']
% ------------------------------------------------------------------------
% Li Ma, April 2015

% defaults & process inputs
% these are here for compatibitily reasons with old codes. could be stripped-out.

    if nargin >= 5
        n = length(varargin);
        for i = 1:2:n-1
            switch upper(varargin{i})
                case 'RETURNFLOW'
                    ReturnFlow = fProcessSwitch(varargin{i+1});
                case 'DTERMS'
                    DTerms = fProcessSwitch(varargin{i+1});
                otherwise
                    error('Invalid option.')
            end
        end
    else
        switch nargin
            case {0,1}
                error('Insufficient input arguments.')
            case 2
                order = 1;
                H = NaN; % does not need H for 1st order calculations
                disp('fDispersionV5: default order = 1 when H not given.')
            case 3
                order = 5;
                disp('fDispersionV5: default order = 5.')
        end
    end
    
    if ~exist('ReturnFlow','var')
        ReturnFlow = 0;
        disp('fDispersionV5: default OFF for return flow.')
    end
    if ~exist('DTerms','var')
        DTerms = 0;
        disp('fDispersionV5: default OFF for d-terms.')
    end

% ----------------------- start calcs -----------------------------------
    omega = 2*pi/T;
    g = 9.81;
    
% initial guess of k
    p = omega^2*d/g;
    q = (tanh(p^0.75))^(-2/3);
    k0 = omega^2*q/g;
    k = fzero(@(k) calcLHS(k,d,omega,H,order,ReturnFlow,DTerms),k0);
end

function [lhs] = calcLHS(k,d,omega,H,order,ReturnFlow,DTerms)
    g = 9.81;
    kd = k*d;
    epsilon = H/2*k;
    S = 1/cosh(2*kd);
    
    C0 = sqrt(tanh(kd));
    C2 = ((sqrt(tanh(kd)))*(2+7*S^2))/(4*(1-S)^2);
    C4 = ((sqrt(tanh(kd)))*(4+32*S - 116*S^2 - 400*S^3- 71*S^4 +146*S^5))/(32*(1-S)^5);
    
    % which dispersion equation to use
    if DTerms % Fenton (1985), calculation of wave number, case 4
        D2 = -0.5*sqrt(1/tanh(kd));
        D4 = (sqrt(1/tanh(kd))*(2 + 4*S + S^2 + 2*S^3))/(8*(1-S)^3);        
    else % Fenton (1985), calculation of wave number, case 3
        D2 = 0;
        D4 = 0;
    end
    
    % incoprating the steepening effect of return flow
    if ReturnFlow
        ur = -((H/2)^2)*k*pi*(1/tanh(kd))/((2*pi/omega)*kd);  
    else
        ur = 0;
    end
    
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

function [out] = fProcessSwitch(in)

    yesList = {'YES','ON','TRUE','Y'};
    noList = {'NO','OFF','FALSE','N'};
    
    if isnumeric(in) | islogical(in)
        out = ~~in;
    else
        if any(strcmpi(in,yesList)) 
            out = true;
        elseif any(strcmpi(in,noList))
            out = false;
        else
            error('Invalid switch value.')
        end
    end
end
