function wp = fStokesFit(t, eta, d, order, fit_type, varargin)

% wp = fFitStokes(t, eta, d, order, fit_type, ...)
% ------------------------------------------------------------------------
% Dispersion equation, with variable order (w.r.t. Hk/2) from 1st to 5th.
% - inputs:
%   t [s], eta [m] - time history of the wave surface profile to be fitted,
%                    will automatically pick the largest crest to fit.
%   d - mean water depth [m]
%   order [-] - order of the Stokes theory to be applied.
%   fit_type - 'H' or 'Cr', wave height (H) or crest elevation (Cr)
% 	... - extra option-value pairs ('option', default or further input):
%       ('plotfit', true) - whether to visualise fitted result
%       ('auto' : [T_choice, H_choice]) - automatically pick a H, T choice: 
%           T: 1) 0-up-x, 2) 0-down-x, 3) trough-to-trough
%           H: 1) trough-to-peak, 2) peak-to-trough, 3) average
%           Use this to embed the function in other scripts.
%       ('ReturnFlow', false), ('DTerms', false) and ('SwlAdjust', true):
%           details included at the beginning of fStokesIn.m
% - outputs:
%   wp - wave properties struct required as the input to:
%        fStokesEta(...), fStokesVel(...) and fStokesAcc(...).
% ------------------------------------------------------------------------
% lm808, 10/2014.
% github.com/lm808, all rights reserved.

%% Defaults
ReturnFlow = 0;
DTerms = 0;
SwlAdjust = 1;
plotfit = 1;
auto = 0;

% Additional options
n = length(varargin);
for i = 1:2:n-1
    switch upper(varargin{i})
        case 'PLOTFIT'
            plotfit = fProcessSwitch(varargin{i+1});
        case 'AUTO'
            auto = 1;
            H_choice = varargin{i+1}(2);
            T_choice = varargin{i+1}(1);
        case 'RETURNFLOW'
            ReturnFlow = fProcessSwitch(varargin{i+1});
        case 'DTERMS'
            DTerms = fProcessSwitch(varargin{i+1});
        case 'SWLADJUST'
            SwlAdjust = fProcessSwitch(varargin{i+1});
        otherwise
            error(['Invalid option: ',varargin{i}])
    end
end

%% extract crest region, and shift crest to t=0
sp = fWavePeakTrough(eta,'peak');
sc = find(eta(sp)==max(eta(sp)));

if length(sp)==1
        its = 1;
        ite = length(t);
elseif length(sp)==2
    if sc == 2
        its = sp(sc-1)+1;
        ite = length(t);
    elseif sc == 1
        its = 1;
        ite = sp(sc+1)-1;
    else
        error('')
    end
else
        its = sp(sc-1)+1;
        ite = sp(sc+1)-1;
end

eta = eta(its:ite);
sp = find(eta == max(eta));

t = t(its:ite);
t = t-t(sp);

% find key points
st = fWavePeakTrough(eta,'trough');

ux = fFind0X(eta,'up');
dx = fFind0X(eta,'down');

%% fit period
T = [diff(t(ux)),diff(t(dx)),diff(t(st))];

if auto
    T = T(T_choice);
else
    fprintf('T choices: \n1) 0-up-x: %f\n2) 0-down-x: %f\n3) trough-to-trough: %f\n',T(1),T(2),T(3))
    T = T(input('T-choice: '));
end

%% fit wave height
switch lower(fit_type)
    case 'h'
        H = [eta(sp) - eta(st(1)),eta(sp) - eta(st(2)),mean([eta(sp) - eta(st(1)),eta(sp) - eta(st(2))])];
        if auto
            H = H(H_choice);
        else
            fprintf('H choices: \n1) trough-to-peak: %f\n2) peak-to-trough: %f\n3) average: %f\n',H(1),H(2),H(3))
            H = H(input('H-choice: '));
        end
    case 'cr'
        eta5 = eta(sp)*2;
        H = eta(sp) - eta(st(1));
        while abs(eta(sp)-eta5)>1e-7
            wp = fStokesIn(d,T,H,order,'ReturnFlow',ReturnFlow,'SwlAdjust',SwlAdjust,'DTerms',DTerms);
            eta5 = fStokesEta(0,0,wp);
            H = H+eta(sp)-eta5;
        end
        disp(['H for specified crest elevation: ',num2str(H)])
    otherwise
        error('Invalid fit type, either ''WaveHeight'' or ''CrestHeight''.')
end

%% generate the fitted wave
wp = fStokesIn(d,T,H,order,'ReturnFlow',ReturnFlow,'SwlAdjust',SwlAdjust,'DTerms',DTerms);
eta5 = fStokesEta(0,t,wp);
if strcmpi(fit_type,'H')
    shift = max(eta)-max(eta5);
else
    shift = 0;
end
wp.shift = shift;

%% plot fitted result
if plotfit
    figure
    plot(t,eta,'b')
    hold on
    plot(t(st),eta(st),'ko');
    plot(t(sp),eta(sp),'mo');
    plot(t(dx),eta(dx),'rx');
    plot(t(ux),eta(ux),'gx');
    rule(0,'h','k--');
    xlabel('t')
    ylabel('\eta');
    plot(t,eta5+shift,'r');
end

end

function [out] = fProcessSwitch(in)
    yesList = {'YES','ON','TRUE','Y'};
    noList = {'NO','OFF','FALSE','N'};
    if isnumeric(in) || islogical(in)
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

