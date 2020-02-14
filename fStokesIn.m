function [wp] = fStokesIn(d, T, H, order, varargin)

% wp = fStokesIn(d, T, H, order, ...)
% ------------------------------------------------------------------------
% Constructs the input wave properties struct for a stokes solution.
% - inputs:
%   d [m] - mean water depth.
%   T [s] - wave period.
%   H [m] - wave height, or the crest elevation. If it is the latter, set
%           'IterateCrest' to true.
%   order [-] - order of the Stokes theory to be applied.
% 	... - extra option-value pairs ('option', default or further input):
%       ('IterateCrest', false) - switch this on if H is specified as a
%                                 target crest elevation.
%       ('ReturnFlow', false) - Eulerian return flow. This affects the wave
%                               number and the horizontal velocity.
%       ('DTerms', false) - D-coefficients in Table 2 of Fenton (195) that
%                           describe the effect of a wave propagating on
%                           top of Stokes drift. These will affect the wave
%                           number calculated.
%       ('SwlAdjust', true) - if order>=3, shift the position of the still
%                             water level. As proposed by "Bowden, K.F. 
%                             (1948) Observations of wave in a tidal
%                             current. Proc. Roy. Soc. A 192:403-25".
% - outputs:
%   wp - wave properties struct required as the input to:
%        fStokesEta(...), fStokesVel(...) and fStokesAcc(...).
% ------------------------------------------------------------------------
% lm808, 10/2019.
% github.com/lm808, all rights reserved.


    % default
    IterateCrest = false;

    % construct  struct
    wp = struct('waveModel','stokes','H',H,'T',T,'d',d,'order',order,'omega',2*pi/T);

    % process input for switches
    if nargin >= 5
        n = length(varargin);
        for i = 1:2:n-1
            switch upper(varargin{i})
                case 'RETURNFLOW'
                    wp.ReturnFlow = fProcessSwitch(varargin{i+1});
                case 'DTERMS'
                    wp.DTerms = fProcessSwitch(varargin{i+1});
                case 'SWLADJUST'
                    wp.SwlAdjust = fProcessSwitch(varargin{i+1});
                case 'ITERATECREST'
                    IterateCrest = fProcessSwitch(varargin{i+1});
                otherwise
                    error('Invalid option.')
            end
        end
    end

    if ~isfield(wp,'ReturnFlow')
        wp.ReturnFlow = 0;
        disp('fStokesIn: default OFF for return flow.')
    end

    if ~isfield(wp,'DTerms')
        wp.DTerms = 0;
        disp('fStokesIn: default OFF for d-terms.')
    end

    if ~isfield(wp,'SwlAdjust')
        wp.SwlAdjust = 1;
        disp('fStokesIn: default ON for SWL adjustment.')
    end

    % compute wave number
    wp.k = fDispersionV5(d,T,H,order,'ReturnFlow',wp.ReturnFlow,'DTerms',wp.DTerms);
    wp.lamda = 2*pi/wp.k;
    
    % if the input H is in fact a target crest for iteration
    if IterateCrest
        c_target = H;
        tol = 1e-12;
        err = 1;
        wp = fStokesIn(d, T, c_target* 2, order, ...
                       'ReturnFlow',wp.ReturnFlow,...
                       'DTerms',wp.DTerms,...
                       'SwlAdjust', wp.SwlAdjust);
        c = fStokesEta(0, 0, wp);
        while err > tol
            wp = fStokesIn(d, T, c_target / c * wp.H, order, ...
                           'ReturnFlow',wp.ReturnFlow,...
                           'DTerms',wp.DTerms,...
                           'SwlAdjust', wp.SwlAdjust);
            c = fStokesEta(0, 0, wp);
            err = abs(c-c_target);
        end
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
