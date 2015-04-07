function [wp] = fStokesIn(d,T,H,order,varargin)
% Constructs the input wave properties struct for stokes solution
% ------------------------------------------------------------------------
% d - mean water depth
% T - wave period
% H - trough-to-crest wave height
% order - order of the stokes theory to be applied
% varargin:
%   'ReturnFlow' : ['on' | 'off']
%   'DTerms' : ['on' | 'off']
%   'SwlAdjust' : ['on' | 'off']
% ------------------------------------------------------------------------
% Li Ma, April 2015

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
        wp.SwlAdjust = 0;
        disp('fStokesIn: default OFF for SWL adjustment.')
    end

    % compute wave number
    wp.k = fDispersionV5(d,T,H,order,'ReturnFlow',wp.ReturnFlow,'DTerms',wp.DTerms);
    wp.lamda = 2*pi/wp.k;
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
