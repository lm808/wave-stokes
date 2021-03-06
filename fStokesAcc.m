function [dudt, dwdt] = fStokesAcc(x, z, t, wp)

% [u, w, ur, uij] = fStokesVel(x,z,t,wp)
% ------------------------------------------------------------------------
% Calculates the water particle accelerations based on Fenton (1985).
% inputs:
%   x, z [m], t [s] - spatial and time points.
%                     Acceptable input formats for x, z and t:
%                        1) Only one is a non-scalar.
%                        2) Arrays of the same size.
%   wp - struct containing the input wave properties. See fStokesIn.m
% outputs:
%   dudt, dwdt [m/s^2] - horizontal and vertical particle accelerations.
% ------------------------------------------------------------------------
% lm808, 08/2013.
% github.com/lm808, all rights reserved.

g = 9.81;

%% Unpack
H = wp.H;
omega = wp.omega;
d = wp.d;
k = wp.k;
order = wp.order;

%% Fenton coefficients
kd = k*d;
S = 1/cosh(2*kd);
epsilon = H/2*k;

C0 = sqrt(tanh(kd));
C = C0*sqrt(g/(k^3));

A = zeros(5,5);
A(1,1) = 1/sinh(kd);
A(2,2) = 3*(S^2)/(2*((1-S)^2));
A(3,1) = (-4-20*S+10*(S^2)-13*(S^3))/(8*sinh(kd)*((1-S)^3));
A(3,3) = (-2*(S^2)+11*(S^3))/(8*sinh(kd)*((1-S)^3));
A(4,2) = (12*S-14*(S^2)-264*(S^3)-45*(S^4)-13*(S^5))/(24*((1-S)^5));
A(4,4) = (10*(S^3)-174*(S^4)+291*(S^5)+278*(S^6))/(48*(3+2*S)*((1-S)^5));
A(5,1) = (-1184+32*S+13232*(S^2)+21712*(S^3)+20940*(S^4)+12554*(S^5)-500*(S^6)-3341*(S^7)-670*(S^8))/(64*sinh(kd)*(3+2*S)*(4+S)*((1-S)^6));
A(5,3) = (4*S+105*(S^2)+198*(S^3)-1376*(S^4)-1302*(S^5)-117*(S^6)+58*(S^7))/(32*sinh(kd)*(3+2*S)*((1-S)^6));
A(5,5) = (-6*(S^3)+272*(S^4)-1552*(S^5)+852*(S^6)+2029*(S^7)+430*(S^8))/(64*sinh(kd)*(3+2*S)*(4+S)*((1-S)^6));

%% adjust SWL (Bowden)
if wp.SwlAdjust && wp.order>1
    swl = ((H/2)^2)*k/(2*sinh(2*kd));
    z = z - swl;
end

%% Components from different orders/harmonics
dudtij = cell(order,order);
dwdtij = cell(order,order);
for i = 1:order
    for j = 1:i
        dudtij{i,j}= - C * (epsilon^i) * (j*j*k*omega) * A(i,j) * cosh(j*k*(z+d)) .* sin(j*(omega*t-k*x));
        dwdtij{i,j} = - C * (epsilon^i) * (j*j*k*omega) * A(i,j) * sinh(j*k*(z+d)) .* cos(j*(omega*t-k*x));
    end
end

%% Sum
dudt = zeros(size(dudtij{1,1}));
dwdt = zeros(size(dwdtij{1,1}));
for i = 1:order
    for j = 1:i
        dudt = dudt + dudtij{i,j};
        dwdt = dwdt + dwdtij{i,j};
    end
end
