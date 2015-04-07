function [eta,etaij,swl] = fStokesEta(x,t,wp)
% Computes surface elevation according to stokes theory (Fenton, 1985)
% Acceptable input space-time coordinate (x,z,t) formats:
% 1) Only one is a non-scalar
% 2) Arrays of same sizes
% 'order'  can be 1~5;
% Li Ma, August 2013

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

B = zeros(5,5);
B(1,1) = 1;
B(2,2) = (1/tanh(kd))*(1+2*S)/(2*(1-S));
B(3,1) = -3*(1+3*S+3*(S^2)+2*(S^3))/(8*((1-S)^3));
B(3,3) = -B(3,1);
B(4,2) = (1/tanh(kd))*(6-26*S-182*(S^2)-204*(S^3)-25*(S^4)+26*(S^5))/(6*(3+2*S)*((1-S)^4));
B(4,4) = (1/tanh(kd))*(24+92*S+122*(S^2)+66*(S^3)+67*(S^4)+34*(S^5))/(24*(3+2*S)*((1-S)^4));
B(5,3) = 9*(132+17*S-2216*(S^2)-5897*(S^3)-6292*(S^4)-2687*(S^5)+194*(S^6)+467*(S^7)+82*(S^8))/(128*(3+2*S)*(4+S)*((1-S)^6));
B(5,5) = 5*(300+1579*S+3176*(S^2)+2949*(S^3)+1188*(S^4)+675*(S^5)+1326*(S^6)+827*(S^7)+130*(S^8))/(384*(3+2*S)*(4+S)*((1-S)^6));
B(5,1) = -(B(5,3)+B(5,5));

B = B/k;

%% Components from different orders/harmonics
etaij = cell(order,order);
for i = 1:order
    for j = 1:i
        etaij{i,j} = (epsilon^i)*B(i,j)*cos(j*(omega*t-k*x));
    end
end

%% Sum
eta = zeros(size(etaij{1,1}));
for i = 1:order
    for j = 1:i
        eta = eta + etaij{i,j};
    end
end

%% Mean water level change (Bowden)
if wp.SwlAdjust && wp.order>1
    swl = ((H/2)^2)*k/(2*sinh(2*kd));
    eta = eta + swl;
end

