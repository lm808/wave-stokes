function [wp] = fStokesIn(H,T,d,order)

wp.H = H;
wp.T = T;
wp.d = d;
wp.order = order;

wp.omega = 2*pi/T;
wp.k = fDispersionV5(d,T,H,order);

wp.lamda = 2*pi/wp.k;
wp.waveModel = 'stokes';
% wp.waveModel = ['stokes-',num2str(order)];

