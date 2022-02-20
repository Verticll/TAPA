Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201;
l = 1e-6;

x =linspace(0,l,nx);
dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

Nd = 1e16 * linspace(1,20,nx); % Const. 1/cm3 (100 cm/m)^3
NetDoping = ones(1,nx).*Nd; % doping

x0 = l/2;
nw = l/20;
npDisturbance = 1e16*1e6*exp(-((x-x0)/nw).^2);

LVbc = 0;
RVbc = 0;

TStop = 14200000*1e-18;
PlDelt = 100000*1e-18;

PlotYAxis = {[-1e-6 2e-6] [-10 10] [-1.5e-2 1.5e-2]...
    [0 2e22] [0 1e22] [0 10e43]...
    [-10e33 10e33] [-1e34 1e34] [-1e8 1e8] ...
    [-1e8 1e8] [-10e-3 10e-3] [0 2e22]};

doPlotImage = 0;
PlotFile = 'Gau2CarRC.gif';
