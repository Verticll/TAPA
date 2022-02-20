% clear all
clearvars
clearvars -GLOBAL
close all
format shorte

set(0, 'DefaultFigureWindowStyle', 'docked')
global C
global Vx x Fx AtomSpacing
global Phi nAtoms time Mass0 Mass1 Pty0in Pty1in
global LJEpsilon LJSigma Phi0 AtomType
global MinX MaxX PhiTot KETot
global nAtoms0 nAtoms1 T T0 T1 MarkerSize
global doPlotImage PlotCount map im PlotSize ScaleV ScaleF

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

MaxX = 0;
MinX = 0;

nAtoms = 0;
MarkerSize = 12;
Limits = [];
doPlot = 1;
doPlotImage = 0;
PlotCount  = 0;
PlotFile = 'image.gif';
PlotSize = [100, 100, 1049, 895];
ScaleV = 0;
ScaleF = 0;
PlotPosOnly = 0;

% Simulation initiallization
% InitThree
% InitBlock
% InitCirc
% InitBlock0
% InitBlock0FD
% InitVStream
% InitHCP
% InitHCPBlob
% InitVStreamHCP
% InitHCPMeltSim
ODMD

MaxX = max(x) * 1.5;
MinX = min(x) * 1.5;



Fx = zeros(1, nAtoms);
Phi = zeros(1, nAtoms);
dx = zeros(1, nAtoms);
dvx = zeros(1, nAtoms);

Pty0in = AtomType == 0;
Pty1in = AtomType == 1;

nAtoms0 = sum(Pty0in);
nAtoms1 = sum(Pty1in);

GetForcesODMD(PhiCutoff,LJEpsilon,LJSigma);
Phi0 = sum(Phi) / nAtoms / 2;

V2 = Vx.*Vx;
% KEc = 1/2*Mass*mean(V2);
% Tc = KEc/C.kb;

t = 0;
c = 1;
time(c) = 0;

PhiTot(c) = sum(Phi) / 2;

V2_0 = (Vx(Pty0in).*Vx(Pty0in));
if nAtoms1
    V2_1 = (Vx(Pty1in).*Vx(Pty1in));
else
    V2_1 = 0;
end

KE0 = mean(V2_0) * Mass0 * 0.5;
KE1 = mean(V2_1) * Mass1 * 0.5;

KETot(c) = (KE0 * nAtoms0 + KE1 * nAtoms1);
T(c) = KETot(c) / nAtoms / C.kb;
T0(c) = KE0 / C.kb;
T1(c) = KE1 / C.kb;

if PlotPosOnly
    PlotOnlyP(c,Limits);
else
    PlotVarsODMD(c, Limits);
end

xp = x - dt * Vx;
xpp = x - 2 * dt * Vx;

Plt0 = PlDelt;

while t < TStop

    %     F = ma
    %     F = m dv/dt

    GetForcesODMD(PhiCutoff,LJEpsilon,LJSigma);

    % Forward difference
    if Method == 'FD'
        %     dv = F/m dt
        %     x = Vx * dt + F/m (dt)^2 / 2

        dvx(Pty0in) = Fx(Pty0in) * dt / Mass0;
        dvx(Pty1in) = Fx(Pty1in) * dt / Mass1;

        Vx = Vx + dvx;
        dx(Pty0in) = Vx(Pty0in) * dt + Fx(Pty0in) * dt^2 / 2 / Mass0;
        dx(Pty1in) = Vx(Pty1in) * dt + Fx(Pty1in) * dt^2 / 2 / Mass1;

        x = xp + dx;

    elseif Method == 'VE'

        x(Pty0in) = -xpp(Pty0in) + 2 * xp(Pty0in) + dt^2 / Mass0 * Fx(Pty0in);
        x(Pty1in) = -xpp(Pty1in) + 2 * xp(Pty1in) + dt^2 / Mass1 * Fx(Pty1in);


        Vx = (x - xpp) / (2 * dt);%+ randn()*sqrt(1.38064852e-23*500/Mass0)
    end

    xpp = xp;

    xp = x;


    c = c + 1;
    t  = t + dt;
    time(c) = t;


    PhiTot(c) = sum(Phi)/2;
    V2_0 = (Vx(Pty0in).*Vx(Pty0in)); 
    if nAtoms1
        V2_1 = (Vx(Pty1in).*Vx(Pty1in)); 
    else
        V2_1 = 0;
    end

    KE0 = mean(V2_0) * Mass0 * 0.5;
    KE1 = mean(V2_1) * Mass1 * 0.5;

    KETot(c) = (KE0 * nAtoms0 + KE1 * nAtoms1);
    T(c) = KETot(c) / nAtoms / C.kb;
    T0(c) = KE0 / C.kb;
    T1(c) = KE1 / C.kb;

    if t > Plt0
        fprintf('time: %g (%5.2g %%)\n', t, t / TStop * 100);

        if PlotPosOnly
            PlotOnlyP(c,Limits);
        else
            PlotVarsODMD(c, Limits);
        end
        

        Plt0 = Plt0 + PlDelt;
        pause(0.00001)
    end

end


if doPlotImage
    imwrite(im, map, PlotFile, 'DelayTime', 0.05, 'LoopCount', inf);
end


