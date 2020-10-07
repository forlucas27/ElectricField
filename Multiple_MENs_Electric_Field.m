%% Approximation of Electric Field at some point from a single MEN
% Rachel Lucas
clear; clc; close all;

% CONSTANTS
ke = 8.987551787e9; % Coulomb's constant in N*m^2*C^-2
el = 1.60217662e-19; % Charge of an electron in Coulombs


% Assumptions
Q = -1.5e-19; % Obtained from Targeted and controlled anticancer drug delivery and release with magnetoelectric nanoparticles
alpha = 100*0.1; % ME coefficient (could be in range from 10-100) in mV/(cm*Oe) *0.1 to convert to V/(m*Oe)

% Variables 
H = [353.6 -353.6]; % Magnetic field in Oe
a = [0 0.5e-9 0.866e-9 1e-9 0.866e-9 0.5e-9 0 -0.5e-9 -0.866e-9 -1e-9 -0.866e-9 -0.5e-9]; % y position of MENs
b = [1e-9 0.866e-9 0.5e-9 0 -0.5e-9 -0.866e-9 -1e-9 -0.866e-9 -0.5e-9 0 0.5e-9 0.866e-9]; % x position of MENs
n = 12; % number of charges

% Calculated values
pol = alpha*H; % polarization of a MEN- given by papers
p = (pol*(3e-8)^3)/ke; % dipole moment of a single MEN calculated from electric field of a dipole assuming distance btwn polarization is 30nm

Vtotal = zeros(51);

[X,Y] = meshgrid(-1.5e-9:0.06e-9:1.5e-9,-1.5e-9:0.06e-9:1.5e-9);
for i=1:1:n
    xhat = X./sqrt((X-b(i)).^2+(Y-a(i)).^2);
    yhat = (Y-a(i))./sqrt((X-b(i)).^2+(Y-a(i)).^2);

    Vpoint = (ke *Q)./sqrt((X-b(i)).^2+(Y-a(i)).^2);

    Vdipole = (ke*(p(1)*xhat + p(2)*yhat))./((X-b(i)).^2+(Y-a(i)).^2);
    
    disp(size(Vtotal));
    disp(size(Vpoint));
    disp(size(Vdipole));
    Vtotal = Vtotal + Vpoint + Vdipole;

end

figure;
contour(X,Y,Vtotal,5000);
[S,T]=gradient(Vtotal);
S(isinf(S)) = nan;
T(isinf(T)) = nan;
hold on; quiver(X,Y,S,T);hold off;






