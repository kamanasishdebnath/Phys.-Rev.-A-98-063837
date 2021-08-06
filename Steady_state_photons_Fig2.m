% The following numerical code can be used to validate
% Fig. 2 reported in Phys. Rev. A 98, 063837 (2018).


rng shuffle
clear
clc

% Change the following parameters as given in Fig.2 (a) and (b)
kappa= 2*pi*160*10^3;    % Cavity linewidth
gamma= 2*pi*0.001;       % Atomic linewidth
g= 2.41;

% Defining X and Y axis of Fig.2  
Ni= logspace(2,7,300);          % number of atoms (Y axis)
del= 2*pi*logspace(-4,3,301);   % Pump rate (X axis)


% This for loop calculates the photon number for different
% values of N and pumping rates
for loop1= 1:length(Ni)
    N1=Ni(loop1);
    disp(loop1)
    for loop2=1:length(del)
        delta= del(loop2);
        photon1(loop1,loop2)= steady_state_photons(N1, gamma, delta, kappa, g); 
    end
end


% Plotting the results
figure(1)
pcolor(del , Ni, log10(abs(photon1)-1)), shading interp
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
colorbar()
colormap('hot')
caxis([-5.5, 2])
xlabel({'Pumping rate (Hz)'})
ylabel({'Number of atoms, N'})
title('Steady state photons')





% This function computes the steady state photon using a
% formula which was computed analytically using the supplied
% Mathematica file.
function callme1= steady_state_photons(N1, gamma, delta, kappa, g)

callme1 = (N1 *(-0.0625* gamma^4 + gamma^2 *delta *(-0.375 *delta - 0.1875 *kappa) + gamma *delta^2 *(-0.25 *delta - 0.1875 *kappa) + ...
    gamma^3* (-0.25 *delta - 0.0625* kappa) + delta^3 *(-0.0625* delta - 0.0625* kappa))*kappa + g^2 *(N1^2 *(-0.25 *gamma^3 - ...
    0.25 *gamma^2 *delta + 0.25 *gamma* delta^2 + 0.25 *delta^3) + N1* ((1. - 0.25* N1)* gamma^2 +  1.5* gamma *delta + ...
    (0.5 + 0.25 *N1)* delta^2)* kappa + (-1. + 1.* N1) *(gamma + delta)*kappa^2) + 0.0625* N1* gamma* sqrt(-4.* g^2 *kappa *(-1. *gamma^3 + ...
    g^2 *(-4. *gamma - 4. *delta) + gamma *delta^2 + gamma^2 *(-1. *delta - 1.* kappa) + delta^2 *(delta + kappa))* (-4.* gamma* kappa - ...
    4. *delta *kappa + N1* (4.* gamma^2 + 8. *gamma *delta + 4.* delta^2 + 4. *gamma* kappa + 4.* delta* kappa)) + ((-1.* gamma^3 + ...
    gamma *delta *(-3. *delta - 2. *kappa) + gamma^2 *(-3. *delta - 1. *kappa) + delta^2 *(-1. *delta - 1. *kappa))* kappa + ...
    g^2 *(-8. *gamma *kappa + N1 *(4. *gamma^2 - 4.* delta^2 + 4. *gamma* kappa - 4. *delta *kappa)))^2) + ...
    0.0625 *N1 *delta* sqrt (-4.* g^2 *kappa *(-1.* gamma^3 + g^2 *(-4. *gamma - 4. *delta) + gamma* delta^2 + ...
    gamma^2 *(-1. *delta - 1. *kappa) + delta^2 *(delta + kappa))* (-4.* gamma* kappa -  4. *delta *kappa + N1* (4.* gamma^2 + 8. *gamma *delta + 4. *delta^2 + ...
    4. *gamma *kappa + 4.* delta* kappa)) + ((-1.* gamma^3 + gamma *delta *(-3. *delta - 2. *kappa) + gamma^2 *(-3. *delta - 1.* kappa) + ...
    delta^2 *(-1.* delta - 1.* kappa))* kappa + g^2 *(-8.* gamma* kappa + N1* (4. *gamma^2 - 4.* delta^2 + 4.* gamma* kappa - ...
    4. *delta *kappa)))^2))/(g^2 *kappa *(-1.* gamma* kappa - 1.* delta* kappa + N1 *(1. *gamma^2 + 2.* gamma* delta + 1.* delta^2 + ...
    1. *gamma* kappa + 1.* delta* kappa)));
end




