
clear
close all

% -------------------------------------------------------------------------
% Averaging AAQ data
% -------------------------------------------------------------------------

% --- Import functions

% --- Configulations

aaq_file = 'test.csv';
max_dep     = 50;    % Maximum depth (m)
dz          = 0.02;   % Depth interval for data binning (m)
grav = 9.81;

% --- Read AAQ data

data = readmatrix(aaq_file);
dep    = data(:,1);   % Water depth (m)
temp   = data(:,2);   % Temperature (degree)
sal    = data(:,3);   % Salinity (psu)
rho    = data(:,6);   % Water density (kg/m3)
sigmat = data(:,7);   % SigmaT
chl    = data(:,9);   % Chl-a (ppb)
tur    = data(:,10);  % Turbidity (FTU)
oxy    = data(:,13);  % DO (mg/L)
par    = data(:,14);  % PAR (umol photon/m2/s)

% --- Binning data

[temp_bin,   dep_bin] = depth_bin(temp, dep, max_dep, dz);
[sal_bin,    dep_bin] = depth_bin(sal, dep, max_dep, dz);
[chl_bin,    dep_bin] = depth_bin(chl, dep, max_dep, dz);
[tur_bin,    dep_bin] = depth_bin(tur, dep, max_dep, dz);
[oxy_bin,    dep_bin] = depth_bin(oxy, dep, max_dep, dz);
[par_bin,    dep_bin] = depth_bin(par, dep, max_dep, dz);
[rho_bin,    dep_bin] = depth_bin(rho, dep, max_dep, dz);
[sigmat_bin, dep_bin] = depth_bin(sigmat, dep, max_dep, dz);

% --- Find surface and bottom layers

for i = 1 : size(dep_bin,1)
    if ~isnan(temp_bin(i))
        sur_layer = i;
        break
    end
end
for i = 1 : size(dep_bin,1)
    z = dep_bin(i);
    if z > 1.0 && isnan(temp_bin(i))
        bot_layer = i-1;
        break
    end
end

% --- Calculate potential energy anomaly (fai: J/m3)

rho_bar = mean(rho_bin(sur_layer:bot_layer));
H = dep_bin(bot_layer)-dep_bin(sur_layer);
fai = 0;
for i = sur_layer : bot_layer
    fai = fai + (rho_bar - rho_bin(i))*grav*(-1)*dep_bin(i)*dz;
end
fai = fai / H;

rho_bar_simple = (rho_bin(sur_layer) + rho_bin(bot_layer))/2;
fai_simple = (rho_bar_simple-rho_bin(sur_layer))*grav*(-1)*dep_bin(sur_layer)*H*0.5 ...
             + (rho_bar_simple-rho_bin(bot_layer))*grav*(-1)*dep_bin(bot_layer)*H*0.5;
fai_simple = fai_simple / H;

% --- Output

output = [H, oxy_bin(sur_layer), oxy_bin(bot_layer), temp_bin(sur_layer), temp_bin(bot_layer),...
    sal_bin(sur_layer), sal_bin(bot_layer), rho_bin(sur_layer), rho_bin(bot_layer), fai, fai_simple];

% --- Plotting

subplot(2,3,1)
plot(temp_bin(1:bot_layer), dep_bin(1:bot_layer))
ax = gca;
ax.YDir = 'reverse';
title('Temp')

subplot(2,3,2)
plot(sal_bin(1:bot_layer), dep_bin(1:bot_layer))
ax = gca;
ax.YDir = 'reverse';
title('Sal')

subplot(2,3,3)
plot(sigmat_bin(1:bot_layer), dep_bin(1:bot_layer))
ax = gca;
ax.YDir = 'reverse';
title('sigmaT')

subplot(2,3,4)
plot(chl_bin(1:bot_layer), dep_bin(1:bot_layer))
ax = gca;
ax.YDir = 'reverse';
title('Chl-a')

subplot(2,3,5)
plot(tur_bin(1:bot_layer), dep_bin(1:bot_layer))
ax = gca;
ax.YDir = 'reverse';
title('Turbidity')

subplot(2,3,6)
plot(oxy_bin(1:bot_layer), dep_bin(1:bot_layer))
ax = gca;
ax.YDir = 'reverse';
title('DO')




















