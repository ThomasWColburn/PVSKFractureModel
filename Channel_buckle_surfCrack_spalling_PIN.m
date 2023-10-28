%--------------%Channel Cracking and Buckling on Soda Lime PIN Devices%--------------%
clear all;

%Thomas Colburn%
%Set the constants%
%Fracture energy release rates for the film layers, in units of J/m^2
GcITO = 5; 
GcNiO = 5; 
GcAgEvap = 5; 
GcAgScreen = 2;
GcPVSKRSPP = 4.4;
GcPVSKhotcast = 0.92;
GcPVSKspincoat = 0.37;
GcPCBM = 0.13;

sigmaITO = -2.23e8; % Thermal stress for ITO on soda lime glass substrate
sigmaNiO = 1.64e8; % Thermal stress for NiOx on soda lime glass substrate
sigmaPVSKRSPP = 1.43e8; sigmaPVSKspincoat = 9.75e7; sigmaPVSKhotcast = 1.10e8;
sigmaPCBM = 2.16e7;

vITO = 0.35; vNiO = 0.416; vAg = 0.33; %Poisson's Ratio

Zchannel = 1.976; % Nondimensional cracking number for channel cracking
Zsurface = 3.951; % Nondimensional cracking number for surface cracking
Zspalling = 0.343; % Nondimensional cracking number for spalling

EprimeITO = 1.75e11; EprimeNiO = 1.64e11; EprimeAgEvap = 1.19e11; EprimeAgScreen = 1.34e10; %Young's modulus divided by 1-v Pa
EprimePVSK = 1.61e10; EprimePCBM = 3.90e9;

EITO = 1.14e11; ENiO = 9.58e10; EAgEvap = 8e10; EAgScreen = 9e9; %Young's modulus


% For steady-state buckling, the energy release rate is:
% G_buckling = 0.5*(sigma^2*h)*(1-v^2)/E
% Calculate hc by setting G_buckling = Gc

% G_channel = Zchannel*(sigma^2*h)/E' , with Zchannel = 1.976
% G_surface = Zsurface*(sigma^2*h)/E' , with Zsurface = 3.951
% Here, E' = E/(1-v) , defined as the biaxial modulus

disp("Channel cracking and buckling hc for sodalimeglass PIN devices:")
%Calculate ITO Buckling hc%
hcITOglass_buckle = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)

%Calculate NiOx channel vs surface crack
hcNiOx_channelglass = GcNiO*EprimeNiO/(Zchannel*sigmaNiO^2)
hcNiOx_surfaceglass = GcNiO*EprimeNiO/(Zsurface*sigmaNiO^2)
hcNiOx_spallingglass = GcNiO*EprimeNiO/(Zspalling*sigmaNiO^2)

% Calculate PVSK and PCBM channel and surface crack on soda lime glass
hcPVSKRSPP_channelglass = GcPVSKRSPP*EprimePVSK/(Zchannel*sigmaPVSKRSPP^2)
hcPVSKRSPP_surfaceglass = GcPVSKRSPP*EprimePVSK/(Zsurface*sigmaPVSKRSPP^2)
hcPVSKRSPP_spallingglass = GcPVSKRSPP*EprimePVSK/(Zspalling*sigmaPVSKRSPP^2)

hcPVSKspincoat_channelglass = GcPVSKspincoat*EprimePVSK/(Zchannel*sigmaPVSKspincoat^2)
hcPVSKspincoat_surfaceglass = GcPVSKspincoat*EprimePVSK/(Zsurface*sigmaPVSKspincoat^2)
hcPVSKspincoat_spallingglass = GcPVSKspincoat*EprimePVSK/(Zspalling*sigmaPVSKspincoat^2)

hcPVSKhotcast_channelglass = GcPVSKhotcast*EprimePVSK/(Zchannel*sigmaPVSKhotcast^2)
hcPVSKhotcast_surfaceglass = GcPVSKhotcast*EprimePVSK/(Zsurface*sigmaPVSKhotcast^2)
hcPVSKhotcast_spallingglass = GcPVSKhotcast*EprimePVSK/(Zspalling*sigmaPVSKhotcast^2)

hcPCBM_channelglass = GcPCBM*EprimePCBM/(Zchannel*sigmaPCBM^2)
hcPCBM_surfaceglass = GcPCBM*EprimePCBM/(Zsurface*sigmaPCBM^2)
hcPCBM_spallingglass = GcPCBM*EprimePCBM/(Zspalling*sigmaPCBM^2)

disp("Channel cracking and buckling hc for polyimide PIN devices:")
%--------------%Buckling on Polyimide for PIN Devices%--------------%
sigmaITO = -9.93e8; sigmaNiO = -3.32E+08; sigmaAgEvap = -5.97E+06; sigmaAgScreen = -1.68e6;
sigmaPVSKRSPP = 1.24e8; sigmaPVSKspincoat = 8.42e7; sigmaPVSKhotcast = 9.54e7;
sigmaPCBM = 1.95e7;

hcITObucklePI = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)
hcNiObucklePI = 2*GcNiO*ENiO/((1-vNiO^2)*sigmaNiO^2)

% Calculate PVSK and PCBM channel and surface crack on soda lime glass
hcPVSKRSPP_channelglass = GcPVSKRSPP*EprimePVSK/(Zchannel*sigmaPVSKRSPP^2)
hcPVSKRSPP_surfaceglass = GcPVSKRSPP*EprimePVSK/(Zsurface*sigmaPVSKRSPP^2)
hcPVSKRSPP_spallingglass = GcPVSKRSPP*EprimePVSK/(Zspalling*sigmaPVSKRSPP^2)

hcPVSKspincoat_channelglass = GcPVSKspincoat*EprimePVSK/(Zchannel*sigmaPVSKspincoat^2)
hcPVSKspincoat_surfaceglass = GcPVSKspincoat*EprimePVSK/(Zsurface*sigmaPVSKspincoat^2)
hcPVSKspincoat_spallingglass = GcPVSKspincoat*EprimePVSK/(Zspalling*sigmaPVSKspincoat^2)

hcPVSKhotcast_channelglass = GcPVSKhotcast*EprimePVSK/(Zchannel*sigmaPVSKhotcast^2)
hcPVSKhotcast_surfaceglass = GcPVSKhotcast*EprimePVSK/(Zsurface*sigmaPVSKhotcast^2)
hcPVSKhotcast_spallingglass = GcPVSKhotcast*EprimePVSK/(Zspalling*sigmaPVSKhotcast^2)

hcPCBM_channelglass = GcPCBM*EprimePCBM/(Zchannel*sigmaPCBM^2)
hcPCBM_surfaceglass = GcPCBM*EprimePCBM/(Zsurface*sigmaPCBM^2)
hcPCBM_spallingglass = GcPCBM*EprimePCBM/(Zspalling*sigmaPCBM^2)

hcAgEvapbucklePI = 2*GcAgEvap*EAgEvap/((1-vAg^2)*sigmaAgEvap^2)
hcAgScreenbucklePI = 2*GcAgScreen*EAgScreen/((1-vAg^2)*sigmaAgScreen^2)

