%--------------%Channel Cracking and Buckling on Soda Lime NIP Devices%--------------%
clear all;
%Thomas Colburn%
%Set the constants%
%Fracture energy release rates, in units of J/m^2
GcITO = 5; 
GcSnOx = 0.7; 
GcAuEvap = 1.2; 
GcAuScreen = 2;
GcPVSKRSPP = 4.4;
GcPVSKhotcast = 0.92;
GcPVSKspincoat = 0.37;
GcSpiro = 0.45;

sigmaITO = -2.23e8; % Thermal stress for ITO on soda lime glass substrate
sigmaSnOx = -1.07e8; % Thermal stress for SnOx on soda lime glass substrate
sigmaPVSKRSPP = 1.43e8; sigmaPVSKspincoat = 9.75e7; sigmaPVSKhotcast = 1.10e8;
sigmaSpiro = 4.80e7;

vITO = 0.35; vSnOx = 0.33; vAu_evap = 0.421; vAu_screen = 0.33; %Poisson's Ratio

Zchannel = 1.976; % Nondimensional cracking number for channel cracking
Zsurface = 3.951; % Nondimensional cracking number for surface cracking
Zspalling = 0.343; % Nondimensional cracking number for spalling

EprimeITO = 1.75e11; EprimeSnOx = 1.31e11; EprimeAuEvap = 1.42e11; EprimeAuScreen = 1.34e10; %Young's modulus divided by 1-v Pa
EprimePVSK = 1.61e10; EprimeSpiro = 2.34e10;
EITO = 1.14e11; ESnOx = 8.8e10; EAuEvap = 8.20e10; EAuScreen = 9e9; %Young's modulus

% For steady-state buckling, the energy release rate is:
% G_buckling = 0.5*(sigma^2*h)*(1-v^2)/E
% Calculate hc by setting G_buckling = Gc

disp("Buckling hc for sodalimeglass NIP devices:")
%Calculate ITO and SnOx Buckling hc%
hcITOglass = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)
hcSnOxglass = 2*GcSnOx*ESnOx/((1-vSnOx^2)*sigmaSnOx^2)

% For channeling and surface cracking, the energy release rate is:
% G_channel = Z_channel * (sigma^2*h)/E'

disp("Channeling and surface cracking for PVSK, and Spiro")
% Calculate channeling, surface cracking for PVSK, Spiro for soda lime
% glass
hcPVSKRSPP_channelglass = GcPVSKRSPP*EprimePVSK/(Zchannel*sigmaPVSKRSPP^2)
hcPVSKRSPP_surfaceglass = GcPVSKRSPP*EprimePVSK/(Zsurface*sigmaPVSKRSPP^2)
hcPVSKRSPP_spallingglass = GcPVSKRSPP*EprimePVSK/(Zspalling*sigmaPVSKRSPP^2)

hcPVSKspincoat_channelglass = GcPVSKspincoat*EprimePVSK/(Zchannel*sigmaPVSKspincoat^2)
hcPVSKspincoat_surfaceglass = GcPVSKspincoat*EprimePVSK/(Zsurface*sigmaPVSKspincoat^2)
hcPVSKspincoat_spallingglass = GcPVSKspincoat*EprimePVSK/(Zspalling*sigmaPVSKspincoat^2)

hcPVSKhotcast_channelglass = GcPVSKhotcast*EprimePVSK/(Zchannel*sigmaPVSKhotcast^2)
hcPVSKhotcast_surfaceglass = GcPVSKhotcast*EprimePVSK/(Zsurface*sigmaPVSKhotcast^2)
hcPVSKhotcast_spallingglass = GcPVSKhotcast*EprimePVSK/(Zspalling*sigmaPVSKhotcast^2)

hcSpiro_channelglass = GcSpiro*EprimeSpiro/(Zchannel*sigmaSpiro^2)
hcSpiro_surfaceglass = GcSpiro*EprimeSpiro/(Zsurface*sigmaSpiro^2)
hcSpiro_spallingglass = GcSpiro*EprimeSpiro/(Zspalling*sigmaSpiro^2)

disp("Buckling hc for polyimide NIP devices:")
%--------------%Buckling on Polyimide for NIP Devices%--------------%
sigmaITO = -9.93e8; sigmaSnOx = -3.31E+08; sigmaAuEvap = -4.53E+07; sigmaAuScreen = -1.68e6;
sigmaPVSKRSPP = 1.24e8; sigmaPVSKspincoat = 8.42e7; sigmaPVSKhotcast = 9.54e7;
sigmaSpiro = 3.52e7;

hcITObucklePI = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)
hcSnOxbucklePI = 2*GcSnOx*ESnOx/((1-vSnOx^2)*sigmaSnOx^2)
hcAuEvapbucklePI = 2*GcAuEvap*EAuEvap/((1-vAu_evap^2)*sigmaAuEvap^2)
hcAuScreenbucklePI = 2*GcAuScreen*EAuScreen/((1-vAu_screen^2)*sigmaAuScreen^2)

% Calculate channeling, surface cracking for PVSK, Spiro for polyimide
hcPVSKRSPP_channelglass = GcPVSKRSPP*EprimePVSK/(Zchannel*sigmaPVSKRSPP^2)
hcPVSKRSPP_surfaceglass = GcPVSKRSPP*EprimePVSK/(Zsurface*sigmaPVSKRSPP^2)
hcPVSKRSPP_spallingglass = GcPVSKRSPP*EprimePVSK/(Zspalling*sigmaPVSKRSPP^2)

hcPVSKspincoat_channelglass = GcPVSKspincoat*EprimePVSK/(Zchannel*sigmaPVSKspincoat^2)
hcPVSKspincoat_surfaceglass = GcPVSKspincoat*EprimePVSK/(Zsurface*sigmaPVSKspincoat^2)
hcPVSKspincoat_spallingglass = GcPVSKspincoat*EprimePVSK/(Zspalling*sigmaPVSKspincoat^2)

hcPVSKhotcast_channelglass = GcPVSKhotcast*EprimePVSK/(Zchannel*sigmaPVSKhotcast^2)
hcPVSKhotcast_surfaceglass = GcPVSKhotcast*EprimePVSK/(Zsurface*sigmaPVSKhotcast^2)
hcPVSKhotcast_spallingglass = GcPVSKhotcast*EprimePVSK/(Zspalling*sigmaPVSKhotcast^2)

hcSpiro_channelglass = GcSpiro*EprimeSpiro/(Zchannel*sigmaSpiro^2)
hcSpiro_surfaceglass = GcSpiro*EprimeSpiro/(Zsurface*sigmaSpiro^2)
hcSpiro_spallingglass = GcSpiro*EprimeSpiro/(Zspalling*sigmaSpiro^2)

