%--------------%Channel Cracking and Buckling on Soda Lime PIN Devices%--------------%
clear all;

%Thomas Colburn%
%Set the constants%
GcITO = 5; GcNiO = 5; GcAgEvap = 5; GcAgScreen = 2; %Fracture energy release rate J/m^2
sigmaITO = -2.23e8; sigmaNiO = 1.64e8; %Thermal stresses
vITO = 0.35; vNiO = 0.416; vAg = 0.33; %Poisson's Ratio
Zchannel = 1.976; Zsurface = 3.951; %Nondimensional geometric parameter
EprimeITO = 1.75e11; EprimeNiO = 1.64e11; EprimeAgEvap = 1.19e11; EprimeAgScreen = 1.34e10; %Young's modulus divided by 1-v Pa
EITO = 1.14e11; ENiO = 9.58e10; EAgEvap = 8e10; EAgScreen = 9e9; %Young's modulus

disp("Channel cracking and buckling hc for sodalimeglass PIN devices:")
%Calculate ITO Buckling hc%
hcITOglass_buckle = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)

%Calculate NiOx channel vs surface crack
hcNiOx_channelglass = GcNiO*EprimeNiO/(Zchannel*sigmaNiO^2)
hcNiOx_surfaceglass = GcNiO*EprimeNiO/(Zsurface*sigmaNiO^2)

disp("Buckling hc for polyimide PIN devices:")
%--------------%Buckling on Polyimide for PIN Devices%--------------%
sigmaITO = -9.93e8; sigmaNiO = -3.32E+08; sigmaAgEvap = -5.97E+06; sigmaAgScreen = -1.68e6;
hcITObucklePI = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)
hcNiObucklePI = 2*GcNiO*ENiO/((1-vNiO^2)*sigmaNiO^2)
hcAgEvapbucklePI = 2*GcAgEvap*EAgEvap/((1-vAg^2)*sigmaAgEvap^2)
hcAgScreenbucklePI = 2*GcAgScreen*EAgScreen/((1-vAg^2)*sigmaAgScreen^2)


%{
%% Thermal analysis for buckling. Relevant layers: ITO, NiO
% Gss_buckling = sigma^2*h / (2*E_eff), where E_eff = effective modulus

% Set heights of ITO and NiO to be the default heights.
h_NiO = 20E-9;
h_ITO = 200E-9;

%Define poisson's ratios, CTE, Modulus, Gc, for all the layers.
v_SodaLimeGlass = 0.23;
v_Polyimide = 0.35;
v_ITO = 0.35;
v_NiO = 0.416;
v_PVSK = 0.31;
v_PCBM = 0.23;
v_Ag = 0.33; % for both screen printed and evap Ag.

% Coefficient of thermal expansions for the materials.
a_SodaLimeGlass = 9.00E-6;
a_Polyimide = 2.00E-5;
a_ITO = 5.81E-6;
a_NiO = 1.26E-5;
a_PVSK = 89.8E-6;
a_PCBM = 1.20E-4;
a_Ag_evap = 1.90E-5;
a_Ag_screen = 1.90E-5;

% Define the Elastic moduli of all the layers including Silver screen
% printed and evaporated Ag.
E_SodaLimeGlass = 6.90E10; % not used in the code
E_Polyimide = 7.50E9; % not used in the code
E_ITO = 1.14E11;
E_NiO = 9.58E10;
E_PVSK = 1.11E10;
E_PCBM = 3.00E9;
E_Ag_evap = 8.00E10;
E_Ag_screen = 9.00E9;

syms delta_T_ITO
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO*E_ITO/(1-v_ITO);
sigma_ITO_Polyimide = (a_ITO - a_Polyimide)*delta_T_ITO*E_ITO/(1-v_ITO);

G_buckle_ITO_SodaLimeGlass = sigma_ITO_SodaLimeGlass^2*h_ITO/(2*E_ITO/(1-v_ITO^2));
G_buckle_ITO_Polyimide = sigma_ITO_Polyimide^2*h_ITO/(2*E_ITO/(1-v_ITO^2));

delta_T_ITO_SodaLimeGlass = vpasolve(G_buckle_ITO_SodaLimeGlass == GcITO, delta_T_ITO, [0,inf]);
delta_T_ITO_SodaLimeGlass = double(delta_T_ITO_SodaLimeGlass)

delta_T_ITO_Polyimide = vpasolve(G_buckle_ITO_Polyimide == GcITO, delta_T_ITO, [0,inf]);
delta_T_ITO_Polyimide = double(delta_T_ITO_Polyimide)


syms delta_T_NiO
sigma_NiO_Polyimide = (a_NiO - a_Polyimide)*delta_T_NiO*E_NiO/(1-v_NiO);

G_buckle_NiO_Polyimide = sigma_NiO_Polyimide^2*h_NiO/(2*E_NiO/(1-v_NiO^2));

delta_T_NiO_Polyimide = vpasolve(G_buckle_NiO_Polyimide == GcNiO, delta_T_NiO, [0,inf]);
delta_T_NiO_Polyimide = double(delta_T_NiO_Polyimide)


figure;
x = [delta_T_ITO_SodaLimeGlass delta_T_ITO_Polyimide delta_T_NiO_Polyimide];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.4940 0.1840 0.5560];
b.CData(2,:) = [0.4940 0.1840 0.5560];
b.CData(3,:) = [0.4660 0.6740 0.1880];
text(1:length(x),x,num2str(x'),'vert','bottom','horiz','center'); 
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'ITO, SodaLimeGlass', 'ITO, Polyimide', 'NiO, Polyimide'});
ylabel('{\Delta}T')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("{\Delta}T_{max} for buckling (P-I-N)");
%}

