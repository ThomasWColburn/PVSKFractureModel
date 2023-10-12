%--------------%Channel Cracking and Buckling on Soda Lime NIP Devices%--------------%
clear all;
%Thomas Colburn%
%Set the constants%
GcITO = 5; GcSnOx = 0.7; GcAuEvap = 1.2; GcAuScreen = 2; %Fracture energy release rate J/m^2
sigmaITO = -2.23e8; sigmaSnOx = -1.07e8; %Thermal stresses
vITO = 0.35; vSnOx = 0.33; vAu_evap = 0.421; vAu_screen = 0.33; %Poisson's Ratio
Zchannel = 1.976; Zsurface = 3.951; %Nondimensional geometric parameter
EprimeITO = 1.75e11; EprimeSnOx = 1.31e11; EprimeAuEvap = 1.42e11; EprimeAuScreen = 1.34e10; %Young's modulus divided by 1-v Pa
EITO = 1.14e11; ESnOx = 8.8e10; EAuEvap = 8.20e10; EAuScreen = 9e9; %Young's modulus

disp("Channel cracking and buckling hc for sodalimeglass NIP devices:")
%Calculate ITO and SnOx Buckling hc%
hcITOglass = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)
hcSnOxglass = 2*GcSnOx*ESnOx/((1-vSnOx^2)*sigmaSnOx^2)

disp("Buckling hc for polyimide NIP devices:")
%--------------%Buckling on Polyimide for PIN Devices%--------------%
sigmaITO = -9.93e8; sigmaSnOx = -3.31E+08; sigmaAuEvap = -4.53E+07; sigmaAuScreen = -1.68e6;
hcITObucklePI = 2*GcITO*EITO/((1-vITO^2)*sigmaITO^2)
hcSnOxbucklePI = 2*GcSnOx*ESnOx/((1-vSnOx^2)*sigmaSnOx^2)
hcAuEvapbucklePI = 2*GcAuEvap*EAuEvap/((1-vAu_evap^2)*sigmaAuEvap^2)
hcAuScreenbucklePI = 2*GcAuScreen*EAuScreen/((1-vAu_screen^2)*sigmaAuScreen^2)

%{
%% Thermal Analysis for buckling. Relevant Layers: ITO, SnOx
% Gss_buckling = sigma^2*h / (2*E_eff), where E_eff = effective modulus
%Define poisson's ratios, CTE, Modulus, Gc, for all the layers.
v_SodaLimeGlass = 0.23;
v_Polyimide = 0.35;
v_ITO = 0.35;
v_SnOx = 0.33;
v_PVSK = 0.31;
v_Spiro = 0.36;
v_Au_evap = 0.421;
v_Au_screen = 0.33;

% Coefficient of thermal expansions for the materials.
a_SodaLimeGlass = 9.00E-6;
a_Polyimide = 2.00E-5;
a_ITO = 5.81E-6;
a_SnOx = 3.75E-6;
a_PVSK = 8.98E-5;
a_Spiro = 5.00E-5;
a_Au_evap = 1.36E-5;
a_Au_screen = 1.90E-5;

% Define the Elastic moduli of all the layers including Silver screen
% printed and evaporated Ag.
E_SodaLimeGlass = 6.90E10; % not used in the code
E_Polyimide = 7.50E9; % not used in the code
E_ITO = 1.14E11;
E_SnOx = 8.80E10;
E_PVSK = 1.11E10;
E_Spiro = 1.50E10;
E_Au_evap = 8.20E10;
E_Au_screen = 9.00E9;

% Define the thicknesses in each layer.
h_SnOx = 50E-9;
h_ITO = 200E-9;


syms delta_T_ITO
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO*E_ITO/(1-v_ITO);
sigma_ITO_Polyimide = (a_ITO - a_Polyimide)*delta_T_ITO*E_ITO/(1-v_ITO);

G_buckle_ITO_SodaLimeGlass = sigma_ITO_SodaLimeGlass^2*h_ITO/(2*E_ITO/(1-v_ITO^2));
G_buckle_ITO_Polyimide = sigma_ITO_Polyimide^2*h_ITO/(2*E_ITO/(1-v_ITO^2));

delta_T_ITO_SodaLimeGlass = vpasolve(G_buckle_ITO_SodaLimeGlass == GcITO, delta_T_ITO, [0,inf]);
delta_T_ITO_SodaLimeGlass = double(delta_T_ITO_SodaLimeGlass)

delta_T_ITO_Polyimide = vpasolve(G_buckle_ITO_Polyimide == GcITO, delta_T_ITO, [0,inf]);
delta_T_ITO_Polyimide = double(delta_T_ITO_Polyimide)


syms delta_T_SnOx
sigma_SnOx_SodaLimeGlass = (a_SnOx - a_SodaLimeGlass)*delta_T_SnOx*E_SnOx/(1-v_SnOx);
sigma_SnOx_Polyimide = (a_SnOx - a_Polyimide)*delta_T_SnOx*E_SnOx/(1-v_SnOx);

G_buckle_SnOx_SodaLimeGlass = sigma_SnOx_SodaLimeGlass^2*h_SnOx/(2*E_SnOx/(1-v_SnOx^2));
G_buckle_SnOx_Polyimide = sigma_SnOx_Polyimide^2*h_SnOx/(2*E_SnOx/(1-v_SnOx^2));

delta_T_SnOx_SodaLimeGlass = vpasolve(G_buckle_SnOx_SodaLimeGlass == GcSnOx, delta_T_SnOx, [0,inf]);
delta_T_SnOx_SodaLimeGlass = double(delta_T_SnOx_SodaLimeGlass)

delta_T_SnOx_Polyimide = vpasolve(G_buckle_SnOx_Polyimide == GcSnOx, delta_T_SnOx, [0,inf]);
delta_T_SnOx_Polyimide = double(delta_T_SnOx_Polyimide)

figure;
x = [delta_T_ITO_SodaLimeGlass delta_T_ITO_Polyimide delta_T_SnOx_SodaLimeGlass delta_T_SnOx_Polyimide];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.4940 0.1840 0.5560];
b.CData(2,:) = [0.4940 0.1840 0.5560];
b.CData(3,:) = [0.4660 0.6740 0.1880];
b.CData(4,:) = [0.4660 0.6740 0.1880];
text(1:length(x),x,num2str(x'),'vert','bottom','horiz','center'); 
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'ITO, SodaLimeGlass', 'ITO, Polyimide', ...
    'SnOx, SodaLimeGlass', 'SnOx, Polyimide'});
ylabel('{\Delta}T')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("{\Delta}T_{max} for buckling (N-I-P)");
%}

