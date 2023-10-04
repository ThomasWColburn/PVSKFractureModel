clear all;

% Define poisson's ratios, CTE, Modulus, Gc, for all the layers.
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
E_SodaLimeGlass = 6.90E10/(1-v_SodaLimeGlass); % not used in the code
E_Polyimide = 7.50E9/(1-v_Polyimide); % not used in the code
E_ITO = 1.14E11/(1-v_ITO);
E_SnOx = 8.80E10/(1-v_SnOx);
E_PVSK = 1.11E10/(1-v_PVSK);
E_Spiro = 1.50E10/(1-v_Spiro);
E_Au_evap = 8.20E10/(1-v_Au_evap);
E_Au_screen = 9.00E9/(1-v_Au_screen);

% Define curvature, "b", and radius of curvature
k = 50E-9;
B = 2e-2;
R = 1/k;

% Define the thicknesses in each layer.
h_Au_evap = 150E-9;
h_Au_screen = 1000E-9;
h_Spiro = 180E-9;
h_PVSK = 400E-9;
h_SnOx = 50E-9;
h_ITO = 200E-9;

% Gc's for all the layers in n-i-p
Gc_ITO = 5;
Gc_SnOx = 0.7;
Gc_PVSKRSPP = 4.4;
Gc_PVSKhotcast = 0.92;
Gc_PVSKspincoat = 0.37;
Gc_Spiro = 0.45;
Gc_Au_evap = 1.2;
Gc_Au_screen = 2;

% Calculate the Gc for each interface (take the minimum of the bulk Gc of
% the layers that compose of the interface)
Gc_ITOSnOx = min(Gc_ITO, Gc_SnOx);
Gc_SnOxPVSKRSPP = min(Gc_SnOx, Gc_PVSKRSPP);
Gc_SnOxPVSKhotcast = min(Gc_SnOx, Gc_PVSKhotcast);
Gc_SnOxPVSKspincoat = min(Gc_SnOx, Gc_PVSKspincoat);

% Define the "default" temperatures felt by the layers during processing
delta_T_ITO_default = 400;
delta_T_SnOx_default = 155;
delta_T_PVSK_spin_default = 75;
delta_T_PVSK_hotcast_default = 85;
delta_T_PVSK_RSPP_default = 110;
delta_T_Spiro_default = 50;
delta_T_Au_evap_default = 50;
delta_T_Au_screen_default = 125;

%% Section 1: Vary delta_T for PVSK, on SodaLimeGlass & Polyimide substr.
% Varying PVSK processing T would impact ITO-SnOx & SnOx-PVSK interfaces
% Calculation of Gss with varying PVSK processing temperatures

%%% PVSK processing temperature variation
delta_T_PVSK = linspace(1, 200, 10000);

% First, let's compute the G_ss for the PVSK-SnOx interface for Glass subs.
% Calculate stresses induced by thermal expansion for Soda Lime Glass.
% Notice that both the ITO and the SnOx are in compression.
% So, no ITO-SnOx delaminatino will occur
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO_default*E_ITO;
sigma_SnOx_SodaLimeGlass = (a_SnOx - a_SodaLimeGlass)*delta_T_SnOx_default*E_SnOx;

sigma_PVSK_SodaLimeGlass = (a_PVSK - a_SodaLimeGlass).*delta_T_PVSK.*E_PVSK;

sigma_Spiro_SodaLimeGlass = (a_Spiro - a_SodaLimeGlass)*delta_T_Spiro_default*E_Spiro;
sigma_Au_evap_SodaLimeGlass = (a_Au_evap - a_SodaLimeGlass)*delta_T_Au_evap_default*E_Au_evap;
sigma_Au_screen_SodaLimeGlass = (a_Au_screen - a_SodaLimeGlass)*delta_T_Au_screen_default*E_Au_screen;

h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);
G_SnOxPVSK_SodaLimeGlass = sigma_Au_screen_SodaLimeGlass^2*h_Au_screen/E_Au_screen + ...
    sigma_Spiro_SodaLimeGlass^2*h_Spiro/E_Spiro + sigma_PVSK_SodaLimeGlass.^2*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);


% Second, let's compute the G_ss for the PVSK-SnOx interface for polyimide
% Calculate stresses induced by thermal expansion for polyimide
% Notice that both the ITO and the SnOx are in compression.
% So, no ITO-SnOx delaminatino will occur
sigma_ITO_Polyimide = (a_ITO - a_Polyimide)*delta_T_ITO_default*E_ITO;
sigma_SnOx_Polyimide = (a_SnOx - a_Polyimide)*delta_T_SnOx_default*E_SnOx;

sigma_PVSK_Polyimide = (a_PVSK - a_Polyimide).*delta_T_PVSK.*E_PVSK;

sigma_Spiro_Polyimide = (a_Spiro - a_Polyimide)*delta_T_Spiro_default*E_Spiro;
sigma_Au_evap_Polyimide = (a_Au_evap - a_Polyimide)*delta_T_Au_evap_default*E_Au_evap;
sigma_Au_screen_Polyimide = (a_Au_screen - a_Polyimide)*delta_T_Au_screen_default*E_Au_screen;

h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);
G_SnOxPVSK_Polyimide = sigma_Au_screen_Polyimide^2*h_Au_screen/E_Au_screen + ...
    sigma_Spiro_Polyimide^2*h_Spiro/E_Spiro + sigma_PVSK_Polyimide.^2*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

figure; hold on
plot(delta_T_PVSK, G_SnOxPVSK_SodaLimeGlass, 'k', "LineWidth", 2);
plot(delta_T_PVSK, G_SnOxPVSK_Polyimide, 'm', "LineWidth", 2);
plot(delta_T_PVSK, Gc_SnOxPVSKhotcast.*ones(1,length(delta_T_PVSK)), '--g', "LineWidth", 1.5);
plot(delta_T_PVSK, Gc_SnOxPVSKspincoat.*ones(1,length(delta_T_PVSK)), '--b', "LineWidth", 1.5);
plot(delta_T_PVSK, Gc_SnOxPVSKRSPP.*ones(1,length(delta_T_PVSK)), '--r', "LineWidth", 1.5);
legend("G_{ss,SodaLimeGlass}", "G_{ss,Polyimide}", "G_{c,SnOx-PVSKhotcast}", "G_{c,SnOx-PVSKspincoat}")
    % "G_{c,SnOx-PVSKRSPP}", )
ylabel("G, (J/m^2)");
xlabel("{\Delta}T_{PVSK}, (K)");
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');

%% Section 2: Computing the delta_T_PVSK,max for SnOx-PVSK delamination
syms delta_T_PVSK

% Compute driving force for SnOx-PVSK delamination for SLG substrate
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO_default*E_ITO;
sigma_SnOx_SodaLimeGlass = (a_SnOx - a_SodaLimeGlass)*delta_T_SnOx_default*E_SnOx;

sigma_PVSK_SodaLimeGlass = (a_PVSK - a_SodaLimeGlass).*delta_T_PVSK.*E_PVSK;

sigma_Spiro_SodaLimeGlass = (a_Spiro - a_SodaLimeGlass)*delta_T_Spiro_default*E_Spiro;
sigma_Au_evap_SodaLimeGlass = (a_Au_evap - a_SodaLimeGlass)*delta_T_Au_evap_default*E_Au_evap;
sigma_Au_screen_SodaLimeGlass = (a_Au_screen - a_SodaLimeGlass)*delta_T_Au_screen_default*E_Au_screen;

h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);
G_SnOxPVSK_SodaLimeGlass = sigma_Au_screen_SodaLimeGlass^2*h_Au_screen/E_Au_screen + ...
    sigma_Spiro_SodaLimeGlass^2*h_Spiro/E_Spiro + sigma_PVSK_SodaLimeGlass.^2*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);


% Compute driving force for SnOx-PVSK delamination for polyimide substrate
sigma_ITO_Polyimide = (a_ITO - a_Polyimide)*delta_T_ITO_default*E_ITO;
sigma_SnOx_Polyimide = (a_SnOx - a_Polyimide)*delta_T_SnOx_default*E_SnOx;

sigma_PVSK_Polyimide = (a_PVSK - a_Polyimide).*delta_T_PVSK.*E_PVSK;

sigma_Spiro_Polyimide = (a_Spiro - a_Polyimide)*delta_T_Spiro_default*E_Spiro;
sigma_Au_evap_Polyimide = (a_Au_evap - a_Polyimide)*delta_T_Au_evap_default*E_Au_evap;
sigma_Au_screen_Polyimide = (a_Au_screen - a_Polyimide)*delta_T_Au_screen_default*E_Au_screen;

h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);
G_SnOxPVSK_Polyimide = sigma_Au_screen_Polyimide^2*h_Au_screen/E_Au_screen + ...
    sigma_Spiro_Polyimide^2*h_Spiro/E_Spiro + sigma_PVSK_Polyimide.^2*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);


% Solve for the processing temperature at which debonding would occur
delta_T_PVSKRSPP_SodaLimeGlass = vpasolve(G_SnOxPVSK_SodaLimeGlass == Gc_SnOxPVSKRSPP, delta_T_PVSK, [0,inf]);
delta_T_PVSKRSPP_SodaLimeGlass = double(delta_T_PVSKRSPP_SodaLimeGlass);

delta_T_PVSKhotcast_SodaLimeGlass = vpasolve(G_SnOxPVSK_SodaLimeGlass == Gc_SnOxPVSKhotcast, delta_T_PVSK, [0,inf]);
delta_T_PVSKhotcast_SodaLimeGlass = double(delta_T_PVSKhotcast_SodaLimeGlass);

delta_T_PVSKspincoat_SodaLimeGlass = vpasolve(G_SnOxPVSK_SodaLimeGlass == Gc_SnOxPVSKspincoat, delta_T_PVSK, [0,inf]);
delta_T_PVSKspincoat_SodaLimeGlass = double(delta_T_PVSKspincoat_SodaLimeGlass);

delta_T_PVSKRSPP_Polyimide = vpasolve(G_SnOxPVSK_Polyimide == Gc_SnOxPVSKRSPP, delta_T_PVSK, [0,inf]);
delta_T_PVSKRSPP_Polyimide = double(delta_T_PVSKRSPP_Polyimide);

delta_T_PVSKhotcast_Polyimide = vpasolve(G_SnOxPVSK_Polyimide == Gc_SnOxPVSKhotcast, delta_T_PVSK, [0,inf]);
delta_T_PVSKhotcast_Polyimide = double(delta_T_PVSKhotcast_Polyimide);

delta_T_PVSKspincoat_Polyimide = vpasolve(G_SnOxPVSK_Polyimide == Gc_SnOxPVSKspincoat, delta_T_PVSK, [0,inf]);
delta_T_PVSKspincoat_Polyimide = double(delta_T_PVSKspincoat_Polyimide);

T_ambient = 25;

figure; hold on
ylim([0 250]);
x = [delta_T_PVSKRSPP_SodaLimeGlass+T_ambient delta_T_PVSKhotcast_SodaLimeGlass+T_ambient ...
    delta_T_PVSKspincoat_SodaLimeGlass+T_ambient delta_T_PVSKRSPP_Polyimide+T_ambient ...
    delta_T_PVSKhotcast_Polyimide+T_ambient delta_T_PVSKspincoat_Polyimide+T_ambient];
b = bar(x, 'FaceColor', 'flat');
text(1:length(x),x,num2str(x'),'vert','bottom','horiz','center'); 
b.CData(1,:) = [0.8500 0.3250 0.0980];
b.CData(2,:) = [0.8500 0.3250 0.0980];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.9290 0.6940 0.1250];
b.CData(5,:) = [0.9290 0.6940 0.1250];
b.CData(6,:) = [0.9290 0.6940 0.1250];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt(2:7), 'XTickLabel', {'RSPP', 'HotCast', 'SpinCoat', 'RSPP', 'HotCast', 'SpinCoat'});
ylabel('T_{PVSK,max}')
box on;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
title("Max Processing T_{PVSK}, PVSK-SnOx Debond, NIP", 'FontWeight','Normal');

% Make the legend by using some manipulation
bh(1) = bar(nan,nan, 'FaceColor', "#D95319");
bh(2) = bar(nan,nan, 'FaceColor', "#EDB120");
legend(bh, ["SodaLimeGlass", "Polyimide"]);

delta_Tm_PVSK = 175;
xlim=get(gca,'xlim');
plot(xlim, [delta_Tm_PVSK+T_ambient delta_Tm_PVSK+T_ambient], "LineWidth", 1.5, 'Color', "k", ...
    "LineStyle", "--", "DisplayName", "T_{decomp,PVSK}");
ylim([0,250]);

