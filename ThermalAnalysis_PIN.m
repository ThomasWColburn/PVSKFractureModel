clear all;

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
E_SodaLimeGlass = 6.90E10/(1-v_SodaLimeGlass); % not used in the code
E_Polyimide = 7.50E9/(1-v_Polyimide); % not used in the code
E_ITO = 1.14E11/(1-v_ITO);
E_NiO = 9.58E10/(1-v_NiO);
E_PVSK = 1.11E10/(1-v_PVSK);
E_PCBM = 3.00E9/(1-v_PCBM);
E_Ag_evap = 8.00E10/(1-v_Ag);
E_Ag_screen = 9.00E9/(1-v_Ag);

% Define curvature, "b", and radius of curvature
k = 50E-9;
B = 2e-2;
R = 1/k;

% Define the thicknesses in each layer.
h_Ag_evap = 150E-9;
h_Ag_screen = 1000E-9;
h_PCBM = 40E-9;
h_PVSK = 400E-9;
h_NiO = 20E-9;
h_ITO = 200E-9;

% Gc for PCBM source
% https://pubs.acs.org/doi/full/10.1021/acsami.5b02202
Gc_ITO = 36;
Gc_NiO = 5;
Gc_PVSKRSPP = 4.4;
Gc_PVSKhotcast = 0.92;
Gc_PVSKspincoat = 0.37;
Gc_PCBM = 0.13;
Gc_Ag_evap = 5;
Gc_Ag_screen = 2;

% Calculate the Gc for each interface (take the minimum of the bulk Gc of
% the layers that compose of the interface)
Gc_ITONiO = min(Gc_NiO, Gc_ITO);
Gc_NiOPVSKRSPP = min(Gc_NiO, Gc_PVSKRSPP);
Gc_NiOPVSKhotcast = min(Gc_NiO, Gc_PVSKhotcast);
Gc_NiOPVSKspincoat = min(Gc_NiO, Gc_PVSKspincoat);
Gc_PVSKPCBM = Gc_PCBM;
Gc_PCBMAgEvap = min(Gc_PCBM, Gc_Ag_evap);
Gc_PCBMAgScreen = min(Gc_PCBM, Gc_Ag_screen);

% Define the "default" temperatures felt by the layers during processing
delta_T_ITO_default = 325;
delta_T_NiO_default = 275;
delta_T_PVSK_spin_default = 75;
delta_T_PVSK_hotcast_default = 85;
delta_T_PVSK_RSPP_default = 110;
delta_T_PCBM_default = 50;
delta_T_Ag_evap_default = 50;
delta_T_Ag_screen_default = 125;

%% Section 1: Vary delta_T for PVSK, SodaLimeGlass & Polyimide.
% Varying PVSK processing T would impact ITO-NiO & NiO-PVSK interfaces

%%% PVSK processing temperature variation
delta_T_PVSK = linspace(1, 300, 10000);

% Calculate stresses induced by thermal expansion for Soda Lime Glass
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO_default*E_ITO;
sigma_NiO_SodaLimeGlass = (a_NiO - a_SodaLimeGlass)*delta_T_NiO_default*E_NiO;

sigma_PVSK_SodaLimeGlass = (a_PVSK - a_SodaLimeGlass).*delta_T_PVSK.*E_PVSK;

sigma_PCBM_SodaLimeGlass = (a_PCBM - a_SodaLimeGlass)*delta_T_PCBM_default*E_PCBM;
sigma_Ag_evap_SodaLimeGlass = (a_Ag_evap - a_SodaLimeGlass)*delta_T_Ag_evap_default*E_Ag_evap;
sigma_Ag_screen_SodaLimeGlass = (a_Ag_screen - a_SodaLimeGlass)*delta_T_Ag_screen_default*E_Ag_screen;

% Calculate Gss for both NiO-ITO delamination and NiO-PVSK delamination
% as a function of delta_T_PVSK for evaporated Ag. (Soda Lime Glass)
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);
G_ITONiO_AgEvap_SodaLimeGlass = sigma_Ag_evap_SodaLimeGlass^2*h_Ag_evap/E_Ag_evap + sigma_PCBM_SodaLimeGlass^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_SodaLimeGlass.^2.*h_PVSK./E_PVSK + sigma_NiO_SodaLimeGlass^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

h_summation = h_Ag_evap + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);
G_NiOPVSK_AgEvap_SodaLimeGlass = sigma_Ag_evap_SodaLimeGlass^2*h_Ag_evap/E_Ag_evap + sigma_PCBM_SodaLimeGlass^2*h_PCBM/E_PCBM ...
    + sigma_PVSK_SodaLimeGlass.^2.*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);


% Calculate stresses induced by thermal expansion for Soda Lime Glass
sigma_ITO_Polyimide = (a_ITO - a_Polyimide)*delta_T_ITO_default*E_ITO;
sigma_NiO_Polyimide = (a_NiO - a_Polyimide)*delta_T_NiO_default*E_NiO;

sigma_PVSK_Polyimide = (a_PVSK - a_Polyimide).*delta_T_PVSK.*E_PVSK;

sigma_PCBM_Polyimide = (a_PCBM - a_Polyimide)*delta_T_PCBM_default*E_PCBM;
sigma_Ag_evap_Polyimide = (a_Ag_evap - a_Polyimide)*delta_T_Ag_evap_default*E_Ag_evap;
sigma_Ag_screen_Polyimide = (a_Ag_evap - a_Polyimide)*delta_T_Ag_screen_default*E_Ag_screen;


% Calculate Gss for NiO-PVSK delamination as a function of delta_T_PVSK for
% evaporated Ag. (Polyimide)
h_summation = h_Ag_evap + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);
G_NiOPVSK_AgEvap_Polyimide = sigma_Ag_evap_Polyimide^2*h_Ag_evap/E_Ag_evap + sigma_PCBM_Polyimide^2*h_PCBM/E_PCBM ...
    + sigma_PVSK_Polyimide.^2.*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

figure; hold on
% plot(delta_T_PVSK, G_ITONiO_AgEvap_SodaLimeGlass, 'k', "LineWidth", 2);
plot(delta_T_PVSK, G_NiOPVSK_AgEvap_SodaLimeGlass, 'k', "LineWidth", 2);
plot(delta_T_PVSK, G_NiOPVSK_AgEvap_Polyimide, 'g', "LineWidth", 2);
plot(delta_T_PVSK, Gc_ITONiO.*ones(1,length(delta_T_PVSK)), '--r', "LineWidth", 1.5);
plot(delta_T_PVSK, Gc_NiOPVSKRSPP.*ones(1,length(delta_T_PVSK)), '--b', "LineWidth", 1.5);
plot(delta_T_PVSK, Gc_NiOPVSKhotcast.*ones(1,length(delta_T_PVSK)), '--', "LineWidth", 1.5);
plot(delta_T_PVSK, Gc_NiOPVSKspincoat.*ones(1,length(delta_T_PVSK)), '--', "LineWidth", 1.5);
legend("G_{ss,SodaLimeGlass}", "G_{ss,Polyimide}")
% "G_{c,NiO-ITO}", ..."G_{c,NiO-PVSKRSPP}", "G_{c,NiO-PVSKhotcast}", "G_{c,NiO-PVSKspin}
ylabel("G, (J/m^2)");
xlabel("{\Delta}T_{PVSK}, (K)");
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');


%% Section for solution for solving the max delta_T_PVSK
% Now, just solve for the delta_T_PVSK solution that would give us
% debonding for each one of the Gc's.
syms delta_T_PVSK

% Calculate stresses induced by thermal expansion for Soda Lime Glass
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO_default*E_ITO;
sigma_NiO_SodaLimeGlass = (a_NiO - a_SodaLimeGlass)*delta_T_NiO_default*E_NiO;

sigma_PVSK_SodaLimeGlass = (a_PVSK - a_SodaLimeGlass).*delta_T_PVSK.*E_PVSK;

sigma_PCBM_SodaLimeGlass = (a_PCBM - a_SodaLimeGlass)*delta_T_PCBM_default*E_PCBM;
sigma_Ag_evap_SodaLimeGlass = (a_Ag_evap - a_SodaLimeGlass)*delta_T_Ag_evap_default*E_Ag_evap;
sigma_Ag_screen_SodaLimeGlass = (a_Ag_evap - a_SodaLimeGlass)*delta_T_Ag_screen_default*E_Ag_screen;

h_summation = h_Ag_evap + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);
G_NiOPVSK_AgEvap_SodaLimeGlass = sigma_Ag_evap_SodaLimeGlass^2*h_Ag_evap/E_Ag_evap + sigma_PCBM_SodaLimeGlass^2*h_PCBM/E_PCBM ...
    + sigma_PVSK_SodaLimeGlass.^2.*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);



% Calculate stresses induced by thermal expansion for Soda Lime Glass
sigma_ITO_Polyimide = (a_ITO - a_Polyimide)*delta_T_ITO_default*E_ITO;
sigma_NiO_Polyimide = (a_NiO - a_Polyimide)*delta_T_NiO_default*E_NiO;

sigma_PVSK_Polyimide = (a_PVSK - a_Polyimide).*delta_T_PVSK.*E_PVSK;

sigma_PCBM_Polyimide = (a_PCBM - a_Polyimide)*delta_T_PCBM_default*E_PCBM;
sigma_Ag_evap_Polyimide = (a_Ag_evap - a_Polyimide)*delta_T_Ag_evap_default*E_Ag_evap;
sigma_Ag_screen_Polyimide = (a_Ag_evap - a_Polyimide)*delta_T_Ag_screen_default*E_Ag_screen;

% Calculate Gss for NiO-PVSK delamination as a function of delta_T_PVSK for
% evaporated Ag. (Polyimide)
h_summation = h_Ag_evap + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);
G_NiOPVSK_AgEvap_Polyimide = sigma_Ag_evap_Polyimide^2*h_Ag_evap/E_Ag_evap + sigma_PCBM_Polyimide^2*h_PCBM/E_PCBM ...
    + sigma_PVSK_Polyimide.^2.*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

delta_T_PVSKRSPP_SodaLimeGlass = vpasolve(G_NiOPVSK_AgEvap_SodaLimeGlass == Gc_NiOPVSKRSPP, delta_T_PVSK, [0,inf]);
delta_T_PVSKRSPP_SodaLimeGlass = double(delta_T_PVSKRSPP_SodaLimeGlass);

delta_T_PVSKhotcast_SodaLimeGlass = vpasolve(G_NiOPVSK_AgEvap_SodaLimeGlass == Gc_NiOPVSKhotcast, delta_T_PVSK, [0,inf]);
delta_T_PVSKhotcast_SodaLimeGlass = double(delta_T_PVSKhotcast_SodaLimeGlass);

delta_T_PVSKspincoat_SodaLimeGlass = vpasolve(G_NiOPVSK_AgEvap_SodaLimeGlass == Gc_NiOPVSKspincoat, delta_T_PVSK, [0,inf]);
delta_T_PVSKspincoat_SodaLimeGlass = double(delta_T_PVSKspincoat_SodaLimeGlass);

delta_T_PVSKRSPP_Polyimide = vpasolve(G_NiOPVSK_AgEvap_Polyimide == Gc_NiOPVSKRSPP, delta_T_PVSK, [0,inf]);
delta_T_PVSKRSPP_Polyimide = double(delta_T_PVSKRSPP_Polyimide);

delta_T_PVSKhotcast_Polyimide = vpasolve(G_NiOPVSK_AgEvap_Polyimide == Gc_NiOPVSKhotcast, delta_T_PVSK, [0,inf]);
delta_T_PVSKhotcast_Polyimide = double(delta_T_PVSKhotcast_Polyimide);

delta_T_PVSKspincoat_Polyimide = vpasolve(G_NiOPVSK_AgEvap_Polyimide == Gc_NiOPVSKspincoat, delta_T_PVSK, [0,inf]);
delta_T_PVSKspincoat_Polyimide = double(delta_T_PVSKspincoat_Polyimide);

T_ambient = 25; % in celsius

figure; hold on
x = [delta_T_PVSKRSPP_SodaLimeGlass+T_ambient delta_T_PVSKhotcast_SodaLimeGlass+T_ambient ...
    delta_T_PVSKspincoat_SodaLimeGlass+T_ambient delta_T_PVSKRSPP_Polyimide+T_ambient...
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
title("Max Processing T_{PVSK}, PVSK-NiO Debond, PIN", 'FontWeight','Normal');

% Make the legend by using some manipulation
bh(1) = bar(nan,nan, 'FaceColor', "#D95319");
bh(2) = bar(nan,nan, 'FaceColor', "#EDB120");
legend(bh, ["SodaLimeGlass", "Polyimide"]);

delta_Tm_PVSK = 175;
xlim=get(gca,'xlim');
plot(xlim, [delta_Tm_PVSK+T_ambient delta_Tm_PVSK+T_ambient], "LineWidth", 1.5, 'Color', "k", "LineStyle", "--", ...
    "DisplayName", "T_{decomp,PVSK}");
ylim([0,450]);

%%
%{
% Calculate Gss for both NiO-ITO delamination and NiO-PVSK delamination
% as a function of delta_T_PVSK_RSPP for screen-printed Ag
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);
G_ITONiO_AgScreen = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK.^2.*h_PVSK./E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

h_summation = h_Ag_screen + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);
G_NiOPVSK_AgScreen = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM ...
    + sigma_PVSK.^2.*h_PVSK./E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

figure; hold on
plot(delta_T_PVSK, G_ITONiO_AgScreen, 'r');
plot(delta_T_PVSK, G_NiOPVSK_AgScreen, 'b');
plot(delta_T_PVSK, Gc_ITONiO.*ones(1,length(delta_T_PVSK)), '--r');
plot(delta_T_PVSK, Gc_NiOPVSKRSPP.*ones(1,length(delta_T_PVSK)), '--b');
plot(delta_T_PVSK, Gc_NiOPVSKhotcast.*ones(1,length(delta_T_PVSK)), '--');
plot(delta_T_PVSK, Gc_NiOPVSKspincoat.*ones(1,length(delta_T_PVSK)), '--');
%}

%% Section 2: Varying the NiO processing temperature
%%% Now, do NiO processing temperature variation
delta_T_NiO = linspace(1, 8000, 10000);

% Calculate stresses induced by thermal expansion for Soda Lime Glass
sigma_ITO_SodaLimeGlass = (a_ITO - a_SodaLimeGlass)*delta_T_ITO_default*E_ITO;

sigma_NiO_SodaLimeGlass = (a_NiO - a_SodaLimeGlass).*delta_T_NiO.*E_NiO;

sigma_PVSK_SodaLimeGlass = (a_PVSK - a_SodaLimeGlass)*delta_T_PVSK_RSPP_default*E_PVSK;
sigma_PCBM_SodaLimeGlass = (a_PCBM - a_SodaLimeGlass)*delta_T_PCBM_default*E_PCBM;
sigma_Ag_evap_SodaLimeGlass = (a_Ag_evap - a_SodaLimeGlass)*delta_T_Ag_evap_default*E_Ag_evap;
sigma_Ag_screen_SodaLimeGlass = (a_Ag_evap - a_SodaLimeGlass)*delta_T_Ag_screen_default*E_Ag_screen;

% Calculate Gss for NiO-ITO delamination as a function of delta_T_NiO for evaporated Ag. (Soda Lime Glass)
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);
G_ITONiO_AgEvap_SodaLimeGlass = sigma_Ag_evap_SodaLimeGlass^2*h_Ag_evap/E_Ag_evap + sigma_PCBM_SodaLimeGlass^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_SodaLimeGlass^2*h_PVSK/E_PVSK + sigma_NiO_SodaLimeGlass.^2.*h_NiO./E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

figure; hold on
% plot(delta_T_PVSK, G_ITONiO_AgEvap_SodaLimeGlass, 'k', "LineWidth", 2);
plot(delta_T_NiO, G_ITONiO_AgEvap_SodaLimeGlass, 'k', "LineWidth", 2);
plot(delta_T_NiO, Gc_ITONiO.*ones(1,length(delta_T_NiO)), '--r', "LineWidth", 2);
legend("G_{ss,SodaLimeGlass}", "G_{c,NiO-ITO}")
ylabel("G, (J/m^2)");
xlabel("{\Delta}T_{NiO}, (K)");
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');

