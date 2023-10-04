%-------------Perovskites on Polyimide-------------%
%%%%%%%%%%%%% Thomas Colburn, Alan Liu, PVSK debond model %%%%%%%%%%%%%%
% PVSK on Polyimide substrate, P-I-N structure %

clear;

% P-I-N architecture for PVSK solar modules on Polyimide substrate
% The pvsk. layers are as follows:
    % Polyimide -> ITO -> NiO -> PVSK -> PCBM -> Ag

% BAR CHART COLOR CODES:
% [0.8500 0.3250 0.0980] = Perovskite RGB code
% [0.5 0.5 0.5] = Ag RGB code
% [0 0.4470 0.7410] = PCBM RGB code
% [0.6350 0.0780 0.1840] = NiO RGB code
% [0.4940 0.1840 0.5560] = ITO RGB code

% Go from top of the stack to the bottom of the stack for delamination
% calculations.

% Define Poisson's ratio for each layer
v_Polyimide = 0.35;
v_ITO = 0.35;
v_NiO = 0.416;
v_PVSK = 0.31;
v_PCBM = 0.23;
v_Ag = 0.33; % for both screen printed and evap Ag.

% Define the coefficient of thermal expansion (CTE) for each layer
a_Polyimide = 2.00E-5;
a_ITO = 5.81E-6;
a_NiO = 1.26E-5;
a_PVSK = 89.8E-6;
a_PCBM = 1.20E-4;
a_Ag = 1.90E-5;

% Define the Elastic moduli of all the layers including Ag screen
% printed and evaporated Ag.
E_Polyimide = 7.50E9/(1-v_Polyimide); % not used in the code
E_ITO = 1.14E11/(1-v_ITO);
E_NiO = 9.58E10/(1-v_NiO);
E_PVSK = 1.11E10/(1-v_PVSK);
E_PCBM = 3.00E9/(1-v_PCBM);
E_Ag_evap = 8.00E10/(1-v_Ag);
E_Ag_screen = 9.00E9/(1-v_Ag);

% Gc for PCBM source
% https://pubs.acs.org/doi/full/10.1021/acsami.5b02202

% Define the fracture toughness for each layer.
Gc_ITO = 36;
Gc_NiO = 5;
Gc_PVSKRSPP = 4.4;
Gc_PVSKhotcast = 0.92;
Gc_PVSKspincoat = 0.37;
Gc_PCBM = 0.13;
Gc_Ag_evap = 5;
Gc_Ag_screen = 2;

% Define curvature, "B" (the in-plane thickness), and radius of curvature
k = 50E-9;
B = 2e-2;
R = 1/k;

% Define delta_T felt by each layer (depends on processing conditions)
T_ITO = 625;
T_ambient = 300;

delta_T_ITO = T_ITO - T_ambient;
delta_T_NiO = 275;

delta_T_PVSK_spin = 75;
delta_T_PVSK_hotcast = 85;
delta_T_PVSK_RSPP = 110;

delta_T_PCBM = 50;
delta_T_Ag_evap = 50;
delta_T_Ag_screen = 125;

% Now, compute stresses induced by thermal expansion constraint of substrate
% on the thin films
sigma_ITO = (a_ITO - a_Polyimide)*delta_T_ITO*E_ITO;
sigma_NiO = (a_NiO - a_Polyimide)*delta_T_NiO*E_NiO;

sigma_PVSK_spin = (a_PVSK - a_Polyimide)*delta_T_PVSK_spin*E_PVSK;
sigma_PVSK_hotcast = (a_PVSK - a_Polyimide)*delta_T_PVSK_hotcast*E_PVSK;
sigma_PVSK_RSPP = (a_PVSK - a_Polyimide)*delta_T_PVSK_RSPP*E_PVSK;

sigma_PCBM = (a_PCBM - a_Polyimide)*delta_T_PCBM*E_PCBM;

sigma_Ag_evap = (a_Ag - a_Polyimide)*delta_T_Ag_evap*E_Ag_evap;
sigma_Ag_screen = (a_Ag - a_Polyimide)*delta_T_Ag_screen*E_Ag_screen;

% Start from the top-most layer and work our way down the stack.


%% Section 1: Ag-PCBM delamination. Ag evap or screen print
%{
% So actually, the Ag layer is undergoing compression on the polyimide
% substrate, and will preferentially BUCKLE instead of delaminate.
% Therefore, the below calculations are NOT VALID. (delamination simply not
% experimentally observed)

syms h_Ag
%Define P symbolically to plug into equation
h_summation = h_Ag;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag^3);

%Solve for Gss
G_AgPCBM = sigma_Ag_evap^2*h_Ag/E_Ag_evap - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12);

% Solve for hc at Gc from literature
hc_AgEvap_AgPCBM = vpasolve(G_AgPCBM == Gc_Ag_evap, h_Ag, [0,inf]);
hc_AgEvap_AgPCBM = double(hc_AgEvap_AgPCBM);

if length(hc_AgEvap_AgPCBM) > 1
    smallest_val = hc_AgEvap_AgPCBM(1);
    for j = 1:length(hc_AgEvap_AgPCBM)
        if hc_AgEvap_AgPCBM(j) < smallest_val
            smallest_val = hc_AgEvap_AgPCBM(j);
        end
    end
    hc_AgEvap_AgPCBM = smallest_val;
end

disp("hc_{Ag} evaporated for Ag-PCBM interface delamination: ")
disp(hc_AgEvap_AgPCBM);


syms h_Ag
%Define P symbolically to plug into equation
h_summation = h_Ag;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3);

%Solve for Gss
G_AgPCBM = sigma_Ag_screen^2*h_Ag/E_Ag_screen - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12);

% Solve for hc at Gc from literature
% Gc for PCBM is used here because Gc_PCBM < Gc_Ag_screen
hc_AgScreen_AgPCBM = vpasolve(G_AgPCBM == Gc_PCBM, h_Ag, [0,inf]);
hc_AgScreen_AgPCBM = double(hc_AgScreen_AgPCBM);

if length(hc_AgScreen_AgPCBM) > 1
    smallest_val = hc_AgScreen_AgPCBM(1);
    for j = 1:length(hc_AgScreen_AgPCBM)
        if hc_AgScreen_AgPCBM(j) < smallest_val
            smallest_val = hc_AgScreen_AgPCBM(j);
        end
    end
    hc_AgScreen_AgPCBM = smallest_val;
end

disp("hc_{Ag} screen for Ag-PCBM interface delamination: ")
disp(hc_AgScreen_AgPCBM);


figure;
x = [hc_AgEvap_AgPCBM*10^9 hc_AgScreen_AgPCBM*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0.5 0.5 0.5];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag evap', 'Ag screen'})
% set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Ag-PCBM delamination");

%}

%% Section 2: PCBM-PVSK delamination
% Section 2.1: PCBM-PVSK delamination, Ag evap. Use Gc_PCBM,varying h_Ag separately
syms h_Ag

% Define the heights of the layers not being varied
h_PCBM=40e-9; % PCBM height is 40 nm by default

% Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag^3 + E_PCBM*h_PCBM^3);

% Solve for Gss
G_PCBMPVSK = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM ...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12);

% Solve for hc at Gc from literature
hc_AgEvap_PCBMPVSK = vpasolve(G_PCBMPVSK == Gc_PCBM, h_Ag, [0,inf]);
hc_AgEvap_PCBMPVSK = double(hc_AgEvap_PCBMPVSK);

if length(hc_AgEvap_PCBMPVSK) > 1
    smallest_val = hc_AgEvap_PCBMPVSK(1);
    for j = 1:length(hc_AgEvap_PCBMPVSK)
        if hc_AgEvap_PCBMPVSK(j) < smallest_val
            smallest_val = hc_AgEvap_PCBMPVSK(j);
        end
    end
    hc_AgEvap_PCBMPVSK = smallest_val;
end

disp("Critical Ag(evap) hc for PCBM-PVSK interface delamination: ")
disp(hc_AgEvap_PCBMPVSK);

% Section 2.2: PCBM-PVSK delamination, Ag screen print 
% Use Gc_PCBM, varying h_Ag
syms h_Ag

% Define the heights of the layers not being varied
h_PCBM=40e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3 + E_PCBM*h_PCBM^3);

% Solve for Gss
G_PCBMPVSK = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM ...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12);

% Solve for hc at Gc from literature
hc_AgScreen_PCBMPVSK = vpasolve(G_PCBMPVSK == Gc_PCBM, h_Ag, [0,inf]);
hc_AgScreen_PCBMPVSK = double(hc_AgScreen_PCBMPVSK);

if length(hc_AgScreen_PCBMPVSK) > 1
    smallest_val = hc_AgScreen_PCBMPVSK(1);
    for j = 1:length(hc_AgScreen_PCBMPVSK)
        if hc_AgScreen_PCBMPVSK(j) < smallest_val
            smallest_val = hc_AgScreen_PCBMPVSK(j);
        end
    end
    hc_AgScreen_PCBMPVSK = smallest_val;
end

disp("Critical Ag(screen) hc for PCBM-PVSK interface delamination: ")
disp(hc_AgScreen_PCBMPVSK);

% Section 2.3: PCBM-PVSK delamination, Ag evap, varying h_PCBM
% Now vary h_PCBM
syms h_PCBM

% Define the heights of the layers not being varied
h_Ag_evap=150e-9; % evaporated silver is 150 nm

% Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3);

% Solve for Gss
G_PCBMPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM ...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12);

% Solve for hc at Gc from literature
hc_PCBM_PCBMPVSK_evap = vpasolve(G_PCBMPVSK == Gc_PCBM, h_PCBM, [0,inf]);
hc_PCBM_PCBMPVSK_evap = double(hc_PCBM_PCBMPVSK_evap);

if length(hc_PCBM_PCBMPVSK_evap) > 1
    smallest_val = hc_PCBM_PCBMPVSK_evap(1);
    for j = 1:length(hc_PCBM_PCBMPVSK_evap)
        if hc_PCBM_PCBMPVSK_evap(j) < smallest_val
            smallest_val = hc_PCBM_PCBMPVSK_evap(j);
        end
    end
    hc_PCBM_PCBMPVSK_evap = smallest_val;
end

disp("PCBM hc for PCBM-PVSK interface delamination (Ag evap): ")
disp(hc_PCBM_PCBMPVSK_evap);


% Section 2.4: PCBM-PVSK delamination, Ag screen, varying h_PCBM
% Now vary h_PCBM
syms h_PCBM

% Define the heights of the layers not being varied
h_Ag_screen=1e-6; % screen silver is 1 micron

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3);

% Solve for Gss
G_PCBMPVSK = sigma_Ag_evap^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM ...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12);

% Solve for hc at Gc from literature
hc_PCBM_PCBMPVSK_screen = vpasolve(G_PCBMPVSK == Gc_PCBM, h_PCBM, [0,inf]);
hc_PCBM_PCBMPVSK_screen = double(hc_PCBM_PCBMPVSK_screen);

if length(hc_PCBM_PCBMPVSK_screen) > 1
    smallest_val = hc_PCBM_PCBMPVSK_screen(1);
    for j = 1:length(hc_PCBM_PCBMPVSK_screen)
        if hc_PCBM_PCBMPVSK_screen(j) < smallest_val
            smallest_val = hc_PCBM_PCBMPVSK_screen(j);
        end
    end
    hc_PCBM_PCBMPVSK_screen = smallest_val;
end

disp("PCBM hc for PCBM-PVSK interface delamination (Ag screen): ")
disp(hc_PCBM_PCBMPVSK_screen);


figure; 
x = [hc_AgEvap_PCBMPVSK*10^9 hc_AgScreen_PCBMPVSK hc_PCBM_PCBMPVSK_evap*10^9 ...
    hc_PCBM_PCBMPVSK_screen*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0.5 0.5 0.5];
b.CData(3,:) = [0 0.4470 0.7410];
b.CData(4,:) = [0 0.4470 0.7410];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag evap', 'Ag screen', 'PCBM (Ag evap)', ...
    'PCBM (Ag screen)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, PCBM-PVSK delamination", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 3.1: NiO-PVSK delamination, varying h_PVSK, h_PCBM, h_Ag. Screen Ag
% Now vary h_Ag
syms h_Ag

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_PCBM=40e-9; 

% Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_NiOPVSK = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_PVSKRSPPNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKRSPP, h_Ag, [0,inf]);
hc_Ag_PVSKRSPPNiO_screen = double(hc_Ag_PVSKRSPPNiO_screen);

if length(hc_Ag_PVSKRSPPNiO_screen) > 1
    smallest_val = hc_Ag_PVSKRSPPNiO_screen(1);
    for j = 1:length(hc_Ag_PVSKRSPPNiO_screen)
        if hc_Ag_PVSKRSPPNiO_screen(j) < smallest_val
            smallest_val = hc_Ag_PVSKRSPPNiO_screen(j);
        end
    end
    hc_Ag_PVSKRSPPNiO_screen = smallest_val;
end

disp("Critical Ag hc for PCBM-PVSK RSPP interface delamination for screen Ag: ")
disp(hc_Ag_PVSKRSPPNiO_screen);


% Solve for hc at Gc from literature for PVSK hot cast
G_NiOPVSK = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

hc_Ag_PVSKhotcastNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKhotcast, h_Ag, [0,inf]);
hc_Ag_PVSKhotcastNiO_screen = double(hc_Ag_PVSKhotcastNiO_screen);

if length(hc_Ag_PVSKhotcastNiO_screen) > 1
    smallest_val = hc_Ag_PVSKhotcastNiO_screen(1);
    for j = 1:length(hc_Ag_PVSKhotcastNiO_screen)
        if hc_Ag_PVSKhotcastNiO_screen(j) < smallest_val
            smallest_val = hc_Ag_PVSKhotcastNiO_screen(j);
        end
    end
    hc_Ag_PVSKhotcastNiO_screen = smallest_val;
end
disp("Critical Ag hc for PCBM-PVSK hot cast interface delamination for screen Ag: ")
disp(hc_Ag_PVSKhotcastNiO_screen);


% Solve for hc at Gc from literature for PVSK spin coat
G_NiOPVSK = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_Ag_PVSKspincoatNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKspincoat, h_Ag, [0,inf]);
hc_Ag_PVSKspincoatNiO_screen = double(hc_Ag_PVSKspincoatNiO_screen);

if length(hc_Ag_PVSKspincoatNiO_screen) > 1
    smallest_val = hc_Ag_PVSKspincoatNiO_screen(1);
    for j = 1:length(hc_Ag_PVSKspincoatNiO_screen)
        if hc_Ag_PVSKspincoatNiO_screen(j) < smallest_val
            smallest_val = hc_Ag_PVSKspincoatNiO_screen(j);
        end
    end
    hc_Ag_PVSKspincoatNiO_screen = smallest_val;
end

disp("Critical Ag hc for PCBM-PVSK spin coat interface delamination for screen Ag: ")
disp(hc_Ag_PVSKspincoatNiO_screen);


% Now vary h_PCBM
syms h_PCBM

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_Ag_screen=1e-6;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_NiOPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_PVSKRSPPNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKRSPP, h_PCBM, [0,inf]);
hc_PCBM_PVSKRSPPNiO_screen = double(hc_PCBM_PVSKRSPPNiO_screen);

if length(hc_PCBM_PVSKRSPPNiO_screen) > 1
    smallest_val = hc_PCBM_PVSKRSPPNiO_screen(1);
    for j = 1:length(hc_PCBM_PVSKRSPPNiO_screen)
        if hc_PCBM_PVSKRSPPNiO_screen(j) < smallest_val
            smallest_val = hc_PCBM_PVSKRSPPNiO_screen(j);
        end
    end
    hc_PCBM_PVSKRSPPNiO_screen = smallest_val;
end

disp("Critical PCBM hc for PCBM-PVSK RSPP interface delamination for screen Ag: ")
disp(hc_PCBM_PVSKRSPPNiO_screen);


% Solve for hc at Gc from literature for PVSK hot cast
G_NiOPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PCBM_PVSKhotcastNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKhotcast, h_PCBM, [0,inf]);
hc_PCBM_PVSKhotcastNiO_screen = double(hc_PCBM_PVSKhotcastNiO_screen);

if length(hc_PCBM_PVSKhotcastNiO_screen) > 1
    smallest_val = hc_PCBM_PVSKhotcastNiO_screen(1);
    for j = 1:length(hc_PCBM_PVSKhotcastNiO_screen)
        if hc_PCBM_PVSKhotcastNiO_screen(j) < smallest_val
            smallest_val = hc_PCBM_PVSKhotcastNiO_screen(j);
        end
    end
    hc_PCBM_PVSKhotcastNiO_screen = smallest_val;
end

disp("Critical PCBM hc for PCBM-PVSK hot cast interface delamination for screen Ag: ")
disp(hc_PCBM_PVSKhotcastNiO_screen);


% Solve for hc at Gc from literature for PVSK spincoat
G_NiOPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PCBM_PVSKspincoatNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKspincoat, h_PCBM, [0,inf]);
hc_PCBM_PVSKspincoatNiO_screen = double(hc_PCBM_PVSKspincoatNiO_screen);

if length(hc_PCBM_PVSKspincoatNiO_screen) > 1
    smallest_val = hc_PCBM_PVSKspincoatNiO_screen(1);
    for j = 1:length(hc_PCBM_PVSKspincoatNiO_screen)
        if hc_PCBM_PVSKspincoatNiO_screen(j) < smallest_val
            smallest_val = hc_PCBM_PVSKspincoatNiO_screen(j);
        end
    end
    hc_PCBM_PVSKspincoatNiO_screen = smallest_val;
end

disp("Critical PCBM hc for PCBM-PVSK spin coat interface delamination for screen Ag: ")
disp(hc_PCBM_PVSKspincoatNiO_screen);


% Now vary h_PVSK
syms h_PVSK

% Define the heights of the layers not being varied
h_PCBM=40e-9; 
h_Ag_screen=1e-6;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_NiOPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKRSPPNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKRSPP, h_PVSK, [0,inf]);
hc_PVSK_PVSKRSPPNiO_screen = double(hc_PVSK_PVSKRSPPNiO_screen);

if length(hc_PVSK_PVSKRSPPNiO_screen) > 1
    smallest_val = hc_PVSK_PVSKRSPPNiO_screen(1);
    for j = 1:length(hc_PVSK_PVSKRSPPNiO_screen)
        if hc_PVSK_PVSKRSPPNiO_screen(j) < smallest_val
            smallest_val = hc_PVSK_PVSKRSPPNiO_screen(j);
        end
    end
    hc_PVSK_PVSKRSPPNiO_screen = smallest_val;
end

disp("Critical PVSK hc for PCBM-PVSK RSPP interface delamination for screen Ag: ")
disp(hc_PVSK_PVSKRSPPNiO_screen);


% Solve for hc at Gc from literature for PVSK hot cast
G_NiOPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PVSK_PVSKhotcastNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKhotcast, h_PVSK, [0,inf]);
hc_PVSK_PVSKhotcastNiO_screen = double(hc_PVSK_PVSKhotcastNiO_screen);

if length(hc_PVSK_PVSKhotcastNiO_screen) > 1
    smallest_val = hc_PVSK_PVSKhotcastNiO_screen(1);
    for j = 1:length(hc_PVSK_PVSKhotcastNiO_screen)
        if hc_PVSK_PVSKhotcastNiO_screen(j) < smallest_val
            smallest_val = hc_PVSK_PVSKhotcastNiO_screen(j);
        end
    end
    hc_PVSK_PVSKhotcastNiO_screen = smallest_val;
end

disp("Critical PVSK hc for PCBM-PVSK hot cast interface delamination: ")
disp(hc_PVSK_PVSKhotcastNiO_screen);


% Solve for hc at Gc from literature for PVSK spin coat
G_NiOPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PVSK_PVSKspincoatNiO_screen = vpasolve(G_NiOPVSK == Gc_PVSKspincoat, h_PVSK, [0,inf]);
hc_PVSK_PVSKspincoatNiO_screen = double(hc_PVSK_PVSKspincoatNiO_screen);

if length(hc_PVSK_PVSKspincoatNiO_screen) > 1
    smallest_val = hc_PVSK_PVSKspincoatNiO_screen(1);
    for j = 1:length(hc_PVSK_PVSKspincoatNiO_screen)
        if hc_PVSK_PVSKspincoatNiO_screen(j) < smallest_val
            smallest_val = hc_PVSK_PVSKspincoatNiO_screen(j);
        end
    end
    hc_PVSK_PVSKspincoatNiO_screen = smallest_val;
end

disp("Critical PVSK hc for PCBM-PVSK spin coat interface delamination: ")
disp(hc_PVSK_PVSKspincoatNiO_screen);

figure;
x = [hc_Ag_PVSKRSPPNiO_screen*10^9 hc_Ag_PVSKhotcastNiO_screen*10^9 hc_Ag_PVSKspincoatNiO_screen*10^9 ... 
    hc_PCBM_PVSKRSPPNiO_screen*10^9 hc_PCBM_PVSKhotcastNiO_screen*10^9 hc_PCBM_PVSKspincoatNiO_screen*10^9 ...
    hc_PVSK_PVSKRSPPNiO_screen*10^9 hc_PVSK_PVSKhotcastNiO_screen*10^9 hc_PVSK_PVSKspincoatNiO_screen*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0.5 0.5 0.5];
b.CData(3,:) = [0.5 0.5 0.5];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0 0.4470 0.7410];
b.CData(6,:) = [0 0.4470 0.7410];
b.CData(7,:) = [0.8500 0.3250 0.0980];
b.CData(8,:) = [0.8500 0.3250 0.0980];
b.CData(9,:) = [0.8500 0.3250 0.0980];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (PVSKRSPP)', 'Ag (PVSKhotcast)', 'Ag (PVSKspincoat)', ...
    'PCBM (PVSKRSPP)', 'PCBM (PVSKhotcast)', 'PCBM (PVSKspincoat)', ...
    'PVSK (PVSKRSPP)', 'PVSK (PVSKhotcast)', 'PVSK (PVSKspincoat)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Screen Ag, PVSK-NiO delamination", "FontWeight", "Normal");
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);



%% Section 3.2: NiO-PVSK delamination, varying h_PVSK, h_PCBM, h_Ag. Evap Ag.
% Now vary h_Ag
syms h_Ag

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_PCBM=40e-9; 

% Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_NiOPVSK = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_PVSKRSPPNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKRSPP, h_Ag, [0,inf]);
hc_Ag_PVSKRSPPNiO_evap = double(hc_Ag_PVSKRSPPNiO_evap);

if length(hc_Ag_PVSKRSPPNiO_evap) > 1
    smallest_val = hc_Ag_PVSKRSPPNiO_evap(1);
    for j = 1:length(hc_Ag_PVSKRSPPNiO_evap)
        if hc_Ag_PVSKRSPPNiO_evap(j) < smallest_val
            smallest_val = hc_Ag_PVSKRSPPNiO_evap(j);
        end
    end
    hc_Ag_PVSKRSPPNiO_evap = smallest_val;
end

disp("Critical Ag hc for PCBM-PVSK RSPP interface delamination for evap Ag: ")
disp(hc_Ag_PVSKRSPPNiO_evap);


% Solve for hc at Gc from literature for PVSK hot cast
G_NiOPVSK = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

hc_Ag_PVSKhotcastNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKhotcast, h_Ag, [0,inf]);
hc_Ag_PVSKhotcastNiO_evap = double(hc_Ag_PVSKhotcastNiO_evap);

if length(hc_Ag_PVSKhotcastNiO_evap) > 1
    smallest_val = hc_Ag_PVSKhotcastNiO_evap(1);
    for j = 1:length(hc_Ag_PVSKhotcastNiO_evap)
        if hc_Ag_PVSKhotcastNiO_evap(j) < smallest_val
            smallest_val = hc_Ag_PVSKhotcastNiO_evap(j);
        end
    end
    hc_Ag_PVSKhotcastNiO_evap = smallest_val;
end
disp("Critical Ag hc for PCBM-PVSK hot cast interface delamination for evap Ag: ")
disp(hc_Ag_PVSKhotcastNiO_evap);


% Solve for hc at Gc from literature for PVSK spin coat
G_NiOPVSK = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_Ag_PVSKspincoatNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKspincoat, h_Ag, [0,inf]);
hc_Ag_PVSKspincoatNiO_evap = double(hc_Ag_PVSKspincoatNiO_evap);

if length(hc_Ag_PVSKspincoatNiO_evap) > 1
    smallest_val = hc_Ag_PVSKspincoatNiO_evap(1);
    for j = 1:length(hc_Ag_PVSKspincoatNiO_evap)
        if hc_Ag_PVSKspincoatNiO_evap(j) < smallest_val
            smallest_val = hc_Ag_PVSKspincoatNiO_evap(j);
        end
    end
    hc_Ag_PVSKspincoatNiO_evap = smallest_val;
end

disp("Critical Ag hc for PCBM-PVSK spin coat interface delamination for evap Ag: ")
disp(hc_Ag_PVSKspincoatNiO_evap);


% Now vary h_PCBM
syms h_PCBM

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_Ag_evap=150e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_NiOPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_PVSKRSPPNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKRSPP, h_PCBM, [0,inf]);
hc_PCBM_PVSKRSPPNiO_evap = double(hc_PCBM_PVSKRSPPNiO_evap);

if length(hc_PCBM_PVSKRSPPNiO_evap) > 1
    smallest_val = hc_PCBM_PVSKRSPPNiO_evap(1);
    for j = 1:length(hc_PCBM_PVSKRSPPNiO_evap)
        if hc_PCBM_PVSKRSPPNiO_evap(j) < smallest_val
            smallest_val = hc_PCBM_PVSKRSPPNiO_evap(j);
        end
    end
    hc_PCBM_PVSKRSPPNiO_evap = smallest_val;
end

disp("Critical PCBM hc for PCBM-PVSK RSPP interface delamination for evap Ag: ")
disp(hc_PCBM_PVSKRSPPNiO_evap);


% Solve for hc at Gc from literature for PVSK hot cast
G_NiOPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PCBM_PVSKhotcastNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKhotcast, h_PCBM, [0,inf]);
hc_PCBM_PVSKhotcastNiO_evap = double(hc_PCBM_PVSKhotcastNiO_evap);

if length(hc_PCBM_PVSKhotcastNiO_evap) > 1
    smallest_val = hc_PCBM_PVSKhotcastNiO_evap(1);
    for j = 1:length(hc_PCBM_PVSKhotcastNiO_evap)
        if hc_PCBM_PVSKhotcastNiO_evap(j) < smallest_val
            smallest_val = hc_PCBM_PVSKhotcastNiO_evap(j);
        end
    end
    hc_PCBM_PVSKhotcastNiO_evap = smallest_val;
end

disp("Critical PCBM hc for PCBM-PVSK hot cast interface delamination for evap Ag: ")
disp(hc_PCBM_PVSKhotcastNiO_evap);


% Solve for hc at Gc from literature for PVSK spincoat
G_NiOPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PCBM_PVSKspincoatNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKspincoat, h_PCBM, [0,inf]);
hc_PCBM_PVSKspincoatNiO_evap = double(hc_PCBM_PVSKspincoatNiO_evap);

if length(hc_PCBM_PVSKspincoatNiO_evap) > 1
    smallest_val = hc_PCBM_PVSKspincoatNiO_evap(1);
    for j = 1:length(hc_PCBM_PVSKspincoatNiO_evap)
        if hc_PCBM_PVSKspincoatNiO_evap(j) < smallest_val
            smallest_val = hc_PCBM_PVSKspincoatNiO_evap(j);
        end
    end
    hc_PCBM_PVSKspincoatNiO_evap = smallest_val;
end

disp("Critical PCBM hc for PCBM-PVSK spin coat interface delamination for evap Ag: ")
disp(hc_PCBM_PVSKspincoatNiO_evap);


% Now vary h_PVSK
syms h_PVSK

% Define the heights of the layers not being varied
h_PCBM=40e-9; 
h_Ag_evap=150e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_NiOPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKRSPPNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKRSPP, h_PVSK, [0,inf]);
hc_PVSK_PVSKRSPPNiO_evap = double(hc_PVSK_PVSKRSPPNiO_evap);

if length(hc_PVSK_PVSKRSPPNiO_evap) > 1
    smallest_val = hc_PVSK_PVSKRSPPNiO_evap(1);
    for j = 1:length(hc_PVSK_PVSKRSPPNiO_evap)
        if hc_PVSK_PVSKRSPPNiO_evap(j) < smallest_val
            smallest_val = hc_PVSK_PVSKRSPPNiO_evap(j);
        end
    end
    hc_PVSK_PVSKRSPPNiO_evap = smallest_val;
end

disp("Critical PVSK hc for PCBM-PVSK RSPP interface delamination for evap Ag: ")
disp(hc_PVSK_PVSKRSPPNiO_evap);


% Solve for hc at Gc from literature for PVSK hot cast
G_NiOPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PVSK_PVSKhotcastNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKhotcast, h_PVSK, [0,inf]);
hc_PVSK_PVSKhotcastNiO_evap = double(hc_PVSK_PVSKhotcastNiO_evap);

if length(hc_PVSK_PVSKhotcastNiO_evap) > 1
    smallest_val = hc_PVSK_PVSKhotcastNiO_evap(1);
    for j = 1:length(hc_PVSK_PVSKhotcastNiO_evap)
        if hc_PVSK_PVSKhotcastNiO_evap(j) < smallest_val
            smallest_val = hc_PVSK_PVSKhotcastNiO_evap(j);
        end
    end
    hc_PVSK_PVSKhotcastNiO_evap = smallest_val;
end

disp("Critical PVSK hc for PCBM-PVSK hot cast interface delamination: ")
disp(hc_PVSK_PVSKhotcastNiO_evap);


% Solve for hc at Gc from literature for PVSK spin coat
G_NiOPVSK = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);
hc_PVSK_PVSKspincoatNiO_evap = vpasolve(G_NiOPVSK == Gc_PVSKspincoat, h_PVSK, [0,inf]);
hc_PVSK_PVSKspincoatNiO_evap = double(hc_PVSK_PVSKspincoatNiO_evap);

if length(hc_PVSK_PVSKspincoatNiO_evap) > 1
    smallest_val = hc_PVSK_PVSKspincoatNiO_evap(1);
    for j = 1:length(hc_PVSK_PVSKspincoatNiO_evap)
        if hc_PVSK_PVSKspincoatNiO_evap(j) < smallest_val
            smallest_val = hc_PVSK_PVSKspincoatNiO_evap(j);
        end
    end
    hc_PVSK_PVSKspincoatNiO_evap = smallest_val;
end

disp("Critical PVSK hc for PCBM-PVSK spin coat interface delamination: ")
disp(hc_PVSK_PVSKspincoatNiO_evap);

figure;
x = [hc_Ag_PVSKRSPPNiO_evap*10^9 hc_Ag_PVSKhotcastNiO_evap*10^9 hc_Ag_PVSKspincoatNiO_evap*10^9 ... 
    hc_PCBM_PVSKRSPPNiO_evap*10^9 hc_PCBM_PVSKhotcastNiO_evap*10^9 hc_PCBM_PVSKspincoatNiO_evap*10^9 ...
    hc_PVSK_PVSKRSPPNiO_evap*10^9 hc_PVSK_PVSKhotcastNiO_evap*10^9 hc_PVSK_PVSKspincoatNiO_evap*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0.5 0.5 0.5];
b.CData(3,:) = [0.5 0.5 0.5];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0 0.4470 0.7410];
b.CData(6,:) = [0 0.4470 0.7410];
b.CData(7,:) = [0.8500 0.3250 0.0980];
b.CData(8,:) = [0.8500 0.3250 0.0980];
b.CData(9,:) = [0.8500 0.3250 0.0980];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (PVSKRSPP)', 'Ag (PVSKhotcast)', 'Ag (PVSKspincoat)', ...
    'PCBM (PVSKRSPP)', 'PCBM (PVSKhotcast)', 'PCBM (PVSKspincoat)', ...
    'PVSK (PVSKRSPP)', 'PVSK (PVSKhotcast)', 'PVSK (PVSKspincoat)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Evap. Ag, PVSK-NiO delamination", "FontWeight", "Normal");
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);



%{
% So actually, the NiO layer is undergoing compression on the polyimide
% substrate, and will preferentially BUCKLE instead of delaminate.
% Therefore, the below calculations are NOT VALID. (delamination simply not
% experimentally observed)

%% Section 4.1: NiO-ITO delamination, varying h_PVSK, h_PCBM, h_Ag, h_NiO. Screen Ag, PVSKRSPP.
% Now vary h_Ag
syms h_Ag

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_PCBM=40e-9;
h_NiO=20e-9; 

% Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

% Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_NiOITO_screen_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_Ag, [0,inf]);
hc_Ag_NiOITO_screen_PVSKRSPP = double(hc_Ag_NiOITO_screen_PVSKRSPP);

if length(hc_Ag_NiOITO_screen_PVSKRSPP) > 1
    smallest_val = hc_Ag_NiOITO_screen_PVSKRSPP(1);
    for j = 1:length(hc_Ag_NiOITO_screen_PVSKRSPP)
        if hc_Ag_NiOITO_screen_PVSKRSPP(j) < smallest_val
            smallest_val = hc_Ag_NiOITO_screen_PVSKRSPP(j);
        end
    end
    hc_Ag_NiOITO_screen_PVSKRSPP = smallest_val;
end

disp("Critical Ag hc for NiO-ITO RSPP interface delamination for screen Ag: ")
disp(hc_Ag_NiOITO_screen_PVSKRSPP);


% Now vary h_PCBM
syms h_PCBM

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_Ag_screen=1e-6;
h_NiO=20e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

% Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_NiOITO_screen_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_PCBM, [0,inf]);
hc_PCBM_NiOITO_screen_PVSKRSPP = double(hc_PCBM_NiOITO_screen_PVSKRSPP);

if length(hc_PCBM_NiOITO_screen_PVSKRSPP) > 1
    smallest_val = hc_PCBM_NiOITO_screen_PVSKRSPP(1);
    for j = 1:length(hc_PCBM_NiOITO_screen_PVSKRSPP)
        if hc_PCBM_NiOITO_screen_PVSKRSPP(j) < smallest_val
            smallest_val = hc_PCBM_NiOITO_screen_PVSKRSPP(j);
        end
    end
    hc_PCBM_NiOITO_screen_PVSKRSPP = smallest_val;
end

disp("Critical PCBM hc for NiO-ITO RSPP interface delamination for screen Ag: ")
disp(hc_PCBM_NiOITO_screen_PVSKRSPP);


% Now vary h_PVSK
syms h_PVSK

% Define the heights of the layers not being varied
h_PCBM=40e-9; 
h_Ag_screen=1e-6;
h_NiO=20e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

% Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_NiOITO_screen_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_PVSK, [0,inf]);
hc_PVSK_NiOITO_screen_PVSKRSPP = double(hc_PVSK_NiOITO_screen_PVSKRSPP);

if length(hc_PVSK_NiOITO_screen_PVSKRSPP) > 1
    smallest_val = hc_PVSK_NiOITO_screen_PVSKRSPP(1);
    for j = 1:length(hc_PVSK_NiOITO_screen_PVSKRSPP)
        if hc_PVSK_NiOITO_screen_PVSKRSPP(j) < smallest_val
            smallest_val = hc_PVSK_NiOITO_screen_PVSKRSPP(j);
        end
    end
    hc_PVSK_NiOITO_screen_PVSKRSPP = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO RSPP interface delamination for screen Ag: ")
disp(hc_PVSK_NiOITO_screen_PVSKRSPP);

% Now vary h_NiO
syms h_NiO

% Define the heights of the layers not being varied
h_PCBM=40e-9; 
h_Ag_screen=1e-6;
h_PVSK=400e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

% Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_NiO_NiOITO_screen_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_NiO, [0,inf]);
hc_NiO_NiOITO_screen_PVSKRSPP = double(hc_NiO_NiOITO_screen_PVSKRSPP);

if length(hc_NiO_NiOITO_screen_PVSKRSPP) > 1
    smallest_val = hc_NiO_NiOITO_screen_PVSKRSPP(1);
    for j = 1:length(hc_NiO_NiOITO_screen_PVSKRSPP)
        if hc_NiO_NiOITO_screen_PVSKRSPP(j) < smallest_val
            smallest_val = hc_NiO_NiOITO_screen_PVSKRSPP(j);
        end
    end
    hc_NiO_NiOITO_screen_PVSKRSPP = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO RSPP interface delamination for screen Ag: ")
disp(hc_NiO_NiOITO_screen_PVSKRSPP);


figure;
x = [hc_Ag_NiOITO_screen_PVSKRSPP*10^9 hc_PCBM_NiOITO_screen_PVSKRSPP*10^9 ...
    hc_PVSK_NiOITO_screen_PVSKRSPP*10^9 hc_NiO_NiOITO_screen_PVSKRSPP*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.6350 0.0780 0.1840];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (NiO-ITO)', 'PCBM (NiO-ITO)', ...
    'PVSK (NiO-ITO)', 'NiO (NiO-ITO)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Screen Ag, NiO-ITO delamination, PVSKRSPP", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 4.2: NiO-ITO delamination, varying h_PVSK, h_PCBM, h_Ag, h_NiO. Screen Ag, PVSKhotcast.
% Now vary h_Ag
syms h_Ag

% Define the heights of the layers not being varied
h_PVSK=400e-9; 
h_PCBM=40e-9;
h_NiO=20e-9; 

%Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_NiOITO_screen_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_Ag, [0,inf]);
hc_Ag_NiOITO_screen_PVSKhotcast = double(hc_Ag_NiOITO_screen_PVSKhotcast);

if length(hc_Ag_NiOITO_screen_PVSKhotcast) > 1
    smallest_val = hc_Ag_NiOITO_screen_PVSKhotcast(1);
    for j = 1:length(hc_Ag_NiOITO_screen_PVSKhotcast)
        if hc_Ag_NiOITO_screen_PVSKhotcast(j) < smallest_val
            smallest_val = hc_Ag_NiOITO_screen_PVSKhotcast(j);
        end
    end
    hc_Ag_NiOITO_screen_PVSKhotcast = smallest_val;
end

disp("Critical Ag hc for NiO-ITO PVSKhotcast interface delamination for screen Ag: ")
disp(hc_Ag_NiOITO_screen_PVSKhotcast);


% Now vary h_PCBM
syms h_PCBM

% Define some variables
h_PVSK=400e-9; 
h_Ag_screen=1e-6;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_NiOITO_screen_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_PCBM, [0,inf]);
hc_PCBM_NiOITO_screen_PVSKhotcast = double(hc_PCBM_NiOITO_screen_PVSKhotcast);

if length(hc_PCBM_NiOITO_screen_PVSKhotcast) > 1
    smallest_val = hc_PCBM_NiOITO_screen_PVSKhotcast(1);
    for j = 1:length(hc_PCBM_NiOITO_screen_PVSKhotcast)
        if hc_PCBM_NiOITO_screen_PVSKhotcast(j) < smallest_val
            smallest_val = hc_PCBM_NiOITO_screen_PVSKhotcast(j);
        end
    end
    hc_PCBM_NiOITO_screen_PVSKhotcast = smallest_val;
end

disp("Critical PCBM hc for NiO-ITO PVSKhotcast interface delamination for screen Ag: ")
disp(hc_PCBM_NiOITO_screen_PVSKhotcast);


% Now vary h_PVSK
syms h_PVSK

% Define some variables
h_PCBM=40e-9; 
h_Ag_screen=1e-6;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_NiOITO_screen_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_PVSK, [0,inf]);
hc_PVSK_NiOITO_screen_PVSKhotcast = double(hc_PVSK_NiOITO_screen_PVSKhotcast);

if length(hc_PVSK_NiOITO_screen_PVSKhotcast) > 1
    smallest_val = hc_PVSK_NiOITO_screen_PVSKhotcast(1);
    for j = 1:length(hc_PVSK_NiOITO_screen_PVSKhotcast)
        if hc_PVSK_NiOITO_screen_PVSKhotcast(j) < smallest_val
            smallest_val = hc_PVSK_NiOITO_screen_PVSKhotcast(j);
        end
    end
    hc_PVSK_NiOITO_screen_PVSKhotcast = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO hotcast interface delamination for screen Ag: ")
disp(hc_PVSK_NiOITO_screen_PVSKhotcast);

% Now vary h_NiO
syms h_NiO

% Define some variables
h_PCBM=40e-9; 
h_Ag_screen=1e-6;
h_PVSK=400e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_NiO_NiOITO_screen_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_NiO, [0,inf]);
hc_NiO_NiOITO_screen_PVSKhotcast = double(hc_NiO_NiOITO_screen_PVSKhotcast);

if length(hc_NiO_NiOITO_screen_PVSKhotcast) > 1
    smallest_val = hc_NiO_NiOITO_screen_PVSKhotcast(1);
    for j = 1:length(hc_NiO_NiOITO_screen_PVSKhotcast)
        if hc_NiO_NiOITO_screen_PVSKhotcast(j) < smallest_val
            smallest_val = hc_NiO_NiOITO_screen_PVSKhotcast(j);
        end
    end
    hc_NiO_NiOITO_screen_PVSKhotcast = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO hotcast interface delamination for screen Ag: ")
disp(hc_NiO_NiOITO_screen_PVSKhotcast);


figure;
x = [hc_Ag_NiOITO_screen_PVSKhotcast*10^9 hc_PCBM_NiOITO_screen_PVSKhotcast*10^9 ...
    hc_PVSK_NiOITO_screen_PVSKhotcast*10^9 hc_NiO_NiOITO_screen_PVSKhotcast*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.6350 0.0780 0.1840];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (NiO-ITO)', 'PCBM (NiO-ITO)', ...
    'PVSK (NiO-ITO)', 'NiO (NiO-ITO)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Screen Ag, NiO-ITO delamination, hotcast")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 4.3: NiO-ITO delamination, varying h_PVSK, h_PCBM, h_Ag, h_NiO. Screen Ag, PVSKspincoat.
% Now vary h_Ag
syms h_Ag

% Define some variables
h_PVSK=400e-9; 
h_PCBM=40e-9;
h_NiO=20e-9; 

%Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_NiOITO_screen_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_Ag, [0,inf]);
hc_Ag_NiOITO_screen_PVSKspincoat = double(hc_Ag_NiOITO_screen_PVSKspincoat);

if length(hc_Ag_NiOITO_screen_PVSKspincoat) > 1
    smallest_val = hc_Ag_NiOITO_screen_PVSKspincoat(1);
    for j = 1:length(hc_Ag_NiOITO_screen_PVSKspincoat)
        if hc_Ag_NiOITO_screen_PVSKspincoat(j) < smallest_val
            smallest_val = hc_Ag_NiOITO_screen_PVSKspincoat(j);
        end
    end
    hc_Ag_NiOITO_screen_PVSKspincoat = smallest_val;
end

disp("Critical Ag hc for NiO-ITO spincoat interface delamination for screen Ag: ")
disp(hc_Ag_NiOITO_screen_PVSKspincoat);


% Now vary h_PCBM
syms h_PCBM

% Define some variables
h_PVSK=400e-9; 
h_Ag_screen=1e-6;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_NiOITO_screen_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_PCBM, [0,inf]);
hc_PCBM_NiOITO_screen_PVSKspincoat = double(hc_PCBM_NiOITO_screen_PVSKspincoat);

if length(hc_PCBM_NiOITO_screen_PVSKspincoat) > 1
    smallest_val = hc_PCBM_NiOITO_screen_PVSKspincoat(1);
    for j = 1:length(hc_PCBM_NiOITO_screen_PVSKspincoat)
        if hc_PCBM_NiOITO_screen_PVSKspincoat(j) < smallest_val
            smallest_val = hc_PCBM_NiOITO_screen_PVSKspincoat(j);
        end
    end
    hc_PCBM_NiOITO_screen_PVSKspincoat = smallest_val;
end

disp("Critical PCBM hc for NiO-ITO spincoat interface delamination for screen Ag: ")
disp(hc_PCBM_NiOITO_screen_PVSKspincoat);


% Now vary h_PVSK
syms h_PVSK

% Define some variables
h_PCBM=40e-9; 
h_Ag_screen=1e-6;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_NiOITO_screen_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_PVSK, [0,inf]);
hc_PVSK_NiOITO_screen_PVSKspincoat = double(hc_PVSK_NiOITO_screen_PVSKspincoat);

if length(hc_PVSK_NiOITO_screen_PVSKspincoat) > 1
    smallest_val = hc_PVSK_NiOITO_screen_PVSKspincoat(1);
    for j = 1:length(hc_PVSK_NiOITO_screen_PVSKspincoat)
        if hc_PVSK_NiOITO_screen_PVSKspincoat(j) < smallest_val
            smallest_val = hc_PVSK_NiOITO_screen_PVSKspincoat(j);
        end
    end
    hc_PVSK_NiOITO_screen_PVSKspincoat = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO spincoat interface delamination for screen Ag: ")
disp(hc_PVSK_NiOITO_screen_PVSKspincoat);

% Now vary h_NiO
syms h_NiO

% Define some variables
h_PCBM=40e-9; 
h_Ag_screen=1e-6;
h_PVSK=400e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_NiO_NiOITO_screen_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_NiO, [0,inf]);
hc_NiO_NiOITO_screen_PVSKspincoat = double(hc_NiO_NiOITO_screen_PVSKspincoat);

if length(hc_NiO_NiOITO_screen_PVSKspincoat) > 1
    smallest_val = hc_NiO_NiOITO_screen_PVSKspincoat(1);
    for j = 1:length(hc_NiO_NiOITO_screen_PVSKspincoat)
        if hc_NiO_NiOITO_screen_PVSKspincoat(j) < smallest_val
            smallest_val = hc_NiO_NiOITO_screen_PVSKspincoat(j);
        end
    end
    hc_NiO_NiOITO_screen_PVSKspincoat = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO spincoat interface delamination for screen Ag: ")
disp(hc_NiO_NiOITO_screen_PVSKspincoat);


figure;
x = [hc_Ag_NiOITO_screen_PVSKspincoat*10^9 hc_PCBM_NiOITO_screen_PVSKspincoat*10^9 ...
    hc_PVSK_NiOITO_screen_PVSKspincoat*10^9 hc_NiO_NiOITO_screen_PVSKspincoat*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.6350 0.0780 0.1840];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (NiO-ITO)', 'PCBM (NiO-ITO)', ...
    'PVSK (NiO-ITO)', 'NiO (NiO-ITO)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Screen Ag, NiO-ITO delamination, spincoat")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 4.4: NiO-ITO delamination, varying h_PVSK, h_PCBM, h_Ag, h_NiO. Evap Ag, PVSKRSPP
% Now vary h_Ag
syms h_Ag

% Define some variables
h_PVSK=400e-9;
h_PCBM=40e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_NiOITO_evap_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_Ag, [0,inf]);
hc_Ag_NiOITO_evap_PVSKRSPP = double(hc_Ag_NiOITO_evap_PVSKRSPP);

if length(hc_Ag_NiOITO_evap_PVSKRSPP) > 1
    smallest_val = hc_Ag_NiOITO_evap_PVSKRSPP(1);
    for j = 1:length(hc_Ag_NiOITO_evap_PVSKRSPP)
        if hc_Ag_NiOITO_evap_PVSKRSPP(j) < smallest_val
            smallest_val = hc_Ag_NiOITO_evap_PVSKRSPP(j);
        end
    end
    hc_Ag_NiOITO_evap_PVSKRSPP = smallest_val;
end

disp("Critical Ag hc for NiO-ITO RSPP interface delamination for evap Ag: ")
disp(hc_Ag_NiOITO_evap_PVSKRSPP);


% Now vary h_PCBM
syms h_PCBM

% Define some variables
h_PVSK=400e-9; 
h_Ag_evap=150e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_NiOITO_evap_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_PCBM, [0,inf]);
hc_PCBM_NiOITO_evap_PVSKRSPP = double(hc_PCBM_NiOITO_evap_PVSKRSPP);

if length(hc_PCBM_NiOITO_evap_PVSKRSPP) > 1
    smallest_val = hc_PCBM_NiOITO_evap_PVSKRSPP(1);
    for j = 1:length(hc_PCBM_NiOITO_evap_PVSKRSPP)
        if hc_PCBM_NiOITO_evap_PVSKRSPP(j) < smallest_val
            smallest_val = hc_PCBM_NiOITO_evap_PVSKRSPP(j);
        end
    end
    hc_PCBM_NiOITO_evap_PVSKRSPP = smallest_val;
end

disp("Critical PCBM hc for NiO-ITO RSPP interface delamination for evap Ag: ")
disp(hc_PCBM_NiOITO_evap_PVSKRSPP);


% Now vary h_PVSK
syms h_PVSK

% Define some variables
h_PCBM=40e-9; 
h_Ag_evap=150e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_NiOITO_evap_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_PVSK, [0,inf]);
hc_PVSK_NiOITO_evap_PVSKRSPP = double(hc_PVSK_NiOITO_evap_PVSKRSPP);

if length(hc_PVSK_NiOITO_evap_PVSKRSPP) > 1
    smallest_val = hc_PVSK_NiOITO_evap_PVSKRSPP(1);
    for j = 1:length(hc_PVSK_NiOITO_evap_PVSKRSPP)
        if hc_PVSK_NiOITO_evap_PVSKRSPP(j) < smallest_val
            smallest_val = hc_PVSK_NiOITO_evap_PVSKRSPP(j);
        end
    end
    hc_PVSK_NiOITO_evap_PVSKRSPP = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO RSPP interface delamination for evap Ag: ")
disp(hc_PVSK_NiOITO_evap_PVSKRSPP);

% Now vary h_NiO
syms h_NiO

% Define some variables
h_PCBM=40e-9; 
h_Ag_evap=150e-9;
h_PVSK=400e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_RSPP^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_NiO_NiOITO_evap_PVSKRSPP = vpasolve(G_NiOITO == Gc_NiO, h_NiO, [0,inf]);
hc_NiO_NiOITO_evap_PVSKRSPP = double(hc_NiO_NiOITO_evap_PVSKRSPP);

if length(hc_NiO_NiOITO_evap_PVSKRSPP) > 1
    smallest_val = hc_NiO_NiOITO_evap_PVSKRSPP(1);
    for j = 1:length(hc_NiO_NiOITO_evap_PVSKRSPP)
        if hc_NiO_NiOITO_evap_PVSKRSPP(j) < smallest_val
            smallest_val = hc_NiO_NiOITO_evap_PVSKRSPP(j);
        end
    end
    hc_NiO_NiOITO_evap_PVSKRSPP = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO RSPP interface delamination for evap Ag: ")
disp(hc_NiO_NiOITO_evap_PVSKRSPP);


% Now, plot all the critical hc, but check that they exist first.
% If no solution, set the value to 0
if isempty(hc_Ag_NiOITO_evap_PVSKRSPP)
    hc_Ag_NiOITO_evap_PVSKRSPP = 0;
end

if isempty(hc_PCBM_NiOITO_evap_PVSKRSPP)
    hc_PCBM_NiOITO_evap_PVSKRSPP = 0;
end

if isempty(hc_PVSK_NiOITO_evap_PVSKRSPP)
    hc_PVSK_NiOITO_evap_PVSKRSPP = 0;
end

if isempty(hc_NiO_NiOITO_evap_PVSKRSPP)
    hc_NiO_NiOITO_evap_PVSKRSPP = 0;
end


figure;
x = [hc_Ag_NiOITO_evap_PVSKRSPP*10^9 hc_PCBM_NiOITO_evap_PVSKRSPP*10^9 hc_PVSK_NiOITO_evap_PVSKRSPP*10^9 ...
    hc_NiO_NiOITO_evap_PVSKRSPP*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.6350 0.0780 0.1840];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (NiO-ITO)', 'PCBM (NiO-ITO)', ...
    'PVSK (NiO-ITO)', 'NiO (NiO-ITO)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Evap. Ag, NiO-ITO delamination, PVSKRSPP")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(x)+1)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 4.5: NiO-ITO delamination, varying h_PVSK, h_PCBM, h_Ag, h_NiO. Evap Ag, PVSKhotcast
% Now vary h_Ag
syms h_Ag

% Define some variables
h_PVSK=400e-9;
h_PCBM=40e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_NiOITO_evap_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_Ag, [0,inf]);
hc_Ag_NiOITO_evap_PVSKhotcast = double(hc_Ag_NiOITO_evap_PVSKhotcast);

if length(hc_Ag_NiOITO_evap_PVSKhotcast) > 1
    smallest_val = hc_Ag_NiOITO_evap_PVSKhotcast(1);
    for j = 1:length(hc_Ag_NiOITO_evap_PVSKhotcast)
        if hc_Ag_NiOITO_evap_PVSKhotcast(j) < smallest_val
            smallest_val = hc_Ag_NiOITO_evap_PVSKhotcast(j);
        end
    end
    hc_Ag_NiOITO_evap_PVSKhotcast = smallest_val;
end

disp("Critical Ag hc for NiO-ITO hotcast interface delamination for evap Ag: ")
disp(hc_Ag_NiOITO_evap_PVSKhotcast);


% Now vary h_PCBM
syms h_PCBM

% Define some variables
h_PVSK=400e-9; 
h_Ag_evap=150e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_NiOITO_evap_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_PCBM, [0,inf]);
hc_PCBM_NiOITO_evap_PVSKhotcast = double(hc_PCBM_NiOITO_evap_PVSKhotcast);

if length(hc_PCBM_NiOITO_evap_PVSKhotcast) > 1
    smallest_val = hc_PCBM_NiOITO_evap_PVSKhotcast(1);
    for j = 1:length(hc_PCBM_NiOITO_evap_PVSKhotcast)
        if hc_PCBM_NiOITO_evap_PVSKhotcast(j) < smallest_val
            smallest_val = hc_PCBM_NiOITO_evap_PVSKhotcast(j);
        end
    end
    hc_PCBM_NiOITO_evap_PVSKhotcast = smallest_val;
end

disp("Critical PCBM hc for NiO-ITO hotcast interface delamination for evap Ag: ")
disp(hc_PCBM_NiOITO_evap_PVSKhotcast);


% Now vary h_PVSK
syms h_PVSK

% Define some variables
h_PCBM=40e-9; 
h_Ag_evap=150e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_NiOITO_evap_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_PVSK, [0,inf]);
hc_PVSK_NiOITO_evap_PVSKhotcast = double(hc_PVSK_NiOITO_evap_PVSKhotcast);

if length(hc_PVSK_NiOITO_evap_PVSKhotcast) > 1
    smallest_val = hc_PVSK_NiOITO_evap_PVSKhotcast(1);
    for j = 1:length(hc_PVSK_NiOITO_evap_PVSKhotcast)
        if hc_PVSK_NiOITO_evap_PVSKhotcast(j) < smallest_val
            smallest_val = hc_PVSK_NiOITO_evap_PVSKhotcast(j);
        end
    end
    hc_PVSK_NiOITO_evap_PVSKhotcast = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO hotcast interface delamination for evap Ag: ")
disp(hc_PVSK_NiOITO_evap_PVSKhotcast);

% Now vary h_NiO
syms h_NiO

% Define some variables
h_PCBM=40e-9; 
h_Ag_evap=150e-9;
h_PVSK=400e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_hotcast^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_NiO_NiOITO_evap_PVSKhotcast = vpasolve(G_NiOITO == Gc_NiO, h_NiO, [0,inf]);
hc_NiO_NiOITO_evap_PVSKhotcast = double(hc_NiO_NiOITO_evap_PVSKhotcast);

if length(hc_NiO_NiOITO_evap_PVSKhotcast) > 1
    smallest_val = hc_NiO_NiOITO_evap_PVSKhotcast(1);
    for j = 1:length(hc_NiO_NiOITO_evap_PVSKhotcast)
        if hc_NiO_NiOITO_evap_PVSKhotcast(j) < smallest_val
            smallest_val = hc_NiO_NiOITO_evap_PVSKhotcast(j);
        end
    end
    hc_NiO_NiOITO_evap_PVSKhotcast = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO hotcast interface delamination for evap Ag: ")
disp(hc_NiO_NiOITO_evap_PVSKhotcast);


% Now, plot all the critical hc, but check that they exist first.
% If no solution, set the value to 0
if isempty(hc_Ag_NiOITO_evap_PVSKhotcast)
    hc_Ag_NiOITO_evap_PVSKhotcast = 0;
end

if isempty(hc_PCBM_NiOITO_evap_PVSKhotcast)
    hc_PCBM_NiOITO_evap_PVSKhotcast = 0;
end

if isempty(hc_PVSK_NiOITO_evap_PVSKhotcast)
    hc_PVSK_NiOITO_evap_PVSKhotcast = 0;
end

if isempty(hc_NiO_NiOITO_evap_PVSKhotcast)
    hc_NiO_NiOITO_evap_PVSKhotcast = 0;
end


figure;
x = [hc_Ag_NiOITO_evap_PVSKhotcast*10^9 hc_PCBM_NiOITO_evap_PVSKhotcast*10^9 hc_PVSK_NiOITO_evap_PVSKhotcast*10^9 ...
    hc_NiO_NiOITO_evap_PVSKhotcast*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.6350 0.0780 0.1840];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (NiO-ITO)', 'PCBM (NiO-ITO)', ...
    'PVSK (NiO-ITO)', 'NiO (NiO-ITO)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Evap. Ag, NiO-ITO delamination, hotcast")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(x)+1)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 4.6: NiO-ITO delamination, varying h_PVSK, h_PCBM, h_Ag, h_NiO. Evap Ag, PVSKspincoat
% Now vary h_Ag
syms h_Ag

% Define some variables
h_PVSK=400e-9;
h_PCBM=40e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag) + E_Ag_evap*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Ag_NiOITO_evap_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_Ag, [0,inf]);
hc_Ag_NiOITO_evap_PVSKspincoat = double(hc_Ag_NiOITO_evap_PVSKspincoat);

if length(hc_Ag_NiOITO_evap_PVSKspincoat) > 1
    smallest_val = hc_Ag_NiOITO_evap_PVSKspincoat(1);
    for j = 1:length(hc_Ag_NiOITO_evap_PVSKspincoat)
        if hc_Ag_NiOITO_evap_PVSKspincoat(j) < smallest_val
            smallest_val = hc_Ag_NiOITO_evap_PVSKspincoat(j);
        end
    end
    hc_Ag_NiOITO_evap_PVSKspincoat = smallest_val;
end

disp("Critical Ag hc for NiO-ITO spincoat interface delamination for evap Ag: ")
disp(hc_Ag_NiOITO_evap_PVSKspincoat);


% Now vary h_PCBM
syms h_PCBM

% Define some variables
h_PVSK=400e-9; 
h_Ag_evap=150e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PCBM_NiOITO_evap_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_PCBM, [0,inf]);
hc_PCBM_NiOITO_evap_PVSKspincoat = double(hc_PCBM_NiOITO_evap_PVSKspincoat);

if length(hc_PCBM_NiOITO_evap_PVSKspincoat) > 1
    smallest_val = hc_PCBM_NiOITO_evap_PVSKspincoat(1);
    for j = 1:length(hc_PCBM_NiOITO_evap_PVSKspincoat)
        if hc_PCBM_NiOITO_evap_PVSKspincoat(j) < smallest_val
            smallest_val = hc_PCBM_NiOITO_evap_PVSKspincoat(j);
        end
    end
    hc_PCBM_NiOITO_evap_PVSKspincoat = smallest_val;
end

disp("Critical PCBM hc for NiO-ITO spincoat interface delamination for evap Ag: ")
disp(hc_PCBM_NiOITO_evap_PVSKspincoat);


% Now vary h_PVSK
syms h_PVSK

% Define some variables
h_PCBM=40e-9; 
h_Ag_evap=150e-9;
h_NiO=20e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_NiOITO_evap_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_PVSK, [0,inf]);
hc_PVSK_NiOITO_evap_PVSKspincoat = double(hc_PVSK_NiOITO_evap_PVSKspincoat);

if length(hc_PVSK_NiOITO_evap_PVSKspincoat) > 1
    smallest_val = hc_PVSK_NiOITO_evap_PVSKspincoat(1);
    for j = 1:length(hc_PVSK_NiOITO_evap_PVSKspincoat)
        if hc_PVSK_NiOITO_evap_PVSKspincoat(j) < smallest_val
            smallest_val = hc_PVSK_NiOITO_evap_PVSKspincoat(j);
        end
    end
    hc_PVSK_NiOITO_evap_PVSKspincoat = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO spincoat interface delamination for evap Ag: ")
disp(hc_PVSK_NiOITO_evap_PVSKspincoat);

% Now vary h_NiO
syms h_NiO

% Define some variables
h_PCBM=40e-9; 
h_Ag_evap=150e-9;
h_PVSK=400e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_evap + h_PCBM + h_PVSK + h_NiO;
P = (k*B/(6*h_summation))*(E_Ag_evap*h_Ag_evap^3 + E_PCBM*h_PCBM^3 + E_PVSK*h_PVSK^3 + E_NiO*h_NiO^3);

%Solve for Gss
G_NiOITO = sigma_Ag_evap^2*h_Ag_evap/E_Ag_evap + sigma_PCBM^2*h_PCBM/E_PCBM + ...
    sigma_PVSK_spin^2*h_PVSK/E_PVSK + sigma_NiO^2*h_NiO/E_NiO...
    - (P^2/(B^2*E_Ag_evap*h_Ag_evap) + E_Ag_evap*h_Ag_evap^3*k^2/12) ...
    - (P^2/(B^2*E_PCBM*h_PCBM) + E_PCBM*h_PCBM^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12) ...
    - (P^2/(B^2*E_NiO*h_NiO) + E_NiO*h_NiO^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_NiO_NiOITO_evap_PVSKspincoat = vpasolve(G_NiOITO == Gc_NiO, h_NiO, [0,inf]);
hc_NiO_NiOITO_evap_PVSKspincoat = double(hc_NiO_NiOITO_evap_PVSKspincoat);

if length(hc_NiO_NiOITO_evap_PVSKspincoat) > 1
    smallest_val = hc_NiO_NiOITO_evap_PVSKspincoat(1);
    for j = 1:length(hc_NiO_NiOITO_evap_PVSKspincoat)
        if hc_NiO_NiOITO_evap_PVSKspincoat(j) < smallest_val
            smallest_val = hc_NiO_NiOITO_evap_PVSKspincoat(j);
        end
    end
    hc_NiO_NiOITO_evap_PVSKspincoat = smallest_val;
end

disp("Critical PVSK hc for NiO-ITO spincoat interface delamination for evap Ag: ")
disp(hc_NiO_NiOITO_evap_PVSKspincoat);


% Now, plot all the critical hc, but check that they exist first.
% If no solution, set the value to 0
if isempty(hc_Ag_NiOITO_evap_PVSKspincoat)
    hc_Ag_NiOITO_evap_PVSKspincoat = 0;
end

if isempty(hc_PCBM_NiOITO_evap_PVSKspincoat)
    hc_PCBM_NiOITO_evap_PVSKspincoat = 0;
end

if isempty(hc_PVSK_NiOITO_evap_PVSKspincoat)
    hc_PVSK_NiOITO_evap_PVSKspincoat = 0;
end

if isempty(hc_NiO_NiOITO_evap_PVSKspincoat)
    hc_NiO_NiOITO_evap_PVSKspincoat = 0;
end


figure;
x = [hc_Ag_NiOITO_evap_PVSKspincoat*10^9 hc_PCBM_NiOITO_evap_PVSKspincoat*10^9 hc_PVSK_NiOITO_evap_PVSKspincoat*10^9 ...
    hc_NiO_NiOITO_evap_PVSKspincoat*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0 0.4470 0.7410];
b.CData(3,:) = [0.8500 0.3250 0.0980];
b.CData(4,:) = [0.6350 0.0780 0.1840];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Ag (NiO-ITO)', 'PCBM (NiO-ITO)', ...
    'PVSK (NiO-ITO)', 'NiO (NiO-ITO)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Evap. Ag, NiO-ITO delamination, spincoat")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(x)+1)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);
%}

%% Find the minimum hc for each layer stack (supplementary figure)

% The value for Ag actually comes from buckling (both screen-Ag and
% evap-Ag). This is because delamination does not occur in Ag

% Find the values for PCBM (both screen-Ag and evap-Ag)
hc_PCBM_screen_PVSKRSPP_values = [hc_PCBM_PCBMPVSK_screen; hc_PCBM_PVSKRSPPNiO_screen];
hc_PCBM_screen_PVSKRSPP = min(hc_PCBM_screen_PVSKRSPP_values);

hc_PCBM_evap_PVSKRSPP_values = [hc_PCBM_PCBMPVSK_evap; hc_PCBM_PVSKRSPPNiO_evap];
hc_PCBM_evap_PVSKRSPP = min(hc_PCBM_evap_PVSKRSPP_values);

hc_PCBM_screen_PVSKspincoat_values = [hc_PCBM_PCBMPVSK_screen; hc_PCBM_PVSKspincoatNiO_screen];
hc_PCBM_screen_PVSKspincoat = min(hc_PCBM_screen_PVSKspincoat_values);

hc_PCBM_evap_PVSKspincoat_values = [hc_PCBM_PCBMPVSK_evap; hc_PCBM_PVSKspincoatNiO_evap];
hc_PCBM_evap_PVSKspincoat = min(hc_PCBM_evap_PVSKspincoat_values);

hc_PCBM_screen_PVSKhotcast_values = [hc_PCBM_PCBMPVSK_screen; hc_PCBM_PVSKhotcastNiO_screen];
hc_PCBM_screen_PVSKhotcast = min(hc_PCBM_screen_PVSKhotcast_values);

hc_PCBM_evap_PVSKhotcast_values = [hc_PCBM_PCBMPVSK_evap; hc_PCBM_PVSKhotcastNiO_evap];
hc_PCBM_evap_PVSKhotcast = min(hc_PCBM_evap_PVSKhotcast_values);

% Find the values for PVSK (all 3 combo, both screen-Ag and evap-Ag)
hc_PVSKRSPP_screen_values = [hc_PVSK_PVSKRSPPNiO_screen];
hc_PVSKRSPP_screen = min(hc_PVSKRSPP_screen_values);

hc_PVSKRSPP_evap_values = [hc_PVSK_PVSKRSPPNiO_evap];
hc_PVSKRSPP_evap = min(hc_PVSKRSPP_evap_values);

hc_PVSKhotcast_screen_values = [hc_PVSK_PVSKhotcastNiO_screen];
hc_PVSKhotcast_screen = min(hc_PVSKhotcast_screen_values);

hc_PVSKhotcast_evap_values = [hc_PVSK_PVSKhotcastNiO_evap];
hc_PVSKhotcast_evap = min(hc_PVSKhotcast_evap_values);

hc_PVSKspincoat_screen_values = [hc_PVSK_PVSKspincoatNiO_screen];
hc_PVSKspincoat_screen = min(hc_PVSKspincoat_screen_values);

hc_PVSKspincoat_evap_values = [hc_PVSK_PVSKspincoatNiO_evap];
hc_PVSKspincoat_evap = min(hc_PVSKspincoat_evap_values);

% We don't need to "find" the minimum value for NiOx because there's only
% one possible value for NiOx (NiOx buckles under PI substrate)

hc_NiO_buckle_polyimide = 1.0510*10^-5;

figure;
x = [hc_PCBM_screen_PVSKRSPP*10^9 hc_PCBM_evap_PVSKRSPP*10^9 hc_PVSKRSPP_screen*10^9 ...
    hc_PVSKRSPP_evap*10^9 hc_NiO_buckle_polyimide*10^9 hc_NiO_buckle_polyimide*10^9; ...
    hc_PCBM_screen_PVSKhotcast*10^9 hc_PCBM_evap_PVSKhotcast*10^9 hc_PVSKhotcast_screen*10^9 ...
    hc_PVSKhotcast_evap*10^9 hc_NiO_buckle_polyimide*10^9 hc_NiO_buckle_polyimide*10^9;...
    hc_PCBM_screen_PVSKspincoat*10^9 hc_PCBM_evap_PVSKspincoat*10^9 hc_PVSKspincoat_screen*10^9 ...
    hc_PVSKspincoat_evap*10^9 hc_NiO_buckle_polyimide*10^9 hc_NiO_buckle_polyimide*10^9;];
b = bar(x, 'FaceColor', 'flat');

for i = 1:size(b,2)
    
    b(1).FaceColor = "#0072BD";
    b(2).FaceColor = "#0000FF";
    b(3).FaceColor = "#EDB120";
    b(4).FaceColor = "#FFFF00";
    b(5).FaceColor = "#A2142F";
    b(6).FaceColor = "#FF0000";

end

%b.CData(1,:) = "#0072BD"; % for Ag
%b.CData(2,:) = "#0072BD"; % for PCBM
% b.CData(3,:) = [0.8500 0.3250 0.0980]; % for PVSK
% b.CData(4,:) = [0.6350 0.0780 0.1840]; % for NiO

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP', 'Hot Cast', 'Spin Coat'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, h_c for each layer")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);
lgd = legend("PCBM (screen)", "PCBM (evap)", "PVSK (screen)", "PVSK (evap)", "NiO (screen)", "NiO (evap)");
lgd.FontSize = 10;



% Make a figure for the paper (grouping by layer, just use evaporated vals)
figure;
x = [hc_PCBM_evap_PVSKRSPP*10^9 hc_PCBM_evap_PVSKhotcast*10^9 hc_PCBM_evap_PVSKspincoat*10^9; ...
    hc_PVSKRSPP_evap*10^9 hc_PVSKhotcast_evap*10^9 hc_PVSKspincoat_evap*10^9;...
    hc_NiO_buckle_polyimide*10^9 hc_NiO_buckle_polyimide*10^9 hc_NiO_buckle_polyimide*10^9];
b = bar(x, 'FaceColor', 'flat');


b(1).CData(1,:) = [0 0.4470 0.7410]; % group 1 1st bar
b(1).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 1st bar
b(1).CData(3,:) = [0.6350 0.0780 0.1840]; % group 3 1st bar

b(2).CData(1,:) = [0 0.4470 0.7410]; % group 1 2nd bar
b(2).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 2nd bar
b(2).CData(3,:) = [0.6350 0.0780 0.1840]; % group 3 2nd bar

b(3).CData(1,:) = [0 0.4470 0.7410]; % group 1 3rd bar
b(3).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 3rd bar
b(3).CData(3,:) = [0.6350 0.0780 0.1840]; % group 3 3rd bar

% b.CData(3,:) = [0.8500 0.3250 0.0980]; % for PVSK
% b.CData(4,:) = [0.6350 0.0780 0.1840]; % for NiO

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'PCBM', 'PVSK', 'NiO'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, h_c for each layer", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Make figures of hc for ITO, NiOx, PVSK, PCBM, and Ag for PIN, evap and screen-printed
%%% THESE FIGURES ARE IN THE PAPER SUPPLEMENTARY

% ITO and NiOx buckles with PI substrate, so we define the hc value for ITO and NiOx buckling
hc_ITOglass_buckle_polyimide = 1.3175*10^-6;
hc_NiO_buckle_polyimide = 1.0510*10^-5;
hc_Ag_screen_buckle_polyimide = 0.0143;
hc_Ag_evap_buckle_polyimide = 0.0252;

figure;
x = [hc_Ag_evap_buckle_polyimide*10^9 hc_PCBM_evap_PVSKRSPP*10^9 hc_PVSKRSPP_evap*10^9 hc_NiO_buckle_polyimide*10^9 hc_ITOglass_buckle_polyimide*10^9; ...
    hc_Ag_evap_buckle_polyimide*10^9 hc_PCBM_evap_PVSKhotcast*10^9 hc_PVSKhotcast_evap*10^9 hc_NiO_buckle_polyimide*10^9 hc_ITOglass_buckle_polyimide*10^9;...
    hc_Ag_evap_buckle_polyimide*10^9 hc_PCBM_evap_PVSKspincoat*10^9 hc_PVSKspincoat_evap*10^9 hc_NiO_buckle_polyimide*10^9 hc_ITOglass_buckle_polyimide*10^9];
b = bar(x, 'FaceColor', 'flat');
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP PVSK', 'Hotcast PVSK', 'Spincoat PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
title("Polyimide substrate, PIN, evap. Ag", 'FontWeight','Normal')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);

legend("Ag", "PCBM", "PVSK", "NiOx", "ITO");

b(1).CData(1,:) = [0.5 0.5 0.5]; % group 1 1st bar
b(1).CData(2,:) = [0.5 0.5 0.5]; % group 2 1st bar
b(1).CData(3,:) = [0.5 0.5 0.5]; % group 3 1st bar

b(2).CData(1,:) = [0 0.4470 0.7410]; % group 1 2nd bar
b(2).CData(2,:) = [0 0.4470 0.7410]; % group 2 2nd bar
b(2).CData(3,:) = [0 0.4470 0.7410]; % group 3 2nd bar

b(3).CData(1,:) = [0.8500 0.3250 0.0980]; % group 1 3rd bar
b(3).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 3rd bar
b(3).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 3rd bar

b(4).CData(1,:) = [0.4660 0.6740 0.1880]; % group 1 4th bar
b(4).CData(2,:) = [0.4660 0.6740 0.1880]; % group 2 4th bar
b(4).CData(3,:) = [0.4660 0.6740 0.1880]; % group 3 4th bar

b(5).CData(1,:) = [0.4940 0.1840 0.5560]; % group 1 5th bar
b(5).CData(2,:) = [0.4940 0.1840 0.5560]; % group 2 5th bar
b(5).CData(3,:) = [0.4940 0.1840 0.5560]; % group 3 5th bar


figure;
x = [hc_Ag_screen_buckle_polyimide*10^9 hc_PCBM_screen_PVSKRSPP*10^9 hc_PVSKRSPP_screen*10^9 hc_NiO_buckle_polyimide*10^9 hc_ITOglass_buckle_polyimide*10^9; ...
    hc_Ag_screen_buckle_polyimide*10^9 hc_PCBM_screen_PVSKhotcast*10^9 hc_PVSKhotcast_screen*10^9 hc_NiO_buckle_polyimide*10^9 hc_ITOglass_buckle_polyimide*10^9;...
    hc_Ag_screen_buckle_polyimide*10^9 hc_PCBM_screen_PVSKspincoat*10^9 hc_PVSKspincoat_screen*10^9 hc_NiO_buckle_polyimide*10^9 hc_ITOglass_buckle_polyimide*10^9];
b = bar(x, 'FaceColor', 'flat');
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP PVSK', 'Hotcast PVSK', 'Spincoat PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
title("Polyimide substrate, PIN, screen-printed Ag", 'FontWeight','Normal')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);

legend("Ag", "PCBM", "PVSK", "NiOx", "ITO");

b(1).CData(1,:) = [0.5 0.5 0.5]; % group 1 1st bar
b(1).CData(2,:) = [0.5 0.5 0.5]; % group 2 1st bar
b(1).CData(3,:) = [0.5 0.5 0.5]; % group 3 1st bar

b(2).CData(1,:) = [0 0.4470 0.7410]; % group 1 2nd bar
b(2).CData(2,:) = [0 0.4470 0.7410]; % group 2 2nd bar
b(2).CData(3,:) = [0 0.4470 0.7410]; % group 3 2nd bar

b(3).CData(1,:) = [0.8500 0.3250 0.0980]; % group 1 3rd bar
b(3).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 3rd bar
b(3).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 3rd bar

b(4).CData(1,:) = [0.4660 0.6740 0.1880]; % group 1 4th bar
b(4).CData(2,:) = [0.4660 0.6740 0.1880]; % group 2 4th bar
b(4).CData(3,:) = [0.4660 0.6740 0.1880]; % group 3 4th bar

b(5).CData(1,:) = [0.4940 0.1840 0.5560]; % group 1 5th bar
b(5).CData(2,:) = [0.4940 0.1840 0.5560]; % group 2 5th bar
b(5).CData(3,:) = [0.4940 0.1840 0.5560]; % group 3 5th bar

