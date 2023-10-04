%-------------Perovskites on Soda Lime Glass-------------%
clear;

%%%%%%%%%%%%%Thomas Colburn, Alan Liu, MSE 358 Project%%%%%%%%%%%%%%
% PVSK on Soda Lime Glass%

% N-I-P architecture for PVSK solar modules.
% The pvsk. layers are as follows:
    % Soda Lime Glass -> ITO -> SnOx nanoparticles -> PVSK. -> Spiro -> Au

% The Au can be evaporated or screen printed in this case.

% BAR CHART COLOR CODES:
% [0.8500 0.3250 0.0980] = Perovskite RGB code
% [0.5 0.5 0.5] = Au RGB code
% [0 0.4470 0.7410] = PCBM RGB code
% [0.6350 0.0780 0.1840] = NiO RGB code
% [0.4940 0.1840 0.5560] = ITO RGB code

% Go from the top of stack to bottom of the stack for delamination calcs.

% Define Poisson's ratio for each layer
v_SodaLimeGlass = 0.23;
v_ITO = 0.35;
v_SnOx = 0.33;
v_PVSK = 0.31;
v_Spiro = 0.36;
v_Au_evap = 0.421;
v_Au_screen = 0.33;

% Define the coefficient of thermal expansion (CTE) for each layer
a_SodaLimeGlass = 9.00E-6;
a_ITO = 5.81E-6;
a_SnOx = 3.75E-6;
a_PVSK = 8.98E-5;
a_Spiro = 5.00E-5;
a_Au_evap = 1.36E-5;
a_Au_screen = 1.90E-5;

% Define the modulus for each layer
E_SodaLimeGlass = 6.90E10/(1-v_SodaLimeGlass);
E_ITO = 1.14E11/(1-v_ITO);
E_SnOx = 8.80E10/(1-v_SnOx);
E_PVSK = 1.11E10/(1-v_PVSK);
E_Spiro = 1.50E10/(1-v_Spiro);
E_Au_evap = 8.20E10/(1-v_Au_evap);
E_Au_screen = 9.00E9/(1-v_Au_screen);

% Define the fracture toughness for each layer.
Gc_ITO = 5;
Gc_SnOx = 0.7;
Gc_PVSKRSPP = 4.4;
Gc_PVSKhotcast = 0.92;
Gc_PVSKspincoat = 0.37;
Gc_Spiro = 0.45;
Gc_Au_evap = 1.2;
Gc_Au_screen = 2;

% Define curvature, "B" (the in-plane thickness), and radius of curvature
k = 50E-9;
B = 2e-2;
R = 1/k;

% Define delta_T felt by each layer (depends on processing conditions)
delta_T_ITO = 325;
delta_T_SnOx = 155;

delta_T_PVSK_spin = 75;
delta_T_PVSK_hotcast = 85;
delta_T_PVSK_RSPP = 110;

delta_T_Spiro = 50;
delta_T_Au_evap = 50;
delta_T_Au_screen = 125;

% Now, compute stresses induced by thermal expansion constraint of substrate
% on the thin films
sigma_ITO = (a_ITO - a_SodaLimeGlass)*delta_T_ITO*E_ITO;
sigma_SnOx = (a_SnOx - a_SodaLimeGlass)*delta_T_SnOx*E_SnOx;

sigma_PVSK_spin = (a_PVSK - a_SodaLimeGlass)*delta_T_PVSK_spin*E_PVSK;
sigma_PVSK_hotcast = (a_PVSK - a_SodaLimeGlass)*delta_T_PVSK_hotcast*E_PVSK;
sigma_PVSK_RSPP = (a_PVSK - a_SodaLimeGlass)*delta_T_PVSK_RSPP*E_PVSK;

sigma_Spiro = (a_Spiro - a_SodaLimeGlass)*delta_T_Spiro*E_Spiro;

sigma_Au_evap = (a_Au_evap - a_SodaLimeGlass)*delta_T_Au_evap*E_Au_evap;
sigma_Au_screen = (a_Au_screen - a_SodaLimeGlass)*delta_T_Au_screen*E_Au_screen;

%% Section 1: Au-Spiro delamination. Au evap or screen print
syms h_Au

% Define P symbolically to plug into equation
h_summation = h_Au;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au^3);

% Solve for Gss
G_AuSpiro = sigma_Au_evap^2*h_Au/E_Au_evap - (P^2/(B^2*E_Au_evap*h_Au) + ...
    E_Au_evap*h_Au^3*k^2/12);

% Solve for hc at Gc from literature
% Gc for evap Au is used here because Gc_Au_evap < Gc_Spiro
hc_AuEvap_AuSpiro = vpasolve(G_AuSpiro == Gc_Au_evap, h_Au, [0,inf]);
hc_AuEvap_AuSpiro = double(hc_AuEvap_AuSpiro);

if length(hc_AuEvap_AuSpiro) > 1
    smallest_val = hc_AuEvap_AuSpiro(1);
    for j = 1:length(hc_AuEvap_AuSpiro)
        if hc_AuEvap_AuSpiro(j) < smallest_val
            smallest_val = hc_AuEvap_AuSpiro(j);
        end
    end
    hc_AuEvap_AuSpiro = smallest_val;
end

disp("hc_{Au} evaporated for Au-Spiro interface delamination: ")
disp(hc_AuEvap_AuSpiro);


syms h_Au
% Define P symbolically to plug into equation
h_summation = h_Au;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au^3);

% Solve for Gss
G_AuSpiro = sigma_Au_screen^2*h_Au/E_Au_screen - (P^2/(B^2*E_Au_screen*h_Au) ...
    + E_Au_screen*h_Au^3*k^2/12);

% Solve for hc at Gc from literature
% Gc for PCBM is used here because Gc_Spiro < Gc_Au_screen
hc_AuScreen_AuSpiro = vpasolve(G_AuSpiro == Gc_Spiro, h_Au, [0,inf]);
hc_AuScreen_AuSpiro = double(hc_AuScreen_AuSpiro);

if length(hc_AuScreen_AuSpiro) > 1
    smallest_val = hc_AuScreen_AuSpiro(1);
    for j = 1:length(hc_AuScreen_AuSpiro)
        if hc_AuScreen_AuSpiro(j) < smallest_val
            smallest_val = hc_AuScreen_AuSpiro(j);
        end
    end
    hc_AuScreen_AuSpiro = smallest_val;
end

disp("hc_{Au} screen for Au-Spiro interface delamination: ")
disp(hc_AuScreen_AuSpiro);


figure;
x = [hc_AuEvap_AuSpiro*10^9 hc_AuScreen_AuSpiro*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [1 1 0];
b.CData(2,:) =[1 1 0];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au evap', 'Au screen'})
% set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Soda Lime Glass, Au-Spiro delamination", "FontWeight", "Normal");


%% Section 2: Spiro-PVSK delamination.

% Section 2.1: Spiro-PVSK delamination, Au evap. Use lowest Gc btwn PVSK
% and Spiro. Varying h_Au separately, with h_Spiro set as default.
syms h_Au

% Define height for the Spiro so that we may vary h_Au separately.
h_Spiro = 180e-9; % in meters

% Define P symbolically to plug into equation.
h_summation = h_Au + h_Spiro;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au^3 + E_Spiro*h_Spiro^3);

% Solve for Gss
G_SpiroPVSK = sigma_Au_evap^2*h_Au/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro ...
    - (P^2/(B^2*E_Au_evap*h_Au) + E_Au_evap*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12);

% Solve for hc when G = Gc from literature
hc_AuEvap_SpiroPVSK = vpasolve(G_SpiroPVSK == Gc_PVSKspincoat);
hc_AuEvap_SpiroPVSK = double(hc_AuEvap_SpiroPVSK);

% Keep the values that are positive, and then check for the minimum.
hc_AuEvap_SpiroPVSK = (hc_AuEvap_SpiroPVSK(hc_AuEvap_SpiroPVSK > 0));

if length(hc_AuEvap_SpiroPVSK) > 1
    smallest_val = hc_AuEvap_SpiroPVSK(1);
    for j = 1:length(hc_AuEvap_SpiroPVSK)
        if hc_AuEvap_SpiroPVSK(j) < smallest_val
            smallest_val = hc_AuEvap_SpiroPVSK(j);
        end
    end
    hc_AuEvap_SpiroPVSK = smallest_val;
end

disp("hc_{Au} evaporated for Spiro-PVSK interface delamination: ")
disp(hc_AuEvap_SpiroPVSK);


% Section 2.2: Spiro-PVSK delamination, Au screen. Use lowest Gc btwn PVSK
% and Spiro. Varying h_Au separately, with h_Spiro set as default.
syms h_Au

% Define height for the Spiro so that we may vary h_Au separately.
h_Spiro = 180e-9; % in meters

% Define P symbolically to plug into equation.
h_summation = h_Au + h_Spiro;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au^3 + E_Spiro*h_Spiro^3);

% Solve for Gss
G_SpiroPVSK = sigma_Au_screen^2*h_Au/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro ...
    - (P^2/(B^2*E_Au_screen*h_Au) + E_Au_screen*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12);

% Solve for hc when G = Gc from literature
hc_AuScreen_SpiroPVSK = vpasolve(G_SpiroPVSK == Gc_PVSKspincoat);
hc_AuScreen_SpiroPVSK = double(hc_AuScreen_SpiroPVSK);

% Keep the values that are positive, and then check for the minimum.
hc_AuScreen_SpiroPVSK = (hc_AuScreen_SpiroPVSK(hc_AuScreen_SpiroPVSK > 0));

if length(hc_AuScreen_SpiroPVSK) > 1
    smallest_val = hc_AuScreen_SpiroPVSK(1);
    for j = 1:length(hc_AuScreen_SpiroPVSK)
        if hc_AuScreen_SpiroPVSK(j) < smallest_val
            smallest_val = hc_AuScreen_SpiroPVSK(j);
        end
    end
    hc_AuScreen_SpiroPVSK = smallest_val;
end

disp("hc_{Au} screen for Spiro-PVSK interface delamination: ")
disp(hc_AuScreen_SpiroPVSK);


% Section 2.3: Spiro-PVSK delamination, Au evap, varying h_Spiro
syms h_Spiro

% Define height for the evap Au so we can vary h_Spiro
h_Au_evap = 150e-9;

% Define P symbolically to plug into equation.
h_summation = h_Au_evap + h_Spiro;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3);

% Solve for Gss
G_SpiroPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro ...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12);

% Solve for hc when G = Gc from literature
hc_Spiro_SpiroPVSK_evap = vpasolve(G_SpiroPVSK == Gc_PVSKspincoat);
hc_Spiro_SpiroPVSK_evap = double(hc_Spiro_SpiroPVSK_evap);

% Keep the values that are positive, and then check for the minimum.
hc_Spiro_SpiroPVSK_evap = (hc_Spiro_SpiroPVSK_evap(hc_Spiro_SpiroPVSK_evap > 0));

if length(hc_Spiro_SpiroPVSK_evap) > 1
    smallest_val = hc_Spiro_SpiroPVSK_evap(1);
    for j = 1:length(hc_Spiro_SpiroPVSK_evap)
        if hc_Spiro_SpiroPVSK_evap(j) < smallest_val
            smallest_val = hc_Spiro_SpiroPVSK_evap(j);
        end
    end
    hc_Spiro_SpiroPVSK_evap = smallest_val;
end

disp("hc_{Spiro} for evap Au Spiro-PVSK interface delamination: ")
disp(hc_Spiro_SpiroPVSK_evap);


% Section 2.4: Spiro-PVSK delamination, Au screen, varying h_Spiro
h_Au_screen = 1000e-9;

% Define P symbolically to plug into equation.
h_summation = h_Au_screen + h_Spiro;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3);

% Solve for Gss
G_SpiroPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro ...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12);

% Solve for hc when G = Gc from literature
hc_Spiro_SpiroPVSK_screen = vpasolve(G_SpiroPVSK == Gc_PVSKspincoat);
hc_Spiro_SpiroPVSK_screen = double(hc_Spiro_SpiroPVSK_screen);

% Keep the values that are positive, and then check for the minimum.
hc_Spiro_SpiroPVSK_screen = (hc_Spiro_SpiroPVSK_screen(hc_Spiro_SpiroPVSK_screen > 0));

if length(hc_Spiro_SpiroPVSK_screen) > 1
    smallest_val = hc_Spiro_SpiroPVSK_screen(1);
    for j = 1:length(hc_Spiro_SpiroPVSK_screen)
        if hc_Spiro_SpiroPVSK_screen(j) < smallest_val
            smallest_val = hc_Spiro_SpiroPVSK_screen(j);
        end
    end
    hc_Spiro_SpiroPVSK_screen = smallest_val;
end

disp("hc_{Spiro} for screen Au Spiro-PVSK interface delamination: ")
disp(hc_Spiro_SpiroPVSK_screen);

figure; 
x = [hc_AuEvap_SpiroPVSK*10^9 hc_AuScreen_SpiroPVSK*10^9 hc_Spiro_SpiroPVSK_evap*10^9 ...
    hc_Spiro_SpiroPVSK_screen*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [1 1 0];
b.CData(2,:) = [1 1 0];
b.CData(3,:) = [0 0.4470 0.7410];
b.CData(4,:) = [0 0.4470 0.7410];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au evap', 'Au screen', 'Spiro (Au evap)', ...
    'Spiro (Au screen)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Soda Lime Glass, Spiro-PVSK delamination", "FontWeight", "Normal");
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 3: SnOx-PVSK delamination. RSPP, hotcast, and spincoat.
% Section 3.1: SnOx-PVSK delamination, varying h_PVSK, h_Spiro, h_Au. screen Au
% Vary h_Au
% with PVSK RSPP
syms h_Au

% Set the heights of the layers not being varied (default values)
h_PVSK = 400e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au) + E_Au_screen*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Au_PVSKRSPPSnOx_screen = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Au, [0,inf]);
hc_Au_PVSKRSPPSnOx_screen = double(hc_Au_PVSKRSPPSnOx_screen);

hc_Au_PVSKRSPPSnOx_screen = (hc_Au_PVSKRSPPSnOx_screen(hc_Au_PVSKRSPPSnOx_screen > 0));

if length(hc_Au_PVSKRSPPSnOx_screen) > 1
    smallest_val = hc_Au_PVSKRSPPSnOx_screen(1);
    for j = 1:length(hc_Au_PVSKRSPPSnOx_screen)
        if hc_Au_PVSKRSPPSnOx_screen(j) < smallest_val
            smallest_val = hc_Au_PVSKRSPPSnOx_screen(j);
        end
    end
    hc_Au_PVSKRSPPSnOx_screen = smallest_val;
end

disp("hc{Au} screen for SnOx-PVSKRSPP interface delamination: ")
disp(hc_Au_PVSKRSPPSnOx_screen);


% Vary h_Au
% with PVSK hotcast
syms h_Au

% Set the "default" heights of PVSK and Spiro, the layers not being varied
h_PVSK = 400e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au) + E_Au_screen*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Au_PVSKhotcastSnOx_screen = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Au, [0,inf]);
hc_Au_PVSKhotcastSnOx_screen = double(hc_Au_PVSKhotcastSnOx_screen);

hc_Au_PVSKhotcastSnOx_screen = (hc_Au_PVSKhotcastSnOx_screen(hc_Au_PVSKhotcastSnOx_screen > 0));

if length(hc_Au_PVSKhotcastSnOx_screen) > 1
    smallest_val = hc_Au_PVSKhotcastSnOx_screen(1);
    for j = 1:length(hc_Au_PVSKhotcastSnOx_screen)
        if hc_Au_PVSKhotcastSnOx_screen(j) < smallest_val
            smallest_val = hc_Au_PVSKhotcastSnOx_screen(j);
        end
    end
    hc_Au_PVSKhotcastSnOx_screen = smallest_val;
end

disp("hc{Au} screen for SnOx-PVSKhotcast interface delamination: ")
disp(hc_Au_PVSKhotcastSnOx_screen);


% Vary h_Au
% with PVSK spincoat
syms h_Au

% Set the "default" heights of PVSK and Spiro
h_PVSK = 400e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au) + E_Au_screen*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Au_PVSKspincoatSnOx_screen = vpasolve(G_SnOxPVSK == Gc_PVSKspincoat, h_Au, [0,inf]);
hc_Au_PVSKspincoatSnOx_screen = double(hc_Au_PVSKspincoatSnOx_screen);

hc_Au_PVSKspincoatSnOx_screen = (hc_Au_PVSKspincoatSnOx_screen(hc_Au_PVSKspincoatSnOx_screen > 0));

if length(hc_Au_PVSKspincoatSnOx_screen) > 1
    smallest_val = hc_Au_PVSKspincoatSnOx_screen(1);
    for j = 1:length(hc_Au_PVSKspincoatSnOx_screen)
        if hc_Au_PVSKspincoatSnOx_screen(j) < smallest_val
            smallest_val = hc_Au_PVSKspincoatSnOx_screen(j);
        end
    end
    hc_Au_PVSKspincoatSnOx_screen = smallest_val;
end

disp("hc{Au} screen for SnOx-PVSKspincoat interface delamination: ")
disp(hc_Au_PVSKspincoatSnOx_screen);



% Now we're going to vary h_Spiro, keeping the others on default values
% PVSK RSPP now
syms h_Spiro

% Set the "default" heights of PVSK and screen-printed Au, the layers not
% being varied in the analysis
h_PVSK = 400e-9;
h_Au_screen = 1000e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Spiro_PVSKRSPPSnOx_screen = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Spiro, [0,inf]);
hc_Spiro_PVSKRSPPSnOx_screen = double(hc_Spiro_PVSKRSPPSnOx_screen);

hc_Spiro_PVSKRSPPSnOx_screen = (hc_Spiro_PVSKRSPPSnOx_screen(hc_Spiro_PVSKRSPPSnOx_screen > 0));

if length(hc_Spiro_PVSKRSPPSnOx_screen) > 1
    smallest_val = hc_Spiro_PVSKRSPPSnOx_screen(1);
    for j = 1:length(hc_Spiro_PVSKRSPPSnOx_screen)
        if hc_Spiro_PVSKRSPPSnOx_screen(j) < smallest_val
            smallest_val = hc_Spiro_PVSKRSPPSnOx_screen(j);
        end
    end
    hc_Spiro_PVSKRSPPSnOx_screen = smallest_val;
end

disp("hc{Spiro} for SnOx-PVSKRSPP interface delamination (screen): ")
disp(hc_Spiro_PVSKRSPPSnOx_screen);

% Now we're going to vary h_Spiro, keeping the others on default values
% PVSK hotcast now
syms h_Spiro

% Set the "default" heights of PVSK and Au screen
h_PVSK = 400e-9;
h_Au_screen = 1000e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Spiro_PVSKhotcastSnOx_screen = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Spiro, [0,inf]);
hc_Spiro_PVSKhotcastSnOx_screen = double(hc_Spiro_PVSKhotcastSnOx_screen);

hc_Spiro_PVSKhotcastSnOx_screen = (hc_Spiro_PVSKhotcastSnOx_screen(hc_Spiro_PVSKhotcastSnOx_screen > 0));

if length(hc_Spiro_PVSKhotcastSnOx_screen) > 1
    smallest_val = hc_Spiro_PVSKhotcastSnOx_screen(1);
    for j = 1:length(hc_Spiro_PVSKhotcastSnOx_screen)
        if hc_Spiro_PVSKhotcastSnOx_screen(j) < smallest_val
            smallest_val = hc_Spiro_PVSKhotcastSnOx_screen(j);
        end
    end
    hc_Spiro_PVSKhotcastSnOx_screen = smallest_val;
end

disp("hc{Spiro} for SnOx-PVSKhotcast interface delamination (screen): ")
disp(hc_Spiro_PVSKhotcastSnOx_screen);



% Now we're going to vary h_Spiro, keeping the others on default values.
% PVSK spincoat now.
syms h_Spiro

% Set the "default" heights of PVSK and screen-printed Au, the layers that
% are not being varied in the analysis.
h_PVSK = 400e-9;
h_Au_screen = 1000e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Spiro_PVSKspincoatSnOx_screen = vpasolve(G_SnOxPVSK == Gc_PVSKspincoat, h_Spiro, [0,inf]);
hc_Spiro_PVSKspincoatSnOx_screen = double(hc_Spiro_PVSKspincoatSnOx_screen);

hc_Spiro_PVSKspincoatSnOx_screen = (hc_Spiro_PVSKspincoatSnOx_screen(hc_Spiro_PVSKspincoatSnOx_screen > 0));

if length(hc_Spiro_PVSKspincoatSnOx_screen) > 1
    smallest_val = hc_Spiro_PVSKspincoatSnOx_screen(1);
    for j = 1:length(hc_Spiro_PVSKspincoatSnOx_screen)
        if hc_Spiro_PVSKspincoatSnOx_screen(j) < smallest_val
            smallest_val = hc_Spiro_PVSKspincoatSnOx_screen(j);
        end
    end
    hc_Spiro_PVSKspincoatSnOx_screen = smallest_val;
end

disp("hc{Spiro} for SnOx-PVSKspincoat interface delamination (screen): ")
disp(hc_Spiro_PVSKspincoatSnOx_screen);



% Vary h_PVSK.
% PVSK RSPP
syms h_PVSK

% Set the "default" values of screen Au and Spiro, the layers that are not
% being varied in the analysis
h_Au_screen = 1000e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKRSPPSnOx_screen = vpasolve(G_SnOxPVSK == Gc_SnOx, h_PVSK, [0,inf]);
hc_PVSK_PVSKRSPPSnOx_screen = double(hc_PVSK_PVSKRSPPSnOx_screen);

hc_PVSK_PVSKRSPPSnOx_screen = (hc_PVSK_PVSKRSPPSnOx_screen(hc_PVSK_PVSKRSPPSnOx_screen > 0));

if length(hc_PVSK_PVSKRSPPSnOx_screen) > 1
    smallest_val = hc_PVSK_PVSKRSPPSnOx_screen(1);
    for j = 1:length(hc_PVSK_PVSKRSPPSnOx_screen)
        if hc_PVSK_PVSKRSPPSnOx_screen(j) < smallest_val
            smallest_val = hc_PVSK_PVSKRSPPSnOx_screen(j);
        end
    end
    hc_PVSK_PVSKRSPPSnOx_screen = smallest_val;
end

disp("hc{PVSKRSPP} for SnOx-PVSKRSPP interface delamination (screen): ")
disp(hc_PVSK_PVSKRSPPSnOx_screen);


% Vary h_PVSK.
% PVSK RSPP
syms h_PVSK

% Set the "default" values of screen Au and Spiro
h_Au_screen = 1000e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKhotcastSnOx_screen = vpasolve(G_SnOxPVSK == Gc_SnOx, h_PVSK, [0,inf]);
hc_PVSK_PVSKhotcastSnOx_screen = double(hc_PVSK_PVSKhotcastSnOx_screen);

hc_PVSK_PVSKhotcastSnOx_screen = (hc_PVSK_PVSKhotcastSnOx_screen(hc_PVSK_PVSKhotcastSnOx_screen > 0));

if length(hc_PVSK_PVSKhotcastSnOx_screen) > 1
    smallest_val = hc_PVSK_PVSKhotcastSnOx_screen(1);
    for j = 1:length(hc_PVSK_PVSKhotcastSnOx_screen)
        if hc_PVSK_PVSKhotcastSnOx_screen(j) < smallest_val
            smallest_val = hc_PVSK_PVSKhotcastSnOx_screen(j);
        end
    end
    hc_PVSK_PVSKhotcastSnOx_screen = smallest_val;
end

disp("hc{PVSKhotcast} for SnOx-PVSKhotcast interface delamination (screen): ")
disp(hc_PVSK_PVSKhotcastSnOx_screen);


% Vary h_PVSK.
% PVSK RSPP
syms h_PVSK

% Set the "default" values of screen Au and Spiro, the layers that are not
% being varied in the analysis
h_Au_screen = 1000e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_screen*h_Au_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_screen^2*h_Au_screen/E_Au_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_screen*h_Au_screen) + E_Au_screen*h_Au_screen^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKspincoatSnOx_screen = vpasolve(G_SnOxPVSK == Gc_PVSKspincoat, h_PVSK, [0,inf]);
hc_PVSK_PVSKspincoatSnOx_screen = double(hc_PVSK_PVSKspincoatSnOx_screen);

hc_PVSK_PVSKspincoatSnOx_screen = (hc_PVSK_PVSKspincoatSnOx_screen(hc_PVSK_PVSKspincoatSnOx_screen > 0));

if length(hc_PVSK_PVSKspincoatSnOx_screen) > 1
    smallest_val = hc_PVSK_PVSKspincoatSnOx_screen(1);
    for j = 1:length(hc_PVSK_PVSKspincoatSnOx_screen)
        if hc_PVSK_PVSKspincoatSnOx_screen(j) < smallest_val
            smallest_val = hc_PVSK_PVSKspincoatSnOx_screen(j);
        end
    end
    hc_PVSK_PVSKspincoatSnOx_screen = smallest_val;
end

disp("hc{PVSKspincoat} for SnOx-PVSKspincoat interface delamination (screen): ")
disp(hc_PVSK_PVSKspincoatSnOx_screen);


% Plot the results for Section 3.1
figure;
x = [hc_Au_PVSKRSPPSnOx_screen*10^9 hc_Au_PVSKhotcastSnOx_screen*10^9 hc_Au_PVSKspincoatSnOx_screen*10^9 ... 
    hc_Spiro_PVSKRSPPSnOx_screen*10^9 hc_Spiro_PVSKhotcastSnOx_screen*10^9 hc_Spiro_PVSKspincoatSnOx_screen*10^9 ...
    hc_PVSK_PVSKRSPPSnOx_screen*10^9 hc_PVSK_PVSKhotcastSnOx_screen*10^9 hc_PVSK_PVSKspincoatSnOx_screen*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [1 1 0];
b.CData(2,:) = [1 1 0];
b.CData(3,:) = [1 1 0];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0 0.4470 0.7410];
b.CData(6,:) = [0 0.4470 0.7410];
b.CData(7,:) = [0.8500 0.3250 0.0980];
b.CData(8,:) = [0.8500 0.3250 0.0980];
b.CData(9,:) = [0.8500 0.3250 0.0980];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au (PVSKRSPP)', 'Au (PVSKhotcast)', 'Au (PVSKspincoat)', ...
    'Spiro (PVSKRSPP)', 'Spiro (PVSKhotcast)', 'Spiro (PVSKspincoat)', ...
    'PVSK (PVSKRSPP)', 'PVSK (PVSKhotcast)', 'PVSK (PVSKspincoat)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Soda Lime Glass, Screen Au, PVSK-SnOx delam.", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


% Section 3.2: SnOx-PVSK delamination, vary h_Au, h_Spiro, h_PVSK. evap Au
% Vary h_Au
% PVSK RSPP
syms h_Au

% Set the "default" heights of PVSK and Spiro
h_PVSK = 400e-9;
h_Spiro = 180e-9;

%Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au) + E_Au_evap*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Au_PVSKRSPPSnOx_evap = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Au, [0,inf]);
hc_Au_PVSKRSPPSnOx_evap = double(hc_Au_PVSKRSPPSnOx_evap);

hc_Au_PVSKRSPPSnOx_evap = (hc_Au_PVSKRSPPSnOx_evap(hc_Au_PVSKRSPPSnOx_evap > 0));

if length(hc_Au_PVSKRSPPSnOx_evap) > 1
    smallest_val = hc_Au_PVSKRSPPSnOx_evap(1);
    for j = 1:length(hc_Au_PVSKRSPPSnOx_evap)
        if hc_Au_PVSKRSPPSnOx_evap(j) < smallest_val
            smallest_val = hc_Au_PVSKRSPPSnOx_evap(j);
        end
    end
    hc_Au_PVSKRSPPSnOx_evap = smallest_val;
end

disp("hc{Au} evap for SnOx-PVSKRSPP interface delamination: ")
disp(hc_Au_PVSKRSPPSnOx_evap);


% Vary h_Au
% PVSK hotcast now.
syms h_Au

% Set the "default" heights of PVSK and Spiro
h_PVSK = 400e-9;
h_Spiro = 180e-9;

%Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au) + E_Au_evap*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Au_PVSKhotcastSnOx_evap = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Au, [0,inf]);
hc_Au_PVSKhotcastSnOx_evap = double(hc_Au_PVSKhotcastSnOx_evap);

hc_Au_PVSKhotcastSnOx_evap = (hc_Au_PVSKhotcastSnOx_evap(hc_Au_PVSKhotcastSnOx_evap > 0));

if length(hc_Au_PVSKhotcastSnOx_evap) > 1
    smallest_val = hc_Au_PVSKhotcastSnOx_evap(1);
    for j = 1:length(hc_Au_PVSKhotcastSnOx_evap)
        if hc_Au_PVSKhotcastSnOx_evap(j) < smallest_val
            smallest_val = hc_Au_PVSKhotcastSnOx_evap(j);
        end
    end
    hc_Au_PVSKhotcastSnOx_evap = smallest_val;
end

disp("hc{Au} evap for SnOx-PVSKhotcast interface delamination: ")
disp(hc_Au_PVSKhotcastSnOx_evap);


% Vary h_Au
% PVSK spincoat now.
syms h_Au

% Set the "default" heights of PVSK and Spiro
h_PVSK = 400e-9;
h_Spiro = 180e-9;

%Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au) + E_Au_evap*h_Au^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Au_PVSKspincoatSnOx_evap = vpasolve(G_SnOxPVSK == Gc_PVSKspincoat, h_Au, [0,inf]);
hc_Au_PVSKspincoatSnOx_evap = double(hc_Au_PVSKspincoatSnOx_evap);

hc_Au_PVSKspincoatSnOx_evap = (hc_Au_PVSKspincoatSnOx_evap(hc_Au_PVSKspincoatSnOx_evap > 0));

if length(hc_Au_PVSKspincoatSnOx_evap) > 1
    smallest_val = hc_Au_PVSKspincoatSnOx_evap(1);
    for j = 1:length(hc_Au_PVSKspincoatSnOx_evap)
        if hc_Au_PVSKspincoatSnOx_evap(j) < smallest_val
            smallest_val = hc_Au_PVSKspincoatSnOx_evap(j);
        end
    end
    hc_Au_PVSKspincoatSnOx_evap = smallest_val;
end

disp("hc{Au} evap for SnOx-PVSKspincoat interface delamination: ")
disp(hc_Au_PVSKspincoatSnOx_evap);



% Now we're going to vary h_Spiro, keeping the others on default values
% PVSK RSPP now
syms h_Spiro

% Set the "default" heights of PVSK and Au screen
h_PVSK = 400e-9;
h_Au_evap = 150e-9;

%Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Spiro_PVSKRSPPSnOx_evap = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Spiro, [0,inf]);
hc_Spiro_PVSKRSPPSnOx_evap = double(hc_Spiro_PVSKRSPPSnOx_evap);

hc_Spiro_PVSKRSPPSnOx_evap = (hc_Spiro_PVSKRSPPSnOx_evap(hc_Spiro_PVSKRSPPSnOx_evap > 0));

if length(hc_Spiro_PVSKRSPPSnOx_evap) > 1
    smallest_val = hc_Spiro_PVSKRSPPSnOx_evap(1);
    for j = 1:length(hc_Spiro_PVSKRSPPSnOx_evap)
        if hc_Spiro_PVSKRSPPSnOx_evap(j) < smallest_val
            smallest_val = hc_Spiro_PVSKRSPPSnOx_evap(j);
        end
    end
    hc_Spiro_PVSKRSPPSnOx_evap = smallest_val;
end

disp("hc{Spiro} for SnOx-PVSKRSPP interface delamination (evap): ")
disp(hc_Spiro_PVSKRSPPSnOx_evap);



% Now we're going to vary h_Spiro, keeping the others on default values
% PVSK hotcast now
syms h_Spiro

% Set the "default" heights of PVSK and Au screen
h_PVSK = 400e-9;
h_Au_evap = 150e-9;

%Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Spiro_PVSKhotcastSnOx_evap = vpasolve(G_SnOxPVSK == Gc_SnOx, h_Spiro, [0,inf]);
hc_Spiro_PVSKhotcastSnOx_evap = double(hc_Spiro_PVSKhotcastSnOx_evap);

hc_Spiro_PVSKhotcastSnOx_evap = (hc_Spiro_PVSKhotcastSnOx_evap(hc_Spiro_PVSKhotcastSnOx_evap > 0));

if length(hc_Spiro_PVSKhotcastSnOx_evap) > 1
    smallest_val = hc_Spiro_PVSKhotcastSnOx_evap(1);
    for j = 1:length(hc_Spiro_PVSKhotcastSnOx_evap)
        if hc_Spiro_PVSKhotcastSnOx_evap(j) < smallest_val
            smallest_val = hc_Spiro_PVSKhotcastSnOx_evap(j);
        end
    end
    hc_Spiro_PVSKhotcastSnOx_evap = smallest_val;
end

disp("hc{Spiro} for SnOx-PVSKhotcast interface delamination (evap): ")
disp(hc_Spiro_PVSKhotcastSnOx_evap);



% Now we're going to vary h_Spiro, keeping the others on default values
% PVSK spincoat now
syms h_Spiro

% Set the "default" heights of PVSK and Au screen
h_PVSK = 400e-9;
h_Au_evap = 150e-9;

%Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_Spiro_PVSKspincoatSnOx_evap = vpasolve(G_SnOxPVSK == Gc_PVSKspincoat, h_Spiro, [0,inf]);
hc_Spiro_PVSKspincoatSnOx_evap = double(hc_Spiro_PVSKspincoatSnOx_evap);

hc_Spiro_PVSKspincoatSnOx_evap = (hc_Spiro_PVSKspincoatSnOx_evap(hc_Spiro_PVSKspincoatSnOx_evap > 0));

if length(hc_Spiro_PVSKspincoatSnOx_evap) > 1
    smallest_val = hc_Spiro_PVSKspincoatSnOx_evap(1);
    for j = 1:length(hc_Spiro_PVSKspincoatSnOx_evap)
        if hc_Spiro_PVSKspincoatSnOx_evap(j) < smallest_val
            smallest_val = hc_Spiro_PVSKspincoatSnOx_evap(j);
        end
    end
    hc_Spiro_PVSKspincoatSnOx_evap = smallest_val;
end

disp("hc{Spiro} for SnOx-PVSKspincoat interface delamination (evap): ")
disp(hc_Spiro_PVSKspincoatSnOx_evap);



% Vary h_PVSK.
% PVSK RSPP
syms h_PVSK

% Set the "default" values of evap Au and Spiro
h_Au_evap = 150e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKRSPPSnOx_evap = vpasolve(G_SnOxPVSK == Gc_SnOx, h_PVSK, [0,inf]);
hc_PVSK_PVSKRSPPSnOx_evap = double(hc_PVSK_PVSKRSPPSnOx_evap);

hc_PVSK_PVSKRSPPSnOx_evap = (hc_PVSK_PVSKRSPPSnOx_evap(hc_PVSK_PVSKRSPPSnOx_evap > 0));

if length(hc_PVSK_PVSKRSPPSnOx_evap) > 1
    smallest_val = hc_PVSK_PVSKRSPPSnOx_evap(1);
    for j = 1:length(hc_PVSK_PVSKRSPPSnOx_evap)
        if hc_PVSK_PVSKRSPPSnOx_evap(j) < smallest_val
            smallest_val = hc_PVSK_PVSKRSPPSnOx_evap(j);
        end
    end
    hc_PVSK_PVSKRSPPSnOx_evap = smallest_val;
end

disp("hc{PVSKRSPP} for SnOx-PVSKRSPP interface delamination (evap): ")
disp(hc_PVSK_PVSKRSPPSnOx_evap);


% Vary h_PVSK.
% PVSK RSPP
syms h_PVSK

% Set the "default" values of screen Au and Spiro
h_Au_evap = 150e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKhotcastSnOx_evap = vpasolve(G_SnOxPVSK == Gc_SnOx, h_PVSK, [0,inf]);
hc_PVSK_PVSKhotcastSnOx_evap = double(hc_PVSK_PVSKhotcastSnOx_evap);

hc_PVSK_PVSKhotcastSnOx_evap = (hc_PVSK_PVSKhotcastSnOx_evap(hc_PVSK_PVSKhotcastSnOx_evap > 0));

if length(hc_PVSK_PVSKhotcastSnOx_evap) > 1
    smallest_val = hc_PVSK_PVSKhotcastSnOx_evap(1);
    for j = 1:length(hc_PVSK_PVSKhotcastSnOx_evap)
        if hc_PVSK_PVSKhotcastSnOx_evap(j) < smallest_val
            smallest_val = hc_PVSK_PVSKhotcastSnOx_evap(j);
        end
    end
    hc_PVSK_PVSKhotcastSnOx_evap = smallest_val;
end

disp("hc{PVSKhotcast} for SnOx-PVSKhotcast interface delamination (evap): ")
disp(hc_PVSK_PVSKhotcastSnOx_evap);


% Vary h_PVSK.
% PVSK RSPP
syms h_PVSK

% Set the "default" values of screen Au and Spiro
h_Au_evap = 150e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Au_evap^2*h_Au_evap/E_Au_evap + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Au_evap*h_Au_evap) + E_Au_evap*h_Au_evap^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12) ...
    - (P^2/(B^2*E_PVSK*h_PVSK) + E_PVSK*h_PVSK^3*k^2/12);

% Solve for hc at Gc from literature for PVSK RSPP
hc_PVSK_PVSKspincoatSnOx_evap = vpasolve(G_SnOxPVSK == Gc_PVSKspincoat, h_PVSK, [0,inf]);
hc_PVSK_PVSKspincoatSnOx_evap = double(hc_PVSK_PVSKspincoatSnOx_evap);

hc_PVSK_PVSKspincoatSnOx_evap = (hc_PVSK_PVSKspincoatSnOx_evap(hc_PVSK_PVSKspincoatSnOx_evap > 0));

if length(hc_PVSK_PVSKspincoatSnOx_evap) > 1
    smallest_val = hc_PVSK_PVSKspincoatSnOx_evap(1);
    for j = 1:length(hc_PVSK_PVSKspincoatSnOx_evap)
        if hc_PVSK_PVSKspincoatSnOx_evap(j) < smallest_val
            smallest_val = hc_PVSK_PVSKspincoatSnOx_evap(j);
        end
    end
    hc_PVSK_PVSKspincoatSnOx_evap = smallest_val;
end

disp("hc{PVSKspincoat} for SnOx-PVSKspincoat interface delamination (evap): ")
disp(hc_PVSK_PVSKspincoatSnOx_evap);



% Plot the results for Section 3.2
figure;
x = [hc_Au_PVSKRSPPSnOx_evap*10^9 hc_Au_PVSKhotcastSnOx_evap*10^9 hc_Au_PVSKspincoatSnOx_evap*10^9 ... 
    hc_Spiro_PVSKRSPPSnOx_evap*10^9 hc_Spiro_PVSKhotcastSnOx_evap*10^9 hc_Spiro_PVSKspincoatSnOx_evap*10^9 ...
    hc_PVSK_PVSKRSPPSnOx_evap*10^9 hc_PVSK_PVSKhotcastSnOx_evap*10^9 hc_PVSK_PVSKspincoatSnOx_evap*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [1 1 0];
b.CData(2,:) = [1 1 0];
b.CData(3,:) = [1 1 0];
b.CData(4,:) = [0 0.4470 0.7410];
b.CData(5,:) = [0 0.4470 0.7410];
b.CData(6,:) = [0 0.4470 0.7410];
b.CData(7,:) = [0.8500 0.3250 0.0980];
b.CData(8,:) = [0.8500 0.3250 0.0980];
b.CData(9,:) = [0.8500 0.3250 0.0980];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au (PVSKRSPP)', 'Au (PVSKhotcast)', 'Au (PVSKspincoat)', ...
    'Spiro (PVSKRSPP)', 'Spiro (PVSKhotcast)', 'Spiro (PVSKspincoat)', ...
    'PVSK (PVSKRSPP)', 'PVSK (PVSKhotcast)', 'PVSK (PVSKspincoat)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Soda Lime Glass, evap. Au, PVSK-SnOx delam.", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Section 4: SnOx-ITO delamination
% Well...SnOx and ITO are in compression actually, so it would
% preferentially buckle, not delaminate.
% So, we don't do calculations for delamination at the SnOx-ITO interface.

% Both the SnOx layer and the ITO layer will buckle.
% (alpha_layer - alpha_substrate) = (-)


%% Find the minimum hc for each layer stack (supplementary figure)

% Find the minimum value for Au (both screen-Au and evap-Au)
hc_Au_screen_PVSKRSPP_values = [hc_AuScreen_AuSpiro; hc_AuScreen_SpiroPVSK; ...
    hc_Au_PVSKRSPPSnOx_screen];
hc_Au_screen_PVSKRSPP = min(hc_Au_screen_PVSKRSPP_values);

hc_Au_evap_PVSKRSPP_values = [hc_AuEvap_AuSpiro; hc_AuEvap_SpiroPVSK; ...
    hc_Au_PVSKRSPPSnOx_evap];
hc_Au_evap_PVSKRSPP = min(hc_Au_evap_PVSKRSPP_values);

hc_Au_screen_PVSKhotcast_values = [hc_AuScreen_AuSpiro; hc_AuScreen_SpiroPVSK;...
    hc_Au_PVSKhotcastSnOx_screen];
hc_Au_screen_PVSKhotcast = min(hc_Au_screen_PVSKhotcast_values);

hc_Au_evap_PVSKhotcast_values = [hc_AuEvap_AuSpiro; hc_AuEvap_SpiroPVSK;...
    hc_Au_PVSKhotcastSnOx_evap];
hc_Au_evap_PVSKhotcast = min(hc_Au_evap_PVSKhotcast_values);

hc_Au_screen_PVSKspincoat_values = [hc_AuScreen_AuSpiro; hc_AuScreen_SpiroPVSK;...
    hc_Au_PVSKspincoatSnOx_screen];
hc_Au_screen_PVSKspincoat = min(hc_Au_screen_PVSKspincoat_values);

hc_Au_evap_PVSKspincoat_values = [hc_AuEvap_AuSpiro; hc_AuEvap_SpiroPVSK;...
    hc_Au_PVSKspincoatSnOx_evap];
hc_Au_evap_PVSKspincoat = min(hc_Au_evap_PVSKspincoat_values);


% Find the minimum value for Spiro (both screen-Au and evap-Au)
hc_Spiro_screen_PVSKRSPP_values = [hc_Spiro_SpiroPVSK_screen; hc_Spiro_PVSKRSPPSnOx_screen];
hc_Spiro_screen_PVSKRSPP = min(hc_Spiro_screen_PVSKRSPP_values);

hc_Spiro_evap_PVSKRSPP_values = [hc_Spiro_SpiroPVSK_evap; hc_Spiro_PVSKRSPPSnOx_evap];
hc_Spiro_evap_PVSKRSPP = min(hc_Spiro_screen_PVSKRSPP_values);

hc_Spiro_screen_PVSKhotcast_values = [hc_Spiro_SpiroPVSK_screen; hc_Spiro_PVSKhotcastSnOx_screen];
hc_Spiro_screen_PVSKhotcast = min(hc_Spiro_screen_PVSKhotcast_values);

hc_Spiro_evap_PVSKhotcast_values = [hc_Spiro_SpiroPVSK_evap; hc_Spiro_PVSKhotcastSnOx_evap];
hc_Spiro_evap_PVSKhotcast = min(hc_Spiro_evap_PVSKhotcast_values);

hc_Spiro_screen_PVSKspincoat_values = [hc_Spiro_SpiroPVSK_screen; hc_Spiro_PVSKspincoatSnOx_screen];
hc_Spiro_screen_PVSKspincoat = min(hc_Spiro_screen_PVSKspincoat_values);

hc_Spiro_evap_PVSKspincoat_values = [hc_Spiro_SpiroPVSK_evap; hc_Spiro_PVSKspincoatSnOx_evap];
hc_Spiro_evap_PVSKspincoat = min(hc_Spiro_evap_PVSKspincoat_values);


% Find the minimum value for PVSK (both screen-Au and evap-Au)
hc_PVSK_screen_PVSKRSPP_values = [hc_PVSK_PVSKRSPPSnOx_screen];
hc_PVSKRSPP_screen = min(hc_PVSK_screen_PVSKRSPP_values);

hc_PVSK_evap_PVSKRSPP_values = [hc_PVSK_PVSKRSPPSnOx_evap];
hc_PVSKRSPP_evap = min(hc_PVSK_evap_PVSKRSPP_values);

hc_PVSK_screen_PVSKhotcast_values = [hc_PVSK_PVSKhotcastSnOx_screen];
hc_PVSKhotcast_screen = min(hc_PVSK_screen_PVSKhotcast_values);

hc_PVSK_evap_PVSKhotcast_values = [hc_PVSK_PVSKhotcastSnOx_evap];
hc_PVSKhotcast_evap = min(hc_PVSK_evap_PVSKhotcast_values);

hc_PVSK_screen_PVSKspincoat_values = [hc_PVSK_PVSKspincoatSnOx_screen];
hc_PVSKspincoat_screen = min(hc_PVSK_screen_PVSKspincoat_values);

hc_PVSK_evap_PVSKspincoat_values = [hc_PVSK_PVSKspincoatSnOx_evap];
hc_PVSKspincoat_evap = min(hc_PVSK_evap_PVSKspincoat_values);

figure;
x = [hc_Au_screen_PVSKRSPP*10^9 hc_Au_evap_PVSKRSPP*10^9 ...
    hc_Spiro_screen_PVSKRSPP*10^9 hc_Spiro_evap_PVSKRSPP*10^9 ...
    hc_PVSKRSPP_screen*10^9 hc_PVSKRSPP_evap*10^9;...
    hc_Au_screen_PVSKhotcast*10^9 hc_Au_evap_PVSKhotcast*10^9 ...
    hc_Spiro_screen_PVSKhotcast*10^9 hc_Spiro_evap_PVSKhotcast*10^9 ...
    hc_PVSKhotcast_screen*10^9 hc_PVSKhotcast_evap*10^9;...
    hc_Au_screen_PVSKspincoat*10^9 hc_Au_evap_PVSKspincoat*10^9 ...
    hc_Spiro_screen_PVSKspincoat*10^9 hc_Spiro_evap_PVSKspincoat*10^9 ...
    hc_PVSKspincoat_screen*10^9 hc_PVSKspincoat_evap*10^9];
b = bar(x, 'FaceColor', 'flat');

for i = 1:size(b,2)

    b(1).FaceColor = "#FFFF00";
    b(2).FaceColor = "#FFFF00";
    b(3).FaceColor = "#0072BD";
    b(4).FaceColor = "#0072BD";
    b(5).FaceColor = "#D95319";
    b(6).FaceColor = "#D95319";

end

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP', 'Hot Cast', 'Spin Coat'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Soda Lime Glass, h_c for each layer", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);
lgd = legend("Au (screen)", "Au (evap)", "Spiro (screen)", "Spiro (evap)", ...
    "PVSK (screen)", "PVSK (evap)");
lgd.FontSize = 10;


%{
%% 
% Make a figure for the paper (grouping by layer, just use evaporated vals)
figure;
x = [hc_Au_evap_PVSKRSPP*10^9 hc_Au_evap_PVSKhotcast*10^9 hc_Au_evap_PVSKspincoat*10^9; ...
    hc_Spiro_evap_PVSKRSPP*10^9 hc_Spiro_evap_PVSKhotcast*10^9 hc_Spiro_evap_PVSKspincoat*10^9; ...
    hc_PVSKRSPP_evap*10^9 hc_PVSKhotcast_evap*10^9 hc_PVSKspincoat_evap*10^9];
b = bar(x, 'FaceColor', 'flat');


b(1).CData(1,:) = [1 1 0]; % group 1 1st bar
b(1).CData(2,:) = [0 0.4470 0.7410]; % group 2 1st bar
b(1).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 1st bar

b(2).CData(1,:) = [1 1 0]; % group 1 2nd bar
b(2).CData(2,:) = [0 0.4470 0.7410]; % group 2 2nd bar
b(2).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 2nd bar

b(3).CData(1,:) = [1 1 0]; % group 1 3rd bar
b(3).CData(2,:) = [0 0.4470 0.7410]; % group 1 3rd bar
b(3).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 3rd bar


xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au', 'Spiro', 'PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Soda Lime Glass, h_c for each layer", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);
%}

%% Make figures of hc for ITO, SnOx, PVSK for NIP, evap & screen Au
%%% THESE FIGURES ARE USED IN THE PAPER SUPPLEMENTARY
% ITO and SnOx buckles under soda lime glass, 
% so we define the hc value for ITO and SnOxbuckling
hc_ITOglass_buckle_sodaLimeGlass = 2.6125*10^-5;
hc_SnOx_buckle_sodaLimeGlass = 1.2076*10^-5;

figure;
x = [hc_PVSKRSPP_evap*10^9 hc_SnOx_buckle_sodaLimeGlass*10^9 hc_ITOglass_buckle_sodaLimeGlass*10^9; ...
    hc_PVSKhotcast_evap*10^9 hc_SnOx_buckle_sodaLimeGlass*10^9 hc_ITOglass_buckle_sodaLimeGlass*10^9;...
    hc_PVSKspincoat_evap*10^9 hc_SnOx_buckle_sodaLimeGlass*10^9 hc_ITOglass_buckle_sodaLimeGlass*10^9];
b = bar(x, 'FaceColor', 'flat');
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP PVSK', 'Hotcast PVSK', 'Spincoat PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
title("Soda lime glass substrate, NIP, evap. Au", 'FontWeight','Normal')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);

legend("PVSK", "SnOx", "ITO");

b(1).CData(1,:) = [0.8500 0.3250 0.0980]; % group 1 1st bar
b(1).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 1st bar
b(1).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 1st bar

b(2).CData(1,:) = [0.9290 0.6940 0.1250]; % group 1 2nd bar
b(2).CData(2,:) = [0.9290 0.6940 0.1250]; % group 2 2nd bar
b(2).CData(3,:) = [0.9290 0.6940 0.1250]; % group 3 2nd bar

b(3).CData(1,:) = [0.4940 0.1840 0.5560]; % group 1 3rd bar
b(3).CData(2,:) = [0.4940 0.1840 0.5560]; % group 2 3rd bar
b(3).CData(3,:) = [0.4940 0.1840 0.5560]; % group 3 3rd bar


figure;
x = [hc_PVSKRSPP_screen*10^9 hc_SnOx_buckle_sodaLimeGlass*10^9 hc_ITOglass_buckle_sodaLimeGlass*10^9; ...
    hc_PVSKhotcast_screen*10^9 hc_SnOx_buckle_sodaLimeGlass*10^9 hc_ITOglass_buckle_sodaLimeGlass*10^9;...
    hc_PVSKspincoat_screen*10^9 hc_SnOx_buckle_sodaLimeGlass*10^9 hc_ITOglass_buckle_sodaLimeGlass*10^9];
b = bar(x, 'FaceColor', 'flat');
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP PVSK', 'Hotcast PVSK', 'Spincoat PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
title("Soda lime glass substrate, NIP, screen-printed Au", 'FontWeight','Normal')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);

legend("PVSK", "SnOx", "ITO");

b(1).CData(1,:) = [0.8500 0.3250 0.0980]; % group 1 1st bar
b(1).CData(2,:) = [0.8500 0.3250 0.0980]; % group 2 1st bar
b(1).CData(3,:) = [0.8500 0.3250 0.0980]; % group 3 1st bar

b(2).CData(1,:) = [0.9290 0.6940 0.1250]; % group 1 2nd bar
b(2).CData(2,:) = [0.9290 0.6940 0.1250]; % group 2 2nd bar
b(2).CData(3,:) = [0.9290 0.6940 0.1250]; % group 3 2nd bar

b(3).CData(1,:) = [0.4940 0.1840 0.5560]; % group 1 3rd bar
b(3).CData(2,:) = [0.4940 0.1840 0.5560]; % group 2 3rd bar
b(3).CData(3,:) = [0.4940 0.1840 0.5560]; % group 3 3rd bar

