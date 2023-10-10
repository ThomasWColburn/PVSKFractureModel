%-------------Perovskites on Polyimide-------------%
clear;

%%%%%%%%%%%%% Thomas Colburn, Alan Liu, PVSK debond model %%%%%%%%%%%%%%
% PVSK on Polyimide, N-I-P structure %

% N-I-P architecture for PVSK solar modules.
% The pvsk. layers are as follows:
    % Polyimide -> ITO -> SnOx nanoparticles -> PVSK. -> Spiro -> Au

% The Au can be evaporated, and we use Ag screen printed because there was
% not a Au screen printed case

% BAR CHART COLOR CODES:
% [0.2 0.094 0.090] = Perovskite RGB code
% [0.980 0.866 0.505] = Au RGB code
% [0.75 0.75 0.75] = Ag-screen RGB code
% [0.192 0.192 0.447] = Spiro RGB code
% [0.694 0.411 0.321] = SnOx RGB code
% [0.533 0.796 0.858] = ITO RGB code

% Go from the top of stack to bottom of the stack for delamination calcs.

% Define Poisson's ratio for each layer
v_Polyimide = 0.35;
v_ITO = 0.35;
v_SnOx = 0.33;
v_PVSK = 0.31;
v_Spiro = 0.36;
v_Au_evap = 0.421;
v_Ag_screen = 0.33;

% Define the coefficient of thermal expansion (CTE) for each layer
a_Polyimide = 2.00E-5;
a_ITO = 5.81E-6;
a_SnOx = 3.75E-6;
a_PVSK = 8.98E-5;
a_Spiro = 5.00E-5;
a_Au_evap = 1.36E-5;
a_Ag_screen = 1.90E-5;

% Define the modulus for each layer
E_Polyimide = 7.50E9/(1-v_Polyimide);
E_ITO = 1.14E11/(1-v_ITO);
E_SnOx = 8.80E10/(1-v_SnOx);
E_PVSK = 1.11E10/(1-v_PVSK);
E_Spiro = 1.50E10/(1-v_Spiro);
E_Au_evap = 8.20E10/(1-v_Au_evap);
E_Ag_screen = 9.00E9/(1-v_Ag_screen);

% Define the fracture toughness for each layer.
Gc_ITO = 5;
Gc_SnOx = 0.7;
Gc_PVSKRSPP = 4.4;
Gc_PVSKhotcast = 0.92;
Gc_PVSKspincoat = 0.37;
Gc_Spiro = 0.45;
Gc_Au_evap = 1.2;
Gc_Ag_screen = 2;

% Define curvature, "B" (the in-plane thickness), and radius of curvature
k = 50E-9;
B = 2e-2;
R = 1/k;

% Define delta_T felt by each layer (default values), depends on the
% processing conditions.
delta_T_ITO = 400;
delta_T_SnOx = 155;

delta_T_PVSK_spin = 75;
delta_T_PVSK_hotcast = 85;
delta_T_PVSK_RSPP = 110;

delta_T_Spiro = 50;
delta_T_Au_evap = 50;
delta_T_Ag_screen = 125;

% Now, define stresses induced by thermal expansion constraint of substrate
% on the thin films.
sigma_ITO = (a_ITO - a_Polyimide)*delta_T_ITO*E_ITO;
sigma_SnOx = (a_SnOx - a_Polyimide)*delta_T_SnOx*E_SnOx;

sigma_PVSK_spin = (a_PVSK - a_Polyimide)*delta_T_PVSK_spin*E_PVSK;
sigma_PVSK_hotcast = (a_PVSK - a_Polyimide)*delta_T_PVSK_hotcast*E_PVSK;
sigma_PVSK_RSPP = (a_PVSK - a_Polyimide)*delta_T_PVSK_RSPP*E_PVSK;

sigma_Spiro = (a_Spiro - a_Polyimide)*delta_T_Spiro*E_Spiro;

sigma_Au_evap = (a_Au_evap - a_Polyimide)*delta_T_Au_evap*E_Au_evap;
sigma_Ag_screen = (a_Ag_screen - a_Polyimide)*delta_T_Ag_screen*E_Ag_screen;

%% Section 1: Au-Spiro delamination. Au evap or screen print
% Since Au buckles under the polyimide substrate, we put the buckling
% values for screen-printed Au and evaporated Au here

% These values are calculated from the script
% "Channel_buckle_fracture_NIP.m"
hc_AuEvapBucklePI = 1.1656*10^-4;
hc_AgScreenBucklePI = 0.0143;

%% Section 2: Spiro-PVSK delamination. Varying h_Spiro and h_Au

% We have to keep in mind that the gold might be thick enough to pull the
% Spiro along with it. It's an issue for screen-printed Au which might be
% upwards of 1 micron thick.

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


% Section 2.2: Spiro-PVSK delamination, Ag screen. Use lowest Gc btwn PVSK
% and Spiro. Varying h_Ag separately, with h_Spiro set as default.
syms h_Ag

% Define height for the Spiro so that we may vary h_Ag separately.
h_Spiro = 180e-9; % in meters

% Define P symbolically to plug into equation.
h_summation = h_Ag + h_Spiro;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag^3 + E_Spiro*h_Spiro^3);

% Solve for Gss
G_SpiroPVSK = sigma_Ag_screen^2*h_Ag/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro ...
    - (P^2/(B^2*E_Ag_screen*h_Ag) + E_Ag_screen*h_Ag^3*k^2/12) ...
    - (P^2/(B^2*E_Spiro*h_Spiro) + E_Spiro*h_Spiro^3*k^2/12);

% Solve for hc when G = Gc from literature
hc_AgScreen_SpiroPVSK = vpasolve(G_SpiroPVSK == Gc_PVSKspincoat);
hc_AgScreen_SpiroPVSK = double(hc_AgScreen_SpiroPVSK);

% Keep the values that are positive, and then check for the minimum.
hc_AgScreen_SpiroPVSK = (hc_AgScreen_SpiroPVSK(hc_AgScreen_SpiroPVSK > 0));

if length(hc_AgScreen_SpiroPVSK) > 1
    smallest_val = hc_AgScreen_SpiroPVSK(1);
    for j = 1:length(hc_AgScreen_SpiroPVSK)
        if hc_AgScreen_SpiroPVSK(j) < smallest_val
            smallest_val = hc_AgScreen_SpiroPVSK(j);
        end
    end
    hc_AgScreen_SpiroPVSK = smallest_val;
end

disp("hc_{Ag} screen for Spiro-PVSK interface delamination: ")
disp(hc_AgScreen_SpiroPVSK);


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


% Section 2.4: Spiro-PVSK delamination, Ag screen, varying h_Spiro
% Define the height of the screen-printed Ag (layer not being varied)
h_Ag_screen = 1000e-9;

% Define P symbolically to plug into equation.
h_summation = h_Ag_screen + h_Spiro;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3);

% Solve for Gss
G_SpiroPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro ...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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
x = [hc_AuEvap_SpiroPVSK*10^9 hc_AgScreen_SpiroPVSK*10^9 hc_Spiro_SpiroPVSK_evap*10^9 ...
    hc_Spiro_SpiroPVSK_screen*10^9];
b = bar(x, 'FaceColor', 'flat');
b.CData(1,:) = [0.980 0.866 0.505];
b.CData(2,:) = [0.75 0.75 0.75];
b.CData(3,:) = [0.192 0.192 0.447];
b.CData(4,:) = [0.192 0.192 0.447];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au evap', 'Ag screen', 'Spiro (Au evap)', ...
    'Spiro (Ag screen)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Spiro-PVSK delamination", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);



%% Section 3: SnOx-PVSK delamination, 
% Section 3.1: SnOx-PVSK delamination, varying h_PVSK, h_Spiro, h_Au. screen Au
% Vary h_Au
% with PVSK RSPP
syms h_Au

% Set the "default" heights of PVSK and Spiro
h_PVSK = 400e-9;
h_Spiro = 180e-9;

%Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Au/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Au) + E_Ag_screen*h_Au^3*k^2/12) ...
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

% Set the "default" heights of PVSK and Spiro
h_PVSK = 400e-9;
h_Spiro = 180e-9;

%Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Au/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Au) + E_Ag_screen*h_Au^3*k^2/12) ...
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

%Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Au/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Au) + E_Ag_screen*h_Au^3*k^2/12) ...
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

% Set the "default" heights of PVSK and Ag screen
h_PVSK = 400e-9;
h_Ag_screen = 1000e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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

% Set the "default" heights of PVSK and Ag screen
h_PVSK = 400e-9;
h_Ag_screen = 1000e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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

% Set the "default" heights of PVSK and Ag screen
h_PVSK = 400e-9;
h_Ag_screen = 1000e-9;

%Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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

% Set the "default" values of screen Au and Spiro
h_Ag_screen = 1000e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_RSPP^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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
h_Ag_screen = 1000e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_hotcast^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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

% Set the "default" values of screen Au and Spiro
h_Ag_screen = 1000e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Ag_screen + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Ag_screen*h_Ag_screen^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

%Solve for Gss
G_SnOxPVSK = sigma_Ag_screen^2*h_Ag_screen/E_Ag_screen + sigma_Spiro^2*h_Spiro/E_Spiro...
    + sigma_PVSK_spin^2*h_PVSK/E_PVSK...
    - (P^2/(B^2*E_Ag_screen*h_Ag_screen) + E_Ag_screen*h_Ag_screen^3*k^2/12) ...
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
b.CData(1,:) = [0.980 0.866 0.505];
b.CData(2,:) = [0.980 0.866 0.505];
b.CData(3,:) = [0.980 0.866 0.505];
b.CData(4,:) = [0.192 0.192 0.447];
b.CData(5,:) = [0.192 0.192 0.447];
b.CData(6,:) = [0.192 0.192 0.447];
b.CData(7,:) = [0.2 0.094 0.090];
b.CData(8,:) = [0.2 0.094 0.090];
b.CData(9,:) = [0.2 0.094 0.090];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au (PVSKRSPP)', 'Au (PVSKhotcast)', 'Au (PVSKspincoat)', ...
    'Spiro (PVSKRSPP)', 'Spiro (PVSKhotcast)', 'Spiro (PVSKspincoat)', ...
    'PVSK (PVSKRSPP)', 'PVSK (PVSKhotcast)', 'PVSK (PVSKspincoat)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, Screen Au, PVSK-SnOx delam.")
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

% Set the "default" heights of PVSK and Spiro, the layers that are not
% being varied
h_PVSK = 400e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
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

% Set the "default" heights of PVSK and Spiro, the heights of the layers
% that are not being varied
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

% Set the "default" heights of PVSK and Au evap, the heights of the
% layers that are not being varied
h_PVSK = 400e-9;
h_Au_evap = 150e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
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

% Set the "default" heights of PVSK and Au evap, the heights of the
% layers not being varied
h_PVSK = 400e-9;
h_Au_evap = 150e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
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

% Set the "default" heights of PVSK and Au evap, the heights of the
% layers not being varied.
h_PVSK = 400e-9;
h_Au_evap = 150e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
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

% Set the "default" values of evap Au and Spiro, the heights of the layers
% not being valued
h_Au_evap = 150e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
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

% Set the "default" values of screen Au and Spiro, the heights of the
% layers not being varied
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

% Set the "default" values of screen Au and Spiro, the heights of the
% layers not being varied
h_Au_evap = 150e-9;
h_Spiro = 180e-9;

% Define P symbolically to plug into equation
h_summation = h_Au_evap + h_Spiro + h_PVSK;
P = (k*B/(6*h_summation))*(E_Au_evap*h_Au_evap^3 + E_Spiro*h_Spiro^3 + E_PVSK*h_PVSK^3);

% Solve for Gss
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
b.CData(1,:) = [0.980 0.866 0.505];
b.CData(2,:) = [0.980 0.866 0.505];
b.CData(3,:) = [0.980 0.866 0.505];
b.CData(4,:) = [0.192 0.192 0.447];
b.CData(5,:) = [0.192 0.192 0.447];
b.CData(6,:) = [0.192 0.192 0.447];
b.CData(7,:) = [0.2 0.094 0.090];
b.CData(8,:) = [0.2 0.094 0.090];
b.CData(9,:) = [0.2 0.094 0.090];
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au (PVSKRSPP)', 'Au (PVSKhotcast)', 'Au (PVSKspincoat)', ...
    'Spiro (PVSKRSPP)', 'Spiro (PVSKhotcast)', 'Spiro (PVSKspincoat)', ...
    'PVSK (PVSKRSPP)', 'PVSK (PVSKhotcast)', 'PVSK (PVSKspincoat)'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, evap. Au, PVSK-SnOx delam.", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = min(x)/10;
y_upperlimit = max(x)*10;
ylim([y_lowerlimit, y_upperlimit]);



%% Section 4: SnOx-ITO delamination
% Well...SnOx and ITO are in compression actually, so it would
% preferentially buckle, not delaminate.
% So, we don't do calculations for delamination at the SnOx-ITO interface.


%% Find the minimum hc for each layer stack (supplementary figure)

% Find the minimum value for the electrode (both screen-Ag and evap-Au)
hc_Ag_screen_PVSKRSPP_values = [hc_AgScreenBucklePI; hc_AgScreen_SpiroPVSK; ...
    hc_Au_PVSKRSPPSnOx_screen];
hc_Ag_screen_PVSKRSPP = min(hc_Ag_screen_PVSKRSPP_values);

hc_Au_evap_PVSKRSPP_values = [hc_AuEvapBucklePI; hc_AuEvap_SpiroPVSK; ...
    hc_Au_PVSKRSPPSnOx_evap];
hc_Au_evap_PVSKRSPP = min(hc_Au_evap_PVSKRSPP_values);

hc_Ag_screen_PVSKhotcast_values = [hc_AgScreenBucklePI; hc_AgScreen_SpiroPVSK;...
    hc_Au_PVSKhotcastSnOx_screen];
hc_Ag_screen_PVSKhotcast = min(hc_Ag_screen_PVSKhotcast_values);

hc_Au_evap_PVSKhotcast_values = [hc_AuEvapBucklePI; hc_AuEvap_SpiroPVSK;...
    hc_Au_PVSKhotcastSnOx_evap];
hc_Au_evap_PVSKhotcast = min(hc_Au_evap_PVSKhotcast_values);

hc_Ag_screen_PVSKspincoat_values = [hc_AgScreenBucklePI; hc_AgScreen_SpiroPVSK;...
    hc_Au_PVSKspincoatSnOx_screen];
hc_Ag_screen_PVSKspincoat = min(hc_Ag_screen_PVSKspincoat_values);

hc_Au_evap_PVSKspincoat_values = [hc_AuEvapBucklePI; hc_AuEvap_SpiroPVSK;...
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
x = [hc_Ag_screen_PVSKRSPP*10^9 hc_Au_evap_PVSKRSPP*10^9 ...
    hc_Spiro_screen_PVSKRSPP*10^9 hc_Spiro_evap_PVSKRSPP*10^9 ...
    hc_PVSKRSPP_screen*10^9 hc_PVSKRSPP_evap*10^9;...
    hc_Ag_screen_PVSKhotcast*10^9 hc_Au_evap_PVSKhotcast*10^9 ...
    hc_Spiro_screen_PVSKhotcast*10^9 hc_Spiro_evap_PVSKhotcast*10^9 ...
    hc_PVSKhotcast_screen*10^9 hc_PVSKhotcast_evap*10^9;...
    hc_Ag_screen_PVSKspincoat*10^9 hc_Au_evap_PVSKspincoat*10^9 ...
    hc_Spiro_screen_PVSKspincoat*10^9 hc_Spiro_evap_PVSKspincoat*10^9 ...
    hc_PVSKspincoat_screen*10^9 hc_PVSKspincoat_evap*10^9];
b = bar(x, 'FaceColor', 'flat');

for i = 1:size(b,2)

    b(1).FaceColor = "#0072BD";
    b(2).FaceColor = "#0000FF";
    b(3).FaceColor = "#EDB120";
    b(4).FaceColor = "#FFFF00";
    b(5).FaceColor = "#A2142F";
    b(6).FaceColor = "#FF0000";

end

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP', 'Hot Cast', 'Spin Coat'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, h_c for each layer", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);
lgd = legend("Au (screen)", "Au (evap)", "Spiro (screen)", "Spiro (evap)", ...
    "PVSK (screen)", "PVSK (evap)");
lgd.FontSize = 10;



%% 
% Make a figure for the paper (grouping by layer, just use screen-prnt vals)
figure;
x = [hc_Ag_screen_PVSKRSPP*10^9 hc_Ag_screen_PVSKhotcast*10^9 hc_Ag_screen_PVSKspincoat*10^9; ...
    hc_Spiro_screen_PVSKRSPP*10^9 hc_Spiro_screen_PVSKhotcast*10^9 hc_Spiro_screen_PVSKspincoat*10^9; ...
    hc_PVSKRSPP_screen*10^9 hc_PVSKhotcast_screen*10^9 hc_PVSKspincoat_screen*10^9];
b = bar(x, 'FaceColor', 'flat');

b(1).CData(1,:) = [0.980 0.866 0.505]; % group 1 1st bar
b(1).CData(2,:) = [0.192 0.192 0.447]; % group 2 1st bar
b(1).CData(3,:) = [0.2 0.094 0.090]; % group 3 1st bar

b(2).CData(1,:) = [0.980 0.866 0.505]; % group 1 2nd bar
b(2).CData(2,:) = [0.192 0.192 0.447]; % group 2 2nd bar
b(2).CData(3,:) = [0.2 0.094 0.090]; % group 3 2nd bar

b(3).CData(1,:) = [0.980 0.866 0.505]; % group 1 3rd bar
b(3).CData(2,:) = [0.192 0.192 0.447]; % group 2 3rd bar
b(3).CData(3,:) = [0.2 0.094 0.090]; % group 3 3rd bar

xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'Au', 'Spiro', 'PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
title("Polyimide, h_c for each layer", "FontWeight", "Normal")
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);


%% Make figures of hc for ITO, SnOx, PVSK for NIP, evap & screen Au
%%% THESE FIGURES ARE USED IN THE PAPER SUPPLEMENTARY
% ITO and SnOx buckles under polyimide, 
% so we define the hc value for ITO and SnOxbuckling
hc_ITOglass_buckle_PI = 1.3175*10^-6;
hc_SnOx_buckle_PI = 1.2619*10^-6;

figure;
x = [hc_Au_evap_PVSKRSPP*10^9 hc_Spiro_evap_PVSKRSPP*10^9 hc_PVSKRSPP_evap*10^9 hc_SnOx_buckle_PI*10^9 hc_ITOglass_buckle_PI*10^9; ...
    hc_Au_evap_PVSKhotcast*10^9 hc_Spiro_evap_PVSKhotcast*10^9 hc_PVSKhotcast_evap*10^9 hc_SnOx_buckle_PI*10^9 hc_ITOglass_buckle_PI*10^9;...
    hc_Au_evap_PVSKspincoat*10^9 hc_Spiro_evap_PVSKspincoat*10^9 hc_PVSKspincoat_evap*10^9 hc_SnOx_buckle_PI*10^9 hc_ITOglass_buckle_PI*10^9];
b = bar(x, 'FaceColor', 'flat');
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP PVSK', 'Hotcast PVSK', 'Spincoat PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
title("Polyimide substrate, NIP, evap. Au", 'FontWeight','Normal')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);

legend("Au", "Spiro", "PVSK", "SnOx", "ITO");

b(1).CData(1,:) = [0.980 0.866 0.505]; % group 1 1st bar
b(1).CData(2,:) = [0.980 0.866 0.505]; % group 2 1st bar
b(1).CData(3,:) = [0.980 0.866 0.505]; % group 3 1st bar

b(2).CData(1,:) = [0.192 0.192 0.447]; % group 1 2nd bar
b(2).CData(2,:) = [0.192 0.192 0.447]; % group 2 2nd bar
b(2).CData(3,:) = [0.192 0.192 0.447]; % group 3 2nd bar

b(3).CData(1,:) = [0.2 0.094 0.090]; % group 1 3rd bar
b(3).CData(2,:) = [0.2 0.094 0.090]; % group 2 3rd bar
b(3).CData(3,:) = [0.2 0.094 0.090]; % group 3 3rd bar

b(4).CData(1,:) = [0.694 0.411 0.321]; % group 1 4th bar
b(4).CData(2,:) = [0.694 0.411 0.321]; % group 2 4th bar
b(4).CData(3,:) = [0.694 0.411 0.321]; % group 3 4th bar

b(5).CData(1,:) = [0.533 0.796 0.858]; % group 1 5th bar
b(5).CData(2,:) = [0.533 0.796 0.858]; % group 2 5th bar
b(5).CData(3,:) = [0.533 0.796 0.858]; % group 3 5th bar


figure;
x = [hc_Ag_screen_PVSKRSPP*10^9 hc_Spiro_screen_PVSKRSPP*10^9 hc_PVSKRSPP_screen*10^9 hc_SnOx_buckle_PI*10^9 hc_ITOglass_buckle_PI*10^9; ...
    hc_Ag_screen_PVSKhotcast*10^9 hc_Spiro_screen_PVSKhotcast*10^9 hc_PVSKhotcast_screen*10^9 hc_SnOx_buckle_PI*10^9 hc_ITOglass_buckle_PI*10^9;...
    hc_Ag_screen_PVSKspincoat*10^9 hc_Spiro_screen_PVSKspincoat*10^9 hc_PVSKspincoat_screen*10^9 hc_SnOx_buckle_PI*10^9 hc_ITOglass_buckle_PI*10^9];
b = bar(x, 'FaceColor', 'flat');
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', {'RSPP PVSK', 'Hotcast PVSK', 'Spincoat PVSK'})
set(gca,'YScale','log')
ylabel('Critical Debond Thickness h_c, (nm)')
title("Polyimide substrate, NIP, screen-printed Ag", 'FontWeight','Normal')
box on;
set(gca, 'FontSize', 14, 'FontName', 'Arial');
% Determine the y-axis limits to use in the bar plot.
y_lowerlimit = (min(min(x))+1)/10;
y_upperlimit = max(max(x))*10;
ylim([y_lowerlimit, y_upperlimit]);

legend("Ag", "Spiro", "PVSK", "SnOx", "ITO");

b(1).CData(1,:) = [0.75 0.75 0.75]; % group 1 1st bar
b(1).CData(2,:) = [0.75 0.75 0.75]; % group 2 1st bar
b(1).CData(3,:) = [0.75 0.75 0.75]; % group 3 1st bar

b(2).CData(1,:) = [0.192 0.192 0.447]; % group 1 2nd bar
b(2).CData(2,:) = [0.192 0.192 0.447]; % group 2 2nd bar
b(2).CData(3,:) = [0.192 0.192 0.447]; % group 3 2nd bar

b(3).CData(1,:) = [0.2 0.094 0.090]; % group 1 3rd bar
b(3).CData(2,:) = [0.2 0.094 0.090]; % group 2 3rd bar
b(3).CData(3,:) = [0.2 0.094 0.090]; % group 3 3rd bar

b(4).CData(1,:) = [0.694 0.411 0.321]; % group 1 4th bar
b(4).CData(2,:) = [0.694 0.411 0.321]; % group 2 4th bar
b(4).CData(3,:) = [0.694 0.411 0.321]; % group 3 4th bar

b(5).CData(1,:) = [0.533 0.796 0.858]; % group 1 5th bar
b(5).CData(2,:) = [0.533 0.796 0.858]; % group 2 5th bar
b(5).CData(3,:) = [0.533 0.796 0.858]; % group 3 5th bar

