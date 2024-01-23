import symulacja_obiektu8y_p4.*
clear;

% Inicjalizacja danych
T_p = 0.5;
k_konc = 400;
y_zad = 1.2;
Y_zad = {[y_zad 1.8 1 2], [y_zad 1 0.3 0.8], [y_zad 0.8 1.9 0.0]};
u_konc = 1;

du_min = -2;
du_max = 2;

u_min = -20;
u_max = 20;

wejscia = 4;
wyjscia = 3;

u = zeros(wejscia, k_konc); % Sterowania obiektu
y = zeros(wyjscia, k_konc); % Wyjścia obiektu

% Punkt pracy
upp = 0; ypp = 0;

% Parametry regulatora PID
K_r = [0.1 0.02 0.5];  % [0.1 0.02 0.002], [0.1 0.02 0.006]
T_i = [0.3 0.1 2];   % [0.3 0.05 0.02], [0.3 0.05 0.02]
T_d = [0.2 0.1 1e-6];    % [0.2 0.1 0.02], [0.2 0.1 0.015]

% Sparowanie dla regulatora PID wejść obiektu do jego wyjść
u_dla_y = [4 1 3];

% Parametery regulatora DMC
D = 400; % Horyzont dynamiki
N = D;   % Horyzont predykcji 98
N_u = D;   % Horyzont sterowania 18
psi = [1 1 1];
lambda = [5.4 8.7 6.1 25]; % 12 9 6 20

for j=1:wejscia
    u(j, :) = zeros(k_konc, 1);
end

for i=1:wyjscia
    y(i, :) = zeros(k_konc, 1);
end

zad = ['Y', 'Y', 'Y', 'N', 'Y']; % Zadania, które będą wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy_4;
end

if strcmp(zad(2), 'Y')
    odp_skokowa_4;
end

if strcmp(zad(3), 'Y')
    k_konc = 2000;
    PID;
end

if strcmp(zad(4), 'Y')
    k_konc = 2000;
    DMC_oszczedny;
end

if strcmp(zad(5), 'Y')
    k_konc = 2000;
    PID_optymalizacja;
end