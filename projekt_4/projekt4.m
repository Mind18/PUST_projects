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

u_min = -10;
u_max = 10;

wejscia = 4;
wyjscia = 3;

u = cell(1, wejscia); % Sterowania obiektu
y = cell(1, wyjscia); % Wyjścia obiektu

% Punkt pracy
upp = 0; ypp = 0;

% Parametry regulatora PID
K_r = [0.1 0.02 0.002];  
T_i = [0.3 0.01 0.02];
T_d = [0.2 0.1 0.02];

% Sparowanie dla regulatora PID wejść obiektu do jego wyjść
u_dla_y = [4 1 3];

% Parametery regulatora DMC
D = 400; % Horyzont dynamiki
N = 10;   % Horyzont predykcji 98
N_u = 15;   % Horyzont sterowania 18
psi = [1 2 3];
lambda = [4 5 6];

for j=1:wejscia
    u{j} = zeros(k_konc, 1);
end

for i=1:wyjscia
    y{i} = zeros(k_konc, 1);
end

zad = ['N', 'N', 'N', 'Y']; % Zadania, które będą wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy_4;
end

if strcmp(zad(2), 'Y')
    odp_skokowa_4;
end

if strcmp(zad(3), 'Y')
    k_konc = 1200;
    PID;
end

if strcmp(zad(4), 'Y')
    DMC_oszczedny;
end