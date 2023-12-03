import symulacja_obiektu8y_p3.*
clear;

% Inicjalizacja danych
T_p = 0.5;
k_konc = 400;
y_zad = -3;
Y_zad = [y_zad 0.1 -2.2 -0.2];

u_min = -1;
u_max = 1;
du_min = -2;
du_max = 2;

n_regulatorow = 5; % Liczba regulatorów
% Strefy rozmycia regulatorów
reg_part = {[-1 -1 -0.6 -0.4], [-0.5 -0.4 -0.3 -0.2], ...
    [-0.3 -0.2 0.05 0.15], [0.05 0.15 0.4 0.6], [0.4 0.6 1 1]}; 
u_konc = -1:(2/(n_regulatorow-1)):1;

% Parametry regulatora PID
K_r = 0.5; T_i = 5; T_d = 0.2; 
% Metoda Zieglera-Nicholsa - K_u=1.4 T_u=7.5

% Parametery regulatora DMC
D = 100; % Horyzont dynamiki
N = 90;   % Horyzont predykcji 98
N_u = 80;   % Horyzont sterowania 18
lambda = 0.5;
Lambda = lambda.*eye(N_u, N_u);
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);
e = zeros(1, k_konc);
e_pid(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu PID
e_dmc(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu DMC

s = {}; % Zestaw dostępnych odp.skokowych
odp_skok = 2; % Odpowiedź skokowa dla nierozmytego regulatora

zad = ['Y' 'Y' 'Y' 'Y' 'N']; % Zadania, które będą wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy_3;
end

if strcmp(zad(2), 'Y')
    odp_skokowa;
end

if strcmp(zad(3), 'Y')
    PID_bez_rozmycia;
end

if strcmp(zad(4), 'Y')
    DMC_bez_rozmycia;
end

if strcmp(zad(5), 'Y')
    fuzzy_division;
end
