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

% Punkt pracy
upp = 0; ypp = 0;

n_regulatorow = 4; % Liczba regulatorów
kryterium = 'u'; % Wybieramy między u lub y - warunek do ustalenia
                % wartości funkcji przynależności
% Strefy rozmycia regulatorów
% reg_part = {[-1 -1 -0.6 -0.4], [-0.5 -0.4 -0.3 -0.2], ...
%     [-0.3 -0.2 0.05 0.15], [0.05 0.15 0.4 0.6], [0.4 0.6 1 1]};

if n_regulatorow == 2 % Dla n_regulatorow=2
    reg_part = {[-1 -1 -0.2 0.1], [-0.1 0.2 1 1]};
    u_konc = [-0.8 0.71];
elseif n_regulatorow == 3 % Dla n_regulatorow=3
    reg_part = {[-1 -1 -0.5 -0.4], [-0.45 -0.3 -0.05 0.1], [0 0.2 1 1]};
    u_konc = [-0.8 -0.2 0.71];
elseif n_regulatorow == 4
    reg_part = {[-1 -1 -0.5 -0.4], [-0.5 -0.4 -0.2 -0.1], ...
        [-0.25 -0.15 0 0.05], [0 0.05 1 1]};
    u_konc = [-0.8 -0.2 -0.1 0.71];
elseif n_regulatorow == 5
    reg_part = {[-1 -1 -0.6 -0.5], [-0.6 -0.5 -0.3 -0.2], ...
        [-0.3 -0.25 -0.15 -0.1], [-0.25 -0.15 0 0.2], [0 0.2 1 1]};
    u_konc = [-0.8 -0.4 -0.2 -0.1 0.71];
elseif n_regulatorow == 6
    reg_part = {[-1 -1 -0.7 -0.6], [-0.7 -0.6 -0.5 -0.4], ...
        [-0.55 -0.45 -0.35 -0.25], [-0.35 -0.25 -0.15 -0.1], ...
        [-0.2 -0.15 0 0.05], [0 0.2 1 1]};
    u_konc = [-0.8 -0.55 -0.4 -0.2 -0.1 0.71];
end

% u_konc = -1:(2/(n_regulatorow-1)):1;

% Wyeliminowanie przypadku w którym mamy zerową odp. skokową
for step=1:size(u_konc, 2)
    if u_konc(step) == 0
        u_konc(step) = 0.01;
    end
end

% Parametry regulatora PID
K_r = 0.22; T_i = 4.75; T_d = 0.45;
% Metoda Zieglera-Nicholsa - K_u=1.4 T_u=7.5

% Uzyskane nastawy dla regulatora nierozmytego
% K_r = 0.22; T_i = 4.75; T_d = 0.45; 

% Parametry rozmytego regulatora PID
if n_regulatorow == 2 
    K_r_lok = [0.13 1.75];
    T_i_lok = [4.75 3];
    T_d_lok = [0.45 0.4];

    D_fuz = [88 88]; % Horyzonty dynamiki
    N_fuz = [30 30];   % Horyzonty predykcji
    N_u_fuz = [5 5];   % Horyzonty sterowania
    lambda_fuz = [1 1]; % Parametry lambda lokalnych 
                                    % regulatorów
elseif n_regulatorow == 3
    K_r_lok = [0.13 0.6 1.75];
    T_i_lok = [4.75 2.75 3];
    T_d_lok = [0.45 0.42 0.4];

    D_fuz = [88 88 88]; % Horyzonty dynamiki
    N_fuz = [30 30 30];   % Horyzonty predykcji
    N_u_fuz = [5 10 5];   % Horyzonty sterowania
    lambda_fuz = [1 1 1]; % Parametry lambda lokalnych 
                                    % regulatorów
elseif n_regulatorow ==4
    K_r_lok = [0.13 0.6 0.6 1.75];
    T_i_lok = [4.75 2.75 3.75 3];
    T_d_lok = [0.45 0.42 0.42 0.4];

    D_fuz = [88 88 88 88]; % Horyzonty dynamiki
    N_fuz = [30 30 30 30];   % Horyzonty predykcji
    N_u_fuz = [5 5 15 5];   % Horyzonty sterowania
    lambda_fuz = [1 1 1 1]; % Parametry lambda lokalnych 
                                    % regulatorów
elseif n_regulatorow == 5
    K_r_lok = [0.13 0.3 0.6 0.6 1.75];
    T_i_lok = [4.75 4.75 2.75 3.75 3];
    T_d_lok = [0.45 0.42 0.42 0.42 0.4];

    D_fuz = [88 88 88 88 88]; % Horyzonty dynamiki
    N_fuz = [30 30 30 30 30];   % Horyzonty predykcji
    N_u_fuz = [5 5 10 15 5];   % Horyzonty sterowania
    lambda_fuz = [1 1 1 1 1]; % Parametry lambda lokalnych 
                                    % regulatorów
elseif n_regulatorow == 6
    K_r_lok = [0.13 0.22 0.3 0.6 0.6 1.75];
    T_i_lok = [4.75 4.75 4.75 2.75 3.75 3];
    T_d_lok = [0.45 0.42 0.42 0.42 0.42 0.4];

    D_fuz = [88 88 88 88 88 88]; % Horyzonty dynamiki
    N_fuz = [30 30 30 30 30 30];   % Horyzonty predykcji
    N_u_fuz = [5 5 5 5 5 5];   % Horyzonty sterowania
    lambda_fuz = [1 1 1 1 1 1]; % Parametry lambda lokalnych 
                                    % regulatorów
end

% Parametery regulatora DMC
D = 88; % Horyzont dynamiki
N = 30;   % Horyzont predykcji 98
N_u = 15;   % Horyzont sterowania 18
lambda = 1;
Lambda = lambda.*eye(N_u, N_u);
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);

% Parametry rozmytego regulatora DMC

Lambda_fuz = {zeros(N_u, N_u), zeros(N_u, N_u), zeros(N_u, N_u), ...
    zeros(N_u, N_u), zeros(N_u, N_u)}; % Alokacja pamięci na macierze
                                       % lambda lokalnych regulatorów
M_fuz = {zeros(N, N_u), zeros(N, N_u), zeros(N, N_u), ...
    zeros(N, N_u), zeros(N, N_u)}; % Macierze M lokalnych regulatorów
M_p_fuz = {zeros(N, D-1), zeros(N, D-1), zeros(N, D-1), ...
    zeros(N, D-1), zeros(N, D-1)}; % Macierze M_p lokalnych regulatorów
U_p_fuz = {zeros(D-1, 1), zeros(D-1, 1), zeros(D-1, 1), ...
    zeros(D-1, 1), zeros(D-1, 1)}; % Wektory U_p lokalnych regulatorów

e = zeros(1, k_konc);
e_pid(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu PID
e_pid_fuz(1:k_konc) = 0; % Błąd średniokwadratowy dla rozmytego PID
e_dmc(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu DMC
e_dmc_fuz(1:k_konc) = 0; % Błąd średniokwadratowy dla rozmytego DMC

w = zeros(1, n_regulatorow);

s = {}; % Zestaw dostępnych odp.skokowych
odp_skok = 4; % Odpowiedź skokowa dla nierozmytego regulatora

zad = ['N' 'Y' 'N' 'N' 'N' 'N' 'Y']; % Zadania, które będą wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy_3;
end

if strcmp(zad(2), 'Y')
    odp_skokowa;
end

k_konc = 400; % Zmiana k_konc na potrzebę wydłużenia testu regulacji

if strcmp(zad(3), 'Y')
    PID_bez_rozmycia;
end

if strcmp(zad(4), 'Y')
    DMC_bez_rozmycia;
end

if strcmp(zad(5), 'Y')
    fuzzy_division;
end

if strcmp(zad(6), 'Y')
    rozmyty_PID;
end

if strcmp(zad(7), 'Y')
    rozmyty_DMC;
end
