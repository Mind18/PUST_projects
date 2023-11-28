% Inicjalizacja danych
T_p = 0.5;
k_konc = 400;

u_min = -1;
u_max = 1;

n_regulatorow = 5; 
u_konc = -1:(2/(n_regulatorow-1)):1;

D = k_konc; % Horyzont dynamiki sterowania D

s = zeros(1, D);

zad = ['Y' 'Y']; % Zadania, które będą wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy_3;
end

if strcmp(zad(2), 'Y')
    odp_skokowa
end
