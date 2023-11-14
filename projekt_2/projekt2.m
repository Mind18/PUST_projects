
% Inicjalizacja danych
T_p = 0.5;
k_konc = 400;
print = 'Y'; % Czy generować wykresy

% Wektory zmiennych procesowych
u = zeros(k_konc, 1); % Sterowanie
z = zeros(k_konc, 1); % Zakłócenie
y = zeros(k_konc, 1); % Wyjście procesu
y_z = zeros(k_konc, 1); % Wyjście procesu

% Ograniczenia regulatorów
du_max = 1;
du_min = -du_max;
u_min = 0.5; % 0.5 z racji na uwarunkowanie obiektu
u_max = 1.5; % 1.5 z racji na uwarunkowanie obiektu

% Punkt pracy
upp = 0; zpp = 0; ypp = 0;
u_konc = [0.3 0.7 1 1.2 1.5]; % u dla których wyznaczana jest odp.skokowa

% Inicjalizacja zmiennych algorytmu DMC z zakłóceniem
D = 229; % Horyzont dynamiki sterowania D
D_z = 56; % Horyzont dynamiki zakłócenia D_z
s = zeros(1, D);
s_z = zeros(1, D_z);
skok_u = 5; % Wybrana odpowiedź skokowa do wyznaczenia wektora s;
skok_z = 3;

zad = ['Y' 'Y']; % Zadanie, które będzie wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy;
end

if strcmp(zad(2), 'Y')
    odpowiedz_skokowa;
end
