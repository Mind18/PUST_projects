import symulacja_obiektu8y_p4.*
clear;

% Inicjalizacja danych
T_p = 0.5;
k_konc = 400;

wejscia = 4;
wyjscia = 3;

u = cell(1, wejscia); % Sterowania obiektu
y = cell(1, wyjscia); % Wyjścia obiektu

for j=1:wejscia
    u{j} = zeros(k_konc, 1);
end

for i=1:wyjscia
    y{i} = zeros(k_konc, 1);
end

zad = ['Y', 'Y']; % Zadania, które będą wykonywane

if strcmp(zad(1), 'Y')
    punkt_pracy_4;
end

if strcmp(zad(2), 'Y')
    odp_skokowa_4;
end