function [sqrd_err] = DMC_SE(lambda, N, N_u)

% Inicjalizacja zmiennych
k_konc = 400;
u_konc = 1.5;
du_max = 1; % 1
du_min = -du_max;
u_min = 0.5; % 0.5
u_max = 1.5; % 1.5
% Wartość zadana dla regulatora DMC
y_zad = 1.2;

% Punkt pracy według funkcji symulacji obiektu
upp = 1;
ypp = 1.7;

D = 200;   % Horyzont dynamiki
% N = dmc_param(1);   % Horyzont predykcji
% N_u = dmc_param(2);   % Horyzont sterowania
% lambda = dmc_param(3); % Parametr kary sterowania lambda
Lambda = lambda.*eye(N_u, N_u);
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);
e_dmc(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu DMC
s = zeros(1, D); % Wektor odpowiedzi skokowej

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
yzad(1:11)=ypp; yzad(12:k_konc)=y_zad;

% Generacja odpoweidzi skokowej

u(1, 1:200) = upp; % U_pp=1
u(1, 201:k_konc) = u_konc;
for k=12:k_konc
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
end

for k=1:D
    s(k) = (y(k+200) - ypp) / (u_konc - upp);
end

% Generacja macierzy M
for j=1:N_u % dla każdej kolumny macierzy M
    for i=j:N % Dla każdego wiersza kolumny j począwszy od przekątnej
        M(i, j) = s(i-j+1);
    end
end

% Generacja macierzy M_p
for j=1:D-1 % dla każdej kolumny macierzy M_p
    for i=1:N % dla każdego wiersza macierzy M_p
        if j+i > D
            p = D;
        else
            p = j+i;
        end
        M_p(i, j) = s(p) - s(j);
    end
end

% Wyznaczenie wektora współczynników K
K = ((M'*M+Lambda)^(-1))*M';

% Wyznaczenie współczynnika k_e
k_e = sum(K(1,:));

% warunki początkowe dla głównej pętli symulacyjnej
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
yzad(1:11)=ypp; yzad(12:k_konc)=y_zad;

for k=12:k_konc
    % symulacja obiektu
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    % wyznaczenie zmiany sterowania
    k_j = 0;
    for j=1:D-1
        k_j = k_j + ((K(1, :)*M_p(:, j))*U_p(j));
    end
    delta_u = k_e*(yzad(k) - y(k)) - k_j;
    % Ograniczenia zmiany sterowania
    if delta_u < du_min
        delta_u = du_min;
    elseif delta_u > du_max
        delta_u = du_max;
    end
    % Zapamiętanie zmiany sterowania do kolejnych iteracji
    for n=D-1:-1:2
        U_p(n,1) = U_p(n-1,1);
    end
    U_p(1) = delta_u;
    % Dokonanie zmiany sterowania
    u(k) = u(k-1) + delta_u;
    % Ograniczenia wartości sterowania
    if u(k) < u_min
        u(k) = u_min;
    elseif u(k) > u_max
        u(k) = u_max;
    end

    e_dmc(k) = e_dmc(k-1) + (yzad(k) - y(k))^2;
end

sqrd_err = e_dmc(k_konc);
