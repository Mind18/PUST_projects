%% TO DO
% Wyregulować algorytmy
% Upewnić się czy generowane grafy są ładne
% 

import symulacja_obiektu8y_p1.*
import PID_SE.*
import DMC_SE.*

% Inicjalizacja 

clear;
zad_2_target = 'u'; % 'u' lub 'y'
u_konc = [0.7 1.3 1.1 0.5 1.5]; % Sygnały u(k) użye do odpowiedzi skokowej
k_konc = 400;
u(1, 1:11) = 0.5; % Sygnał początkowy do zad.1
u(1, 12:k_konc) = 1; % Sygnał końcowy do zad.2
y = zeros(1, k_konc);
% Ograniczenia regulatorów
du_max = 1;
du_min = -du_max;
u_min = 0.5; % 0.5 z racji na uwarunkowanie obiektu
u_max = 1.5; % 1.5 z racji na uwarunkowanie obiektu
% Wartość zadana dla regulatora DMC
y_zad = 1.2;

% Inicjalizacja danych algorytmu DMC
D = 200; % Horyzont dynamiki D
s = zeros(1, D);
skok_u = 5; % Wybrana odpowiedź skokowa do wyznaczenia wektora s;

% Punkt pracy według funkcji symulacji obiektu
upp = 1;
ypp = 1.7;

%% Realizacja zadania 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1

% Symulacja obiektu
for k=12:k_konc
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
end

% Narysowanie wykresu u(k)
figure;
plot(1:k_konc, u);
hold on;
xlabel('k');
ylabel('u(k)');
ylim([0 1.2]);
title('Wykres u(k)');
hold off;
export_fig('./pliki_wynikowe/test_punktu_pracy_u(k).pdf');

% Narysowanie wykresu y(k)
figure;
plot(1:k_konc, y);
hold on;
xlabel('k');
ylabel('y(k)');
title('Wykres y(k)');
hold off;
export_fig('./pliki_wynikowe/test_punktu_pracy_y(k).pdf');

%% Realizacja zadania 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wektory dla przechowywania danych statycznych
u_stat = zeros(1, size(u_konc, 2)+1);
u_stat(1) = upp;
u_stat(2:size(u_konc, 2)+1) = u_konc;
y_stat = zeros(1, size(u_konc, 2)+1);
y_stat(1) = ypp;
% wektory dla przechowywania danych statycznych

% Narysowanie wykresów sygnałów u(k) odp. skokowej
figure;
for i=1:size(u_konc, 2)
    u(1, 1:k_konc) = upp; % U_pp=1
    u(1, 15:k_konc) = u_konc(i);
    plot(1:k_konc, u);
    ylim([0.4 1.6]);
    hold on;
end

% Dodanie informacji do wygenerowanego wykresu
legend_list = strings([1 size(u_konc, 2)]);
xlabel('k');
ylabel('u(k)');
ylim([0.4 1.6]);
title('Wykres u(k)');
for i=1:size(u_konc, 2)
    legend_text = 'u(k)=' + string(u_konc(i));
    legend_list(i) = legend_text;
end
legend(legend_list, 'Location', 'best');
hold off;
export_fig('./pliki_wynikowe/zad2_u(k).pdf');

% Narysowanie wykresów y(k) odp. skokowej
figure;
for i=1:size(u_konc, 2)
    u(1, 1:k_konc) = upp; % U_pp=1
    u(1, 15:k_konc) = u_konc(i);
    y(1, 1:k_konc) = ypp;
    for k=12:k_konc
        y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    end

    if skok_u == i
        for k=1:D
            % Realizacja zadania 3
            s(k) = (y(k+15) - ypp) / 0.5; % (u_konc(skok_u) - upp);
        end
    end

    plot(1:k_konc, y);
    hold on;
  
    y_stat(i+1) = y(k_konc); % zapis danych dla wektora statycznego y
end

% Dodanie informacji do wygenerowanego wykresu
legend_list = strings([1 size(u_konc, 2)]);
xlabel('k');
ylabel('y(k)');
title('Wykres y(k)');
for i=1:size(u_konc, 2)
    legend_text = 'y(k)=' + string(y_stat(i+1));
    legend_list(i) = legend_text;
end
legend(legend_list, 'Location', 'southeast');
hold off;
export_fig('./pliki_wynikowe/zad2_y(k).pdf');

% ch-ka stat%
% sortowanie i rysowanie danych
figure;
[u_stat, sortIndex] = sort(u_stat);
y_stat = y_stat(sortIndex);

plot(u_stat, y_stat);
hold on;
plot(u_stat, y_stat, '.', 'MarkerSize',12);
legend("Interpolacja danych statycznych", "Dane statyczne", 'Location', 'southeast');
k_stat = rdivide(y_stat, u_stat);
xlabel('u');
ylabel('y(u)');
title('Dane statyczne y(u)');
hold off;
export_fig('./pliki_wynikowe/zad2_y(u)_stat.pdf');

% Narysowanie odpowiedzi skokowej
figure;
plot(1:200, s);
title_name = "Odpowiedź skokowa układu dla u_{konc}=" + ...
    string(u_konc(skok_u));
title(title_name);
xlabel('k');
ylabel('y(k)');
export_fig('./pliki_wynikowe/odpowiedz_skokowa.pdf');

%% Realizacja zadania 4

%% Algorytm regulacji PID
k_konc = 1200;
Y_zad = [y_zad 1.8 1.3 2.3];

% Inicjalizacja zmiennych
K_r = 0.32; T_p = 0.5; T_i = 7.3; T_d = 0.3; 
% K_r - 0.5 -> 0.35 -> 0.32
% T_i - 9 -> 8 -> 7.3
% T_d - 2 -> 0.2
% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
% Generacja zmiennej trajektori
yzad(1:11)=ypp;
yzad(12:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

e(1:k_konc)=0; e_pid(1:k_konc) = 0;

% Współczynniki algorytmu
r2 = (K_r * T_d) / T_p;
r1 = K_r * (T_p/(2*T_i) - 2*(T_d / T_p) - 1);
r0 = K_r * (1 + (T_p / (2*T_i)) + (T_d/T_p));

for k=12:k_konc % główna pętla symulacyjna
    % symulacja obiektu
    y(k)=symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    % uchyb regulacji
    e(k)=yzad(k) - y(k);
    % sygnał sterujący regulatora PID
    u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    % Ograniczenia zmiany sterowania
    du = u(k) - u(k-1);
    if du < du_min
        u(k) = u(k-1) + du_min;
    elseif du > du_max
        u(k) = u(k-1) + du_max;
    end
    % Ograniczenia wartości sterowania
    if u(k) < u_min
        u(k) = u_min;
    elseif u(k) > u_max
        u(k) = u_max;
    end

    % Błąd średniokwadratowy dla algorytmu PID
    e_pid(k) = e_pid(k-1) + (yzad(k) - y(k))^2;
end

E = e_pid(k_konc);

% Narysowanie wykresów
figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([0.4 1.6]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm PID u(k): K_r=" + string(K_r) ...
    + " T_i=" + string(T_i) + " T_D=" + string(T_d) + " E=" ...
    + string(E);
title(title_str);
export_fig('./pliki_wynikowe/regulator_pid_u(k).pdf');

figure;
stairs(y); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([0.9 2.5]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm PID y(k): K_r=" + string(K_r) ...
    + " T_i=" + string(T_i) + " T_D=" + string(T_d) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
export_fig('./pliki_wynikowe/regulator_pid_y(k).pdf');

%% Algorytm regulacji DMC
D = 200;   % Horyzont dynamiki
N = 80;   % Horyzont predykcji
N_u = N;   % Horyzont sterowania
lambda = 15;
Lambda = lambda.*eye(N_u, N_u);
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);
e = zeros(k_konc, 1);
e_dmc(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu DMC

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:k_konc)=upp; y(1:11)=ypp;
% Generacja zmiennej trajektori
yzad(1:11)=ypp;
yzad(12:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

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

for k=12:k_konc
    % symulacja obiektu
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    % wyznaczenie zmiany sterowania
    k_j = 0;
    e(k) = yzad(k) - y(k);
    for j=1:D-1
        k_j_inc = K(1, :)*M_p(:, j);
        k_j = k_j + k_j_inc*U_p(j);
    end
    delta_u = k_e*e(k) - k_j;
    % Ograniczenia zmiany sterowania
    if delta_u < du_min
        delta_u = du_min;
    end
    if delta_u > du_max
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
        U_p(1) = u(k)-u(k-1);
    end
    if u(k) > u_max
        u(k) = u_max;
        U_p(1) = u(k)-u(k-1);
    end

    e_dmc(k) = e_dmc(k-1) + e(k)^2;
end

E = e_dmc(k_konc); % Błąd średniokwadratowy algorytmu

% Narysowanie wykresów
figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([0.4 1.6]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm DMC u(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
export_fig('./pliki_wynikowe/regulator_dmc_u(k).pdf');

figure;
stairs(y); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([0.9 2.5]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm DMC y(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
export_fig('./pliki_wynikowe/regulator_dmc_y(k).pdf');

%% Zadanie 6

% Optymalizacja parametrów PID

x0 = [K_r, T_i, T_d];
A = [1, 1, 1];
b = 100;
pid_params = fmincon(@PID_SE, x0, A, b);

% Inicjalizacja zmiennych
K_r = pid_params(1); T_p = 0.5; T_i = pid_params(2); T_d = pid_params(3); 

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
% Generacja zmiennej trajektori
yzad(1:11)=ypp;
yzad(12:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

e(1:k_konc)=0; e_pid_fmincon(1:k_konc) = 0;

% Współczynniki algorytmu
r2 = (K_r * T_d) / T_p;
r1 = K_r * (T_p/(2*T_i) - 2*(T_d / T_p) - 1);
r0 = K_r * (1 + (T_p / (2*T_i)) + (T_d/T_p));

for k=12:k_konc % główna pętla symulacyjna
    % symulacja obiektu
    y(k)=symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    % uchyb regulacji
    e(k)=yzad(k) - y(k);
    % sygnał sterujący regulatora PID
    u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    % Ograniczenia zmiany sterowania
    du = u(k) - u(k-1);
    if du < du_min
        u(k) = u(k-1) + du_min;
    elseif du > du_max
        u(k) = u(k-1) + du_max;
    end
    % Ograniczenia wartości sterowania
    if u(k) < u_min
        u(k) = u_min;
    elseif u(k) > u_max
        u(k) = u_max;
    end

    % Błąd średniokwadratowy dla algorytmu PID
    e_pid_fmincon(k) = e_pid_fmincon(k-1) + (yzad(k) - y(k))^2;
end

E = e_pid_fmincon(k_konc);

% Narysowanie wykresów
figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([0.4 1.6]);
xlabel('k');
ylabel('u(k)');
title_str = "PID u(k) - optymalizacja: K_r=" + string(K_r) ...
    + " T_i=" + string(T_i) + " T_D=" + string(T_d) + " E=" ...
    + string(E);
title(title_str);
export_fig('./pliki_wynikowe/regulator_pid_u(k)_optymalizacja.pdf');

figure;
stairs(y); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([0.9 2.5]);
xlabel('k');
ylabel('y(k)');
title_str = "PID y(k) - optymalizacja: K_r=" + string(K_r) ...
    + " T_i=" + string(T_i) + " T_D=" + string(T_d) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
export_fig('./pliki_wynikowe/regulator_pid_y(k)_optymalizacja.pdf');

% Optymalizacja parametrów DMC

% Inicjalizacja parametrów
e_min = Inf;

x0 = 15;
A = 0.5;
b = 100;
N_test = 80;
N_u_test = 80;

% Wyznaczenie optymalnej wartości lambda dla N_test=200 i N_u_test=200
f = @(lambda_test)DMC_SE(lambda_test, N_test, N_u_test);
dmc_param = fmincon(f, x0, A, b);

% Zrealizuj symulację algorytmu DMC
D = 200;   % Horyzont dynamiki
N = N_test;   % Horyzont predykcji
N_u = N_u_test;   % Horyzont sterowania
lambda = dmc_param(1);
Lambda = lambda.*eye(N_u, N_u);
k_j = 0;
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);
e = zeros(k_konc, 1);
e_dmc(1:k_konc) = 0; % Błąd średniokwadratowy dla algorytmu DMC

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
% Generacja zmiennej trajektori
yzad(1:11)=ypp;
yzad(12:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

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

for k=12:k_konc
    % symulacja obiektu
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    % wyznaczenie zmiany sterowania
    k_j = 0;
    e(k) = yzad(k) - y(k);
    for j=1:D-1
        k_j_inc = K(1, :)*M_p(:, j);
        k_j = k_j + k_j_inc*U_p(j);
    end
    delta_u = k_e*e(k) - k_j;
    % Ograniczenia zmiany sterowania
    if delta_u < du_min
        delta_u = du_min;
    end
    if delta_u > du_max
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
        U_p(1) = u(k)-u(k-1);
    end
    if u(k) > u_max
        u(k) = u_max;
        U_p(1) = u(k)-u(k-1);
    end

    e_dmc(k) = e_dmc(k-1) + e(k)^2;
end

E = e_dmc(k_konc);

% Narysowanie wykresów
figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([0.4 1.6]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm DMC u(k) - optymalizacja: N=" + string(N_test) ...
    + " N_u=" + string(N_u_test) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
export_fig('./pliki_wynikowe/regulator_dmc_u(k)_optymalizacja.pdf');

figure;
stairs(y); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([0.9 2.5]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm DMC y(k) - optymalizacja: N=" + string(N_test) ...
    + " N_u=" + string(N_u_test) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
export_fig('./pliki_wynikowe/regulator_dmc_y(k)_optymalizacja.pdf');
