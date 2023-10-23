import symulacja_obiektu8y_p1.*

% Inicjalizacja

clear;
zad_2_target = 'u'; % 'u' lub 'y'
u_konc = [0.7 1.3 1.1 0.5 1.5]; % Sygnały u(k) użye do odpowiedzi skokowej
k_konc = 400;
u(1, 1:11) = 0.5; % Sygnał początkowy do zad.1
u(1, 12:k_konc) = 1; % Sygnał końcowy do zad.2
y = zeros(1, k_konc);
% Ograniczenia regulatora DMC
du_max = 1;
du_min = -du_max;
u_min = 0.5;
u_max = 1.5;
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
    u(1, 1:200) = upp; % U_pp=1
    u(1, 201:k_konc) = u_konc(i);
    plot(1:k_konc, u);
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
    u(1, 1:200) = upp; % U_pp=1
    u(1, 201:k_konc) = u_konc(i);
    for k=12:k_konc
        y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    end

    if skok_u == i
        for k=1:D
            % Realizacja zadania 3
            s(k) = (y(k+200) - ypp) / (u_konc(skok_u) - upp);
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
export_fig('./pliki_wynikowe/odpowiedź_skokowa.pdf');

%% Realizacja zadania 4

%% Algorytm regulacji PID

% Inicjalizacja zmiennych
K_r = 0.5; T_p = 0.5; T_i = 9; T_d = 2; e_pid = 0;

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
yzad(1:11)=0; yzad(12:k_konc)=y_zad;
e(1:k_konc)=0;

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
    e_pid = e_pid + (yzad(k) - y(k))^2;
end

% Narysowanie wykresów
figure;
stairs(u);
xlabel('k');
ylabel('u(k)');
title("Sygnał sterujący u(k) algorytmu PID");

figure;
stairs(y);
hold on;
stairs(yzad, ':');
xlabel('k');
ylabel('y(k)');
title("Sygnał wyjściowy y(K) algorytmu PID");
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
%% Algorytm regulacji DMC

D = 200;   % Horyzont dynamiki
N = 200;   % Horyzont predykcji
N_u = 200;   % Horyzont sterowania
lambda = 170;
Lambda = lambda.*eye(N_u, N_u);
k_j = 0;
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);
e_dmc = 0; %Błąd średniokwadratowy dla algorytmu DMC

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
yzad(1:11)=ypp; yzad(12:k_konc)=y_zad;

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

    e_dmc = e_dmc + (yzad(k) - y(k))^2;
end

% Narysowanie wykresów
figure;
stairs(u);
xlabel('k');
ylabel('u(k)');
title("Sygnał sterujący u(k) algorytmu DMC");

figure;
stairs(y);
hold on;
stairs(yzad, ':');
xlabel('k');
ylabel('y(k)');
title("Sygnał wyjściowy y(K) algorytmu DMC");
legend('y(k)', 'y^{zad}', 'Location', 'southeast');