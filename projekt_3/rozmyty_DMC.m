% Rozmyty regulator DMC
% 1. Zasymulowanie obiektu
% 2. Obliczenie wartości funkcji przynależności każdej reguły
% 3. Obliczenie wartości sterowania dla każdego regulatora
% 4. Obliczenie wartości sterowania jako średnia ważona sterowań
%    z wartościami funkcji przynależności jako wagi

import trapezoid_assign.*
k_konc = 1200;

% Inicjalizacja zmiennych
o = '2_5'; % Nazwa wykresu

% warunki początkowe
u = zeros(1, k_konc);
u_w = {zeros(1, k_konc), zeros(1, k_konc), zeros(1, k_konc), ...
    zeros(1, k_konc), zeros(1, k_konc)};
u_part = zeros(1, n_regulatorow); % Wektor do trzymania częściowych
                                  % Obliczeń sterowania
y_dmc_fuz = zeros(1, k_konc);
u(1:14)=upp; y_dmc_fuz(1:14)=ypp;

e = zeros(1, k_konc);

% Generacja zmiennej trajektori
yzad(1:14)=ypp;
yzad(15:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

% yzad(1:14)=ypp;
% yzad(15:k_konc)=Y_zad(1);

% Generacje macierzy Lambda
for r=1:n_regulatorow
    Lambda_fuz{r} = lambda_fuz(r).*eye(N_u, N_u);
end

% Generacja macierzy M
for r=1:n_regulatorow
    for j=1:N_u % dla każdej kolumny macierzy M
        for i=j:N % Dla każdego wiersza kolumny j począwszy od przekątnej
            M_fuz{r}(i, j) = s{r}(i-j+1);
        end
    end
end

% Generacja macierzy M_p
for r=1:n_regulatorow
    for j=1:D-1 % dla każdej kolumny macierzy M_p
        for i=1:N % dla każdego wiersza macierzy M_p
            if j+i > D
                p = D;
            else
                p = j+i;
            end
            M_p_fuz{r}(i, j) = s{r}(p) - s{r}(j);
        end
    end
end

% Wyznaczenie lokalnych wektorów współczynników K
for r=1:n_regulatorow
    K_fuz{r} = ((M_fuz{r}'*M_fuz{r}+Lambda_fuz{r})^(-1))*M_fuz{r}';
end

% Wyznaczenie współczynnika k_e
for r=1:n_regulatorow
    k_e_fuz(r) = sum(K_fuz{r}(1,:));
end

for k=15:k_konc % Główna pętla symulacji
    
    if k == 17
        u;
    end

    y_dmc_fuz(k) = symulacja_obiektu8y_p3(u(k-5), u(k-6), ...
        y_dmc_fuz(k-1), y_dmc_fuz(k-2));

    % uchyb regulacji
    e(k)=yzad(k) - y_dmc_fuz(k);
    k_j_fuz = zeros(1, n_regulatorow);

    % Wybranie kryterium dla funkcji przynależności
    if strcmp(kryterium, 'u') % Jeżeli kryterium jest sterowanie
        x_w = u(k-1); % Wykorzystaj poprzednią wartość sterowania
    elseif strcmp(kryterium, 'y') % Jeżeli kryterium jest sygnał wyjściowy
        x_w = y(k); % Wykorzystaj aktualną wartość sygnału wyjściowego
    end

    % Logika rozmyta
    for r=1:n_regulatorow % dla każdego z regulatorów lokalnych
        % Wyznaczenie wartości funkcji przynależności
        w(r) = trapezoid_assign(x_w, reg_part{r});
        
        for j=1:D-1
            k_j_inc = K_fuz{r}(1, :)*M_p_fuz{r}(:, j);
            k_j_fuz(r) = k_j_fuz(r) + k_j_inc*U_p_fuz{r}(j);
        end
        % sygnał sterujący regulatora lokalnego PID
        delta_u = k_e_fuz(r)*e(k) - k_j_fuz(r);
        % Ograniczenia zmiany sterowania
        if delta_u < du_min
            delta_u = du_min;
        elseif delta_u > du_max
            delta_u = du_max;
        end
        % Zapamiętanie zmiany sterowania do kolejnych iteracji
        for n=D-1:-1:2
            U_p_fuz{r}(n,1) = U_p_fuz{r}(n-1,1);
        end
        U_p_fuz{r}(1) = delta_u;
        % Dokonanie zmiany sterowania
        u_w{r}(k) = u(k-1) + delta_u;
        % Ograniczenia wartości sterowania
        if u_w{r}(k) < u_min
            u_w{r}(k) = u_min;
            U_p_fuz{r}(1) = u_w{r}(k)-u(k-1);
        elseif u_w{r}(k) > u_max
            u_w{r}(k) = u_max;
            U_p_fuz{r}(1) = u_w{r}(k)-u(k-1);
        end

        u_part(r) = w(r)*u_w{r}(k);
    end
    
    % Wyznaczenie sterowania regulatora
    u(k) = sum(u_part) / sum(w);

    % Błąd średniokwadratowy dla rozmytego algorytmu PID
    e_dmc_fuz(k) = e_dmc_fuz(k-1) + (yzad(k) - y_dmc_fuz(k))^2;
end

E = e_dmc_fuz(k_konc); % Błąd średniokwadratowy algorytmu

% Narysowanie wykresów
figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([min(u)-0.1 max(u)+0.1]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm rozmyty DMC u(k): Liczba regulatorów: " + ...
    string(n_regulatorow) + " E=" + string(E);
title(title_str);
filenameu = "./pliki_wynikowe/"+"regulator_dmc_rozmyty_u(k)"+ ...
    string(o)+".pdf";
export_fig(filenameu);

figure;
stairs(y_dmc_fuz); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([min(y_dmc_fuz)-0.1 max(y_dmc_fuz)+0.1]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm rozmyty DMC y(k): Liczba regulatorów: " + ...
    string(n_regulatorow) + " E=" + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
filenamey = "./pliki_wynikowe/"+"regulator_dmc_rozmyty_y(k)"+ ...
    string(o)+".pdf";
export_fig(filenamey);