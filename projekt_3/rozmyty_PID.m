% Rozmyty regulator PID
% 1. Zasymulowanie obiektu
% 2. Obliczenie wartości funkcji przynależności każdej reguły
% 3. Obliczenie wartości sterowania dla każdego regulatora
% 4. Obliczenie wartości sterowania jako średnia ważona sterowań
%    z wartościami funkcji przynależności jako wagi

import trapezoid_assign.*
k_konc = 1200;

% Inicjalizacja zmiennych
o = 2; % Nazwa wykresu

% warunki początkowe
u = zeros(1, k_konc);
u_w = {zeros(1, k_konc), zeros(1, k_konc), zeros(1, k_konc), ...
    zeros(1, k_konc), zeros(1, k_konc)};
u_part = zeros(1, n_regulatorow); % Wektor do trzymania częściowych
                                  % Obliczeń sterowania
e = zeros(1, k_konc);

y_pid_fuz = zeros(1, k_konc);
u(1:14)=upp; y_pid_fuz(1:14)=ypp;

% Generacja zmiennej trajektori
yzad(1:14)=ypp;
yzad(15:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

% yzad(1:14)=ypp;
% yzad(15:k_konc)=Y_zad(1);

% Współczynniki algorytmu
r2_fuz = zeros(1, n_regulatorow);
r1_fuz = zeros(1, n_regulatorow);
r0_fuz = zeros(1, n_regulatorow);
for i=1:n_regulatorow
    r2_fuz(i) = (K_r_lok(i) * T_d_lok(i)) / T_p;
    r1_fuz(i) = K_r_lok(i) * (T_p/(2*T_i_lok(i)) - 2*(T_d_lok(i) / T_p) - 1);
    r0_fuz(i) = K_r_lok(i) * (1 + (T_p / (2*T_i_lok(i))) + (T_d_lok(i)/T_p));
end

for k=15:k_konc % Główna pętla symulacji
    y_pid_fuz(k) = symulacja_obiektu8y_p3(u(k-5), u(k-6), ...
        y_pid_fuz(k-1), y_pid_fuz(k-2));

    % uchyb regulacji
    e(k)=yzad(k) - y_pid_fuz(k);

    % Wybranie kryterium dla funkcji przynależności
    if strcmp(kryterium, 'u') % Jeżeli kryterium jest sterowanie
        x_w = u(k-1); % Wykorzystaj poprzednią wartość sterowania
    elseif strcmp(kryterium, 'y') % Jeżeli kryterium jest sygnał wyjściowy
        x_w = y(k); % Wykorzystaj aktualną wartość sygnału wyjściowego
    end

    % Logika rozmyta
    for i=1:n_regulatorow % dla każdego z regulatorów lokalnych
        % Wyznaczenie wartości funkcji przynależności
        w(i) = trapezoid_assign(x_w, reg_part{i});

        % sygnał sterujący regulatora lokalnego PID
        u_w{i}(k)=r2_fuz(i)*e(k-2)+r1_fuz(i)*e(k-1)+r0_fuz(i)*e(k)+u(k-1);
        % Ograniczenia zmiany sterowania
        du = u_w{i}(k) - u(k-1);
        if du < du_min
            u_w{i}(k) = u(k-1) + du_min;
        elseif du > du_max
            u_w{i}(k) = u(k-1) + du_max;
        end
        % Ograniczenia wartości sterowania
        if u_w{i}(k) < u_min
            u_w{i}(k) = u_min;
        elseif u_w{i}(k) > u_max
            u_w{i}(k) = u_max;
        end

        u_part(i) = w(i)*u_w{i}(k);
    end
    
    % Wyznaczenie sterowania regulatora
    u(k) = sum(u_part) / sum(w);

    % Błąd średniokwadratowy dla rozmytego algorytmu PID
    e_pid_fuz(k) = e_pid_fuz(k-1) + (yzad(k) - y_pid_fuz(k))^2;
end

E = e_pid_fuz(k_konc);

figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([min(u)-0.1 max(u)+0.1]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm rozmyty PID u(k): Liczba regulatorów: " + ...
    string(n_regulatorow) + " E=" + string(E);
title(title_str);
filenameu = "./pliki_wynikowe/"+"regulator_pid_rozmyty_u(k)"+ ...
    string(o)+".pdf";
export_fig(filenameu);

figure;
stairs(y_pid_fuz); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([min(y_pid_fuz)-0.1 max(y_pid_fuz)+0.1]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm rozmyty PID y(k): Liczba regulatorów: " + ...
    string(n_regulatorow) + " E=" + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
filenamey = "./pliki_wynikowe/"+"regulator_pid_rozmyty_y(k)"+ ...
    string(o)+".pdf";
export_fig(filenamey);