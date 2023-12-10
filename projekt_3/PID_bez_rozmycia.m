%% Algorytm regulacji PID

% Inicjalizacja zmiennych
o = 2; % Nazwa wykresu
% K_r - 0.5 -> 0.35 -> 0.32
% T_i - 9 -> 8 -> 7.3
% T_d - 2 -> 0.2 -> 0.7 -> 0.74
% warunki początkowe
u = zeros(1, k_konc); y_pid = zeros(1, k_konc);
u(1:14)=upp; y_pid(1:14)=ypp;

% Generacja zmiennej trajektori
yzad(1:14)=ypp;
yzad(15:300)=Y_zad(1);
yzad(301:600)=Y_zad(2);
yzad(601:900)=Y_zad(3);
yzad(901:k_konc)=Y_zad(4);

% Generacja trajektorii do testów

% yzad(1:14)=ypp;
% yzad(15:k_konc)=Y_zad(1);

% Generacja trajektorii do testu regulatorów lokalnych

% yzad(1:14)=ypp;
% yzad(15:k_konc)=-1;

% Współczynniki algorytmu
r2 = (K_r * T_d) / T_p;
r1 = K_r * (T_p/(2*T_i) - 2*(T_d / T_p) - 1);
r0 = K_r * (1 + (T_p / (2*T_i)) + (T_d/T_p));

for k=15:k_konc % główna pętla symulacyjna
    % symulacja obiektu
    y_pid(k)=symulacja_obiektu8y_p3(u(k-5), u(k-6), y_pid(k-1), ...
        y_pid(k-2));
    % uchyb regulacji
    e(k)=yzad(k) - y_pid(k);
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
    e_pid(k) = e_pid(k-1) + (yzad(k) - y_pid(k))^2;
end

E = e_pid(k_konc);

% Narysowanie wykresów
figure;
stairs(u); % Dodać wartość błędu średniokwadratowego do tytułu
ylim([min(u)-0.1 max(u)+0.1]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm PID u(k): K_r=" + string(K_r) ...
    + " T_i=" + string(T_i) + " T_D=" + string(T_d) + " E=" ...
    + string(E);
title(title_str);
filenameu = "./pliki_wynikowe/"+"regulator_pid_u(k)"+string(o)+".pdf";
export_fig(filenameu);

figure;
stairs(y_pid); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([min(y_pid)-0.1 max(y_pid)+0.1]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm PID y(k): K_r=" + string(K_r) ...
    + " T_i=" + string(T_i) + " T_D=" + string(T_d) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'best');
filenamey = "./pliki_wynikowe/"+"regulator_pid_y(k)"+string(o)+".pdf";
export_fig(filenamey);