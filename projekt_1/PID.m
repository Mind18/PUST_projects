%% Algorytm regulacji PID
k_konc = 1200;
Y_zad = [y_zad 1.8 1.3 2.3];

% Inicjalizacja zmiennych
o = 2; % Nazwa wykresu
K_r = 0.35; T_p = 0.5; T_i = 7.3; T_d = 0.74; 
% K_r - 0.5 -> 0.35 -> 0.32
% T_i - 9 -> 8 -> 7.3
% T_d - 2 -> 0.2 -> 0.7 -> 0.74
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
filenameu = "./PIDDMC/"+"regulator_pid_u(k)"+string(o)+".pdf";
export_fig(filenameu);

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
filenamey = "./PIDDMC/"+"regulator_pid_y(k)"+string(o)+".pdf";
export_fig(filenamey);