import symulacja_obiektu8y_p4.*

% Alokacja wektora sterowań w odpowiednim rozmiarze
u = zeros(wejscia, k_konc);

% Inicjacja macierzy trajektorii zadanej
yzad = zeros(wyjscia, k_konc);
y = zeros(wyjscia, k_konc);

e = cell(1, wyjscia);
e_tmp = 0;
e_pid(1:k_konc) = 0;

K_r = [0.1 0.02 0.002];
T_i = [0.3 0.01 0.02];
T_d = [0.2 0.1 0.02];

x0 = [K_r(1), T_i(1), T_d(1), K_r(2), T_i(2), T_d(2), K_r(3), T_i(3), T_d(3)];
A = [10, 7, 3, 5, 1, 1, 2, 0.5, 1]; % 22.974, 1, 9.81
b = 10;
lb = [0, 0, 0, 0, 0, 0, 0, 0, 0];
ub = [2, 2, 2, 2, 2, 2, 2, 2, 2];
pid_params = fmincon(@PID_SE, x0, A, b, [], [], lb, ub);

% Warunki początkowe
for j=1:wejscia
    u(j, :) = zeros(k_konc, 1);
end

for i=1:wyjscia
    e{i}(1:k_konc) = 0;
    y(i, :) = zeros(k_konc, 1);

    % Współczynniki algorytmu
    r2(i) = (pid_params(1+(i-1)*3) * pid_params(3+(i-1)*3)) / T_p;
    r1(i) = pid_params(1+(i-1)*3) * (T_p/(2*pid_params(2+(i-1)*3)) - ...
        2*(pid_params(3+(i-1)*3) / T_p) - 1);
    r0(i) = pid_params(1+(i-1)*3) * ...
        (1 + (T_p / (2*pid_params(2+(i-1)*3))) + (pid_params(3+(i-1)*3)/T_p));

    % Generacja zmiennej trajektori
    yzad(i, 1:9)=ypp;
    yzad(i, 10:500)=Y_zad{i}(1);
    yzad(i, 501:1000)=Y_zad{i}(2);
    yzad(i, 1001:1500)=Y_zad{i}(3);
    yzad(i, 1501:k_konc)=Y_zad{i}(4);
end

for k=10:k_konc
    % symulacja obiektu
    [y(1, k), y(2, k), y(3, k)] = symulacja_obiektu8y_p4(u(1, k-1), ...
        u(1, k-2), u(1, k-3), u(1, k-4), u(2, k-1), u(2, k-2), u(2, k-3), ...
        u(2, k-4), u(3, k-1), u(3, k-2), u(3, k-3), u(3, k-4), u(4, k-1), ...
        u(4, k-2), u(4, k-3), u(4, k-4), y(1, k-1), y(1, k-2), y(1, k-3), ...
        y(1, k-4), y(2, k-1), y(2, k-2), y(2, k-3), y(2, k-4), y(3, k-1), ...
        y(3, k-2), y(3, k-3), y(3, k-4));
    
    for i=1:wyjscia
        % Dobierane sterowanie
        ster_u = u_dla_y(i);

        % uchyb regulacji
        e{i}(k)=yzad(i, k) - y(i, k);

        % sygnał sterujący regulatora PID
        u(ster_u, k)=r2(i)*e{i}(k-2)+r1(i)*e{i}(k-1)+ ...
            r0(i)*e{i}(k)+u(ster_u, k-1);

        % Ograniczenia zmiany sterowania
        du = u(ster_u, k) - u(ster_u, k-1);
        if du < du_min
            u(ster_u, k) = u(ster_u, k-1) + du_min;
        elseif du > du_max
            u(ster_u, k) = u(ster_u, k-1) + du_max;
        end
        % Ograniczenia wartości sterowania
        if u(ster_u, k) < u_min
            u(ster_u, k) = u_min;
        elseif u(ster_u, k) > u_max
            u(ster_u, k) = u_max;
        end

        e_tmp = e_tmp + (yzad(i, k) - y(i, k))^2;
    end

    % Błąd średniokwadratowy dla algorytmu PID
    e_pid(k) = e_pid(k-1) + e_tmp;
    e_tmp = 0;
end

E = e_pid(k_konc); % Zapamiętanie błędu średniokwadratowego E symulacji

figure;
hold on;
for i=1:wejscia
    plot(1:k_konc, u(i, :));
end
title('u(k) - PID optymalizacja');
legend('u_1', 'u_2', 'u_3', 'u_4');
hold off;
export_fig("./pliki_wynikowe/uzad.pdf")

figure;
hold on;
for i=1:wyjscia
    plot(1:k_konc, y(i, :));
    plot(1:k_konc, yzad(i, 1:k_konc));
end
title('y(k) - DMC optymalizacja E=' + string(E));
legend('y_1', 'y^{zad}_1', 'y_2', 'y^{zad}_2', 'y_3', 'y^{zad}_3');
hold off;
export_fig("./pliki_wynikowe/yzad.pdf")
    