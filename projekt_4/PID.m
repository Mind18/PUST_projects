import symulacja_obiektu8y_p4.*

e = cell(1, wyjscia);
e_tmp = 0;
e_pid(1:k_konc) = 0;

% Warunki początkowe
for j=1:wejscia
    u{j} = zeros(k_konc, 1);
end

for i=1:wyjscia
    e{i}(1:k_konc) = 0;
    y{i} = zeros(k_konc, 1);

    % Współczynniki algorytmu
    r2(i) = (K_r(i) * T_d(i)) / T_p;
    r1(i) = K_r(i) * (T_p/(2*T_i(i)) - 2*(T_d(i) / T_p) - 1);
    r0(i) = K_r(i) * (1 + (T_p / (2*T_i(i))) + (T_d(i)/T_p));

    % Generacja zmiennej trajektori
    yzad{i}(1:9)=ypp;
    yzad{i}(10:300)=Y_zad{i}(1);
    yzad{i}(301:600)=Y_zad{i}(2);
    yzad{i}(601:900)=Y_zad{i}(3);
    yzad{i}(901:k_konc)=Y_zad{i}(4);
end

for k=10:k_konc
    % symulacja obiektu
    [y{1}(k), y{2}(k), y{3}(k)] = symulacja_obiektu8y_p4(u{1}(k-1), ...
        u{1}(k-2), u{1}(k-3), u{1}(k-4), u{2}(k-1), u{2}(k-2), u{2}(k-3), ...
        u{2}(k-4), u{3}(k-1), u{3}(k-2), u{3}(k-3), u{3}(k-4), u{4}(k-1), ...
        u{4}(k-2), u{4}(k-3), u{4}(k-4), y{1}(k-1), y{1}(k-2), y{1}(k-3), ...
        y{1}(k-4), y{2}(k-1), y{2}(k-2), y{2}(k-3), y{2}(k-4), y{3}(k-1), ...
        y{3}(k-2), y{3}(k-3), y{3}(k-4));
    
    for i=1:wyjscia
        % Dobierane sterowanie
        ster_u = u_dla_y(i);

        % uchyb regulacji
        e{i}(k)=yzad{i}(k) - y{i}(k);

        % sygnał sterujący regulatora PID
        u{ster_u}(k)=r2(i)*e{i}(k-2)+r1(i)*e{i}(k-1)+ ...
            r0(i)*e{i}(k)+u{ster_u}(k-1);

        % Ograniczenia zmiany sterowania
        du = u{ster_u}(k) - u{ster_u}(k-1);
        if du < du_min
            u{ster_u}(k) = u{ster_u}(k-1) + du_min;
        elseif du > du_max
            u{ster_u}(k) = u{ster_u}(k-1) + du_max;
        end
        % Ograniczenia wartości sterowania
        if u{ster_u}(k) < u_min
            u{ster_u}(k) = u_min;
        elseif u{ster_u}(k) > u_max
            u{ster_u}(k) = u_max;
        end

        e_tmp = e_tmp + (yzad{i}(k) - y{i}(k))^2;
    end

    % Błąd średniokwadratowy dla algorytmu PID
    e_pid(k) = e_pid(k-1) + e_tmp;
    e_tmp = 0;
end

E = e_pid(k_konc);

figure;
hold on;
for i=1:wejscia
    plot(1:k_konc, u{i});
end
title('u(k)');
legend('u_1', 'u_2', 'u_3', 'u_4');
hold off;

figure;
hold on;
for i=1:wyjscia
    plot(1:k_konc, y{i});
    plot(1:k_konc, yzad{i});
end
title('y(k)');
legend('y_1', 'y^{zad}_1', 'y_2', 'y^{zad}_2', 'y_3', 'y^{zad}_3');
hold off;
    