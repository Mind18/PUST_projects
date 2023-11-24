%% Algorytm regulacji DMC

%tryb zakłócenia
% 3 oznacza brak zakłócenia
% 2 oznacza skok
% 1 oznacza zakłócenie z sinusem
% 0 oznacza zakłócenie szumem
% zaklocenie = 2; 

% warunki początkowe
u = zeros(1, k_konc); z = zeros(1,k_konc); y = zeros(1, k_konc);
u(1:k_konc)=upp; y(1:11)=ypp;
% Generacja zmiennej trajektori
yzad(1:11)=ypp;
yzad(12:k_konc)=Y_zad(1);

%% Zakłócenia
for i=90:k_konc
    if (zaklocenie==0) %szum
        % Do uzyskania zakresu szumu <0.85, 1>
        % z(i) = 0.85 + (0.15)*rand(1,1);
        % Do uzyskania następnych
        z(i) = -1 + (2)*rand(1, 1);
    elseif (zaklocenie==1) %sinusoidalne
        zakres_sin = 0:0.1:13*pi;
        a = 0.05*sin(zakres_sin);
        z(i) = a(i);
    elseif (zaklocenie==2) %skok
        z(i)=1;
    elseif (zaklocenie==3) %brak zakłócenia
        z(i)=0;
    end
end

% yzad(301:600)=Y_zad(2);
% yzad(601:900)=Y_zad(3);
% yzad(901:k_konc)=Y_zad(4);

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
    y(k) = symulacja_obiektu8y_p2(u(k-6), u(k-7), z(k-1), z(k-2), ...
        y(k-1), y(k-2));
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
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm DMC u(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda)+ " E=" ...
    + string(E);
title(title_str);
export_fig('./pliki_wynikowe/regulator_dmc_u(k)_zak.pdf');

figure;
stairs(y); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([0 2.4]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm DMC y(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
export_fig('./pliki_wynikowe/regulator_dmc_y(k)_zak.pdf');