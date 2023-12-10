%% Algorytm regulacji DMC
o = 20; % Nazwa wykresu

% warunki początkowe
u = zeros(1, k_konc); y_dmc = zeros(1, k_konc);
u(1:k_konc)=upp; y_dmc(1:11)=ypp;
s_i = s{odp_skok};
e = zeros(1, k_konc);

% Generacja zmiennej trajektori
% yzad(1:11)=ypp;
% yzad(15:300)=Y_zad(1);
% yzad(301:600)=Y_zad(2);
% yzad(601:900)=Y_zad(3);
% yzad(901:k_konc)=Y_zad(4);

% yzad(1:14)=ypp;
% yzad(15:k_konc)=Y_zad(1);

% Generacja trajektorii do testu regulatorów lokalnych

yzad(1:14)=ypp;
yzad(15:k_konc)=-0.06;

% Generacja macierzy M
for j=1:N_u % dla każdej kolumny macierzy M
    for i=j:N % Dla każdego wiersza kolumny j począwszy od przekątnej
        M(i, j) = s_i(i-j+1);
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
        M_p(i, j) = s_i(p) - s_i(j);
    end
end

% Wyznaczenie wektora współczynników K
K = ((M'*M+Lambda)^(-1))*M';

% Wyznaczenie współczynnika k_e
k_e = sum(K(1,:));

for k=15:k_konc

    if k == 17
        u;
    end

    % symulacja obiektu
    y_dmc(k) = symulacja_obiektu8y_p3(u(k-5), u(k-6), y_dmc(k-1), ...
        y_dmc(k-2));
    % wyznaczenie zmiany sterowania
    k_j = 0;
    e(k) = yzad(k) - y_dmc(k);
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
ylim([min(u)-0.1 max(u)+0.1]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm DMC u(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
filenameu = "./pliki_wynikowe/"+"regulator_dmc_bez_rozmycia_u(k)"+ ...
    string(o)+".pdf";
export_fig(filenameu);

figure;
stairs(y_dmc); % Dodać wartość błędu średniokwadratowego do tytułu
hold on;
stairs(yzad, ':');
ylim([min(y_dmc)-0.1 max(y_dmc)+0.1]);
xlabel('k');
ylabel('y(k)');
title_str = "Algorytm DMC y(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
legend('y(k)', 'y^{zad}', 'Location', 'southeast');
filenamey = "./pliki_wynikowe/"+"regulator_dmc_bez_rozmycia_y(k)"+ ...
    string(o)+".pdf";
export_fig(filenamey);