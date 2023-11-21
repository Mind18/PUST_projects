%% Algorytm regulacji DMC

% warunki początkowe
u = zeros(1, k_konc); z = zeros(1,k_konc); y = zeros(1, k_konc);
u(1:k_konc)=upp; y(1:11)=ypp;
% Generacja zmiennej trajektori
yzad(1:11)=ypp;
yzad(12:k_konc)=Y_zad(1);
z(91:k_konc) = 1;
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

% % Generacja macierzy M_p_z
% for j=1:D_z % dla każdej kolumny macierzy M_p
%     for i=1:N % dla każdego wiersza macierzy M_p
%         if j+i > D_z
%             p = D_z;
%         else
%             p = j+i-1;
%         end
%         if j == 1
%             M_p_z(i, j) = s_z(p);
%         else
%             M_p_z(i, j) = s_z(p) - s_z(j);
%         end
%     end
% end

M_p_z=zeros(N,D_z);
for i=1:(D_z)
    M_p_z(i:N,i)=s(1:N-i+1);
end

% Wyznaczenie wektora współczynników K
K = ((M'*M+Lambda)^(-1))*M';

% Wyznaczenie współczynnika k_e
k_e = sum(K(1,:));

for k=12:k_konc
    % symulacja obiektu
    if k==91
        k;
    end
    y(k) = symulacja_obiektu8y_p2(u(k-6), u(k-7), z(k-1), z(k-2), ...
        y(k-1), y(k-2));
    % Aktualizacja wektora zmian zakłócenia
    for n=D_z:-1:2
        Z_p(n,1) = Z_p(n-1, 1);
    end
    Z_p(1) = z(k) - z(k-1);
    % wyznaczenie zmiany sterowania
    k_j = 0;
    k_j_z = 0;
    e(k) = yzad(k) - y(k);
    for j=1:D-1
        k_j_inc = K(1, :)*M_p(:, j);
        k_j = k_j + k_j_inc*U_p(j);
    end
    for j=0:D_z-1
        k_j_z_inc = K(1, :)*M_p_z(:, j+1);
        k_j_z = k_j_z + k_j_z_inc*Z_p(j+1);
    end
    delta_u = k_e*e(k) - k_j - k_j_z;
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
ylim([0 1.6]);
xlabel('k');
ylabel('u(k)');
title_str = "Algorytm DMC u(k): N=" + string(N) ...
    + " N_u=" + string(N_u) + " λ=" + string(lambda) + " E=" ...
    + string(E);
title(title_str);
export_fig('./pliki_wynikowe/regulator_dmc_komp_u(k)_zak.pdf');

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
export_fig('./pliki_wynikowe/regulator_dmc_komp_y(k)_zak.pdf');