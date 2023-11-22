import sendControlsToG1AndDisturbance.*
clear x y k Yzad Z;

W1 = 50;

initial_size = 1e4;

% Ograniczenia
du_max = 100;
du_min = -du_max;
u_max = 100;
u_min = 0;

% Punkty pracy
upp = 28;
ypp = 32;
yzad = ypp;
z = 0;

D = 300;   % Horyzont dynamiki
D_z = 325; % Horyzont dynamiki zakłóceń
N = 200;   % Horyzont predykcji
N_u = 170;   % Horyzont sterowania
lambda = 0.002;
Lambda = lambda.*eye(N_u, N_u);
M = zeros(N, N_u);
M_p = zeros(N, D-1);
M_p_z = zeros(N, D_z-1);
Z_p = zeros(D_z, 1);
U_p = zeros(D-1, 1);
e = zeros(initial_size, 1);
e_dmc(1:initial_size) = 0; %    Błąd średniokwadratowy dla algorytmu DMC

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

% Generacja macierzy M_p_z
for j=1:D_z % dla każdej kolumny macierzy M_p
    for i=1:N % dla każdego wiersza macierzy M_p
        if j+i > D_z
            p = D_z;
        else
            p = j+i-1;
        end
        if j == 1
            M_p_z(i, j) = s_z(p);
        else
            M_p_z(i, j) = s_z(p) - s_z(j);
        end
    end
end

% Wyznaczenie wektora współczynników K
K = ((M'*M+Lambda)^(-1))*M';

% Wyznaczenie współczynnika k_e
k_e = sum(K(1,:));

u = zeros(initial_size,1);
k = 1;
while (1)
    addpath ('D:\SerialCommunication') ; % add a path
    initSerialControl COM7 % initialise com port
    % Dokonanie pomiaru sygnału wyjściowego
    measurements = readMeasurements (1);
    y(k) = measurements;
    x(k) = k;
    Z(k) = z;
    % Aktualizacja wektora zmian zakłócenia
    for n=D_z:-1:2
        Z_p(n,1) = Z_p(n-1, 1);
    end
    
    if k == 1
        Z_p(1) = Z(k);
    else
        Z_p(1) = Z(k) - Z(k-1);
    end
    % Wyznaczenie zmiany sterowania
    e(k) = yzad - y(k);
    k_j = 0;
    k_j_z = 0;
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
    if(k == 1)
       u(k) = upp + delta_u; 
    else
       u(k) = u(k-1) + delta_u;
    end
    % Ograniczenia wartości sterowania
   
    if u(k) < u_min
        u(k) = u_min;
        U_p(1) = u(k)-u(k-1);
    end
    if u(k) > u_max
        u(k) = u_max;
        U_p(1) = u(k)-u(k-1);
    end
    
    if(k == 1)
       e_dmc(k) = e(k)^2; 
    else
       e_dmc(k) = e_dmc(k-1) + e(k)^2;
    end
    
    sendControls ([ 1 , 2, 3, 4, 5, 6] , [W1, 0 , 0 , 0 ,u(k), 0]) ;
    sendControlsToG1AndDisturbance(u(k), Z);
    
    %% synchronising with the control process
%    plot(x, u);
    subplot(2,1,1);
    plot(x,y , 'Color', [0 0.4470 0.7410]);
    hold on;
    Yzad(k) = yzad;
    Z(k) = z;
    E = e_dmc(k);
    plot(x, Yzad, 'Color', [0.8500 0.3250 0.0980]);
    title('y(k)/z(k) E=' + string(E));
    legend('y(k)', 'y^{zad}(k)', 'Location', 'southeast');
    subplot(2,1,2);
    plot(x, Z, 'Color', [0.4660 0.6740 0.1880]);
    legend('z(k)');
    %hold off;

    drawnow;
    
    waitForNewIteration () ; % wait for new iteration
 
    disp(u(k))
    k = k+1;
    
    if k == 20
        yzad = 35;
    end
    
    if k == 470
        z = 25;
    end
    
    if k == 920
       z = 15; 
    end
%     if k == 600
%         yzad = 40;
%     end
    
 
end