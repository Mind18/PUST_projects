import MinimalWorkingExample.*
clear x y k Yzad;

W1 = 50;

initial_size = 1e4;

% Ograniczenia
du_max = 100;
du_min = -du_max;
u_max = 100;
u_min = 0;

% Punkty pracy
upp = 28;
ypp = 31.7;
y_s = y35;
yzad = ypp;

D = 280;   % Horyzont dynamiki
N = 200;   % Horyzont predykcji
N_u = 200;   % Horyzont sterowania
lambda = 0.1;
Lambda = lambda.*eye(N_u, N_u);
M = zeros(N, N_u);
M_p = zeros(N, D-1);
U_p = zeros(D-1, 1);
e = zeros(initial_size, 1);
e_dmc(1:initial_size) = 0; % Błąd średniokwadratowy dla algorytmu DMC
s = zeros(D, 1);

for k=1:D
    % Realizacja zadania 3
    s(k) = (y_s(k) - ypp) / 7; % (u_konc(skok_u) - upp);
end

for j=1:N_u % dla każdej kolumny macierzy M
    for i=j:N % Dla każdego wiersza kolumny j począwszy od przekątnej
        M(i, j) = y1(i-j+1);
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
        M_p(i, j) = y1(p) - y1(j);
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
    
    % Wyznaczenie zmiany sterowania
    e(k) = yzad - y(k);
    k_j = 0;
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
    
    %% synchronising with the control process
%    plot(x, u);
    plot(x,y , 'b');
    hold on;
    Yzad(k) = yzad;
    E = e_dmc(k);
    plot(x, Yzad, 'r');
    title('E=' + string(E));
    %hold off;

    drawnow;
    
    waitForNewIteration () ; % wait for new iteration
 
    disp(u(k))
    k = k+1;
    
    if k == 200
        yzad = 35;
    end
    
    if k == 600
        yzad = 40;
    end
    
 
end
