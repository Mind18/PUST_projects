clear x y k Yzad;

kryterium = 'u';
n_regulatorow = 3;

k = 1;
initial_size = 1e4;

reg_part = {[0, 0, 40, 45], [40, 45, 60, 65], ...
    [60, 65, 100, 100]};
u_w = {};
s = {y30, y40, y50};

u_part = [];

% Ograniczenia
du_max = 100;
du_min = -du_max;
u_max = 100;
u_min = 0;

upp = 28;
ypp = 31.7;
yzad = ypp;

% Inicjalizacja zmiennych
D_fuz = [566 760 431]; % Horyzonty dynamiki
N_fuz = [566 760 431];   % Horyzonty predykcji
N_u_fuz = [566 76 431];   % Horyzonty sterowania
lambda_fuz = [1 1 1]; % Parametry lambda lokalnych 
                                    % regulatorów
k_j_fuz = zeros(1 ,n_regulatorow);

for r=1:n_regulatorow
    U_p_fuz{r} = zeros(D_fuz(r)-1, 1);;
end

du = 0;
e = zeros(initial_size, 1);
e_dmc(1:initial_size) = 0; % Błąd średniokwadratowy dla algorytmu DMC

% Generacje macierzy Lambda
for r=1:n_regulatorow
    Lambda_fuz{r} = lambda_fuz(r).*eye(N_u_fuz(r), N_u_fuz(r));
end

% Generacja macierzy M
for r=1:n_regulatorow
    for j=1:N_u_fuz(r) % dla każdej kolumny macierzy M
        for i=j:N_fuz(r) % Dla każdego wiersza kolumny j począwszy od przekątnej
            M_fuz{r}(i, j) = s{r}(i-j+1);
        end
    end
end

% Generacja macierzy M_p
for r=1:n_regulatorow
    for j=1:D_fuz(r)-1 % dla każdej kolumny macierzy M_p
        for i=1:N_fuz(r) % dla każdego wiersza macierzy M_p
            if j+i > D_fuz(r)
                p = D_fuz(r);
            else
                p = j+i;
            end
            M_p_fuz{r}(i, j) = s{r}(p) - s{r}(j);
        end
    end
end

% Wyznaczenie lokalnych wektorów współczynników K
for r=1:n_regulatorow
    K_fuz{r} = ((M_fuz{r}'*M_fuz{r}+Lambda_fuz{r})^(-1))*M_fuz{r}';
end

% Wyznaczenie współczynnika k_e
for r=1:n_regulatorow
    k_e_fuz(r) = sum(K_fuz{r}(1,:));
end

while (1)
    addpath ('D:\SerialCommunication') ; % add a path
    initSerialControl COM7 % initialise com port
    % Uzyskanie sygnału wyjściowego
    measurements = readMeasurements (1);
    y(k) = measurements;
    x(k) = k;
    e(k) = yzad - y(k);
    % sygnał sterujący regulatora PID
    
    % Wybranie kryterium dla funkcji przynależności
    if strcmp(kryterium, 'u') % Jeżeli kryterium jest sterowanie
        if k~=1
            x_w = u(k-1); % Wykorzystaj poprzednią wartość sterowania
        else
            x_w = upp;
        end
    elseif strcmp(kryterium, 'y') % Jeżeli kryterium jest sygnał wyjściowy
        x_w = y(k); % Wykorzystaj aktualną wartość sygnału wyjściowego
    end
    
    for i=1:n_regulatorow % dla każdego z regulatorów lokalnych
        % Wyznaczenie wartości funkcji przynależności
        w(i) = trapezoid_assign(x_w, reg_part{i});

        for j=1:D_fuz(i)-1
            k_j_inc = K_fuz{i}(1, :)*M_p_fuz{i}(:, j);
            k_j_fuz(i) = k_j_fuz(i) + k_j_inc*U_p_fuz{i}(j);
        end
        % sygnał sterujący regulatora lokalnego PID
        delta_u = k_e_fuz(i)*e(k) - k_j_fuz(i);
        % Ograniczenia zmiany sterowania
        if delta_u < du_min
            delta_u = du_min;
        elseif delta_u > du_max
            delta_u = du_max;
        end
        % Zapamiętanie zmiany sterowania do kolejnych iteracji
        for n=D_fuz(i)-1:-1:2
            U_p_fuz{i}(n,1) = U_p_fuz{i}(n-1,1);
        end
        U_p_fuz{r}(1) = delta_u;
        % Dokonanie zmiany sterowania
        if k ~= 1
            u_w{i}(k) = u(k-1) + delta_u;
        else
            u_w{i}(k) = upp + delta_u;
        end
        % Ograniczenia wartości sterowania
        if u_w{i}(k) < u_min
            u_w{i}(k) = u_min;
            U_p_fuz{i}(1) = u_w{r}(k)-u(k-1);
        elseif u_w{i}(k) > u_max
            u_w{i}(k) = u_max;
            U_p_fuz{i}(1) = u_w{r}(k)-u(k-1);
        end

        u_part(i) = w(i)*u_w{i}(k);
    end
    
    % Wyznaczenie sterowania regulatora
    u(k) = sum(u_part) / sum(w);

    % Błąd średniokwadratowy dla algorytmu PID
    if k == 1
        e_dmc(k) = (e(k))^2;
    else
        e_dmc(k) = e_dmc(k-1) + (e(k))^2;
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
        yzad = 38.81;
    end
    
    if k == 600
        yzad = 48.81;
    end
    
    if k == 1000
       yzad = 33.81;
    end
end