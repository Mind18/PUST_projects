clear x y k Yzad;
k = 1;
initial_size = 1e4;

% Ograniczenia
du_max = 100;
du_min = -du_max;
u_max = 100;
u_min = 0;

upp = 28;
ypp = 31.7;
yzad = ypp;

% Inicjalizacja zmiennych
K_r = 0.28173; T_p = 1; T_i = 18.8296; T_d = -3.5362;

du = 0;
e = zeros(initial_size, 1);
e_pid(1:initial_size) = 0; % Błąd średniokwadratowy dla algorytmu DMC

% Współczynniki algorytmu
r2 = (K_r * T_d) / T_p;
r1 = K_r * (T_p/(2*T_i) - 2*(T_d / T_p) - 1);
r0 = K_r * (1 + (T_p / (2*T_i)) + (T_d/T_p));

while (1)
    addpath ('D:\SerialCommunication') ; % add a path
    initSerialControl COM7 % initialise com port
    % Uzyskanie sygnału wyjściowego
    measurements = readMeasurements (1);
    y(k) = measurements;
    x(k) = k;
    e(k) = yzad - y(k);
    % sygnał sterujący regulatora PID
    if k==1
        u(k)=r0*e(k)+upp;
    elseif k==2
        u(k)=r1*e(k-1)+r0*e(k)+u(k-1);
    else
        u(k)=r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
    end
    % Ograniczenia zmiany sterowania
    if k > 1
       du = u(k) - u(k-1);
    end
    if (du < du_min) && (k > 1)
        u(k) = u(k-1) + du_min;
    elseif (du > du_max) && (k > 1)
        u(k) = u(k-1) + du_max;
    end
    % Ograniczenia wartości sterowania
    if u(k) < u_min
        u(k) = u_min;
    elseif u(k) > u_max
        u(k) = u_max;
    end

    % Błąd średniokwadratowy dla algorytmu PID
    if k == 1
        e_pid(k) = (e(k))^2;
    else
        e_pid(k) = e_pid(k-1) + (e(k))^2;
    end
    sendControls ([ 1 , 2, 3, 4, 5, 6] , [W1, 0 , 0 , 0 ,u(k), 0]) ;
    
    %% synchronising with the control process
%    plot(x, u);
    plot(x,y , 'b');
    hold on;
    Yzad(k) = yzad;
    E = e_pid(k);
    plot(x, Yzad, 'r');
    title('E=' + string(E));
    %hold off;

    drawnow;
    
    waitForNewIteration () ; % wait for new iteration
 
    disp(u(k))
    k = k+1;
    
    if k == 300
        yzad = 35;
    end
    
%     if k == 400
%         yzad = 40;
%     end
end