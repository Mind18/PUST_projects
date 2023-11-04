function [sqrd_err] = PID_SE(pid_param)

% Stałe funkcji
k_konc = 400;
T_p = 0.5;
upp = 1;
ypp = 1.7;
y_zad = 1.2;

% Ograniczenia regulatora
du_max = 1; % 1
du_min = -du_max;
u_min = 0.5; % 0.5
u_max = 1.5; % 1.5

% warunki początkowe
u = zeros(1, k_konc); y = zeros(1, k_konc);
u(1:11)=upp; y(1:11)=ypp;
yzad(1:11)=0; yzad(12:k_konc)=y_zad;
e(1:k_konc)=0; e_pid(1:k_konc) = 0;

% Współczynniki algorytmu PID
r2 = (pid_param(1) * pid_param(3)) / T_p;
r1 = pid_param(1) * (T_p/(2*pid_param(2)) - 2*(pid_param(3) / T_p) - 1);
r0 = pid_param(1) * (1 + (T_p / (2*pid_param(2))) + (pid_param(3)/T_p));

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

sqrd_err = e_pid(k_konc);
