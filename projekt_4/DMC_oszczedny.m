import symulacja_obiektu8y_p4.*

e = cell(1, wyjscia);
e_tmp = 0;
e_dmc(1:k_konc) = 0;

% Generowanie macierzy wagowych
Psi = zeros(wyjscia*N, wyjscia*N);
psi_inputed = 1;
for i=1:wyjscia*N
    Psi(i, i) = psi(psi_inputed);
    if psi_inputed == wyjscia
        psi_inputed = 1;
    else
        psi_inputed = psi_inputed + 1;
    end
end

Lambda = zeros(wejscia*N_u, wejscia*N_u);
lambda_inputed = 1;
for i=1:wejscia*N_u
    Lambda(i, i) = lambda(lambda_inputed);
    if lambda_inputed == wyjscia
        lambda_inputed = 1;
    else
        lambda_inputed = lambda_inputed + 1;
    end
end

M_cell = cell(N, N_u);
for j=1:N_u
    for i=1:N
        M_cell{i, j} = zeros(wyjscia, wejscia);
    end
end

for n=1:N_u
    for m=n:N
        M_cell{m, n} = S{m-n+1};
    end
end
M = cell2mat(M_cell);

M_p_cell = cell(N, (D-1));
for j=1:(D-1)
    for i=1:N
        M_p_cell{i, j} = zeros(wyjscia, wejscia);
    end
end

for j=1:(D-1)
    for i=1:N
        if j+i > D
            p = D;
        else
            p = j+i;
        end
        M_p_cell{i, j} = S{p} - S{j};
    end
end
M_p = cell2mat(M_p_cell);

for i=1:wyjscia
    e{i}(1:k_konc) = 0;
    y{i} = zeros(k_konc, 1);

    % Generacja zmiennej trajektori
    yzad{i}(1:9)=ypp;
    yzad{i}(10:300)=Y_zad{i}(1);
    yzad{i}(301:600)=Y_zad{i}(2);
    yzad{i}(601:900)=Y_zad{i}(3);
    yzad{i}(901:k_konc)=Y_zad{i}(4);
end