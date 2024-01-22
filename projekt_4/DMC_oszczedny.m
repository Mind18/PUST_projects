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

M = zeros(wyjscia*N, wejscia*N_u);
i = 0;
for n=1:size(S{1}, 2):size(M, 2)
    k = 1;
    for m=(size(S{1}, 1)*i)+1:size(S{1}, 1):size(M, 1)
        M(m:m+wyjscia-1, n:n+wejscia-1) = S{k};
        k = k + 1;
    end
    i = i + 1;
end

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