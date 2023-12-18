%% Symulacja odpowiedzi skokowych - zadanie 2

% Inicjalizacja macierzy odpowiedzi skokowych

s = cell(3, 4);
for i=1:3
    for j=1:4
        s{i,j} = zeros(k_konc, 1);
    end
end

