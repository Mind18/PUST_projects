%% Symulacja odpowiedzi skokowych - zadanie 2

% Inicjalizacja macierzy odpowiedzi skokowych

s = cell(wyjscia, wejscia);
for i=1:wyjscia
    for j=1:wejscia
        s{i,j} = zeros(k_konc, 1);
    end
end

figure;
for i=1:wejscia
    subplot(3, 1, i);
end
