%% Symulacja odpowiedzi skokowych - zadanie 2

% Inicjalizacja macierzy odpowiedzi skokowych

s = cell(wyjscia, wejscia);
for i=1:wyjscia
    for j=1:wejscia
        s{i,j} = zeros(k_konc, 1);
    end
end

figure;
for j=1:wejscia
    if mod(wejscia, 2) ~= 0 && j == wejscia
        subplot(round(wejscia / 2), 2, [j, j+1]);
    else
        subplot(round(wejscia / 2), 2, j)
    end
    u{j}(1:k_konc, 1) = upp;
    u{j}(15:k_konc, 1) = u_konc;
    plot(1:k_konc, u{j});
    xlabel('k');
    ylabel('u_' + string(j) + '(k)');
    ylim([-0.1 1.1]);
end
export_fig("./pliki_wynikowe/odp_skokowe_u(k).pdf")

figure;
for j=1:wejscia
    for i=1:wyjscia
        y{i}(1:k_konc) = ypp;
    end

    for k=15:k_konc
        
    end
end

