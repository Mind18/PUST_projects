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
    u{j}(1:k_konc, 1) = upp;
end
export_fig("./pliki_wynikowe/odp_skokowe_u(k).pdf")

figure;
for j=1:wejscia
    u{j}(1:k_konc, 1) = upp;
    u{j}(15:k_konc, 1) = u_konc;

    for i=1:wyjscia
        y{i}(1:k_konc) = ypp;
    end

    for k=15:k_konc
        [y{1}(k), y{2}(k), y{3}(k)] = symulacja_obiektu8y_p4(u{1}(k-1), ...
        u{1}(k-2), u{1}(k-3), u{1}(k-4), u{2}(k-1), u{2}(k-2), u{2}(k-3), ...
        u{2}(k-4), u{3}(k-1), u{3}(k-2), u{3}(k-3), u{3}(k-4), u{4}(k-1), ...
        u{4}(k-2), u{4}(k-3), u{4}(k-4), y{1}(k-1), y{1}(k-2), y{1}(k-3), ...
        y{1}(k-4), y{2}(k-1), y{2}(k-2), y{2}(k-3), y{2}(k-4), y{3}(k-1), ...
        y{3}(k-2), y{3}(k-3), y{3}(k-4));
    end

    for i=1:wyjscia
        s{i, j} = y{i};
        subplot(wejscia, wyjscia, wyjscia*(j-1) + i);
        plot(1:k_konc, s{i, j});
        xlabel('k');
        ylabel('y_{' + string(wyjscia*(j-1) + i) + '}^{' + string(i) + ...
            ',' + string(j) + '}');
        ylim padded;
    end

    u{j}(1:k_konc, 1) = upp;
end
export_fig("./pliki_wynikowe/odp_skokowe_s(k).pdf")

