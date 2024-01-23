%% Symulacja odpowiedzi skokowych - zadanie 2

% Inicjalizacja macierzy odpowiedzi skokowych

s = cell(wyjscia, wejscia); % Odp. skokowe do rysowania wykresu
S = cell(1, D);             % Odp. skokowe do DMC
for i=1:wyjscia
    for j=1:wejscia
        s{i,j} = zeros(k_konc, 1);
    end
end

for k=1:D
    S{k} = zeros(wyjscia, wejscia);
end

figure;
for j=1:wejscia
    if mod(wejscia, 2) ~= 0 && j == wejscia
        subplot(round(wejscia / 2), 2, [j, j+1]);
    else
        subplot(round(wejscia / 2), 2, j)
    end
    u(j, 1:k_konc) = upp;
    u(j, 15:k_konc) = u_konc;

    plot(1:k_konc, u(j, :));
    xlabel('k');
    ylabel('u_' + string(j) + '(k)');
    ylim([-0.1 1.1]);
    u(j, 1:k_konc) = upp;
end
export_fig("./pliki_wynikowe/odp_skokowe_u(k).pdf")

figure;
for n_u=1:wejscia
    u(n_u, 1:k_konc) = upp;
    u(n_u, 15:k_konc) = u_konc;

    for i=1:wyjscia
        y(i, 1:k_konc) = ypp;
    end

    for k=15:k_konc
        [y(1, k), y(2, k), y(3, k)] = symulacja_obiektu8y_p4(u(1, k-1), ...
        u(1, k-2), u(1, k-3), u(1, k-4), u(2, k-1), u(2, k-2), u(2, k-3), ...
        u(2, k-4), u(3, k-1), u(3, k-2), u(3, k-3), u(3, k-4), u(4, k-1), ...
        u(4, k-2), u(4, k-3), u(4, k-4), y(1, k-1), y(1, k-2), y(1, k-3), ...
        y(1, k-4), y(2, k-1), y(2, k-2), y(2, k-3), y(2, k-4), y(3, k-1), ...
        y(3, k-2), y(3, k-3), y(3, k-4));
        
        if k > 15
            for n_y=1:wyjscia
                S{k-15}(n_y, n_u) = y(n_y, k);
            end
        end
    end

    for n_y=1:wyjscia
        s{n_y, n_u} = y(n_y, :);
        subplot(wejscia, wyjscia, wyjscia*(n_u-1) + n_y);
        plot(1:k_konc, s{n_y, n_u});
        xlabel('k');
        ylabel('y_{' + string(wyjscia*(n_u-1) + n_y) + '}^{' + string(n_y) + ...
            ',' + string(n_u) + '}');
        ylim padded;
    end

    u(n_u, 1:k_konc) = upp;
end
export_fig("./pliki_wynikowe/odp_skokowe_s(k).pdf")

