import symulacja_obiektu8y_p4.*

% Inicjalizacja danych
T_p = 0.5;
print = 'Y'; % Czy generować wykresy

% Wektory zmiennych procesowych

%% Sprawdzanie punktu pracy - zadanie 1

for k=8:k_konc
    [y(1, k), y(2, k), y(3, k)] = symulacja_obiektu8y_p4(u(1, k-1), ...
        u(1, k-2), u(1, k-3), u(1, k-4), u(2, k-1), u(2, k-2), u(2, k-3), ...
        u(2, k-4), u(3, k-1), u(3, k-2), u(3, k-3), u(3, k-4), u(4, k-1), ...
        u(4, k-2), u(4, k-3), u(4, k-4), y(1, k-1), y(1, k-2), y(1, k-3), ...
        y(1, k-4), y(2, k-1), y(2, k-2), y(2, k-3), y(2, k-4), y(3, k-1), ...
        y(3, k-2), y(3, k-3), y(3, k-4));
end

% Narysowanie wykresów

if strcmp(print, 'Y')
    figure;
    subplot(2, 2, 1);
    plot(1:k_konc, u(1, :));
    hold on;
    title("Wykres u_1(k)");
    xlabel('k');
    ylabel('u_1(k)');
    ylim([-0.2 0.2]);

    subplot(2, 2, 2);
    plot(1:k_konc, u(2, :));
    hold on;
    title("Wykres u_2(k)");
    xlabel('k');
    ylabel('u_2(k)');
    ylim([-0.2 0.2]);

    subplot(2, 2, 3);
    plot(1:k_konc, u(3, :));
    hold on;
    title("Wykres u_3(k)");
    xlabel('k');
    ylabel('u_3(k)');
    ylim([-0.2 0.2]);

    subplot(2, 2, 4);
    plot(1:k_konc, u(4, :));
    hold on;
    title("Wykres u_4(k)");
    xlabel('k');
    ylabel('u_4(k)');
    ylim([-0.2 0.2]);

    export_fig("./pliki_wynikowe/test_punktu_pracy_u.pdf");

    figure;
    subplot(3, 1, 1);
    plot(1:k_konc, y(1, :), 'Color', "#D95319");
    hold on;
    title("Wykres y_1(k)");
    xlabel('k');
    ylabel('y_1(k)');
    ylim([-0.2 0.2]);

    subplot(3, 1, 2);
    plot(1:k_konc, y(2, :), 'Color', "#D95319");
    hold on;
    title("Wykres y_2(k)");
    xlabel('k');
    ylabel('y_2(k)');
    ylim([-0.2 0.2]);

    subplot(3, 1, 3);
    plot(1:k_konc, y(3, :), 'Color', "#D95319");
    hold on;
    title("Wykres y_3(k)");
    xlabel('k');
    ylabel('y_3(k)');
    ylim([-0.2 0.2]);

    export_fig("./pliki_wynikowe/test_punktu_pracy_y.pdf");
end