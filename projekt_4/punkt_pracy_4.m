import symulacja_obiektu8y_p4.*

% Inicjalizacja danych
T_p = 0.5;
print = 'Y'; % Czy generować wykresy

% Wektory zmiennych procesowych
u_1 = zeros(k_konc, 1); % Sterowanie
u_2 = zeros(k_konc, 1);
u_3 = zeros(k_konc, 1);
u_4 = zeros(k_konc, 1);
y_1 = zeros(k_konc, 1); % Wyjście procesu
y_3 = zeros(k_konc, 1);
y_2 = zeros(k_konc, 1);

%% Sprawdzanie punktu pracy - zadanie 1

for k=8:k_konc
    [y_1(k), y_2(k), y_3(k)] = symulacja_obiektu8y_p4(u_1(k-1), ...
        u_1(k-2), u_1(k-3), u_1(k-4), u_2(k-1), u_2(k-2), u_2(k-3), ...
        u_2(k-4), u_3(k-1), u_3(k-2), u_3(k-3), u_3(k-4), u_4(k-1), ...
        u_4(k-2), u_4(k-3), u_4(k-4), y_1(k-1), y_1(k-2), y_1(k-3), ...
        y_1(k-4), y_2(k-1), y_2(k-2), y_2(k-3), y_2(k-4), y_3(k-1), ...
        y_3(k-2), y_3(k-3), y_3(k-4));
end

% Narysowanie wykresów

if strcmp(print, 'Y')
    figure;
    subplot(2, 2, 1);
    plot(1:k_konc, u_1);
    hold on;
    title("Wykres u_1(k)");
    xlabel('k');
    ylabel('u_1(k)');
    ylim([-0.2 0.2]);

    subplot(2, 2, 2);
    plot(1:k_konc, u_2);
    hold on;
    title("Wykres u_2(k)");
    xlabel('k');
    ylabel('u_2(k)');
    ylim([-0.2 0.2]);

    subplot(2, 2, 3);
    plot(1:k_konc, u_3);
    hold on;
    title("Wykres u_3(k)");
    xlabel('k');
    ylabel('u_3(k)');
    ylim([-0.2 0.2]);

    subplot(2, 2, 4);
    plot(1:k_konc, u_4);
    hold on;
    title("Wykres u_4(k)");
    xlabel('k');
    ylabel('u_4(k)');
    ylim([-0.2 0.2]);

    export_fig("./pliki_wynikowe/test_punktu_pracy_u.pdf");

    figure;
    subplot(3, 1, 1);
    plot(1:k_konc, y_1, 'Color', "#D95319");
    hold on;
    title("Wykres y_1(k)");
    xlabel('k');
    ylabel('y_1(k)');
    ylim([-0.2 0.2]);

    subplot(3, 1, 2);
    plot(1:k_konc, y_2, 'Color', "#D95319");
    hold on;
    title("Wykres y_2(k)");
    xlabel('k');
    ylabel('y_2(k)');
    ylim([-0.2 0.2]);

    subplot(3, 1, 3);
    plot(1:k_konc, y_3, 'Color', "#D95319");
    hold on;
    title("Wykres y_3(k)");
    xlabel('k');
    ylabel('y_3(k)');
    ylim([-0.2 0.2]);

    export_fig("./pliki_wynikowe/test_punktu_pracy_y.pdf");
end