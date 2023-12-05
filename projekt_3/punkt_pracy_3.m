import symulacja_obiektu8y_p2.*

% Inicjalizacja danych
T_p = 0.5;
print = 'Y'; % Czy generować wykresy

% Wektory zmiennych procesowych
u = zeros(k_konc, 1); % Sterowanie
y = zeros(k_konc, 1); % Wyjście procesu

%% Sprawdzenie punktu pracy - zadanie 1

for k=8:k_konc
    y(k) = symulacja_obiektu8y_p3(u(k-5), u(k-6), y(k-1), y(k-2));
end

% Narysowanie wykresu

if strcmp(print, 'Y')
    figure;
    subplot(2, 1, 1)
    plot(1:k_konc, u);
    hold on;
    title("Wykres u(k)");
    xlabel('k');
    ylabel('u(k)');
    ylim([-0.2 0.2]);
    
    subplot(2, 1, 2)
    plot(1:k_konc, y, Color="#EDB120");
    title("Wykres y(k)");
    xlabel('k');
    ylabel('y(k)');
    ylim([-0.2 0.2]);
    export_fig("./pliki_wynikowe/test_punktu_pracy.pdf");
end
