%% Symulacja odpowiedzi skokowej - zadanie 2

% Wektory dla przechowywania danych statycznych
u_stat = zeros(1, size(u_konc, 2)+1);
u_stat(1) = upp;
u_stat(2:size(u_konc, 2)+1) = u_konc;
y_stat = zeros(1, size(u_konc, 2)+1);
y_stat(1) = ypp;
s = {};

% Narysowanie wykresów sygnałów
figure;
for i=1:size(u_konc, 2)
    u(1:k_konc, 1) = upp; % U_pp=0
    u(15:k_konc, 1) = u_konc(i);
    plot(1:k_konc, u);
    ylim([-0.1 2.5]);
    hold on;
end

% Dodanie informacji do wygenerowanego wykresu
legend_list = strings([1 size(u_konc, 2)]);
xlabel('k');
ylabel('u(k)');
ylim([-1.1 1.1]);
title('Wykres u(k)');
for i=1:size(u_konc, 2)
    legend_text = 'u(k)=' + string(u_konc(i));
    legend_list(i) = legend_text;
end
legend(legend_list, 'Location', 'best');
hold off;
export_fig('./pliki_wynikowe/zad2_u(k).pdf');

% Narysowanie wykresu y(k) odp. skokowej - tor wejście-wyjście

figure;
for i=1:size(u_konc, 2)
    u(1:14) = upp;
    u(15:k_konc) = u_konc(i);
    y(1:k_konc) = ypp;

    for k=12:k_konc
        y(k) = symulacja_obiektu8y_p3(u(k-5), u(k-6), y(k-1), y(k-2));
    end

    y_stat(i+1) = y(k_konc); % zapis danych dla wektora statycznego y
    for k=1:D
        if u_konc(i) ~= upp
            s{i}(k) = (y(k+15) - ypp) / u_konc(i) - upp;
        else
            s{i}(k) = 0;
        end
    end

    plot(1:k_konc, y);
    hold on;
end

% Dodanie informacji do wygenerowanego wykresu
legend_list = strings([1 size(u_konc, 2)]);
xlabel('k');
ylabel('y(k)');
title('Wykres y(k)');
for i=1:size(u_konc, 2)
    legend_text = 'y(k)=' + string(y_stat(i+1));
    legend_list(i) = legend_text;
end
legend(legend_list, 'Location', 'southeast');
hold off;
export_fig('./pliki_wynikowe/zad2_y(k).pdf');

% Ch-ka statyczna

figure;
[u_stat, sortIndex] = sort(u_stat);
y_stat = y_stat(sortIndex);

plot(u_stat, y_stat);
hold on;
plot(u_stat, y_stat, '.', 'MarkerSize',12);
legend("Interpolacja danych statycznych", "Dane statyczne", 'Location', 'southeast');
k_stat = rdivide(y_stat, u_stat);
xlabel('u');
ylabel('y(u)');
title('Dane statyczne y(u)');
hold off;
export_fig('./pliki_wynikowe/zad2_y(u)_stat.pdf');
