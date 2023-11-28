%% Symulacja odpowiedzi skokowej - zadanie 2

u_konc = [0.3 0.7 1 1.2 1.5];
z_konc = [0.3 0.5 1 1.2 1.5];

% Wektory dla przechowywania danych statycznych
u_stat = zeros(1, size(u_konc, 2)+1);
u_stat(1) = upp;
u_stat(2:size(u_konc, 2)+1) = u_konc;
y_stat = zeros(1, size(u_konc, 2)+1);
y_stat(1) = ypp;
y_z_stat = zeros(1, size(u_konc, 2)+1);
y_z_stat(1) = ypp;

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
ylim([-0.1 1.6]);
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
    z(1:k_konc) = 0;
    u(1:14) = upp;
    u(15:k_konc) = u_konc(i);
    y(1:k_konc) = ypp;

    for k=12:k_konc
        y(k) = symulacja_obiektu8y_p2(u(k-6), u(k-7), z(k-1), z(k-2), ...
        y(k-1), y(k-2));
    end

    plot(1:k_konc, y);
    hold on;
  
    y_stat(i+1) = y(k_konc); % zapis danych dla wektora statycznego y
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
export_fig('./pliki_wynikowe/zad2_y(k)_od_u.pdf');

figure;
for i=1:size(z_konc, 2)
    z(1:k_konc, 1) = 0; % U_pp=0
    z(15:k_konc, 1) = z_konc(i);
    plot(1:k_konc, z);
    ylim([-0.1 2.5]);
    hold on;
end

% Dodanie informacji do wygenerowanego wykresu
legend_list = strings([1 size(z_konc, 2)]);
xlabel('k');
ylabel('z(k)');
ylim([-0.1 1.6]);
title('Wykres z(k)');
for i=1:size(z_konc, 2)
    legend_text = 'z(k)=' + string(z_konc(i));
    legend_list(i) = legend_text;
end
legend(legend_list, 'Location', 'best');
hold off;
export_fig('./pliki_wynikowe/zad2_z(k).pdf');

figure;
for i=1:size(z_konc, 2)
    u(1:k_konc) = 0;
    z(1:14) = 0;
    z(15:k_konc) = z_konc(i);
    y_z(1:k_konc) = ypp;
    
    for k=16:k_konc
        y_z(k) = symulacja_obiektu8y_p2(u(k-6), u(k-7), z(k-1), z(k-2), ...
        y_z(k-1), y_z(k-2));
    end

    plot(1:k_konc, y_z);
    hold on;
  
    y_z_stat(i+1) = y_z(k_konc); % zapis danych dla wektora statycznego y
end

% Dodanie informacji do wygenerowanego wykresu
legend_list = strings([1 size(u_konc, 2)]);
xlabel('k');
ylabel('y^z(k)');
title('Wykres y^z(k)');
for i=1:size(u_konc, 2)
    legend_text = 'y^z(k)=' + string(y_z_stat(i+1));
    legend_list(i) = legend_text;
end
legend(legend_list, 'Location', 'southeast');
hold off;
export_fig('./pliki_wynikowe/zad2_y(k)_od_z.pdf');

%% Realizacja zadania 3
% Inicjalizacja danych dla skoku sterowania
z(1:k_konc) = 0;
u(1:14) = 0;
u(15:k_konc) = 1;
y(1:k_konc) = ypp;

for k=16:k_konc
    y(k) = symulacja_obiektu8y_p2(u(k-6), u(k-7), z(k-1), z(k-2), ...
        y(k-1), y(k-2));
end

for k=1:D
    s(k) = y(k+15); % ypp=0 i upp=0 - pomijamy
end

figure;
stairs(1:D, s);
title("Odp. skokowa sygnału wyjściowego s(k)")
xlabel('k');
ylabel('s(k)');
ylim([0 2.6])
export_fig('./pliki_wynikowe/zad3_s(k)_od_u.pdf');

% Inicjacja danych dla skoku zakłócenia
u(1:k_konc) = 0;
z(1:14) = 0;
z(15:k_konc) = 1;
y_z(1:k_konc) = ypp;

for k=16:k_konc
    y_z(k) = symulacja_obiektu8y_p2(u(k-6), u(k-7), z(k-1), z(k-2), ...
        y_z(k-1), y_z(k-2));
end

for k=1:D_z
    s_z(k) = y_z(k+15); % ypp=0 i upp=0 - pomijamy
end

figure;
stairs(1:D_z, s_z);
title("Odp. skokowa sygnału wyjściowego s^z(k)")
xlabel('k');
ylabel('s^z(k)');
export_fig('./pliki_wynikowe/zad3_s_z(k)_od_z.pdf');
