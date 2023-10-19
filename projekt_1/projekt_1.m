import symulacja_obiektu8y_p1.*

% Inicjalizacja

clear;
zad_2_target = 'u'; % 'u' lub 'y'
u_konc = [0.7 1.3 1.1 0.5 1.5];
k_konc = 400;
u(1, 1:11) = 0.5;
u(1, 12:k_konc) = 1;
y = zeros(1, k_konc);

% Realizacja zadania 1

for k=12:k_konc
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
end

figure;
plot(1:k_konc, u);
hold on;
xlabel('k');
ylabel('u(k)');
ylim([0 1.2]);
title('Wykres u(k)');
hold off;
export_fig('./pliki wynikowe/test_punktu_pracy_u(k).pdf');

figure;
plot(1:k_konc, y);
hold on;
xlabel('k');
ylabel('y(k)');
title('Wykres y(k)');
hold off;
export_fig('./pliki wynikowe/test_punktu_pracy_y(k).pdf');

% Realizacja zadania 2

figure;
for i=1:size(u_konc, 2)
    u(1, 1:200) = 1; % U_pp=1
    u(1, 201:k_konc) = u_konc(i);
    y = zeros(1, k_konc);
    
    for k=12:k_konc
        y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
    end

    if zad_2_target == 'u'
        plot(1:k_konc, u);
    elseif zad_2_target == 'y'
        plot(1:k_konc, y);
    end
    hold on;
end

if zad_2_target == 'u'
    legend_list = strings([1 size(u_konc, 2)]);
    xlabel('k');
    ylabel('u(k)');
    ylim([0.4 1.6]);
    title('Wykres u(k)');
    for i=1:size(u_konc, 2)
        legend_text = 'u(k)=' + string(u_konc(i));
        legend_list(i) = legend_text;
    end
    legend(legend_list, 'Location', 'best');
    hold off;
    export_fig('./pliki wynikowe/zad2_u(k).pdf');
elseif zad_2_target == 'y'
    legend_list = strings([1 size(u_konc, 2)]);
    xlabel('k');
    ylabel('y(k)');
    title('Wykres y(k)');
    for i=1:size(u_konc, 2)
        legend_text = 'y(k)=' + string(u_konc(i));
        legend_list(i) = legend_text;
    end
    legend(legend_list, 'Location', 'southeast');
    hold off;
    export_fig('./pliki wynikowe/zad2_y(k).pdf');
end
