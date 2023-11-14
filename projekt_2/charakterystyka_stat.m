figure;
[u_stat, sortIndex] = sort(u_stat);
y_stat = y_stat(sortIndex);
y_z_stat = y_z_stat(sortIndex);

subplot(2, 1, 1)
hold on;
xlabel('u');
ylabel('y(u)');
plot(u_stat, y_stat, '.', 'MarkerSize',12);
plot(u_stat, y_stat);
legend("Interpolacja danych statycznych", "Dane statyczne y(u)", 'Location', 'southeast');
title('Dane statyczne y(u)');
hold off

subplot(2, 1, 2)
hold on
plot(u_stat, y_z_stat, '.', 'MarkerSize',12);
plot(u_stat, y_z_stat);
xlabel('u');
ylabel('y(z)');
legend("Interpolacja danych statycznych", "Dane statyczne y(z)", 'Location', 'southeast');
k_stat = rdivide(y_stat, u_stat);
xlabel('u');
ylabel('y(z)');
title('Dane statyczne y(z)');
hold off;
export_fig('./pliki_wynikowe/zad2_y(u,z)_char_stat.pdf');

