

figure;
plot(1:size(y20, 2), y20);
hold on;
plot(1:size(y30, 2), y30);
plot(1:size(y40, 2), y40);
plot(1:size(y50, 2), y50);
plot(1:size(y60, 2), y60);
plot(1:size(y70, 2), y70);
plot(1:size(y80, 2), y80);
xlabel('G1');
ylabel('T1(G1)');
legend('y20', 'y30', 'y40', 'y50', 'y60', 'y70', ...
    'y80', 'Loction', 'best');
title('Odp. skokowe T1(G1) dla różych skoków G1');
print('./ss3/odp_skokowe_razem_lab.png', '-dpng');
