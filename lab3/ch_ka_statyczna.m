% Badane sterowania
u_stat = [20; 30; 40; 50; 60; 70; 80];
% Punkty stabilizacji dla konkretnych sterowań
y_stat = [31; 35.06; 39.87; 43.81; 46.56; 48.62; 51.06];

figure;
plot(u_stat, y_stat);
xlabel('G1');
ylabel('T1(G1)');
title('Charakterystyka statyczna obiektu T1(G1)');
print('./ss3/ch_ka_statyczna.png', '-dpng');
