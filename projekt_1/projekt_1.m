import symulacja_obiektu8y_p1.*

clear;
k_konc = 200;
u(1, 1:11) = 0;
u(1, 12:k_konc) = 1;
y = zeros(1, k_konc);

for k=12:k_konc
    y(k) = symulacja_obiektu8y_p1(u(k-10), u(k-11), y(k-1), y(k-2));
end

figure;
plot(1:k_konc, u);
hold on;
xlabel('k');
ylabel('u(k)');
title('Wykres u(k)');
hold off;

figure;
plot(1:k_konc, y);
hold on;
xlabel('k');
ylabel('y(k)');
title('Wykres y(k)');
hold off;
