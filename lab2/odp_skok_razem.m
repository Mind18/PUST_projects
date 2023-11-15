figure;
subplot(2,1,1);
plot(x, z10);
hold on;
plot(x, z30);
plot(x, z50);
xlabel('k');
ylabel('z(k)');
ylim([0 55]);
legend('z(k)=10', 'z(k)=30', 'z(k)=50', 'Location', 'best');

subplot(2,1,2)
plot(x, y10);
hold on;
plot(x, y30);
plot(x, y50);
xlabel('k');
ylabel('y(k)');
ylim([30 45]);
legend('y(k)=10', 'y(k)=30', 'y(k)=50', 'Location', 'best');
