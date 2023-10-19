clear all;

n=400; % Ilość sprawdzeń
U(1:n)=1; %Punkty pracy wejściowe
Y(1:11)=1.7; %Punkty pracy wyjściowe
%Y(12:n)=0;
PP=0;
% Punkty pracy
Upp = 1.000000e+00;
Ypp = 1.700000e+00;

% Wywołanie funkcji dla punktów pracy
for k=12:n
    Y(k)=symulacja_obiektu8y_p1(U(k-10),U(k-11),Y(k-1),Y(k-2));
end

for i=1:n
    if U(i) == 1 && Y(i) == 1.7
        PP = PP+1;
    end
end
    if PP == 400
        disp('Punkty pracy są poprawne')
    else
        disp('Punkty pracy nie są poprawne')
    end

    subplot(2,1,1);
    stairs(U);
    xlabel('k');
    ylabel('U');
    hold on;
    subplot(2,1,2);
    stairs(Y);
    xlabel('k');
    ylabel('Y');
    hold on;