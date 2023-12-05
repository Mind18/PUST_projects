import trapezoid_assign.*

figure;
% Sprawdzenie poprawności danych
if n_regulatorow == size(reg_part, 2)
    lin = -1:0.01:1;
    y_trap = zeros(1, size(lin, 2));
    for i=1:n_regulatorow
        for k=1:size(lin, 2)
            y_trap(k) = trapezoid_assign(lin(k), reg_part{i});
        end
        % y = trapezoid_assign(lin, reg_part{i});
        plot(lin, y_trap);
        hold on;
    end
else
    disp('Niepoprawna liczba regulatorów')
end