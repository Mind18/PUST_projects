function [E] = Approximation_lab(params, y50, Td)

E = 0;

alfa_1 = exp(-1/params(1));
alfa_2 = exp(-1/params(2));
a1 = -alfa_1 - alfa_2;
a2 = alfa_1 * alfa_2;
b1 = (params(3) / (params(1) - params(2)))*(params(1)*(1 - alfa_1)- ...
    params(2)*(1-alfa_2));
b2 = (params(3) / (params(1) - params(2)))*(alfa_1*params(2)*(1 - alfa_2)- ...
    alfa_2*params(1)*(1-alfa_1));

dim = size(y50, 2);
u_param = zeros(1, dim); y_param = zeros(1, dim);
u_param(1, 1:Td+2) = 28; y_param(1, 1:Td+2) = 31.75;

for i=Td+3:dim
    y_param(i) = b1*u_param(i-1-Td) + b2*u_param(i-2-Td) - ...
        a1*y_param(i-1) - a2*y_param(i-2);
    E = E + (y50(i) - y_param(i))^2;
end
end