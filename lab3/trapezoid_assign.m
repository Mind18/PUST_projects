function [w] = trapezoid_assign(x, trap_params)
%trapezoid_assign Function assignes trapezoid fuzzy logic assignment
%   Function reprezents trapezoid type of fuuzy logic assingment function -
%   if value is in range of <high_begin, high_end> function returns 1
%   If value is in either range <rise_begin, high_begin> or 
%   <high_end, fall_end> it returns a value of a linear function connecting
%   points (rise_begin, 0), (high_begin, 1) or (high_end, 1), (fall_end, 0)
%   In other cases function returns 0.
    if (x >= trap_params(2)) && (x <= trap_params(3))
        w = 1;
    elseif (x >= trap_params(1)) && (x <= trap_params(2))
        a = 1/(trap_params(2)-trap_params(1));
        b = -(trap_params(1)/(trap_params(2)-trap_params(1)));
        w = a*x+b;
    elseif (x >= trap_params(3)) && (x <= trap_params(4))
        a = 1/(trap_params(3)-trap_params(4));
        b = -(trap_params(4)/(trap_params(3)-trap_params(4)));
        w = a*x+b;
    else
        w = 0;
    end
end

