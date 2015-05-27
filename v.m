function out = v(x, flag)

global n;
a = 16;
b = 45.92;
c = 4;

% In fact, most interesting thing there - argument FLAG.
% Depending on its value function change its specialisation,
% so anywhere in program we can use it now in such terms, 
% in which it will be more usefull for us.
if (flag == 1)
    out = sym('out', [n, 1]);
end
if (flag == 0)
    out = zeros(n, 1);
end

% Now we can write right part of original DS there, 
% and we have not to write it elsewhere in program
 out(1) = a*(x(2) - x(1));
 out(2) = x(1)*(b - x(3)) - x(2);
 out(3) = x(1)*x(2) - c*x(3);

end