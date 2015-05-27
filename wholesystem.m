function out = wholesystem(y)
global n n1 n2 DK DRho Q_I k JJ;

out = zeros(n2, 1);
out(1:n) = v(y(1:n), 0);

% shows the value of Jacobi matrix for current moment of time
curJ = (Q_I.')*J(y(1:n))*Q_I;

% will change value of Jacobi matrix in derivative of matrix K notation to 
% a current one
TMP = subs(DK, JJ(:, :), curJ(:, :));

% will substitute initial conditions and will get a second system of DEs
out(n + 1 : n + n1) = double(subs(TMP, k(:), y(n + 1 : n + n1)));


TMP2 = subs(DRho, JJ(:, :), curJ(:, :));
% will substitute initial conditions and will get a third system of DEs
out(n + n1 + 1 : n2) = double(subs(TMP2, k(:), y(n + 1 : n + n1)));

end	