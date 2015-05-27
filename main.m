clear all;
global n n1 n2 h DK DRho Jacobian Q_I x k JJ;

% Store start time for running time calculation
time_start = fix(clock);
time_start = time_start(4)*3600 + time_start(5)*60 + time_start(6);


n = 3;
n1 = n*(n-1)/2;
n2 = n*(n+3)/2;

T_end = 100;
h = 0.001;
y0 = [0 1 0];


% calculate J
x = sym('x', [n, 1]);
Jacobian = sym('Jacobian', [n, n]);
for i = 1:n
    Jacobian(:, i) = diff(v(x, 1), x(i));
end

% set lower triangle of matrix K
k = sym('k', [n1, 1]);
% define K as a skew-symmetric matrix
K = sym('K', [n, n]);
counter = 1;
for i = 1:n
    for j = 1:n
        if (i < j)
            K(i, j) = -k(counter);
            counter = counter + 1;
        end
        if (i == j)
            K(i, j) = 0;
        end
        if (i > j)
            K(i, j) = -K(j, i);
        end
    end
end

% intermediary matrices
I = eye(n);
G = I+K;
H = inv(G);
JJ = sym('JJ', [n, n]);
J2 = G * JJ * (G.');

 
S = sym('S', [n, n]);
TMP = (H.') * J2 * H;
for i = 1:n
    for j = 1:n
        if (i < j)
            S(i, j) = 1/2*TMP(j, i);
        end
        if (i == j)
            S(i, j) = 0;
        end
        if (i > j)
            S(i, j) = -1/2*TMP(i, j);
        end
    end
end

% matrix of derivative of K
TMP = (G.') * S * G;


% Lower triangle of derivative of matrix K. It will be another part
% of differential equations, which have to be solved, and will be 
% an argument of WHOLE_SYSTEM function
DK = sym('DK', [n1, 1]);
counter = 1;
for i = 1:n
    for j = 1:n
        if (i < j)
            DK(counter) = -TMP(i, j);
            counter = counter + 1;
        end
    end
end



DRho = sym('DRho', [n, 1]);
for i = 1:n
    DRho(i) = (H(:, i).') * J2 * H(:, i);
end


y_start = zeros(n2, 1);
y_start(1:n) = y0;



% Calculating with Runge-Kutta manually
num_of_iter = fix(T_end/h);
t0 = 0;
y0 = y_start;
KK = zeros(n);
T = zeros(num_of_iter, 1);
Y = zeros(num_of_iter, n);
T(1) = 0;
Y(1, :) = zeros(n, 1);
rho_I = zeros(n, 1);
Q_I = I;
for iter = 2:num_of_iter
    y0 = RK(y0);
    t0 = t0 + h;
    
	% change matrix K for current moment of time
    counter = 1;
    for i = 1:n
        for j = 1:n
            if (i < j)
                KK(i, j) = -y0(n + counter);
                KK(j, i) = y0(n + counter);
                counter = counter + 1;
            end
            if (i == j)
                KK(i, j) = 0;
            end
        end
    end
    
    Q = Q_I*(2*I/(I + KK) - I);
    
    if (norm(KK, 1) >= 0.2)
        rho_I = rho_I + y0(n + n1 + 1 : n2);
        Q_I = Q;
        y0(n+1 : n2) = 0;
    end

    T(iter) = t0;
    Y(iter, :) = (rho_I + y0(n + n1 + 1 : n2))/t0;
end
plot(T, Y)
xlabel('Time')
ylabel('LCEs')


time_end = fix(clock);
time_end = time_end(4)*3600 + time_end(5)*60 + time_end(6);
time_total = time_end - time_start;
time_total = [fix(time_total/3600) fix(mod(time_total, 3600)/60) mod(time_total, 60)];
array2table(time_total, 'VariableNames',{'Hours' 'Minutes' 'Seconds'})


