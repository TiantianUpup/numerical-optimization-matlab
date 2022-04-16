% 基于不完全LU分解的预处理的共轭梯度法求解方程AX=B
% 
% INPUT
% A ...... 方程系数矩阵
% B
% X ...... 初始迭代点
% eta .... 收敛条件
% i_max .. 最大迭代次数
%
% OUTPUT
% X ...... 方程的解
% iten ... 算法迭代次数
% cputime . 算法耗时
function runhist = icmcg(A, B, X, eta, i_max)
    timect = cputime;
    n = size(A, 1); %获取方程维度
    i = 0;  
    R0 = A*X-B;
    %求解MY=R0 其中M=L*U [L, U]=ilu(A)
    A = sparse(A); %IC只能应用在稀疏矩阵中
    [L, U] = ilu(A); %A = L * U; 
    temp = L \ R0;
    Y0  = U \ temp;
    P = -Y0;
    R_norm = norm(reshape(R0, [1, n^2]));
 
    while R_norm > eta && i < i_max
        alpha = (reshape(R0, [1, n^2]) * reshape(Y0, [n^2, 1])) /  (reshape(P, [1, n^2]) * reshape(A*P, [n^2, 1]));
        X = X + alpha * P;
        R1 = R0 + alpha * A * P;
        R_norm = norm(reshape(R1, [1, n^2]));
        temp = L \ R1;
        Y1  = U \ temp;

        beta = (reshape(R1, [1, n^2]) * reshape(Y1, [n^2, 1])) / (reshape(R0, [1, n^2]) * reshape(Y0, [n^2, 1]));
        P = -Y1 + beta * P;
        R0 = R1;
        Y0 = Y1;
       
        i = i + 1;
 
    end
    
    runhist.X = X;
    runhist.iten = i;
    runhist.cputime = cputime - timect;