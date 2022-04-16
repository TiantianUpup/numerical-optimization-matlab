% ���ڲ���ȫLU�ֽ��Ԥ����Ĺ����ݶȷ���ⷽ��AX=B
% 
% INPUT
% A ...... ����ϵ������
% B
% X ...... ��ʼ������
% eta .... ��������
% i_max .. ����������
%
% OUTPUT
% X ...... ���̵Ľ�
% iten ... �㷨��������
% cputime . �㷨��ʱ
function runhist = icmcg(A, B, X, eta, i_max)
    timect = cputime;
    n = size(A, 1); %��ȡ����ά��
    i = 0;  
    R0 = A*X-B;
    %���MY=R0 ����M=L*U [L, U]=ilu(A)
    A = sparse(A); %ICֻ��Ӧ����ϡ�������
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