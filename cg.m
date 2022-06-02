%���Ժ���
function main
    % ��mat�ļ��ж�ȡ����
    % bundle1.mat gr_30_30.mat
    S1 = load('bundle1.mat');
    BP1 = struct2cell(S1);
    Mymat1 = cell2mat(BP1);
    A = Mymat1.A;
    n = size(A, 1); %��ȡ�����ά��
    x_opt = randn(n, 1);
    b = A * x_opt;
    tol = 1e-6;
    maxit = 1000;
    %Ԥ����ķ���
    %0 ������Ԥ���� 
    %1 �Խ�Ԥ���� 
    %3 incomplete Cholesky����Ԥ����
    pcond = 0; 
    x0 = zeros(n, 1);
    start_time = cputime;
    [x, flag, relres, iter, resvec] = mycg(A, b, tol, maxit, pcond, x0);
    
    fprintf('the dimension of A is %2d, the pcond is %2d \n', n, pcond);
    fprintf('flag    iter    relres       cost \n');
    fprintf('%2d     %2d     %1.6e   %1.6f \n', flag, iter, relres, cputime-start_time);
end

function [x, flag, relres, iter, resvec] = mycg(A, b, tol, maxit, pcond, x0)
    iter = 0; %��������
    r0 = A * x0 - b;
    x = x0;

    r_norm = norm(r0);
    flag = 0;

    % ʹ�öԽǽ���Ԥ����
    if (pcond == 1)
        tempA = diag(A); % ��ȡ���Խ����ϵ�Ԫ��
        dA = diag(tempA); % ��ɶԽǾ���
        y0 = dA \ r0;
        p = -y0;
        % ʹ�� incomplete Cholesky����Ԥ����
    elseif (pcond == 3)
        A = sparse(A); %ICֻ��Ӧ����ϡ�������
        L = ichol(A); %A = L * L';
        temp = L \ r0;
        y0 = L' \ temp;
        p = -y0;
    else
        p = -r0;
    end

    % �Ƿ�������������
    while (r_norm > tol && iter <= maxit)

        if (pcond == 0)
            alpha = (r0' * r0) / (p' * A * p);
        else
            alpha = (r0' * y0) / (p' * A * p);
        end

        x = x + alpha * p;
        r1 = r0 + alpha * A * p;

        if (pcond == 1)
            y1 = A \ r1;
            beta = (r1' * y1) / (r0' * y0);
        elseif (pcond == 3)
            temp = L \ r1;
            y1 = L' \ temp;
            beta = (r1' * y1) / (r0' * y0);
        else
            beta = (r1' * r1) / (r0' * r0);
        end

        if (pcond == 0)
            p = -r1 + beta * p;
        else
            p = -y1 + beta * p;
        end
       
        r0 = r1;
        if (pcond ~=0)
            y0 = y1;
        end    
        r_norm = norm(r0);
        iter = iter + 1;
        relres = norm(b - A * x) / norm(b);
        resvec = b - A * x;

    end

    if (iter > maxit)
        flag = 1;
    end

end
