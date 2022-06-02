% A =[1, 0, 2, 0;
%     0, 4, 0, 1;
%     2, 0, 7, 0;
%     0, 1, 0, 10];
% A = sparse(A);
% L = ichol(A);
% A = L * L';
% disp(A);
% disp("=========================")
% L = chol(A);
% A = L' * L;
% disp(A);
% A = [1, 1;
%      1, 3];
% B = [1, 2;
%      1, 6];
% X = A \ B;
size = 10000;
A = rand(size, size);
A = sparse(A);
timect = cputime;
[L,U] = ilu(A);
fprintf("time cost is %3.5f", cputime - timect)
% disp(L);
% disp(U);
% C = L * U;
% disp(C);

       

