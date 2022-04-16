A =[1, 0, 2, 0;
    0, 4, 0, 1;
    2, 0, 7, 0;
    0, 1, 0, 10];
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
A = sparse(A);
[L,U] = ilu(A);
% disp(L);
% disp(U);
C = L * U;
disp(C);

       

