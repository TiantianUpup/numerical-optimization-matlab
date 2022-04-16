% A =[1, 0, 2, 0;
%     0, 4, 0, 1;
%     2, 0, 7, 0;
%     0, 1, 0, 10];

% B = [1, 3, 1, 1;
%      1, 9, 7, 0;
%      0, 2, 31, 10;
%      8, 3, 1, 1];

size = 1000;
A = rand(size, size); 
A = A*A'; %产生一个半正定随机整数矩阵
B = rand(size, size);
AB = A * B;
X = zeros(size);
eta = 1e-8;
i_max = 5000;
%runcghist = mcg(A, AB, X, eta, i_max);
runiccghist = icmcg(A, AB, X, eta, i_max);
runilucghist = ilumcg(A, AB, X, eta, i_max);
%fprintf("cg iten is %d, cputime is %5.2f, rss is %3.4e\n", runcghist.iten, runcghist.cputime, norm(reshape(runcghist.X-B, [1, size^2])));
fprintf("iccg iten is %d, cputime is %5.2f, rss is %3.4e\n", runiccghist.iten, runiccghist.cputime, norm(reshape(runiccghist.X-B, [1, size^2])));
fprintf("ilucg iten is %d, cputime is %5.2f, rss is %3.4e\n", runilucghist.iten, runilucghist.cputime, norm(reshape(runilucghist.X-B, [1, size^2])));
