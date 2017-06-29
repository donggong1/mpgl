% generate synthetic data with model y = Ax+n
%%
outpath = './data/syn/';
if(~exist(outpath, 'dir'))
    mkdir(outpath);
end

%%
% ls: least square
% flsa: fused lasso signal approximation (FLSA), A is identity matrix.
type_A = 'ls'; % 'flsa'; 
% fl_1d: 1d fused lasso
% gfl_2dtv: generalized fused lasso with 2d tv guided model
type_x = 'fl_1d'; % 'gfl_2dtv';

fprintf('type_A=%s, type_x=%s\n', type_A, type_x);
n = 1000;    
nsig = 0.05; % noise level
%%
randNum = 'default';
if(strcmp(type_x, 'fl_1d'))
    fprintf('Dim=%d\n', n);
    x = zeros(n,1);
    x(round(n/2)-round(n/20):round(n/2)+round(n/20)) = 0.5;
    %
    D.edge_head = [2:n, 1]';
    D.edge_tail = [1:n]';
    d.D = D;
    d.graph_type = {'1d fused lasso', n};
elseif(strcmp(type_x, 'gfl_2dtv'))
    t = floor(sqrt(n));
    n1 = t;
    n2 = t;
    n = n1*n2;
    fprintf('dim=%d\n', n);
    x = zeros(n1,n2);
    x(round(n1/2)-round(n1/20):round(n1/2)+round(n1/20),...
        round(n2/2)-round(n2/20):round(n2/2)+round(n2/20)) = 0.5;
    x = x(:);
    %
    % rand_edge_num = round(n*0.2); % number of the randomly generated edges
    rand_edge_num = 0; % only 2d tv based graph
    [edge_head, edge_tail] = gen_graph_2dtv(t,t);
    rhead = max(floor(rand(rand_edge_num, 1)*n),0);
    rtail = max(floor(rand(rand_edge_num, 1)*n),0);
    D.edge_head = [edge_head; rhead];
    D.edge_tail = [edge_tail; rhead]; % modified at 1016
    fprintf('graph edge number=%d\n', length(D.edge_head));
    d.D = D;
    d.graph_type = {'gen fused lasso, 2dtv', n};
    d.coor2d = [n1, n2];
else
    fprintf('Wrong type_x\n');
end


if(strcmp(type_A, 'ls'))
    m = n;
    rng(randNum);
    A=randn(m,n);
    y = A*x;
elseif(strcmp(type_A, 'flsa'))
    y = x; 
    A = 1;
else
    fprintf('Wrong type_x\n');
end

noise=randn(size(y))*nsig;
y = y+noise;

%%
d.x_gt = x;
d.y = y;
d.A = A;
d.noise = noise;
d.nsig = nsig;

save([outpath, 'syn_1dfl_dim_', num2str(n), '.mat'], 'd');



