% entry_mpgl
addpath('utils');
addpath('solver_ls')
addpath('solver_proj')
addpath('mexsolver')

inpath = './data/syn/';
outpath = './genlasso_result/syn/';

if(~exist(outpath, 'dir'))
    mkdir(outpath);
end
%%
n = 1000;
name = ['syn_1dfl_dim_', num2str(n), '.mat'];
% n = 10000;
% name = ['syn_1dgfl_proj_dim_', num2str(n), '.mat'];
load([inpath, name]);
%%
mpgl_lambda = 0.005;
opts_mpgl.ite_max_mp = 3;
opts_mpgl.ite_max_inner = 3;5;
opts_mpgl.ite_max_inner_final = 5;5;
opts_mpgl.x_diff_tol_inner = 1e-4;
opts_mpgl.x_diff_tol_mp = 1e-3; % rough setting in this implementation, cannot rely on this
opts_mpgl.wa_lambda = 100;
opts_mpgl.gt_x = d.x_gt;
opts_mpgl.rho = 1;
opts_mpgl.remove_rate = 0;0.35;
opts_mpgl.cg_max_ite_inner = 3;2;
opts_mpgl.cg_max_ite_inner_final = 3;3;
opts_mpgl.rho_rate = 1;
opts_mpgl.kappa_rate = 0.5;
opts_mpgl.type_D = 'graphfuse'; D = d.D; D.sI = 1; D.edge_num = length(D.edge_head);  0.6;
%%
edge_num = length(D.edge_head);
Dx = [d.x_gt; mex_graph_diff(d.x_gt, D.edge_head, D.edge_tail, edge_num)];
xx = abs(Dx)>0; sd = sum(xx);
opts_mpgl.kappa = max(floor(sd/(opts_mpgl.ite_max_mp+1))+15*(2*max(1, n/1000)), 1);
D.edge_head = max(D.edge_head, 1);
D.edge_tail = max(D.edge_tail, 1);
opts_mpgl.kappa = -1;
if(d.A==1)
    [x, supp_map, time] = mpgenlasso_proj_mex(d.y, D, D.edge_num, ...
        int32(D.edge_head-1), int32(D.edge_tail-1), D.sI, mpgl_lambda, opts_mpgl);
else
    [x, supp_map, time] = mpgenlasso_mex(d.A, d.y, D, D.edge_num,...
        int32(D.edge_head-1), int32(D.edge_tail-1), D.sI, mpgl_lambda, opts_mpgl);
end
%%
t = time;
r = x-d.x_gt;
nr = norm(r);
mse = mean(abs(r).^2);
re = d.y-d.noise-d.A*x;
nre = norm(re);
fprintf('time:%f, norm error: %f, MSE:%e, recon_error:%f\n',t, nr, mse,nre);
%%
figure;
plot(x, '-.');
hold on
plot(r, 'g');
legend('estimated x', 'error')
title('MPGL');