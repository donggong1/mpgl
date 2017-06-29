function [edge_head, edge_tail] = gen_graph_2dtv(m,n)
% 2dtv graph
coorh = [1:m]';
coorh = repmat(coorh, 1, n);
coorw = [1:n];
coorw = repmat(coorw, m, 1);

dirv_coorh = [coorh(2:end, :); coorh(1,:)];
dirv_coorw = coorw;


dirh_coorh = coorh;
dirh_coorw = [coorw(:, 2:end), coorw(:,1)];

coor1d = (coorw-1)*m + coorh;
coor1d_v = (dirv_coorw-1)*m + dirv_coorh;

coor1d_h = (dirh_coorw-1)*m + dirh_coorh;

edge_head = [coor1d_v; coor1d_h];
edge_head = edge_head(:);
edge_tail = [coor1d; coor1d];
edge_tail = edge_tail(:);
return