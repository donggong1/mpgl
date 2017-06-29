% mex c files
mex('utils/mex_graph_diff.cpp',  '-outdir', 'utils')
mex('utils/mex_graph_diffT.cpp', '-outdir', 'utils')

% your/path/to/openblas/lib/
openblas_path = 'OpenBLAS-0.2.18/lib/';
% mex('-O', 'mexsolver/mex_solve_x_cg.cpp', [openblas_path, 'libopenblas.a'], '-outdir', 'mexsolver');
% mex('-O', 'mexsolver/mex_solve_x_cg_proj.cpp', [openblas_path, 'libopenblas.a'], '-outdir', 'mexsolver');

mex('-O', 'mexsolver/mex_subadmm.cpp', [openblas_path, 'libopenblas.a'], '-outdir', 'mexsolver');
mex('-O', 'mexsolver/mex_subadmm_proj.cpp', [openblas_path, 'libopenblas.a'], '-outdir', 'mexsolver');