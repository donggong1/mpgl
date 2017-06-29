#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	// get data from input
	int num_edge;
	double* edge_head_list;
	double* edge_tail_list;
	double* x_list;
	double* diff_list;
	int idx;

	x_list = (double*)mxGetPr(prhs[0]);
	edge_head_list = (double*)mxGetPr(prhs[1]);
	edge_tail_list = (double*)mxGetPr(prhs[2]);
	num_edge = (int)(*(double*)mxGetPr(prhs[3]));
//	num_edge = (int)mxGetM(prhs[1]);
	// output
	plhs[0] = mxCreateDoubleMatrix(num_edge,1,mxREAL); 
	diff_list = (double *)mxGetPr(plhs[0]);
	for(idx=0; idx<num_edge; idx++)
	{
		diff_list[idx] = x_list[(int)edge_head_list[idx]-1] - x_list[(int)edge_tail_list[idx]-1]; 
	}	
}



