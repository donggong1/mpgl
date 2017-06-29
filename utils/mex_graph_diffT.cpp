#include "mex.h"
#include <cstring>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	// get data from input
	int num_edge, num_elem;
	//int ; // dim of output
	double* edge_head_list, *edge_tail_list , *z_list, * diffT_list;
	double tmpz;
	int idx_edge;

	z_list = (double*)mxGetPr(prhs[0]);
	num_edge = (int)mxGetM(prhs[0]);

	edge_head_list = (double*)mxGetPr(prhs[1]);
	edge_tail_list = (double*)mxGetPr(prhs[2]);
	num_elem = (int)(*(double*)mxGetPr(prhs[3]));

	// output
	plhs[0] = mxCreateDoubleMatrix(num_elem,1,mxREAL); 
	diffT_list = (double *)mxGetPr(plhs[0]);

	memset(diffT_list, 0, sizeof(diffT_list));
	for(idx_edge=0; idx_edge<num_edge; idx_edge++)
	{
		tmpz = z_list[idx_edge];
		diffT_list[(int)edge_head_list[idx_edge]-1] += tmpz;
		diffT_list[(int)edge_tail_list[idx_edge]-1] -= tmpz;
	}	
}



