#include <stdio.h>
#include <stdlib.h>
#include <memory>
#include <cstring>
#include "math.h"
#include "mex.h"
#include "cblas.h"
// #include "mex_solve_x_cg.h"
/*
This implementation is mainly for the graph guided flsa (projection) problem.
*/
/*
struct graph
{
	int *edge_head;
	int *edge_tail;
	int edge_num;
	int node_num;
};*/


void graph_cal_Dx(double *z, double *x, int *edge_head, int *edge_tail, int num_edge)
{
	// z = Dx;
	for (int idx_edge = 0; idx_edge < num_edge; ++idx_edge)
	{
		z[idx_edge] = x[edge_head[idx_edge]] - x[edge_tail[idx_edge]];
	}
}

void graph_cal_DTz(double *dtz, const double *z, int *edge_head, int *edge_tail, int num_edge)
{
	// dtz = dtz + DTz;
	double tmpz;
	for (int idx_edge = 0; idx_edge < num_edge; ++idx_edge)
	{
		tmpz = z[idx_edge];
		dtz[edge_head[idx_edge]] += tmpz;
		dtz[edge_tail[idx_edge]] -= tmpz;
	}
}

void graph_cal_DDTx(double *ddtx, const double *x, int *edge_head, int *edge_tail, int num_edge)
{
	// ddtx = ddtx + DDTx;
/*	int idx_edge, tmp_idx_head, tmp_idx_tail;
	double tmpz;
	for (idx_edge = 0; idx_edge < num_edge; ++idx_edge)
	{
		tmp_idx_head = edge_head[idx_edge];
		tmp_idx_tail = edge_tail[idx_edge];
		tmpz = x[tmp_idx_head] - x[tmp_idx_tail];
		ddtx[tmp_idx_head] += tmpz;
		ddtx[tmp_idx_tail] -= tmpz;
	}
*/
	double tmpz;
	for (int idx_edge = 0; idx_edge < num_edge; ++idx_edge)
	{
		tmpz = x[edge_head[idx_edge]] - x[edge_tail[idx_edge]];
		ddtx[edge_head[idx_edge]] += tmpz;
		ddtx[edge_tail[idx_edge]] -= tmpz;
	}
}

int solve_x_cg(double *x, const double *ATy, const double *zuI, const double *zuG, 
	int *edge_head, int *edge_tail, int edge_num, int node_num, double sI, double rho, int cg_max_ite, double cg_tol2)
{
	// setting for blas
	enum CBLAS_ORDER order;
    enum CBLAS_TRANSPOSE transa;
    int inc = 1;
    order = CblasColMajor;
	transa = CblasNoTrans;
	// graph information
//	int edge_num, node_num;
//	double *edge_head, *edge_tail;
	// temp variables
	double sI2;
	double *r, *tmpDTDx, *p, *Ap;
	double r_prev2, alpha, r2, beta;
	int cg_ite;
	sI2 = sI*sI;

	// r = -b = -(ATy + rho.*( mex_graph_diffT(zuG, edge_head, edge_tail, s)+sI.*zuI));
	r = new double[node_num];
	tmpDTDx = new double[node_num];
	for (int i = 0; i < node_num; ++i)
	{
		r[i] = sI*zuI[i];
		tmpDTDx[i] = sI2*x[i];
	}
/*	memset(r, 0, sizeof(double)*node_num);
	cblas_daxpy(node_num, sI, zuI, inc, r, inc);
	memset(tmpDTDx, 0, sizeof(double)*node_num);
	cblas_daxpy(node_num, sI2, x, inc, tmpDTDx, inc);*/

	graph_cal_DTz(r, zuG, edge_head, edge_tail, edge_num); // -r = DTzu;
	/*for (int i = 0; i < node_num; ++i)
	{
		r[i] = -ATy[i] - rho*r[i];
	}*/
	cblas_dscal(node_num, -rho, r, inc); // r = -b = -rho.*DTzu
	cblas_daxpy(node_num, -1, ATy, inc, r, inc); // ATy + rho.*DTzu,
	
	//=====
	// DTDx = mex_graph_diffT(Dx, edge_head, edge_tail, s) + sI2.*x
/*	tmpDTDx = new double[node_num];
	for (int i = 0; i < node_num; ++i)
	{
		tmpDTDx[i] = sI2*x[i];
	}*/
	graph_cal_DDTx(tmpDTDx, x, edge_head, edge_tail, edge_num);
	// rho.*( mex_graph_diffT(Dx, edge_head, edge_tail, s) + sI2.*x ) - b;
	cblas_daxpy(node_num, rho, tmpDTDx, inc, r, inc);  // note: r is -b: rho.*DTDx-b
	// r = ATA*x + rho.*( mex_graph_diffT(Dx, edge_head, edge_tail, s) + sI2.*x ) - b;

	cblas_daxpy(node_num, 1, x, inc, r, inc);
//	cblas_dgemv(order, transa, node_num, node_num, 1, ATA, node_num, x, inc, 1, r, inc); 
	
	// r_prev2 = r'*r;
	r_prev2 = cblas_ddot(node_num, r, inc, r, inc);
	
	// p = -r;
	p = new double[node_num];
	Ap = new double[node_num];
	/*for (int i = 0; i < node_num; ++i)
	{
		p[i] = -1*r[i];
		Ap[i] = sI2*p[i];
	}*/
	memset(p, 0, sizeof(double)*node_num);
	cblas_daxpy(node_num, -1, r, inc, p, inc);
	memset(Ap, 0, sizeof(double)*node_num);
	cblas_daxpy(node_num, sI2, p, inc, Ap, inc);
	
	// =====
	// Ap = ATA*p + rho.*( mex_graph_diffT(Dp, edge_head, edge_tail, s) + sI2.*p );
/*	Ap = new double[node_num];
	for (int i = 0; i < node_num; ++i)
	{
		Ap[i] = sI2*p[i];
	}*/
	
	graph_cal_DDTx(Ap, p, edge_head, edge_tail, edge_num);
	cblas_dscal(node_num, rho, Ap, inc);
	cblas_daxpy(node_num, 1, p, inc, Ap, inc);
//	cblas_dgemv(order, transa, node_num, node_num, 1, ATA, node_num, p, inc, rho, Ap, inc);	
	// cg 
	for (cg_ite = 0; cg_ite < cg_max_ite; ++cg_ite)
	{
		alpha = r_prev2 / cblas_ddot(node_num, p, inc, Ap, inc);
		cblas_daxpy(node_num, alpha, p, inc, x, inc); // x = alpha.*p + x;
		cblas_daxpy(node_num, alpha,Ap, inc, r, inc); // r = alpha.*Ap + r;
		r2 = cblas_ddot(node_num, r, inc, r, inc);
		beta = r2/r_prev2;
		cblas_dscal(node_num, beta, p, inc); // beta.*p
		cblas_daxpy(node_num, -1, r, inc, p, inc); // p = -r + beta.*p; 
		/*for (int i = 0; i < node_num; ++i)
		{
			Ap[i] = sI2*p[i];
		}*/
		memset(Ap, 0, sizeof(double)*node_num);
		cblas_daxpy(node_num, sI2, p, inc, Ap, inc);
		graph_cal_DDTx(Ap, p, edge_head, edge_tail, edge_num);
		cblas_dscal(node_num, rho, Ap, inc);
		cblas_daxpy(node_num, 1, p, inc, Ap, inc);
//		cblas_dgemv(order, transa, node_num, node_num, 1, ATA, node_num, p, inc, rho, Ap, inc);
		//
		r_prev2 = r2;
		if(r2<cg_tol2)
		{
			break;
		}
	}


	delete[] r;
	delete[] tmpDTDx;
	delete[] p;
	delete[] Ap;
	return 0;
//	return (double)cg_ite+1;
}


// only for generalized graph guided fused lasso
int solve_subadmm(double *x, double *z, double *u, double *dx, double *zu,
	const double *ATy, int *supp_idx_list, int non_zero_num,
	int *edge_head, int *edge_tail, int edge_num, int node_num, double sI, 
	double lambda, double rho, int admm_max_ite, double admm_tol2, 
	int cg_max_ite, double cg_tol2)
{//*dx, *u, *z, *zu,
	double  *z_prev, *x_prev, x_snorm, x_prev_snorm=0, x_diff_snorm, x_diff_rel;
	double lambda_rho, m_rate=1.8, m_rate2;
	int inc=1, s, admm_ite;
	lambda_rho = lambda/rho;
	s = node_num+edge_num;
//	dx = new double[s];
//	zu = new double[s];
//	z = new double[s];
//	u = new double[s];
	z_prev = new double[s];
	x_prev = new double[node_num];
	std::size_t doublesizeS = sizeof(double)*s;
	memset(z_prev, 0, doublesizeS);
	memset(z, 0, doublesizeS);
	memset(u, 0, doublesizeS);
	m_rate2 = m_rate-1;

//	cblas_dcopy(s, x_prev, inc, x, inc);
//	x_prev_snorm = cblas_ddot(node_num, x_prev, inc, x_prev, inc);
	// Dx
	for (int i = 0; i < node_num; ++i)
	{
		dx[i] = sI * x[i];
	}
	graph_cal_Dx(dx+node_num, x, edge_head, edge_tail, edge_num);


	int tmpidx;
//	double x_prev_snorm=0;
//	double tmpval, tmpsng;
	// specific operation for fused lasso (sparse data with l1 reg term)
	///////////
	double* tmpx;
	tmpx = new double[node_num];
//  specific operation for fused lasso
/*    memset(tmpx, 0, sizeof(double)*node_num);
	for (int i = 0; i < non_zero_num; ++i)
	{
		tmpidx = supp_idx_list[i];
		if(tmpidx<node_num)
			tmpx[tmpidx] = x[tmpidx];
	}
	cblas_dcopy(node_num, tmpx, inc, x, inc);
	///////////  */
	for (admm_ite = 0; admm_ite < admm_max_ite; ++admm_ite)
	{
		// z updating
		memset(z, 0, doublesizeS);
		for (int i = 0; i < non_zero_num; ++i)
		{// all nonzero elements
			tmpidx = supp_idx_list[i];
			if(dx[tmpidx]>=0)
				z[tmpidx] = dx[tmpidx] - lambda_rho;
			else
				z[tmpidx] = dx[tmpidx] + lambda_rho;
		}
		cblas_dcopy(s, z, inc, zu, inc);
		if(admm_ite>0)
		{
			cblas_daxpy(s, -1, u, inc, zu, inc);
		}
		// x updating
		solve_x_cg(x, ATy, zu, zu+node_num, edge_head, edge_tail, edge_num, node_num, sI, rho, cg_max_ite, cg_tol2);
		/////////// specific operation for fused lasso
		memset(tmpx, 0, sizeof(double)*node_num);
		for (int i = 0; i < non_zero_num; ++i)
		{
			tmpidx = supp_idx_list[i];
			if(tmpidx <node_num)
				tmpx[tmpidx] = x[tmpidx];
		}
		cblas_dcopy(node_num, tmpx, inc, x, inc);
		///////////  
        // dx
		for (int i = 0; i < node_num; ++i)
		{
			dx[i] = sI * x[i];
		}
		graph_cal_Dx(dx+node_num, x, edge_head, edge_tail, edge_num);

		// checking stopping condition
		x_snorm = cblas_ddot(node_num, x, inc, x, inc);
		if(admm_ite>0)
		{
			x_diff_snorm = x_snorm + x_prev_snorm - 2*cblas_ddot(node_num, x, inc, x_prev, inc);
			x_diff_rel = x_diff_snorm/x_prev_snorm;
			if(x_diff_rel<admm_tol2)
				break;
		}
		cblas_dcopy(node_num, x, inc, x_prev, inc);
		x_prev_snorm = x_snorm;
		// u (dual variable) updating 
		for (int i = 0; i < s; ++i)
		{
			u[i] += (m_rate*dx[i] - m_rate2*z_prev[i] - z[i]);
//			dx[i] += u[i];
//			z_prev[i] = z[i];
		}
		cblas_daxpy(s, 1, u, inc, dx, inc);
		cblas_dcopy(s, z, inc, z_prev, inc);
//		cblas_daxpy(s, 1, u, inc, dx, inc);
	}

	delete[] dx;
	delete[] zu;
	delete[] u;
	delete[] z_prev;
	delete[] x_prev;
	return 0;
}





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	double *x_in, *x_out, *ATy, *zuI, *zuG, *edge_head_tmp, *edge_tail_tmp;
	int *edge_head, *edge_tail, *supp_idx_list;
	int edge_num, node_num, non_zero_num,admm_max_ite, cg_max_ite;
	double sI, lambda, rho,admm_tol2, cg_tol2, cg_ite_out;

	x_in = (double*)mxGetPr(prhs[0]);
//	ATA  = (double*)mxGetPr(prhs[1]);
	ATy  = (double*)mxGetPr(prhs[1]);
	supp_idx_list = (int *)mxGetPr(prhs[2]);
	non_zero_num = (int)(*(double *)mxGetPr(prhs[3]));
	edge_head = (int*)mxGetPr(prhs[4]); 
	edge_tail = (int*)mxGetPr(prhs[5]);
	edge_num = (int)(*(double*)mxGetPr(prhs[6]));
	node_num = (int)(*(double*)mxGetPr(prhs[7]));
	sI  = *(double *)mxGetPr(prhs[8]);
	lambda = *(double *)mxGetPr(prhs[9]);
	rho = *(double *)mxGetPr(prhs[10]);
	admm_max_ite = (int)(*(double *)mxGetPr(prhs[11]));
	admm_tol2 = *(double *)mxGetPr(prhs[12]);
	cg_max_ite = (int)(*(double *)mxGetPr(prhs[13]));
	cg_tol2 = *(double *)mxGetPr(prhs[14]);
	

/*	edge_head_tmp = (double*)mxGetPr(prhs[5]); 
	edge_tail_tmp = (double*)mxGetPr(prhs[6]);
	edge_head = new int[edge_num];
	edge_tail = new int[edge_num];
	for (int i = 0; i < edge_num; ++i)
	{
		edge_head[i] = (int)edge_head_tmp[i];
		edge_tail[i] = (int)edge_tail_tmp[i];
	}*/


	plhs[0] = mxCreateDoubleMatrix(node_num,1,mxREAL); 
	x_out = (double *)mxGetPr(plhs[0]);
	for (int i = 0; i < node_num; ++i)
	{
		x_out[i] = x_in[i];
	}

/*	double *z;
	plhs[1] = mxCreateDoubleMatrix(node_num+edge_num,1,mxREAL); 
	z = (double *)mxGetPr(plhs[1]);
*/

	double *z, *u, *dx, *zu;
	int s = node_num+edge_num;
	dx = new double[s];
	zu = new double[s];
	z = new double[s];
	u = new double[s];

	solve_subadmm(x_out, z, u, dx, zu,
		ATy, supp_idx_list, non_zero_num, edge_head, edge_tail, edge_num, node_num, sI, 
		lambda, rho, admm_max_ite, admm_tol2, cg_max_ite, cg_tol2);

//	solve_x_cg(x_out, ATA, ATy, zuI, zuG, edge_head, edge_tail, edge_num, node_num, sI, rho, cg_max_ite, cg_tol2);
}



/*
int solve_subadmm(double *x, const double *ATA, const double *ATy, int *supp_idx_list, int non_zero_num,
	int *edge_head, int *edge_tail, int edge_num, int node_num, double sI, 
	double lambda, double rho, int admm_max_ite, double admm_tol, 
	int cg_max_ite, double cg_tol2)
{
	double *dxI, *dxG, *zI, *zG, *uI, *uG;
	double lambda_rho;

	lambda_rho = lambda/rho;
	dxI = new double[node_num];
	dxG = new double[edge_num];
	zI = new double[node_num];
	zG = new double[edge_num];
	uI = new double[node_num];
	uG = new double[edge_num];

	for (int admm_ite = 0; admm_ite < admm_max_ite; ++admm_ite)
	{
		// z updating
		for (int i = 0; i < non_zero_num; ++i)
		{
			
		}
	}



	return 0;
}

*/