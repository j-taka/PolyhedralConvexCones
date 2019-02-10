// PCCCalculator.cpp

#include "PCCCalculator.h"
#include "Combination.h"

// #define _PCCDEBUG

#ifdef WIN32
#define NRANSI
#ifdef __cplusplus
extern "C" {
#endif   
#endif
#include <libqhull/qhull_a.h>
qhT qh_qh;
#ifdef WIN32
#ifdef __cplusplus
}
#endif   
#endif




#ifdef _PCCDEBUG
#include <iostream>
#endif

int PCCCalculator::DDOF() const
{
	return cvDim + 1;
}

void PCCCalculator::Dual(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, const Eigen::MatrixXd &src)
{
	esrow = src.rows();
	if (src.rows() >= src.cols()) {
		// main
		Calc(dest, dest2, src);
	}
	else {
		Eigen::MatrixXd tmp(src.cols(), src.cols());
		tmp.block(0, 0, src.rows(), src.cols()) = src;
		for (int i(src.rows()); i < src.cols(); ++i) {
			tmp.row(i) = tmp.row(i - 1);
		}
		// main
		Calc(dest, dest2, tmp);
	}
}

void PCCCalculator::Dual(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, const Eigen::MatrixXd &src, const Eigen::MatrixXd &src2)
{
	assert(src.rows() + src2.rows() > 0); // no data
	assert(src.cols() == src2.cols()); // same size
	// copy
	esrow = src.rows() * 2 + src2.rows();
	const int mat_row = (esrow >= src.cols() ? esrow : src.cols());
	Eigen::MatrixXd v_mat(mat_row, src.cols());
	for (int i(0); i < src.rows(); ++i) {
		v_mat.row(2 * i    ) =  src.row(i);
		v_mat.row(2 * i + 1) = -src.row(i);
	}
	for (int i(0); i < src2.rows(); ++i) {
		v_mat.row(2 * src.rows() + i) = src2.row(i);
	}
	for (int i(esrow); i < v_mat.rows(); ++i) {
		v_mat.row(i) = v_mat.row(i - 1);
	}
	// main 
	Calc(dest, dest2, v_mat);
}

void PCCCalculator::Calc(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, const Eigen::MatrixXd &src)
{
	// clear
	cvDim = -1;
	/* divide zero space */
	Eigen::MatrixXd V1, Hi;
	DivideZeroSpace(dest, V1, Hi, src);
#ifdef _PCCDEBUG
	std::cout << V1 << std::endl << std::endl;
	std::cout << Hi << std::endl << std::endl;
	std::cout << dest << std::endl << std::endl;
#endif 	
	if (src.cols() - dest.rows() != 1) { 
		/* calculate candidates of planes */
		std::vector<Eigen::VectorXd> edges;
		std::vector<bool> calcH;
		CalcIntersection(edges, calcH, Hi);
		if (!edges.empty()) {
#ifdef _PCCDEBUG
			for (size_t i(0); i < edges.size(); ++i) {
				std::cout << edges[i].transpose() << std::endl;
			}
			std::cout << "calculate hyper plane" << std::endl;
#endif
			// calculate hyperplance of intersection to calculate convex hull
			Eigen::VectorXd H;
			CalcHyperPlane(H, calcH, Hi);
#ifdef _PCCDEBUG
			std::cout << H.transpose() << std::endl << std::endl;
			std::cout << "calculate convex hull" << std::endl;
#endif 
			// calculate convex hull
			Eigen::MatrixXd conv;
			CalcConvexHull(conv, H, edges);
#ifdef _PCCDEBUG
			std::cout << conv << std::endl << std::endl;
			std::cout << "transform coordinate system" << std::endl;
#endif 
			CoordTrans(dest2, conv, V1);
		}
	}
	else {
		// spacial operation in the case RANK is one
		int i(1);
		for (; i < Hi.rows(); i++) {
			if (Hi.row(0).dot(Hi.row(i)) < 0) { break; }
		}
		if (i >= Hi.rows()) {
			Eigen::VectorXd temp = Hi(0, 0) * V1.col(0);
			temp /= temp.norm(); /* normalize */
			dest2 = Eigen::MatrixXd(1, temp.size());
			dest2.row(0) = temp.transpose();
			cvDim = 0;
		}
		else {
			dest2.resize(0, 0);
		}
	}
}

void PCCCalculator::DivideZeroSpace(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, Eigen::MatrixXd &dest3, const Eigen::MatrixXd &src) const
{
	// singular value decomposition
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(src, Eigen::ComputeThinU | Eigen::ComputeThinV);
	const double wmax = svd.singularValues().maxCoeff();
	const double wmin = wmax * SVDEPS; // zero space criterion
	// count
	int dim(0);
	for (int i(0); i < svd.singularValues().size(); ++i) {
		if (svd.singularValues()[i] < wmin) { // basis of the zero space
			dim++;
		}	
	}
	// assign
	dest = Eigen::MatrixXd(dim, src.cols());
	for (int i(0), r(0); i < svd.singularValues().size(); ++i) {
		if (svd.singularValues()[i] < wmin) { // basis of the zero space
			dest.row(r) = RZ(svd.matrixV().col(i).transpose());
			r++;
		}
	}
	/* basis of non-zero space */
	dest2 = Eigen::MatrixXd(src.cols(), src.cols() - dim);
	dest3 = Eigen::MatrixXd(esrow, src.cols() - dim);
	for (int i(0), c(0); i < src.cols(); i++) {
		if (svd.singularValues()[i] >= wmin) { // non-zero
			dest2.col(c) = RZ(svd.matrixV().col(i));
			dest3.col(c) = RZ(svd.matrixU().block(0, i, esrow, 1) * svd.singularValues()[i]);
			c++;
		}
	}
	// normalize
	for (int i(0); i < dest3.rows(); ++i) {
		dest3.row(i) /= dest3.row(i).norm();
	}
}

#if 0
/* ì¡àŸílï™âÇµ, É[ÉçãÛä‘Çï™ó£
* @param dest  : É[ÉçãÛä‘ÇÃäÓíÍ(à€éùâ^ìÆ)
* @param src   : î˜è¨ïœà ïsìôéÆ(ílÇ™ïœÇÌÇ¡ÇƒÇµÇ‹Ç§)
*/
void PCCCalculator::
DivideZeroSpace(CVLVectdList &dest, double **src) const
{
	double **v, *w, wmax = 0.0, wmin;
	int i, j;
	CVLVectd temp(col);
	/* ÉÅÉÇÉäämï€ */
	v = dmatrix(1, col, 1, col);
	w = dvector(1, col);
	/* ì¡àŸílï™â */
	dsvdcmp(src, row, col, w, v);
	/* É[ÉçãÛä‘ÇÃíäèo */
	for (i = 1; i <= col; i++) if (w[i]>wmax) wmax = w[i];
	wmin = wmax*MTSSVDEPS; /* Ç«ÇÃîÕàÕÇ≈É[ÉçãÛä‘Ç∆å©ÇÈÇ© */
	for (i = 1; i <= col; i++)
		if (w[i]<wmin) { /* É[ÉçãÛä‘ÇÃäÓíÍ(ìÆÇ©Ç»Ç¢ï™ÇÃÇ‡â¡Ç¶Çƒ) */
			for (j = 0; j<col; j++) temp[j] = v[j + 1][i];
			dest.push_back(RZ(temp));
		}
	/* ÉÅÉÇÉäâï˙ */
	free_dmatrix(v, 1, col, 1, col);
	free_dvector(w, 1, col);
}
#endif

void PCCCalculator::CalcHyperPlane(Eigen::VectorXd &dest, const std::vector<bool> &src2, const Eigen::MatrixXd &src) const
{
	// initialize
	dest = Eigen::VectorXd::Zero(src.cols());
	for (int i(0); i < src.rows(); i++) {
		dest += src.row(i).transpose();
	}
	/* possibility if zero */
	if (fabs(dest.norm()) < _EPS) {
		int i(0);
		for (; i < src.rows(); i++){
			if (src2[i]) {
				dest += src.row(i).transpose();
				break;
			}
		}
		assert(i != src.rows()); // choose at least one of them
	}
	dest /= dest.norm(); // normalize
	dest = RZ(dest);
}

void PCCCalculator::CalcIntersection(std::vector<Eigen::VectorXd> &dest, std::vector<bool> &dest2, const Eigen::MatrixXd &src) const
{
	const int r = src.cols();
	Combination comb(src.rows(), r - 1);
	// initialize
	dest2.resize(src.rows(), false);
	for (;;) {
		Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(r - 1, r);
		int res(0);
		for (; res < r - 1; res++) {
			mat.row(res) = src.row(comb[res]);
			// calculate rank
			Eigen::FullPivLU<Eigen::MatrixXd> lu(mat);
			// should set the appropriate threshold
			// lu.setThreshold();
			if (lu.rank() < res + 1) {
				break;
			}
		}
		if (res == r - 1) {
			// calculate
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(mat.transpose() * mat);
			// sorted increasing order 
			if (eig.eigenvalues()[1] > eig.eigenvalues()[r - 1] * SVDEPS * SVDEPS) {
				// rank is r - 1
				Eigen::VectorXd edge = RZ(eig.eigenvectors().col(0)); // eigen vectors 
				// check direction
				int sign(0), i(0);
				for (; i < src.rows(); ++i) {
					const double cur = src.row(i) * edge;
					if (fabs(cur) > _EPS) {
						if (sign == 0) {
							sign = (cur > 0 ? +1 : -1);
						}
						else if ((double)sign * cur < 0) {
							break;
						}
					}
				}
				if (i == src.rows()) {
					dest.push_back(edge);
					for (int j(0); j < src.rows(); ++j) {
						if (fabs(src.row(j) * edge) > _EPS) {
							dest2[j] = true;
						}
					}
				}
			}
		}
		if (!comb.Next(r - 1 - res)) return;
	}
}

void PCCCalculator::CalcConvexHull(Eigen::MatrixXd &dest, const Eigen::VectorXd &src, const std::vector<Eigen::VectorXd> &src2)
{
	if (src2.empty()) {
		return; // no data
	}
	// calculat intersection
	Eigen::MatrixXd pts(src2.size(), src2.front().size());
	for (size_t i(0); i < src2.size(); ++i){
		pts.row(i) = 1.0 / (src2[i].dot(src)) * src2[i].transpose();
	}
	if (pts.rows() == 1) {
		dest = pts;
		cvDim = 0;
		return;
	}
	// dimension of convex hull
	// average
	Eigen::VectorXd ave = Eigen::VectorXd::Zero(src2.front().size());
	for (int i(0); i < pts.rows(); ++i) {
		ave += pts.row(i).transpose();
	}
	ave /= pts.rows();
	// subtract ave
	pts.rowwise() -= ave.transpose();
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(pts.transpose() * pts);
	// sort increasing order
	const int dim = eig.eigenvalues().size();
	// special case
	if (eig.eigenvalues()[dim - 1] < SVDEPS) {
		dest = Eigen::MatrixXd(1, ave.size());
		dest.row(0) = ave.transpose();
		return;
	}
	for (int i(1); i < dim; ++i) {
		if (eig.eigenvalues()[dim - i - 1] < eig.eigenvalues()[dim - 1] * SVDEPS) {
			cvDim = i;
			break;
		}
	}
	Eigen::MatrixXd R(cvDim, dim);
	for (int i(0); i < cvDim; ++i) {
		R.row(i) = RZ(eig.eigenvectors().col(dim - i - 1).transpose());
	}
	const Eigen::MatrixXd cvPos = pts * R.transpose();
	// convex hull
	Eigen::MatrixXd ch;
	if (cvDim == 1) {
		Calc1DConvexHull(ch, cvPos);
	}
	else {
		CalcNDConvexHull(ch, cvPos);
	}
	/* transform in the original coordinates */
	dest = ch * R;
	dest.rowwise() += ave.transpose();
}

void PCCCalculator::Calc1DConvexHull(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src) const
{
	const double min_v = src.minCoeff();
	const double max_v = src.maxCoeff();
	if (max_v - min_v > _EPS) {
		dest = Eigen::MatrixXd(2, 1);
		dest(0, 0) = max_v;
		dest(1, 0) = min_v;
	}
	else {
		dest = Eigen::MatrixXd(1, 1);
		dest(0, 0) = max_v;
	}
}

void PCCCalculator::CalcNDConvexHull(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src) const
{
	QHULL_LIB_CHECK

	coordT *points, *point;
	boolT ismalloc = False;
	char flags[250];
#ifdef _PCCDEBUG
	FILE *outfile = stdout, *errfile = stderr;
#else
	FILE *outfile = NULL, *errfile = stderr;
#endif
	vertexT *vertex;
	int curlong, totlong, num = src.rows(), dim = src.cols(), exitcode;
	/* memory allocation */
	points = new coordT[dim * num];
	sprintf(flags, "qhull "); /* set option */
	point = points;
	/* set */
	for (int i(0); i < src.rows(); ++i) {
		for (int j(0); j < src.cols(); j++, point++) {
			*point = src(i, j);
		}
	}
	exitcode = qh_new_qhull(dim, num, points, ismalloc, flags, outfile, errfile);
	if (!exitcode) {                  /* if no error */
		int count(0);
		dest = Eigen::MatrixXd(qh num_vertices, dim);
		FORALLvertices{
			for (int i(0); i < qh hull_dim; ++i) {
				dest(count, i) = vertex->point[i];
			}
			count++;
		}
	}
	qh_freeqhull(!qh_ALL);                   /* free long memory  */
	qh_memfreeshort(&curlong, &totlong);    /* free short memory and memory allocator */
	if (curlong || totlong)
		fprintf(errfile, "qhull internal warning (user_eg, #1): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
	delete[]points;
}

void PCCCalculator::CoordTrans(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src, const Eigen::MatrixXd &src2) const
{
	dest = src * src2.transpose();
	dest.rowwise().normalize();
	dest = RZ(dest);
}
