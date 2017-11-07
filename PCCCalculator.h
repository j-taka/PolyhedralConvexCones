// PCCcalculator.h
/* calculate dual of Polyhedral Convex Cones (PCCs)
 *
 * Hirukawa et.al: "‘½–Ê‘ÌŠÔ‚ÌÚG‚É‚æ‚éS‘©ğŒ‚ÌS‘©‰ğ–@‚Æ
 * ‚»‚Ì—£’E“®ìŒv‰æ‚Ö‚Ì‰—p", RSJ1991 ‚ğƒAƒŒƒ“ƒW
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

static const double _EPS = 1.0e-6;

/* class for calculating dual of PCCs */
class PCCCalculator
{
protected:
	int esrow;
	int cvDim; /* dimension */
	double SVDEPS; 

public:
	// constructor
	PCCCalculator() : SVDEPS(1.0e-6) {}
	// dof of detaching directions
	int DDOF() const;

	/*! calcualte dual of pcc
	* @param dest  : ˆÛ‰^“®(+‚à-‚àŠî’ê)
	* @param dest2 : —£’E‰^“®(+‚Ì‚İ‚ªŠî’ê)
	* @param src   : dual‚ğŒvZ‚µ‚½‚¢s—ñ
	*/
	void Dual(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, const Eigen::MatrixXd &src);
	/* dual‚ğŒvZ‚·‚é
	* @param dest  : ˆÛ‰^“®(+‚à-‚àŠî’ê)
	* @param dest2 : —£’E‰^“®(+‚Ì‚İ‚ªŠî’ê)
	* @param src   : “™®•”•ª
	* @param src2  : •s“™®•”•ª
	*/
	void Dual(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, const Eigen::MatrixXd &src, const Eigen::MatrixXd &src2);

protected:
	/* calculate */
	void Calc(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, const Eigen::MatrixXd &src);
	void DivideZeroSpace(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src) const;

private:
	/* divide zero space using singular value decomposition
	 * @param dest  basis of zero space
 @  * @param dest2 the matrix to embed to non-zero space 
	 * @param dest3 candidates of PCC's hyperplanes
	 * @param src   simultaneous inequarities
	 */
	void DivideZeroSpace(Eigen::MatrixXd &dest, Eigen::MatrixXd &dest2, Eigen::MatrixXd &dest3, const Eigen::MatrixXd &src) const;
	/* calculate edges of the PCC
  	 * @param dest  edges
	 * @param dest2 candidates of hyperplane of intersection
	 * @param src   candidates of PCC's hyperplanes
	 */
	void CalcIntersection(std::vector<Eigen::VectorXd> &dest, std::vector<bool> &dest2, const Eigen::MatrixXd &src) const;
	/* calculate hyperplane H_{i} of intersection
	 * @param dest hyperplane of intersection
	 * @param src  edges of the PCC
	 */
	void CalcHyperPlane(Eigen::VectorXd &dest, const std::vector<bool> &src2, const Eigen::MatrixXd &src) const;
	/* calculat convex hull
 	 * @param dest points on convex hull
	 * @param src  surface normal of the hyperplane
	 * @param src2 edges of the PCC
	 */
	void CalcConvexHull(Eigen::MatrixXd &dest, const Eigen::VectorXd &src, const std::vector<Eigen::VectorXd> &src2);
	/* calculate convex hull in 1D (i.e., calculate minimum and maximum)
	 * @param dest points on convex hull
	 * @param src  points
	 */
	void Calc1DConvexHull(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src) const;
	/* calculate convex hull using qhull
	 * @param dest points on convex hull
	 * @param src  points
	 */
	void CalcNDConvexHull(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src) const;
	/* transform in the original coordinates
	 * @param dest dual of PCC
	 * @param src  points on convex hull
	 * @param src2 V_{1} transformation matrix
	 */
	void CoordTrans(Eigen::MatrixXd &dest, const Eigen::MatrixXd &src, const Eigen::MatrixXd &src2) const;
	static double RZ(double src) {
		return (fabs(src) < _EPS ? 0 : src);
	}
	static Eigen::MatrixXd RZ(const Eigen::MatrixXd &src) {
		Eigen::MatrixXd dest(src.rows(), src.cols());
		for (int i(0); i < src.rows(); ++i) {
			for (int j(0); j < src.cols(); ++j) {
				dest(i, j) = RZ(src(i, j));
			}
		}
		return dest;
	}
};
