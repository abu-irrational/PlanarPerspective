/*
**
	Algorithm from https://www.cs.cmu.edu/~ph/texfund/texfund.pdf
		Fundamentals of Texture Mapping and Image Warping, Paul Heckbert
	Original Code:
		https://www.cs.cmu.edu/~ph/src/texfund/
**
*/

#include "quadmap.h"

static
double det2x2(double a, double b, double c, double d) {
	return a * d - b * c;
}

/*
  Compute the adjugate matrix of A
  In homogeneous algebra, the adjugate matrix always exists,
  while the inverse does not if the matrix is singular (i.e if determinant is 0).

    A x Adj(A) = det(A)*I
	then  if det(A) is not zero,
	inv(M) = Adj(M)/det(A).

	== return the adjugate matrix (by argument) AND the determinat of A
*/
static
double adj3x3(const mat3x3 A, mat3x3 B) {
	B[0][0] = det2x2(A[1][1], A[1][2], A[2][1], A[2][2]);
	B[1][0] = det2x2(A[1][2], A[1][0], A[2][2], A[2][0]);
	B[2][0] = det2x2(A[1][0], A[1][1], A[2][0], A[2][1]);

	B[0][1] = det2x2(A[2][1], A[2][2], A[0][1], A[0][2]);
	B[1][1] = det2x2(A[2][2], A[2][0], A[0][2], A[0][0]);
	B[2][1] = det2x2(A[2][0], A[2][1], A[0][0], A[0][1]);

	B[0][2] = det2x2(A[0][1], A[0][2], A[1][1], A[1][2]);
	B[1][2] = det2x2(A[0][2], A[0][0], A[1][2], A[1][0]);
	B[2][2] = det2x2(A[0][0], A[0][1], A[1][0], A[1][1]);
	 // return the A's determinant.
	return A[0][0] * B[0][0] + A[0][1] * B[0][1] + A[0][2] * B[0][2];
}

/*
**  3x3 Matrix multiplication :  C = A x B
*/
static
void mul3x3(const mat3x3 A, const mat3x3 B, mat3x3 C) {
	C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
	C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
	C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

	C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
	C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
	C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

	C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
	C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
	C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}


#define TOLERANCE 1e-13
#define ZERO(x) ((x)<TOLERANCE && (x)>-TOLERANCE)

static
enum pmapType pmap_square_quad(const quad Q, mat3x3 SQ) {
	double px, py;

	px = (Q[0].x - Q[1].x) - (Q[3].x - Q[2].x);
	py = (Q[0].y - Q[1].y) - (Q[3].y - Q[2].y);

	// Note: we presume Q is not a degenere-quad
	if (ZERO(px) && ZERO(py)) {
		SQ[0][0] = Q[1].x - Q[0].x ;
		SQ[1][0] = Q[2].x - Q[1].x;
		SQ[2][0] = Q[0].x;

		SQ[0][1] = Q[1].y - Q[0].y;
		SQ[1][1] = Q[2].y - Q[1].y;
		SQ[2][1] = Q[0].y;

		SQ[0][2] = 0.0;
		SQ[1][2] = 0.0;
		SQ[2][2] = 1.0;
		return PMAP_AFFINE;
	}
	// else ..
	double dx1, dx2, dy1, dy2, del;

	dx1 = Q[1].x - Q[2].x;
	dx2 = Q[3].x - Q[2].x;
	dy1 = Q[1].y - Q[2].y;
	dy2 = Q[3].y - Q[2].y;
	del = det2x2(dx1, dx2, dy1, dy2);
	if (del == 0.0) {
		//  ERROR("pmap_square_quad: bad mapping\n");
		return PMAP_BAD;
	}
	SQ[0][2] = det2x2(px, dx2, py, dy2) / del;
	SQ[1][2] = det2x2(dx1, px, dy1, py) / del;
	SQ[2][2] = 1.0;
	SQ[0][0] = Q[1].x - Q[0].x + SQ[0][2] * Q[1].x;
	SQ[1][0] = Q[3].x - Q[0].x + SQ[1][2] * Q[3].x;
	SQ[2][0] = Q[0].x;
	SQ[0][1] = Q[1].y - Q[0].y + SQ[0][2] * Q[1].y;
	SQ[1][1] = Q[3].y - Q[0].y + SQ[1][2] * Q[3].y;
	SQ[2][1] = Q[0].y;
	return PMAP_PROJECTIVE;
}

/*
 * map_quad_rect: find mapping between quadrilateral and rectangle.
 * The correspondence is:
 *
 *      quad[0] --> (u0,v0)
 *      quad[1] --> (u1,v0)
 *      quad[2] --> (u1,v1)
 *      quad[3] --> (u0,v1)
 *
 * This method of computing the adjugate-matrix  numerically is cheaper than
 * computing it symbolically.
 */
static
enum pmapType pmap_quad_rect(const double u0, const double v0, const double u1, const double v1,
	const quad Q,
	mat3x3 QR) {
	enum pmapType ret;
	double du, dv;
	mat3x3 RQ;   // Rect to Quad transform

	du = u1 - u0;
	dv = v1 - v0;
	if (du == 0.0 || dv == 0.0) {
		// ERROR("pmap_quad_rect: null rectangle\n");
		return PMAP_BAD;
	}
    /* first,  find mapping from unit uv square to quad Q */
	ret = pmap_square_quad(Q, RQ);
	if (ret == PMAP_BAD) return PMAP_BAD;

    /* concatenate transform from uv rectangle (u0,v0,u1,v1) to unit square */
	RQ[0][0] /= du;
	RQ[1][0] /= dv;
	RQ[2][0] -= RQ[0][0] * u0 + RQ[1][0] * v0;

	RQ[0][1] /= du;
	RQ[1][1] /= dv;
	RQ[2][1] -= RQ[0][1] * u0 + RQ[1][1] * v0;

	RQ[0][2] /= du;
	RQ[1][2] /= dv;
	RQ[2][2] -= RQ[0][2] * u0 + RQ[1][2] * v0;

    /* now RQ is transform from uv rectangle to xy quadrilateral */
    /* QR = inverse transform, which maps xy to uv */
	if (adj3x3(RQ, QR) == 0.0) {
		// ERROR("pmap_quad_rect: warning - determinant is 0\n");   //    ..  ?? ma non ritorni errore ???
		return PMAP_BAD;
	}
	return ret;
}

/*
 * map_quad_quad: find the projective mapping ST from quad A to quad B.
 *
 * This method is faster than solving an 8x8 system of equations.
 */
static
enum pmapType pmap_quad_quad(const quad A, const quad B, mat3x3 ST) {
	enum pmapType type1, type2;
	mat3x3 MS;
	mat3x3 SM;
	mat3x3 MT;

	type1 = pmap_square_quad(A, MS);
	if ( type1 == PMAP_BAD)  return PMAP_BAD;

	 // Since MS in NOT PMAP_BAD, the resulting adjugate matrix SM = k * inv(MS)
	if (adj3x3(MS, SM) == 0.0) {
		// ERROR("pmap_quad_quad: warning - determinant is 0\n");
		return PMAP_BAD;
	}
	type2 = pmap_square_quad(B, MT);
	if (type2 == PMAP_BAD)  return PMAP_BAD;

	mul3x3(SM, MT, ST);
	return (type1 == PMAP_AFFINE && type2 == PMAP_AFFINE) ? PMAP_AFFINE : PMAP_PROJECTIVE;
}

/*
** pmap_quads: find the projective mapping ST from quad A to quad B
*/
static
enum pmapType _quadtoquad(const quad A, const quad B, mat3x3 ST) {
	// if B's edges 0-1 and 2.3 are horizontal and 1-2 and 3-0 are vertical
	// then B is rectangle
	if (B[0].y == B[1].y && B[2].y == B[3].y && B[1].x == B[2].x && B[3].x == B[0].x) {
		return pmap_quad_rect(B[0].x, B[0].y, B[2].x, B[2].y, A, ST);
	}

	// if B's edges 0-1 and 2.3 are vertical and 1-2 and 3-0 are horizontal
	// then B is a (rotated) rectangle (-> A should be rotated, too )
	if (B[0].x == B[1].x && B[2].x == B[3].x && B[1].y == B[2].y && B[3].y == B[0].y) {
		quad C;
		C[0] = A[1];
		C[1] = A[2];
		C[2] = A[3];
		C[3] = A[0];
		return pmap_quad_rect(B[1].x, B[1].y, B[3].x, B[3].y, C, ST);
	}
	// else, B is not an orthogonally-oriented rectangle
	return pmap_quad_quad(A, B, ST);
}

enum pmapType quadtoquad(const quad A, const quad B, mat3x3 M) {
	enum pmapType res = _quadtoquad(A, B, M);
	 // Normalize the matrix
	if ( res != PMAP_BAD ) {
		double m22 = M[2][2];
		if (m22==0.0) return PMAP_BAD;
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				M[i][j] = M[i][j]/m22;
			}
		}
	}
	return res;
}

  // prereq: M is a (not PMAP_BAD) Normalized  matrix (m22==1)
int is_affine(const mat3x3 M) {
	return M[0][2] == 0 && M[1][2] == 0;
}

  /*
  ** Code grabbed and adaptet from "Graphics Gems III"
  ** Faster Line Segment Intersection
  **   by Franklin Antonio
  **
  ** https://github.com/erich666/GraphicsGems/blob/master/gemsiii/insectc.c
  **
  ** Changes:
  **  + works with float numbers.
  **  + don't compute the intersection. Only returns 0 1
  */

#define DONT_INTERSECT 0
#define DO_INTERSECT   1
#define PARALLEL       2

static
int segment_intersection(
	double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4
	) {
	double Ax = x2-x1;
	double Bx = x3-x4;
	{
		double x1lo, x1hi;
		if(Ax<0) {						/* X bound box test*/
		  x1lo=x2; x1hi=x1;
		} else {
		  x1hi=x2; x1lo=x1;
		}
		if(Bx>0) {
		  if(x1hi < x4 || x3 < x1lo) return DONT_INTERSECT;
		} else {
		  if(x1hi < x3 || x4 < x1lo) return DONT_INTERSECT;
		}
	}

	double Ay = y2-y1;
	double By = y3-y4;
	{
		double y1lo, y1hi;
		if(Ay<0) {						/* Y bound box test*/
		  y1lo=y2; y1hi=y1;
		} else {
		  y1hi=y2; y1lo=y1;
		}
		if(By>0) {
		  if(y1hi < y4 || y3 < y1lo) return DONT_INTERSECT;
		} else {
		  if(y1hi < y3 || y4 < y1lo) return DONT_INTERSECT;
		}
	}

	double Cx = x1-x3;
	double Cy = y1-y3;

	double f = Ay*Bx - Ax*By;					/* both denominator*/

	if(f==0) return PARALLEL;

	double d = By*Cx - Bx*Cy;					/* alpha numerator*/
	if(f>0) {						/* alpha tests*/
	  if(d<0 || d>f) return DONT_INTERSECT;
	} else {
	  if(d>0 || d<f) return DONT_INTERSECT;
	}

	double e = Ax*Cy - Ay*Cx;					/* beta numerator*/
	if(f>0) {						/* beta tests*/
	  if(e<0 || e>f) return DONT_INTERSECT;
	} else {
	  if(e>0 || e<f) return DONT_INTERSECT;
	}

#if 0    //don't care for the intersection coordinates

#define SAME_SIGNS( a, b ) sign(a)==sign(b)

	num = d*Ax;						/* numerator */
	offset = SAME_SIGNS(num,f) ? f/2 : -f/2;		/* round direction*/
	*x = x1 + (num+offset) / f;				/* intersection x */

	num = d*Ay;
	offset = SAME_SIGNS(num,f) ? f/2 : -f/2;
	*y = y1 + (num+offset) / f;				/* intersection y */

#endif
	return DO_INTERSECT;
}

int is_quad_convex(const quad A) {
	// is_convex <==> quad's diagonals intersect
	return segment_intersection(
		A[0].x, A[0].y, A[2].x, A[2].y,   // first diagonal
		A[1].x, A[1].y, A[3].x, A[3].y   // second diagonal
		);
}


// -  TEST -------------------------------------------------------------------

#if 0

void matPrint(const mat3x3 M) {
	printf("--------------\n");
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			printf("%10.2lf ", M[i][j]);
		}
		printf("\n");
	}
}

int main() {
	quad A = { 0,0,   1,0,   1,1,   0,1 };
	quad B = { 10,10, 20,10, 20,20, 10,18 };

	mat3x3 T;

	quadtoquad(A, B, T);
	matPrint(T);
}

#endif
