/* 
  ==========================================================================
  quadmap.h
  ==========================================================================
*/
#ifndef __QUADMAP_INCLUDED
#define __QUADMAP_INCLUDED

/*
 * For C++ compilers, use extern "C"
 */
#ifdef __cplusplus
extern "C" {
#endif

/*
**  A quadrilater
**  4 points listed in a circular order
**    P0 --- P1
**     |     |
**    P3 --- P2
**  It's not important the clockwise/counter-clockwise direction,
**  nor if P0.y is greater or lesser than P3.y .
**  The only important thing is that P0 and P2 are opposite corners (so, P1 and P3, too)
*/
typedef struct {double x, y;} quad[4];

typedef double mat3x3[3][3];

enum pmapType {
	PMAP_BAD,
	PMAP_AFFINE,
	PMAP_PROJECTIVE
};

enum pmapType quadtoquad(const quad A, const quad B, mat3x3 M);

int is_quad_convex(const quad A);
int is_affine(const mat3x3 M);

/*
 * end block for C++
 */
#ifdef __cplusplus
}
#endif

#endif
