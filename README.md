Planar Perspective Transformation
---------------------------------

Function 'blxPathPerspectiveTransform()' is an extension of the Blend2d library,
which allows you to apply a perspective transformation to the geometry of a BLPath.

	  ** WARNING **
	This perspective transformation is not completely correct from a mathematical
	 point of view, in the sense that the curves that form a BLPath are not
	 totally correctly deformed; only the control points of these curves are deformed.
	In practice, all the straight segments are correctly deformed, and the error
	 in the curves is perceptible only in the case of accentuated perspective deformations.

---
The core of this transformation is the function

	BLResult blxPathPerspectiveTransform(BLPath* path, const mat3x3x M)
where

	M is a 3x3 matrix (mat3x3) that can be obtained by calling the 'quadtoquad()' function.

Here is the quadtoquad() function:

	enum pmapType quadtoquad(const quad Q1, const quad Q2, mat3x3 M)
where

	'quad' (quadrilater) is an array of 4 '2D Points' (i.e an array of 8 double)
	'enum pmapType' is  { PMAP_AFFINE, PMAP_PERSPECTIVE, PMAP_BAD }.

Before calling 'quadtoquad' you must be sure that Q1 and Q2 are convex quadrilaters,
so you can call the function 'is_convex(Q)'

Please pay attention, after calling quadtoquad(), always check the returned pmapType;
if it's equal to PMAP_BAD, then there's no transformation from Q1 to Q2 and M is invalid.

Finally, once you get M, you can apply this transformation to a BlPath p
   blxPathPerspectiveTransform(&p, M)


Full Example
-----------
	#include "blx_path.h"  // Blend2d extension

	BLPath path1;
	 // prepare path1 ....
    path1.addCircle(BLCircle(100, 100, 90));
    path1.addCircle(BLCircle(100, 100, 70));

	 // quad vertices should be set in clockwise (or counter-clockwise) order
	 //  ??  what happens if Q1 and Q2 are listed with different clokwise order ?
	quad Q1 = { {  0,  0}, {200,  0}, {200,200}, {  0,200} };
	quad Q2 = { { 10, 20}, {200,  0}, {200,200}, {  0,200} };
	// optional : be sure Q1 and Q2 are CONVEX quads
	if ( !is_convex(Q1) || !is_convex(Q2) ) ... error ...

	mat3x3 M;
	// ---------------------------------------------
	enum pmapType mapType = quadtoquad( Q1, Q2, M );
	if (mapType == PMAP_BAD) ... error ...
	// ---------------------------------------------

	// apply the transformation
	//  Hint:  maybe you should work on a copy
	//    BLPath path2 = path1;
	blxPathPerspectiveTransform(&path2, M)
	// .. then you can can stroke path2 ...
