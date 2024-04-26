/* 
  ==========================================================================
  Planar Perspective Transformation
  - Blend2d extension
  ---
   2024 - A.Buratti fecit
  ==========================================================================
*/

#include <blend2d.h>
#include "blx_path.h"

BLResult blxPathPerspectiveTransform( BLPathCore* self, const mat3x3 M ) {

	const BLPoint* P = blPathGetVertexData(self);
	size_t N = blPathGetSize(self);
	for (size_t i = 0; i<N; i++, ++P) {
		double x = P->x*M[0][0] + P->y*M[1][0] + M[2][0];
		double y = P->x*M[0][1] + P->y*M[1][1] + M[2][1];
		double w = P->x*M[0][2] + P->y*M[1][2] + M[2][2];

		blPathSetVertexAt(self, i, BL_PATH_CMD_PRESERVE, x/w, y/w);
	}
	return BL_SUCCESS;
}
