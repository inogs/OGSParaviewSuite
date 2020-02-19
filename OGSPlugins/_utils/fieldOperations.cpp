/*=========================================================================

  Program:   Field Operations
  Module:    fieldOperations.cpp

  Useful operations with field arrays.

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cstdlib>
#include <cmath>

#include "fieldOperations.h"

namespace field 
{
	/* COUNTDEPTHLEVELS

		Counts the number of depth levels (unique values in Z direction) and
		returns the values and the field connectivity.

		If uniquevals is not empty, it will return the depth levels at the
		desired uniquevals.

	*/
	field::Field<int> countDepthLevels(v3::V3v &xyz, std::vector<double> &uniquevals, double epsi, bool addlayer) {

		field::Field<int> cId2zId(xyz.len(),1);

		// Loop the mesh using iterators
		v3::V3v::iterator itxyz;
		field::Field<int>::iterator itId;

		if (uniquevals.empty()) {
			// This is the case where the user has not inputed any value of z
			// therefore we must search for all the depth levels
			for (itId = cId2zId.begin(), itxyz = xyz.begin(); itxyz != xyz.end(); 
				++itxyz, ++itId) {
				// Start by assuming that the point is unique
				// We shall compare it with the points stored in uniquevals
				bool isunique = true;
				for (double zunique : uniquevals) { // Range based foor loop here
					// If we have already stored the point, then it is not unique
					if ( std::fabs(zunique - itxyz[2]) < epsi) { isunique = false; break; }
				} 
				// If the point is unique, increase the counter and store it
				if (isunique) uniquevals.push_back(itxyz[2]);
				// Now, set the connectivity
				itId[0] = -1;
				for (int ii = 0; ii < uniquevals.size(); ii++)
					if (std::fabs(uniquevals[ii] - itxyz[2]) < epsi)
						itId[0] = ii;
			}
		} else {
			// Add one extra depth level in uniquevals corresponding to the bottom
			if (addlayer) uniquevals.push_back( 2.*xyz[-1][2] );
			// The user has inputed a range of depth levels, we shall search within
			// these range and set the connectivity matrix appropriately.
			for (itId = cId2zId.begin(), itxyz = xyz.begin(); itxyz != xyz.end(); 
				++itxyz, ++itId) {
				// Loop on the unique values and set the connectivity
				itId[0] = -1;
				for (int ii = 0; ii < uniquevals.size(); ++ii)
					if (itxyz[2] > uniquevals[ii]) { itId[0] = ii; break; }
				// If the point has not been set, default to the last one
				itId[0] = (itId[0] < 0) ? uniquevals.size()-1 : itId[0];
			}
		}

		return cId2zId;
	}

}