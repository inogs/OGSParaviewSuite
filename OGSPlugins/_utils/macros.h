/*=========================================================================

  Module:    Macros

  Define macros for all the ParaView filters. 

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/* ARRAY PRECISION

	The macros here define the array precision for all the
	ParaView filters.
		> FLDARRAY: Field arrays precision
		> VTKARRAY: VTK array precision
		> FLDMASK:  Mask array precision (uint8_t = 8 byte)
		> VTKMASK:  VTK mask array precision
*/
#define FLDARRAY double
#define VTKARRAY vtkDoubleArray
#define FLDMASK uint8_t
#define VTKMASK vtkTypeUInt8Array