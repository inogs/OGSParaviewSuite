/*=========================================================================

  Program:   OGSDensity
  Module:    vtkOGSDensity.cxx

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOGSDensity.h"

#include "vtkDataSet.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include <string>
#include <sstream>
#include <vector>

#ifdef __linux__
// Include OpenMP when working with GCC
#include <omp.h>
#define OMP_NUM_THREADS omp_get_num_threads()
#define OMP_THREAD_NUM  omp_get_thread_num()
#else
#define OMP_NUM_THREADS 1
#define OMP_THREAD_NUM  0
#endif

#define POW2(x) (x)*(x)
#define POW3(x) (x)*(x)*(x)
#define POW4(x) (x)*(x)*(x)*(x)
#define POW5(x) (x)*(x)*(x)*(x)*(x)

#define Z2P(x) -9.80665*997.0474*1e-5*(x) // meter 2 bar (correct the sign as z are negative)

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSDensity, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(vtkOGSDensity);

//----------------------------------------------------------------------------

#include "macros.h"
#include "V3.h"
#include "field.h"
#include "vtkFields.h"
#include "vtkOperations.h"

//----------------------------------------------------------------------------

void rhoLinT(field::Field<FLDARRAY> &rho, field::Field<FLDARRAY> &T, double ralpha, double rau0) {
/*
	linear equation of state function of temperature only
		rdn(t) = ( rho(t) - rau0 ) / rau0 = 0.028 - ralpha * t
		rhop(t,s)  = rho(t,s)

	from OGSTM PHYS/eos.f90
*/
	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<rho.get_n(); ii+=OMP_NUM_THREADS) {
		double rdn = 0.028 - ralpha*T[ii][0]; // Density anomaly
		rho[ii][0] = rau0*rdn + rau0;
	}
	}
}

void rhoLinTS(field::Field<FLDARRAY> &rho, field::Field<FLDARRAY> &T, field::Field<FLDARRAY> &S,
	double ralpha, double rbeta, double rau0) {
/*
	Atmosphere, Ocean and Climate Dynamics, Marshall & Plumb

	Equation of state eq. (9-5)
*/
	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<rho.get_n(); ii+=OMP_NUM_THREADS) {
		double rdn = rbeta*S[ii][0] - ralpha*T[ii][0];
		rho[ii][0] = rau0*rdn + rau0;
	}
	}	
}

void rhoJMD94(field::Field<FLDARRAY> &rho, field::Field<FLDARRAY> &T, 
	field::Field<FLDARRAY> &S, v3::V3v &xyz, bool useDepth) {
/*
	Jackett and McDougall (1994) equation of state.
*/
	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<rho.get_n(); ii+=OMP_NUM_THREADS) {
		// potential temperature and salinity
		double zt = T[ii][0];
		double zs = S[ii][0];

		// depth
		double zh = (useDepth) ? -xyz[ii][2] : 0.; // since depth are negative in z

		// square root salinity
		double zsr = sqrt( fabs( zs ) );

		// compute volumic mass pure water at atm pressure
		double zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt-9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594;

		// seawater volumic mass atm pressure
		double zr2 = ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt -4.0899e-3 ) *zt+0.824493;
		double zr3 = ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3;
		double zr4 = 4.8314e-4;

		// potential volumic mass (reference to the surface)
		double zrhop = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1;

		// add the compression terms
		double ze  = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6;
		double zbw = (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4;
		double zb  = zbw + ze * zs;

		double zd  = -2.042967e-2;
		double zc  =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896;
		double zaw = ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt -4.721788;
		double za  = ( zd*zsr + zc ) *zs + zaw;

		double zb1 =   (-0.1909078*zt+7.390729 ) *zt-55.87545;
		double za1 = ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077;
		double zkw = ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt +2098.925 ) *zt+190925.6;
		double zk0 = ( zb1*zsr + za1 )*zs + zkw;

		// in situ density
		rho[ii][0] = zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  );
	}
	}	
}

void rhoJAOT12(field::Field<FLDARRAY> &rho, field::Field<FLDARRAY> &T, 
	field::Field<FLDARRAY> &S, v3::V3v &xyz, bool useDepth) {
/*
	Density of Sea Water using the Jackett and McDougall 1995 (JAOT 12) polynomial
*/	
	double eosJMDCFw[] = { 999.842594,    6.793952e-02, -9.095290e-03,  1.001685e-04, -1.120083e-06, 6.536332e-09 }; // density of fresh water at p = 0
	double eosJMDCSw[] = { 8.244930e-01, -4.089900e-03,  7.643800e-05, -8.246700e-07,  5.387500e-09, -5.724660e-03, 1.022700e-04, -1.654600e-06, 4.831400e-04 }; // density of sea water at p = 0

	double eosJMDCKFw[] = { 1.965933e+04,  1.444304e+02, -1.706103e+00,  9.648704e-03, -4.190253e-05 }; // secant bulk modulus K of fresh water at p = 0
	double eosJMDCKSw[] = { 5.284855e+01, -3.101089e-01,  6.283263e-03, -5.084188e-05,  3.886640e-01,  9.085835e-03, -4.619924e-04 }; // secant bulk modulus K of sea water at p = 0
	double eosJMDCKP[]  = { 3.186519e+00,  2.212276e-02, -2.984642e-04,  1.956415e-06,  6.704388e-03, -1.847318e-04,  2.059331e-07, 1.480266e-04, 2.102898e-04, -1.202016e-05, 1.394680e-07, -2.040237e-06, 6.128773e-08, 6.207323e-10 }; // secant bulk modulus K of sea water at p

	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<rho.get_n(); ii+=OMP_NUM_THREADS) {
		double s3o2 = S[ii][0]*sqrt(S[ii][0]);

		// density of freshwater at the surface
		rho[ii][0] = eosJMDCFw[0] 
				   + eosJMDCFw[1]*T[ii][0] 
				   + eosJMDCFw[2]*POW2(T[ii][0]) 
				   + eosJMDCFw[3]*POW3(T[ii][0])
				   + eosJMDCFw[4]*POW4(T[ii][0])
				   + eosJMDCFw[5]*POW5(T[ii][0]);

		// density of sea water at the surface
		rho[ii][0] += S[ii][0]*( eosJMDCSw[0] 
							   + eosJMDCSw[1]*T[ii][0] 
							   + eosJMDCSw[2]*POW2(T[ii][0])
							   + eosJMDCSw[3]*POW3(T[ii][0])
							   + eosJMDCSw[4]*POW4(T[ii][0]) )
					    + s3o2*( eosJMDCSw[5]
					    	   + eosJMDCSw[6]*T[ii][0] 
					    	   + eosJMDCSw[7]*POW2(T[ii][0]) )
					    + eosJMDCSw[8]*POW2(S[ii][0]);
		
		// pressure
		double P = (useDepth) ? Z2P(xyz[ii][2]) : 0.; // bar

		// secant bulk modulus of fresh water at the surface
		double bulkmod = eosJMDCKFw[0]
					   + eosJMDCKFw[1]*T[ii][0] 
					   + eosJMDCKFw[2]*POW2(T[ii][0])
					   + eosJMDCKFw[3]*POW3(T[ii][0])
					   + eosJMDCKFw[4]*POW4(T[ii][0]);
		
		// secant bulk modulus of sea water at the surface
		bulkmod += S[ii][0]*( eosJMDCKSw[0]
							+ eosJMDCKSw[1]*T[ii][0] 
							+ eosJMDCKSw[2]*POW2(T[ii][0])
							+ eosJMDCKSw[3]*POW3(T[ii][0]) )
					 + s3o2*( eosJMDCKSw[4]
					 		+ eosJMDCKSw[5]*T[ii][0] 
					 		+ eosJMDCKSw[6]*POW2(T[ii][0]) );

		// secant bulk modulus of sea water at pressure p
		bulkmod += P*( eosJMDCKP[0]
					 + eosJMDCKP[1]*T[ii][0] 
					 + eosJMDCKP[2]*POW2(T[ii][0])
					 + eosJMDCKP[3]*POW3(T[ii][0]) )
		+ P*S[ii][0]*( eosJMDCKP[4]
					 + eosJMDCKP[5]*T[ii][0] 
					 + eosJMDCKP[6]*POW2(T[ii][0]) )
		+ P*s3o2*eosJMDCKP[7]
		+ POW2(P)*( eosJMDCKP[8]
				  + eosJMDCKP[9]*T[ii][0]
				  + eosJMDCKP[10]*POW2(T[ii][0]) )
		+ POW2(P)*S[ii][0]*( eosJMDCKP[11]
						   + eosJMDCKP[12]*T[ii][0] 
						   + eosJMDCKP[13]*POW2(T[ii][0]) );

		// final density
		rho[ii][0] /= 1. - P/bulkmod;
	}
	}
}

void rhoJAOT20(field::Field<FLDARRAY> &rho, field::Field<FLDARRAY> &T, 
	field::Field<FLDARRAY> &S, v3::V3v &xyz, bool useDepth) {
/*
	Density of Sea Water using McDougall et al. 2003 (JAOT 20) polynomial
*/	
	// coefficients nonlinear equation of state in pressure coordinates for
	double eosMDJWFnum[] = { 7.35212840e+00, -5.45928211e-02, 3.98476704e-04, 2.96938239e+00, -7.23268813e-03, 2.12382341e-03, 1.04004591e-02, 1.03970529e-07, 5.18761880e-06, -3.24041825e-08, -1.23869360e-11, 9.99843699e+02 };
	double eosMDJWFden[] = { 7.28606739e-03, -4.60835542e-05, 3.68390573e-07, 1.80809186e-10, 2.14691708e-03, -9.27062484e-06, -1.78343643e-10, 4.76534122e-06, 1.63410736e-09, 5.30848875e-06, -3.03175128e-16, -1.27934137e-17, 1.00000000e+00 };
	double epsi = 0.;

	// Loop the mesh
	#pragma omp parallel
	{
	for (int ii=OMP_THREAD_NUM; ii<rho.get_n(); ii+=OMP_NUM_THREADS) {
		// pressure
		double P = (useDepth) ? 10.*Z2P(xyz[ii][2]) : 0.; // dbar

		double sp5  = sqrt(S[ii][0]);
		double p1t1 = P*T[ii][0];

		double num = eosMDJWFnum[11] + T[ii][0]*(eosMDJWFnum[0] + T[ii][0]*(eosMDJWFnum[1] + eosMDJWFnum[2]*T[ii][0]))
				   + S[ii][0]*(eosMDJWFnum[3] + eosMDJWFnum[4]*T[ii][0]  + eosMDJWFnum[5]*S[ii][0])
				   + P*(eosMDJWFnum[6] + eosMDJWFnum[7]*POW2(T[ii][0]) + eosMDJWFnum[8]*S[ii][0] + P*(eosMDJWFnum[9] + eosMDJWFnum[10]*POW2(T[ii][0])));

		double den =  eosMDJWFden[12] + T[ii][0]*(eosMDJWFden[0] + T[ii][0]*(eosMDJWFden[1] + T[ii][0]*(eosMDJWFden[2] + T[ii][0]*eosMDJWFden[3])))
				   + S[ii][0]*(eosMDJWFden[4] + T[ii][0]*(eosMDJWFden[5] + eosMDJWFden[6]*POW2(T[ii][0])) + sp5*(eosMDJWFden[7] + eosMDJWFden[8]*POW2(T[ii][0])))
				   + P*(eosMDJWFden[9] + p1t1*(eosMDJWFden[10]*POW2(T[ii][0]) + eosMDJWFden[11]*P));

		// final density
		rho[ii][0] = num/(epsi+den);
	}
	}
}

//----------------------------------------------------------------------------
vtkOGSDensity::vtkOGSDensity() {
	this->Tarrname = NULL;
	this->Sarrname = NULL;
	this->method   = 0;
	this->useDepth = true;
	this->rau0     = 1020.;
	this->ralpha   = 2.e-4;
	this->rbeta    = 0.001;
	this->nProcs   = 0;
	this->procId   = 0;

	#ifdef PARAVIEW_USE_MPI
		this->Controller = NULL;
		this->SetController(vtkMultiProcessController::GetGlobalController());
	#endif
}

//----------------------------------------------------------------------------
vtkOGSDensity::~vtkOGSDensity() {
	this->SetTarrname(0);
	this->SetSarrname(0);

	#ifdef PARAVIEW_USE_MPI
		this->SetController(NULL);	
	#endif
}

//----------------------------------------------------------------------------
int vtkOGSDensity::RequestInformation(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector) {
	
	/* SET UP THE PARALLEL CONTROLLER

		The MPI threads come initialized by the ParaView server. Here
		we set up the environment for this filter.

	*/
	#ifdef PARAVIEW_USE_MPI
	if (this->Controller->GetNumberOfProcesses() > 1) {
		this->nProcs = this->Controller->GetNumberOfProcesses();
		this->procId = this->Controller->GetLocalProcessId();
	}

	// Stop all threads except from the master to execute
	if (this->procId > 0) return 1;
	#endif

	return 1;
}

//----------------------------------------------------------------------------
int vtkOGSDensity::RequestData(vtkInformation *vtkNotUsed(request), 
	vtkInformationVector **inputVector, vtkInformationVector *outputVector) {

	// get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Stop all threads except from the master to execute
	#ifdef PARAVIEW_USE_MPI
	if (this->procId > 0) return 1;
	#endif

	// get the input and output
	vtkDataSet *input = vtkDataSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *output = vtkDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	output->ShallowCopy(input);
	this->UpdateProgress(0.);

	// Decide whether we have cell or point data
	int n_cell_vars  = input->GetCellData()->GetNumberOfArrays();
	int n_point_vars = input->GetPointData()->GetNumberOfArrays();

	bool iscelld = (n_cell_vars > n_point_vars) ? true : false;

	// Load temperature and salinity arrays
	VTKARRAY *vtkT, *vtkS;
	if (iscelld) {
		vtkT = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(this->Tarrname) );
		vtkS = VTKARRAY::SafeDownCast( input->GetCellData()->GetArray(this->Sarrname) );
	} else {
		vtkT = VTKARRAY::SafeDownCast( input->GetPointData()->GetArray(this->Tarrname) );
		vtkS = VTKARRAY::SafeDownCast( input->GetPointData()->GetArray(this->Sarrname) );
	}

	if (!vtkT) {
		vtkErrorMacro("Error loading <"<<this->Tarrname<<">! Aborting!");
		return 0;
	}

	if (!vtkS) {
		vtkErrorMacro("Error loading <"<<this->Sarrname<<">! Aborting!");
		return 0;
	}

	field::Field<FLDARRAY> T = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkT);
	field::Field<FLDARRAY> S = VTK::createFieldfromVTK<VTKARRAY,FLDARRAY>(vtkS);
	
	field::Field<FLDARRAY> rho(T.get_n(),1,0.);

	this->UpdateProgress(0.25);

	// Recover Metadata array (depth factor)
	vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
		input->GetFieldData()->GetAbstractArray("Metadata"));
	double dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;
	
	if (vtkmetadata == NULL) 
		vtkWarningMacro("Field array Metadata not found! Depth factor set to 1000. automatically.");

	// Select, according to the user, which density model to use
	switch(this->method) {
		case 0: // Linear depending on temperature
			rhoLinT(rho,T,this->ralpha,this->rau0); // Not used!
			break;
		case 1: // Linear depending on temperature and salinity
			rhoLinTS(rho,T,S,this->ralpha,this->rbeta,this->rau0); // Not used!
			break;
		case 2: // JMD 94
			{
				v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,dfact) : VTK::getVTKCellPoints(input,dfact);
				rhoJMD94(rho,T,S,xyz,this->useDepth);
			}
			break;
		case 3: // JAOT 12
			{
				v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,dfact) : VTK::getVTKCellPoints(input,dfact);
				rhoJAOT12(rho,T,S,xyz,this->useDepth);
			}
			break;
		case 4: // JAOT 20
			{
				v3::V3v xyz = (iscelld) ? VTK::getVTKCellCenters(input,dfact) : VTK::getVTKCellPoints(input,dfact);
				rhoJAOT20(rho,T,S,xyz,this->useDepth);
			}
			break;
		case 5: // TEOS10
			{
				vtkWarningMacro("TEOS10 not implemented yet!");
			}
			break;
		default:
			vtkErrorMacro("Unimplemented model <"<<this->method<<">! Aborting!");
			return 0;
	}
	this->UpdateProgress(0.75);

	// Output field
	VTKARRAY *vtkrho;
	vtkrho = VTK::createVTKfromField<VTKARRAY,FLDARRAY>("Density",rho);

	if (iscelld)
		output->GetCellData()->AddArray(vtkrho);
	else
		output->GetPointData()->AddArray(vtkrho);

	vtkrho->Delete();

	// Return
	this->UpdateProgress(1.0);
	return 1;
}