/*=========================================================================

  Program:   OGS Gradient
  Module:    vtkOGSGradient.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/*----------------------------------------------------------------------------
 Copyright (c) Sandia Corporation
 See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.
----------------------------------------------------------------------------*/

#include "vtkOGSGradient.h"

#include "vtkCell.h"
#include "vtkGenericCell.h"
#include "vtkCellData.h"
#include "vtkCellDataToPointData.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkIdList.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStringArray.h"

#include <limits>
#include <vector>
#include <omp.h>
int omp_get_num_threads();
int omp_get_thread_num();

#include "../_utils/matrixMN.h"

//-----------------------------------------------------------------------------

vtkStandardNewMacro(vtkOGSGradient);

namespace
{
// special template macro for only float and double types since we will never
// have other data types for output arrays and this helps with reducing template expansion
#define vtkFloatingPointTemplateMacro(call)         \
  vtkTemplateMacroCase(VTK_DOUBLE, double, call);   \
  vtkTemplateMacroCase(VTK_FLOAT, float, call);

// helper function to replace the gradient of a vector
// with the vorticity/curl of that vector
//-----------------------------------------------------------------------------
  template<class data_type>
  void ComputeVorticityFromGradient(data_type* gradients, data_type* vorticity)
  {
    vorticity[0] = gradients[7] - gradients[5];
    vorticity[1] = gradients[2] - gradients[6];
    vorticity[2] = gradients[3] - gradients[1];
  }

  template<class data_type>
  void ComputeDivergenceFromGradient(data_type* gradients, data_type* divergence)
  {
    divergence[0] = gradients[0]+gradients[4]+gradients[8];
  }

  template<class data_type>
  void ComputeQCriterionFromGradient(data_type* gradients, data_type* qCriterion)
  {
    // see http://public.kitware.com/pipermail/paraview/2015-May/034233.html for
    // paper citation and formula on Q-criterion.
    qCriterion[0] =
      - (gradients[0]*gradients[0]+gradients[4]*gradients[4]+gradients[8]*gradients[8])/2.
      - (gradients[1]*gradients[3]+gradients[2]*gradients[6]+gradients[5]*gradients[7]);
  }

  template<class data_type>
  void ComputeLambda2CriterionFromGradient(data_type* gradients, data_type* l2Criterion, int n)
  {
    // Compute S, O and R matrices
    matMN::matrixMN<data_type> grad(3,gradients), V(3,(data_type)(0.));
    matMN::matrixMN<data_type> A = 0.5*(grad + grad.t()); //A /= std::sqrt(A.norm2());
    matMN::matrixMN<data_type> B = 0.5*(grad - grad.t()); //B /= std::sqrt(B.norm2());
    matMN::matrixMN<data_type> C = (A^A) + (B^B);         // Dot product

    // Compute the eigenvalues and eigenvectors of R
    data_type e[3] = {0.,0.,0.};
    int n_iter = matMN::eigen(C,e,V,n);

    // Store Lambda2
    l2Criterion[0] = e[1];
  }

  template<class data_type>
  void ComputeOmegaCriterionFromGradient(data_type* gradients, data_type* omegaCriterion, data_type* aux, data_type epsi)
  {
    // Compute A and B matrices
    matMN::matrixMN<data_type> grad(3,gradients);
    matMN::matrixMN<data_type> A = 0.5*(grad + grad.t());
    matMN::matrixMN<data_type> B = 0.5*(grad - grad.t());

    // Compute a and b and store them in aux and Omega
    aux[0]            = A.norm2(); // a (Frobenius norm squared)
    omegaCriterion[0] = B.norm2(); // b (Frobenius norm squared)

    // Store the maximum of (b - a)
    data_type ba_aux = omegaCriterion[0] - aux[0];
    epsi = (ba_aux > epsi) ? ba_aux : epsi;
  }

  // Functions for unstructured grids and polydatas
  template<class data_type>
  void ComputePointGradientsUG(
    vtkDataSet *structure, vtkDataArray *array, data_type *gradients,
    int numberOfInputComponents, data_type* vorticity, data_type* qCriterion, data_type* L2Criterion,
    data_type* omegaCriterion, data_type* divergence, double eps0, double dfact, int highestCellDimension, int contributingCellOption);

  int GetCellParametricData(
    vtkIdType pointId, double pointCoord[3], vtkCell *cell, int & subId,
    double parametricCoord[3]);

  template<class data_type>
  void ComputeCellGradientsUG(
    vtkDataSet *structure, vtkDataArray *array, data_type *gradients,
    int numberOfInputComponents, data_type* vorticity, data_type* qCriterion, 
    data_type* L2Criterion, data_type* omegaCriterion, data_type* divergence, double eps0, double dfact);

  // Functions for image data and structured grids
  template<class Grid, class data_type>
  void ComputeGradientsSG(Grid output, vtkDataArray* array, data_type* gradients,
                          int numberOfInputComponents, int fieldAssociation,
                          data_type* vorticity, data_type* qCriterion, data_type* L2Criterion,
                          data_type* omegaCriterion, data_type* divergence, double eps0, double dfact);

  bool vtkGradientFilterHasArray(vtkFieldData *fieldData,
                                 vtkDataArray *array)
  {
    int numarrays = fieldData->GetNumberOfArrays();
    for (int i = 0; i < numarrays; i++)
    {
      if (array == fieldData->GetArray(i))
      {
        return true;
      }
    }
    return false;
  }

  // generic way to get the coordinate for either a cell (using
  // the parametric center) or a point
  void GetGridEntityCoordinate(vtkDataSet* grid, int fieldAssociation, vtkGenericCell *cell,
                               vtkIdType index, double coords[3])
  {
    if(fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_POINTS)
    {
      grid->GetPoint(index, coords);
    }
    else
    {
      grid->GetCell(index,cell);
      double pcoords[3];
      int subId = cell->GetParametricCenter(pcoords);
      std::vector<double> weights(cell->GetNumberOfPoints()+1);
      cell->EvaluateLocation(subId, pcoords, coords, &weights[0]);
    }
  }

  template<class data_type>
  void Fill(vtkDataArray* array, data_type vtkNotUsed(data), int replacementValueOption)
  {
    switch (replacementValueOption)
    {
    case vtkOGSGradient::Zero:
      array->Fill(0);
      return;
    case vtkOGSGradient::NaN:
      array->Fill(vtkMath::Nan());
      return;
    case vtkOGSGradient::DataTypeMin:
      array->Fill(std::numeric_limits<data_type>::min());
      return;
    case vtkOGSGradient::DataTypeMax:
      array->Fill(std::numeric_limits<data_type>::max());
      return;
    }
  }
  template<class data_type>
  int GetOutputDataType(data_type vtkNotUsed(data))
  {
    if(sizeof(data_type)>=8)
    {
      return VTK_DOUBLE;
    }
    return VTK_FLOAT;
  }
} // end anonymous namespace

//-----------------------------------------------------------------------------
vtkOGSGradient::vtkOGSGradient()
{
  this->ResultArrayName = nullptr;
  this->DivergenceArrayName = nullptr;
  this->VorticityArrayName = nullptr;
  this->QCriterionArrayName = nullptr;
  this->Lambda2CriterionArrayName = nullptr;
  this->OmegaCriterionArrayName = nullptr;
  this->FasterApproximation = 0;
  this->ComputeGradient = 1;
  this->ComputeDivergence = 0;
  this->ComputeVorticity = 0;
  this->ComputeQCriterion = 0;
  this->ComputeLambda2Criterion = 0;
  this->ComputeOmegaCriterion = 0;
  this->epsi = 1.e-3;
  this->ContributingCellOption = vtkOGSGradient::All;
  this->ReplacementValueOption = vtkOGSGradient::Zero;
  this->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
                        vtkDataSetAttributes::SCALARS);
}

//-----------------------------------------------------------------------------
vtkOGSGradient::~vtkOGSGradient()
{
  this->SetResultArrayName(nullptr);
  this->SetDivergenceArrayName(nullptr);
  this->SetVorticityArrayName(nullptr);
  this->SetQCriterionArrayName(nullptr);
  this->SetLambda2CriterionArrayName(nullptr);
  this->SetOmegaCriterionArrayName(nullptr);
}

//-----------------------------------------------------------------------------
void vtkOGSGradient::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "ResultArrayName:"
     << (this->ResultArrayName ? this->ResultArrayName : "Gradients") << endl;
  os << indent << "DivergenceArrayName:"
     << (this->DivergenceArrayName ? this->DivergenceArrayName : "Divergence") << endl;
  os << indent << "VorticityArrayName:"
     << (this->VorticityArrayName ? this->VorticityArrayName : "Vorticity") << endl;
  os << indent << "QCriterionArrayName:"
     << (this->QCriterionArrayName ? this->QCriterionArrayName : "Q-criterion") << endl;
  os << indent << "FasterApproximation:" << this->FasterApproximation << endl;
  os << indent << "ComputeGradient:"  << this->ComputeGradient << endl;
  os << indent << "ComputeDivergence:"  << this->ComputeDivergence << endl;
  os << indent << "ComputeVorticity:" << this->ComputeVorticity << endl;
  os << indent << "ComputeQCriterion:" << this->ComputeQCriterion << endl;
  os << indent << "ComputeLambda2Criterion:" << this->ComputeLambda2Criterion << endl;
  os << indent << "ComputeOmegaCriterion:" << this->ComputeOmegaCriterion << endl;
  os << indent << "ContributingCellOption:" << this->ContributingCellOption << endl;
  os << indent << "ReplacementValueOption:" << this->ReplacementValueOption << endl;
}

//-----------------------------------------------------------------------------
void vtkOGSGradient::SetInputScalars(int fieldAssociation, const char *name)
{
  if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_CELLS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) )
  {
    vtkErrorMacro("Input Array must be associated with points or cells.");
    return;
  }

  this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, name);
}

//-----------------------------------------------------------------------------
void vtkOGSGradient::SetInputScalars(int fieldAssociation,
                                        int fieldAttributeType)
{
  if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_CELLS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) )
  {
    vtkErrorMacro("Input Array must be associated with points or cells.");
    return;
  }

  this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, fieldAttributeType);
}

//-----------------------------------------------------------------------------
int vtkOGSGradient::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Technically, this code is only correct for pieces extent types.  However,
  // since this class is pretty inefficient for data types that use 3D extents,
  // we'll punt on the ghost levels for them, too.
  int piece, numPieces, ghostLevels;

  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(
                   vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevels = outInfo->Get(
             vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  if (numPieces > 1)
  {
    ++ghostLevels;
  }

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
              numPieces);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
              ghostLevels);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}

//-----------------------------------------------------------------------------
int vtkOGSGradient::RequestData(vtkInformation *vtkNotUsed(request),
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector)
{
  vtkDebugMacro("RequestData");

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkDataSet *input
    = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output
    = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataArray *array = this->GetInputArrayToProcess(0, inputVector);

  if (input->GetNumberOfCells() == 0)
  {
    // need cells to compute the gradient so if we don't have cells. we can't compute anything.
    // if we have points and an array with values provide a warning letting the user know that
    // no gradient will be computed because of the lack of cells. otherwise the dataset is
    // assumed empty and we can skip providing a warning message to the user.
    if (input->GetNumberOfPoints() && array && array->GetNumberOfTuples())
    {
      vtkWarningMacro("Cannot compute gradient for datasets without cells");
    }
    output->ShallowCopy(input);
    return 1;
  }

  if (array == nullptr)
  {
    vtkErrorMacro("No input array. If this dataset is part of a composite dataset"
                  << " check to make sure that all non-empty blocks have this array.");
    return 0;
  }
  if (array->GetNumberOfComponents() == 0)
  {
    vtkErrorMacro("Input array must have at least one component.");
    return 0;
  }
    
  // Recover Metadata array (depth factor)
  // This is necessary to be consistent with the stretching of the mesh
  // imposed in the third dimension of the reader.
  vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
    input->GetFieldData()->GetAbstractArray("Metadata"));
  double dfact = (vtkmetadata != NULL) ? std::stod( vtkmetadata->GetValue(2) ) : 1000.;

  // we can only compute vorticity and Q criterion if the input
  // array has 3 components. if we can't compute them because of
  // this we only mark internally the we aren't computing them
  // since we don't want to change the state of the filter.
  bool computeVorticity = this->ComputeVorticity != 0;
  bool computeDivergence = this->ComputeDivergence != 0;
  bool computeQCriterion = this->ComputeQCriterion != 0;
  bool computeLambda2Criterion = this->ComputeLambda2Criterion != 0;
  bool computeOmegaCriterion = this->ComputeOmegaCriterion != 0;
  if( (this->ComputeOmegaCriterion || this->ComputeLambda2Criterion || this->ComputeQCriterion || this->ComputeVorticity ||
       this->ComputeDivergence) && array->GetNumberOfComponents() != 3)
  {
    vtkWarningMacro("Input array must have exactly three components with "
                    << "ComputeDivergence, ComputeVorticity, ComputeQCriterion, computeLambda2Criterion or computeOmegaCriterion flag enabled."
                    << "Skipping divergence, vorticity, Q-criterion, Lambda2-criterion and Omega-criterion computation.");
    computeVorticity = false;
    computeQCriterion = false;
    computeDivergence = false;
    computeLambda2Criterion = false;
    computeOmegaCriterion = false;
  }

  int fieldAssociation;
  if (vtkGradientFilterHasArray(input->GetPointData(), array))
  {
    fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
  }
  else if (vtkGradientFilterHasArray(input->GetCellData(), array))
  {
    fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_CELLS;
  }
  else
  {
    vtkErrorMacro("Input arrays do not seem to be either point or cell arrays.");
    return 0;
  }

  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());

  if(output->IsA("vtkImageData") || output->IsA("vtkStructuredGrid") ||
          output->IsA("vtkRectilinearGrid") )
  {
    this->ComputeRegularGridGradient(
      array, fieldAssociation, computeVorticity, computeQCriterion, 
      computeLambda2Criterion, computeOmegaCriterion, computeDivergence, dfact, output);
  }
  else
  {
    this->ComputeUnstructuredGridGradient(
      array, fieldAssociation, input, computeVorticity, computeQCriterion, 
      computeLambda2Criterion, computeOmegaCriterion, computeDivergence, dfact, output);
  }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkOGSGradient::ComputeUnstructuredGridGradient(
  vtkDataArray* array, int fieldAssociation, vtkDataSet* input,
  bool computeVorticity, bool computeQCriterion, bool computeLambda2Criterion, 
  bool computeOmegaCriterion, bool computeDivergence, double dfact, vtkDataSet* output)
{
  int arrayType = this->GetOutputArrayType(array);
  int numberOfInputComponents = array->GetNumberOfComponents();
  vtkSmartPointer<vtkDataArray> gradients = nullptr;
  if(this->ComputeGradient)
  {
    gradients.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    gradients->SetNumberOfComponents(3*numberOfInputComponents);
    gradients->SetNumberOfTuples(array->GetNumberOfTuples());
    switch (arrayType)
    {
      vtkFloatingPointTemplateMacro(Fill(gradients, static_cast<VTK_TT>(0), this->ReplacementValueOption));
    }
    if (this->ResultArrayName)
    {
      gradients->SetName(this->ResultArrayName);
    }
    else
    {
      gradients->SetName("Gradients");
    }
  }
  vtkSmartPointer<vtkDataArray> divergence = nullptr;
  if(computeDivergence)
  {
    divergence.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    divergence->SetNumberOfTuples(array->GetNumberOfTuples());
    switch (arrayType)
    {
      vtkFloatingPointTemplateMacro(Fill(divergence, static_cast<VTK_TT>(0), this->ReplacementValueOption));
    }
    if (this->DivergenceArrayName)
    {
      divergence->SetName(this->DivergenceArrayName);
    }
    else
    {
      divergence->SetName("Divergence");
    }
  }
  vtkSmartPointer<vtkDataArray> vorticity;
  if(computeVorticity)
  {
    vorticity.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    vorticity->SetNumberOfComponents(3);
    vorticity->SetNumberOfTuples(array->GetNumberOfTuples());
    switch (arrayType)
    {
      vtkFloatingPointTemplateMacro(Fill(vorticity, static_cast<VTK_TT>(0), this->ReplacementValueOption));
    }
    if (this->VorticityArrayName)
    {
      vorticity->SetName(this->VorticityArrayName);
    }
    else
    {
      vorticity->SetName("Vorticity");
    }
  }
  vtkSmartPointer<vtkDataArray> qCriterion;
  if(computeQCriterion)
  {
    qCriterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    qCriterion->SetNumberOfTuples(array->GetNumberOfTuples());
    switch (arrayType)
    {
      vtkFloatingPointTemplateMacro(Fill(qCriterion, static_cast<VTK_TT>(0), this->ReplacementValueOption));
    }
    if (this->QCriterionArrayName)
    {
      qCriterion->SetName(this->QCriterionArrayName);
    }
    else
    {
      qCriterion->SetName("Q-criterion");
    }
  }
  vtkSmartPointer<vtkDataArray> L2Criterion;
  if(computeLambda2Criterion)
  {
    L2Criterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    L2Criterion->SetNumberOfTuples(array->GetNumberOfTuples());
    switch (arrayType)
    {
      vtkFloatingPointTemplateMacro(Fill(L2Criterion, static_cast<VTK_TT>(0), this->ReplacementValueOption));
    }
    if (this->Lambda2CriterionArrayName)
    {
      L2Criterion->SetName(this->Lambda2CriterionArrayName);
    }
    else
    {
      L2Criterion->SetName("L2-criterion");
    }
  }
  vtkSmartPointer<vtkDataArray> OCriterion;
  if(computeOmegaCriterion)
  {
    OCriterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    OCriterion->SetNumberOfTuples(array->GetNumberOfTuples());
    switch (arrayType)
    {
      vtkFloatingPointTemplateMacro(Fill(OCriterion, static_cast<VTK_TT>(0), this->ReplacementValueOption));
    }
    if (this->OmegaCriterionArrayName)
    {
      OCriterion->SetName(this->OmegaCriterionArrayName);
    }
    else
    {
      OCriterion->SetName("Omega-criterion");
    }
  }

  int highestCellDimension = 0;
  if (this->ContributingCellOption == vtkOGSGradient::DataSetMax)
  {
    int maxDimension = input->IsA("vtkPolyData") == 1 ? 2 : 3;
    for (vtkIdType i=0;i<input->GetNumberOfCells();i++)
    {
      int dim = input->GetCell(i)->GetCellDimension();
      if (dim > highestCellDimension)
      {
        highestCellDimension = dim;
        if (highestCellDimension == maxDimension)
        {
          break;
        }
      }
    }
  }

  if (fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_POINTS)
  {
    if (!this->FasterApproximation)
    {
      switch (arrayType)
      { // ok to use template macro here since we made the output arrays ourselves
        vtkFloatingPointTemplateMacro(ComputePointGradientsUG(
                           input, array,
                           (gradients == nullptr ? nullptr :
                            static_cast<VTK_TT *>(gradients->GetVoidPointer(0))),
                           numberOfInputComponents,
                           (vorticity == nullptr ? nullptr :
                            static_cast<VTK_TT *>(vorticity->GetVoidPointer(0))),
                           (qCriterion == nullptr ? nullptr :
                            static_cast<VTK_TT *>(qCriterion->GetVoidPointer(0))),
                           (L2Criterion == nullptr ? nullptr :
                            static_cast<VTK_TT *>(L2Criterion->GetVoidPointer(0))),
                           (OCriterion == nullptr ? nullptr :
                            static_cast<VTK_TT *>(OCriterion->GetVoidPointer(0))),
                           (divergence == nullptr ? nullptr :
                            static_cast<VTK_TT *>(divergence->GetVoidPointer(0))),
                           this->epsi, dfact, highestCellDimension, this->ContributingCellOption));
      }
      if(gradients)
      {
        output->GetPointData()->AddArray(gradients);
      }
      if(divergence)
      {
        output->GetPointData()->AddArray(divergence);
      }
      if(vorticity)
      {
        output->GetPointData()->AddArray(vorticity);
      }
      if(qCriterion)
      {
        output->GetPointData()->AddArray(qCriterion);
      }
      if(L2Criterion)
      {
        output->GetPointData()->AddArray(L2Criterion);
      }
      if(OCriterion)
      {
        output->GetPointData()->AddArray(OCriterion);
      }
    }
    else // this->FasterApproximation
    {
      // The cell computation is faster and works off of point data anyway.  The
      // faster approximation is to use the cell algorithm and then convert the
      // result to point data.
      vtkSmartPointer<vtkDataArray> cellGradients = nullptr;
      if(gradients)
      {
        cellGradients.TakeReference(vtkDataArray::CreateDataArray(arrayType));
        cellGradients->SetName(gradients->GetName());
        cellGradients->SetNumberOfComponents(3*array->GetNumberOfComponents());
        cellGradients->SetNumberOfTuples(input->GetNumberOfCells());
      }
      vtkSmartPointer<vtkDataArray> cellDivergence = nullptr;
      if(divergence)
      {
        cellDivergence.TakeReference(vtkDataArray::CreateDataArray(arrayType));
        cellDivergence->SetName(divergence->GetName());
        cellDivergence->SetNumberOfTuples(input->GetNumberOfCells());
      }
      vtkSmartPointer<vtkDataArray> cellVorticity = nullptr;
      if(vorticity)
      {
        cellVorticity.TakeReference(vtkDataArray::CreateDataArray(arrayType));
        cellVorticity->SetName(vorticity->GetName());
        cellVorticity->SetNumberOfComponents(3);
        cellVorticity->SetNumberOfTuples(input->GetNumberOfCells());
      }
      vtkSmartPointer<vtkDataArray> cellQCriterion = nullptr;
      if(qCriterion)
      {
        cellQCriterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
        cellQCriterion->SetName(qCriterion->GetName());
        cellQCriterion->SetNumberOfTuples(input->GetNumberOfCells());
      }
      vtkSmartPointer<vtkDataArray> cellL2Criterion = nullptr;
      if(L2Criterion)
      {
        cellL2Criterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
        cellL2Criterion->SetName(L2Criterion->GetName());
        cellL2Criterion->SetNumberOfTuples(input->GetNumberOfCells());
      }
      vtkSmartPointer<vtkDataArray> cellOCriterion = nullptr;
      if(L2Criterion)
      {
        cellOCriterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
        cellOCriterion->SetName(OCriterion->GetName());
        cellOCriterion->SetNumberOfTuples(input->GetNumberOfCells());
      }

      switch (arrayType)
      { // ok to use template macro here since we made the output arrays ourselves
        vtkFloatingPointTemplateMacro(
          ComputeCellGradientsUG(
            input, array,
            (cellGradients == nullptr ? nullptr :
             static_cast<VTK_TT *>(cellGradients->GetVoidPointer(0))),
            numberOfInputComponents,
            (vorticity == nullptr ? nullptr :
             static_cast<VTK_TT *>(cellVorticity->GetVoidPointer(0))),
            (qCriterion == nullptr ? nullptr :
             static_cast<VTK_TT *>(cellQCriterion->GetVoidPointer(0))),
            (L2Criterion == nullptr ? nullptr :
             static_cast<VTK_TT *>(cellL2Criterion->GetVoidPointer(0))),
            (OCriterion == nullptr ? nullptr :
             static_cast<VTK_TT *>(cellOCriterion->GetVoidPointer(0))),
            (divergence == nullptr ? nullptr :
             static_cast<VTK_TT *>(cellDivergence->GetVoidPointer(0))),
            this->epsi, dfact));
      }

      // We need to convert cell Array to points Array.
      vtkSmartPointer<vtkDataSet> dummy;
      dummy.TakeReference(input->NewInstance());
      dummy->CopyStructure(input);
      if(cellGradients)
      {
        dummy->GetCellData()->AddArray(cellGradients);
      }
      if(divergence)
      {
        dummy->GetCellData()->AddArray(cellDivergence);
      }
      if(vorticity)
      {
        dummy->GetCellData()->AddArray(cellVorticity);
      }
      if(qCriterion)
      {
        dummy->GetCellData()->AddArray(cellQCriterion);
      }
      if(L2Criterion)
      {
        dummy->GetCellData()->AddArray(cellL2Criterion);
      }
      if(OCriterion)
      {
        dummy->GetCellData()->AddArray(cellOCriterion);
      }

      vtkNew<vtkCellDataToPointData> cd2pd;
      cd2pd->SetInputData(dummy);
      cd2pd->PassCellDataOff();
      cd2pd->SetContributingCellOption(this->ContributingCellOption);
      cd2pd->Update();

      // Set the gradients array in the output and cleanup.
      if(gradients)
      {
        output->GetPointData()->AddArray(
          cd2pd->GetOutput()->GetPointData()->GetArray(gradients->GetName()));
      }
      if(qCriterion)
      {
        output->GetPointData()->AddArray(
          cd2pd->GetOutput()->GetPointData()->GetArray(qCriterion->GetName()));
      }
      if(L2Criterion)
      {
        output->GetPointData()->AddArray(
          cd2pd->GetOutput()->GetPointData()->GetArray(L2Criterion->GetName()));
      }
      if(OCriterion)
      {
        output->GetPointData()->AddArray(
          cd2pd->GetOutput()->GetPointData()->GetArray(OCriterion->GetName()));
      }
      if(divergence)
      {
        output->GetPointData()->AddArray(
          cd2pd->GetOutput()->GetPointData()->GetArray(divergence->GetName()));
      }
      if(vorticity)
      {
        output->GetPointData()->AddArray(
          cd2pd->GetOutput()->GetPointData()->GetArray(vorticity->GetName()));
      }
    }
  }
  else  // fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_CELLS
  {
    // We need to convert cell Array to points Array.
    vtkSmartPointer<vtkDataSet> dummy;
    dummy.TakeReference(input->NewInstance());
    dummy->CopyStructure(input);
    dummy->GetCellData()->SetScalars(array);

    vtkNew<vtkCellDataToPointData> cd2pd;
    cd2pd->SetInputData(dummy);
    cd2pd->PassCellDataOff();
    cd2pd->SetContributingCellOption(this->ContributingCellOption);
    cd2pd->Update();
    vtkDataArray *pointScalars
      = cd2pd->GetOutput()->GetPointData()->GetScalars();
    pointScalars->Register(this);

    switch (arrayType)
    { // ok to use template macro here since we made the output arrays ourselves
      vtkFloatingPointTemplateMacro(ComputeCellGradientsUG(
                         input, pointScalars,
                         (gradients == nullptr ? nullptr :
                          static_cast<VTK_TT *>(gradients->GetVoidPointer(0))),
                         numberOfInputComponents,
                         (vorticity == nullptr ? nullptr :
                          static_cast<VTK_TT *>(vorticity->GetVoidPointer(0))),
                         (qCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(qCriterion->GetVoidPointer(0))),
                         (L2Criterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(L2Criterion->GetVoidPointer(0))),
                         (OCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(OCriterion->GetVoidPointer(0))),
                         (divergence == nullptr ? nullptr :
                          static_cast<VTK_TT *>(divergence->GetVoidPointer(0))),
                         this->epsi, dfact));
    }

    if(gradients)
    {
      output->GetCellData()->AddArray(gradients);
    }
    if(vorticity)
    {
      output->GetCellData()->AddArray(vorticity);
    }
    if(divergence)
    {
      output->GetCellData()->AddArray(divergence);
    }
    if(qCriterion)
    {
      output->GetCellData()->AddArray(qCriterion);
    }
    if(L2Criterion)
    {
      output->GetCellData()->AddArray(L2Criterion);
    }
    if(OCriterion)
    {
      output->GetCellData()->AddArray(OCriterion);
    }
    pointScalars->UnRegister(this);
  }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkOGSGradient::ComputeRegularGridGradient(
  vtkDataArray* array, int fieldAssociation, bool computeVorticity,
  bool computeQCriterion, bool computeLambda2Criterion, bool computeOmegaCriterion, 
  bool computeDivergence, double dfact, vtkDataSet* output)
{
  int arrayType = this->GetOutputArrayType(array);
  int numberOfInputComponents = array->GetNumberOfComponents();
  vtkSmartPointer<vtkDataArray> gradients = nullptr;
  if(this->ComputeGradient)
  {
    gradients.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    gradients->SetNumberOfComponents(3*numberOfInputComponents);
    gradients->SetNumberOfTuples(array->GetNumberOfTuples());
    if (this->ResultArrayName)
    {
      gradients->SetName(this->ResultArrayName);
    }
    else
    {
      gradients->SetName("Gradients");
    }
  }
  vtkSmartPointer<vtkDataArray> divergence = nullptr;
  if(computeDivergence)
  {
    divergence.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    divergence->SetNumberOfTuples(array->GetNumberOfTuples());
    if (this->DivergenceArrayName)
    {
      divergence->SetName(this->DivergenceArrayName);
    }
    else
    {
      divergence->SetName("Divergence");
    }
  }
  vtkSmartPointer<vtkDataArray> vorticity;
  if(computeVorticity)
  {
    vorticity.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    vorticity->SetNumberOfComponents(3);
    vorticity->SetNumberOfTuples(array->GetNumberOfTuples());
    if (this->VorticityArrayName)
    {
      vorticity->SetName(this->VorticityArrayName);
    }
    else
    {
      vorticity->SetName("Vorticity");
    }
  }
  vtkSmartPointer<vtkDataArray> qCriterion;
  if(computeQCriterion)
  {
    qCriterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    qCriterion->SetNumberOfTuples(array->GetNumberOfTuples());
    if (this->QCriterionArrayName)
    {
      qCriterion->SetName(this->QCriterionArrayName);
    }
    else
    {
      qCriterion->SetName("Q-criterion");
    }
  }
  vtkSmartPointer<vtkDataArray> L2Criterion;
  if(computeLambda2Criterion)
  {
    L2Criterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    L2Criterion->SetNumberOfTuples(array->GetNumberOfTuples());
    if (this->Lambda2CriterionArrayName)
    {
      L2Criterion->SetName(this->Lambda2CriterionArrayName);
    }
    else
    {
      L2Criterion->SetName("L2-criterion");
    }
  }
  vtkSmartPointer<vtkDataArray> OCriterion;
  if(computeOmegaCriterion)
  {
    OCriterion.TakeReference(vtkDataArray::CreateDataArray(arrayType));
    OCriterion->SetNumberOfTuples(array->GetNumberOfTuples());
    if (this->OmegaCriterionArrayName)
    {
      OCriterion->SetName(this->OmegaCriterionArrayName);
    }
    else
    {
      OCriterion->SetName("Omega-criterion");
    }
  }

  if(vtkStructuredGrid* structuredGrid = vtkStructuredGrid::SafeDownCast(output))
  {
    switch (arrayType)
    { // ok to use template macro here since we made the output arrays ourselves
      vtkFloatingPointTemplateMacro(ComputeGradientsSG(
                         structuredGrid, array,
                         (gradients == nullptr ? nullptr :
                          static_cast<VTK_TT *>(gradients->GetVoidPointer(0))),
                         numberOfInputComponents, fieldAssociation,
                         (vorticity == nullptr ? nullptr :
                          static_cast<VTK_TT *>(vorticity->GetVoidPointer(0))),
                         (qCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(qCriterion->GetVoidPointer(0))),
                         (L2Criterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(L2Criterion->GetVoidPointer(0))),
                         (OCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(OCriterion->GetVoidPointer(0))),
                         (divergence == nullptr ? nullptr :
                          static_cast<VTK_TT *>(divergence->GetVoidPointer(0))),
                         this->epsi, dfact));

    }
  }
  else if(vtkImageData* imageData = vtkImageData::SafeDownCast(output))
  {
    switch (arrayType)
    { // ok to use template macro here since we made the output arrays ourselves
      vtkFloatingPointTemplateMacro(ComputeGradientsSG(
                         imageData, array,
                         (gradients == nullptr ? nullptr :
                          static_cast<VTK_TT *>(gradients->GetVoidPointer(0))),
                         numberOfInputComponents, fieldAssociation,
                         (vorticity == nullptr ? nullptr :
                          static_cast<VTK_TT *>(vorticity->GetVoidPointer(0))),
                         (qCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(qCriterion->GetVoidPointer(0))),
                         (L2Criterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(L2Criterion->GetVoidPointer(0))),
                         (OCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(OCriterion->GetVoidPointer(0))),
                         (divergence == nullptr ? nullptr :
                          static_cast<VTK_TT *>(divergence->GetVoidPointer(0))),
                         this->epsi, dfact));
    }
  }
  else if(vtkRectilinearGrid* rectilinearGrid = vtkRectilinearGrid::SafeDownCast(output))
  {
    switch (arrayType)
    { // ok to use template macro here since we made the output arrays ourselves
      vtkFloatingPointTemplateMacro(ComputeGradientsSG(
                         rectilinearGrid, array,
                         (gradients == nullptr ? nullptr :
                          static_cast<VTK_TT *>(gradients->GetVoidPointer(0))),
                         numberOfInputComponents, fieldAssociation,
                         (vorticity == nullptr ? nullptr :
                          static_cast<VTK_TT *>(vorticity->GetVoidPointer(0))),
                         (qCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(qCriterion->GetVoidPointer(0))),
                         (L2Criterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(L2Criterion->GetVoidPointer(0))),
                         (OCriterion == nullptr ? nullptr :
                          static_cast<VTK_TT *>(OCriterion->GetVoidPointer(0))),                         
                         (divergence == nullptr ? nullptr :
                          static_cast<VTK_TT *>(divergence->GetVoidPointer(0))),
                         this->epsi, dfact));

    }
  }
  if(fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_POINTS)
  {
    if(gradients)
    {
      output->GetPointData()->AddArray(gradients);
    }
    if(vorticity)
    {
      output->GetPointData()->AddArray(vorticity);
    }
    if(qCriterion)
    {
      output->GetPointData()->AddArray(qCriterion);
    }
    if(L2Criterion)
    {
      output->GetPointData()->AddArray(L2Criterion);
    }
    if(OCriterion)
    {
      output->GetPointData()->AddArray(OCriterion);
    }
    if(divergence)
    {
      output->GetPointData()->AddArray(divergence);
    }
  }
  else if(fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_CELLS)
  {
    if(gradients)
    {
      output->GetCellData()->AddArray(gradients);
    }
    if(vorticity)
    {
      output->GetCellData()->AddArray(vorticity);
    }
    if(qCriterion)
    {
      output->GetCellData()->AddArray(qCriterion);
    }
    if(L2Criterion)
    {
      output->GetCellData()->AddArray(L2Criterion);
    }
    if(OCriterion)
    {
      output->GetCellData()->AddArray(OCriterion);
    }
    if(divergence)
    {
      output->GetCellData()->AddArray(divergence);
    }
  }
  else
  {
    vtkErrorMacro("Bad fieldAssociation value " << fieldAssociation << endl);
  }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkOGSGradient::GetOutputArrayType(vtkDataArray* array)
{
  int retType = VTK_DOUBLE;
  switch (array->GetDataType())
  {
    vtkTemplateMacro(retType = GetOutputDataType(static_cast<VTK_TT>(0)));
  }
  return retType;
}

namespace {
//-----------------------------------------------------------------------------
  template<class data_type>
  void ComputePointGradientsUG(
    vtkDataSet *structure, vtkDataArray *array, data_type *gradients,
    int numberOfInputComponents, data_type* vorticity, data_type* qCriterion, data_type* L2Criterion,
    data_type* OCriterion, data_type* divergence, double eps0, double dfact, int highestCellDimension, int contributingCellOption)
  {

    vtkIdType numpts = structure->GetNumberOfPoints();
    int numberOfOutputComponents = 3*numberOfInputComponents;
    data_type epsi = 0.;
    std::vector<data_type> aux(numpts);

    #pragma omp parallel firstprivate(numpts,contributingCellOption,numberOfOutputComponents) reduction(max:epsi,highestCellDimension)
    {
    std::vector<data_type> g(numberOfOutputComponents);
    vtkNew<vtkIdList> currentPoint; currentPoint->SetNumberOfIds(1);
    vtkNew<vtkIdList> cellsOnPoint;

    // if we are doing patches for contributing cell dimensions we want to keep track of
    // the maximum expected dimension so we can exit out of the check loop quicker
    const int maxCellDimension = structure->IsA("vtkPolyData") ? 2 : 3;

    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

    for (vtkIdType point = omp_get_thread_num(); point < numpts; point+=omp_get_num_threads())
    {
      currentPoint->SetId(0, point);
      double pointcoords[3];
      structure->GetPoint(point, pointcoords);

      // Get all cells touching this point.
      #pragma omp critical
      structure->GetCellNeighbors(-1, currentPoint, cellsOnPoint);
      vtkIdType numCellNeighbors = cellsOnPoint->GetNumberOfIds();

      for(int i=0;i<numberOfOutputComponents;i++)
      {
        g[i] = 0;
      }

      if (contributingCellOption == vtkOGSGradient::Patch)
      {
        highestCellDimension = 0;
        for (vtkIdType neighbor = 0; neighbor < numCellNeighbors; neighbor++)
        {
          structure->GetCell(cellsOnPoint->GetId(neighbor),cell);
          int cellDimension = cell->GetCellDimension();
          if (cellDimension > highestCellDimension)
          {
            highestCellDimension = cellDimension;
            if (highestCellDimension == maxCellDimension)
            {
              break;
            }
          }
        }
      }
      vtkIdType numValidCellNeighbors = 0;

      // Iterate on all cells and find all points connected to current point
      // by an edge.
      for (vtkIdType neighbor = 0; neighbor < numCellNeighbors; neighbor++)
      {
        structure->GetCell(cellsOnPoint->GetId(neighbor),cell);
        if (cell->GetCellDimension() >= highestCellDimension)
        {
          int subId;
          double parametricCoord[3];
          if(GetCellParametricData(point, pointcoords, cell,
                                   subId, parametricCoord))
          {
            numValidCellNeighbors++;
            for(int inputComponent=0;inputComponent<numberOfInputComponents;inputComponent++)
            {
              int numberOfCellPoints = cell->GetNumberOfPoints();
              std::vector<double> values(numberOfCellPoints);
              // Get values of Array at cell points.
              for (int i = 0; i < numberOfCellPoints; i++)
              {
                values[i] = array->GetComponent(cell->GetPointId(i), inputComponent);
              }

              double derivative[3];
              // Get derivative of cell at point.
              cell->Derivatives(subId, parametricCoord, &values[0], 1, derivative);

              g[inputComponent*3]   += static_cast<data_type>(derivative[0]);
              g[inputComponent*3+1] += static_cast<data_type>(derivative[1]);
              g[inputComponent*3+2] += static_cast<data_type>(dfact*derivative[2]);
            } // iterating over Components
          } // if(GetCellParametricData())
        } // if(cell->GetCellDimension () >= highestCellDimension
      } // iterating over neighbors

      if (numValidCellNeighbors > 0)
      {
        for(int i=0;i<3*numberOfInputComponents;i++)
        {
          g[i] /= numValidCellNeighbors;
        }

        if(vorticity)
        {
          ComputeVorticityFromGradient(&g[0], vorticity+3*point);
        }
        if(qCriterion)
        {
          ComputeQCriterionFromGradient(&g[0], qCriterion+point);
        }
        if(L2Criterion)
        {
          ComputeLambda2CriterionFromGradient(&g[0], L2Criterion+point, 100);
        }
        if(OCriterion)
        {
          ComputeOmegaCriterionFromGradient(&g[0], OCriterion+point, &aux[point], epsi);
        }
        if(divergence)
        {
          ComputeDivergenceFromGradient(&g[0], divergence+point);
        }
        if(gradients)
        {
          for(int i=0;i<numberOfOutputComponents;i++)
          {
            gradients[point*numberOfOutputComponents+i] = g[i];
          }
        }
      }
    }  // iterating over points in grid
    }

    if(OCriterion)
    {
      epsi *= eps0;
      #pragma omp parallel firstprivate(epsi)
      {
      for (vtkIdType point = omp_get_thread_num(); point < numpts; point+=omp_get_num_threads())
        OCriterion[point] /= (aux[point] + OCriterion[point] + epsi);
      }
    }
  }

//-----------------------------------------------------------------------------
  int GetCellParametricData(vtkIdType pointId, double pointCoord[3],
                            vtkCell *cell, int &subId, double parametricCoord[3])
  {
    // Watch out for degenerate cells.  They make the derivative calculation
    // fail.
    vtkIdList *pointIds = cell->GetPointIds();
    int timesPointRegistered = 0;
    for (int i = 0; i < pointIds->GetNumberOfIds(); i++)
    {
      if (pointId == pointIds->GetId(i))
      {
        timesPointRegistered++;
      }
    }
    if (timesPointRegistered != 1)
    {
      // The cell should have the point exactly once.  Not good.
      return 0;
    }

    double dummy;
    int numpoints = cell->GetNumberOfPoints();
    std::vector<double> values(numpoints);
    // Get parametric position of point.
    cell->EvaluatePosition(pointCoord, nullptr, subId, parametricCoord,
                           dummy, &values[0]/*Really another dummy.*/);

    return 1;
  }

//-----------------------------------------------------------------------------
  template<class data_type>
    void ComputeCellGradientsUG(
      vtkDataSet *structure, vtkDataArray *array, data_type *gradients,
      int numberOfInputComponents, data_type* vorticity, data_type* qCriterion, 
      data_type* L2Criterion, data_type* OCriterion, data_type* divergence, double eps0, double dfact)
  {
    vtkIdType numcells = structure->GetNumberOfCells();
    data_type epsi = 0.;
    std::vector<data_type> aux(numcells);

    #pragma omp parallel firstprivate(numcells,numberOfInputComponents) reduction(max:epsi)
    {
    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    std::vector<double> values(8);
    std::vector<data_type> cellGradients(3*numberOfInputComponents);

    for (vtkIdType cellid = omp_get_thread_num(); cellid < numcells; cellid+=omp_get_num_threads())
    {
      structure->GetCell(cellid,cell);
      double cellCenter[3];
      int subId = cell->GetParametricCenter(cellCenter);

      int numpoints = cell->GetNumberOfPoints();
      if(static_cast<size_t>(numpoints) > values.size())
      {
        values.resize(numpoints);
      }
      double derivative[3];
      for(int inputComponent=0;inputComponent<numberOfInputComponents;
          inputComponent++)
      {
        for (int i = 0; i < numpoints; i++)
        {
          values[i] = array->GetComponent(cell->GetPointId(i), inputComponent);
        }

        cell->Derivatives(subId, cellCenter, &values[0], 1, derivative);
        cellGradients[inputComponent*3] =
          static_cast<data_type>(derivative[0]);
        cellGradients[inputComponent*3+1] =
          static_cast<data_type>(derivative[1]);
        cellGradients[inputComponent*3+2] =
          static_cast<data_type>(dfact*derivative[2]);
      }

      if(gradients)
      {
        for(int i=0;i<3*numberOfInputComponents;i++)
        {
          gradients[cellid*3*numberOfInputComponents+i] = cellGradients[i];
        }
      }
      if(vorticity)
      {
        ComputeVorticityFromGradient(&cellGradients[0], vorticity+3*cellid);
      }
      if(qCriterion)
      {
        ComputeQCriterionFromGradient(&cellGradients[0], qCriterion+cellid);
      }
      if(L2Criterion)
      {
        ComputeLambda2CriterionFromGradient(&cellGradients[0], L2Criterion+cellid, 100);
      }
      if(OCriterion)
      {
        ComputeOmegaCriterionFromGradient(&cellGradients[0], OCriterion+cellid, &aux[cellid], epsi);
      }
      if(divergence)
      {
        ComputeDivergenceFromGradient(&cellGradients[0], divergence+cellid);
      }
    }

    if(OCriterion)
    {
      epsi *= eps0;

      #pragma omp parallel firstprivate(epsi)
      {
      for (vtkIdType cellid = omp_get_thread_num(); cellid < numcells; cellid+=omp_get_num_threads())
        OCriterion[cellid] /= (aux[cellid] + OCriterion[cellid] + epsi);
      }
    }
    }
  }

//-----------------------------------------------------------------------------
  template<class Grid, class data_type>
  void ComputeGradientsSG(Grid output, vtkDataArray* array, data_type* gradients,
                          int numberOfInputComponents, int fieldAssociation,
                          data_type* vorticity, data_type* qCriterion, data_type* L2Criterion,
                          data_type* OCriterion, data_type* divergence, double eps0, double dfact)
  {
    int dims[3];
    output->GetDimensions(dims);
    if(fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_CELLS)
    {
      // reduce the dimensions by 1 for cells
      for(int i=0;i<3;i++)
      {
        dims[i]--;
      }
    }
    int ijsize = dims[0]*dims[1];

    data_type epsi = 0.;
    std::vector<data_type> aux(dims[0]*dims[1]*dims[2]);

    #pragma omp parallel firstprivate(dims,fieldAssociation,numberOfInputComponents) reduction(max:epsi)
    {

    int idx, idx2, inputComponent;
    double xp[3], xm[3], factor;
    xp[0] = xp[1] = xp[2] = xm[0] = xm[1] = xm[2] = factor = 0;
    double xxi, yxi, zxi, xeta, yeta, zeta, xzeta, yzeta, zzeta;
    yxi = zxi = xeta = yeta = zeta = xzeta = yzeta = zzeta = 0;
    double aj, xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz;
    xix = xiy = xiz = etax = etay = etaz = zetax = zetay = zetaz = 0;
    // for finite differencing -- the values on the "plus" side and
    // "minus" side of the point to be computed at
    std::vector<double> plusvalues(numberOfInputComponents);
    std::vector<double> minusvalues(numberOfInputComponents);

    std::vector<double> dValuesdXi(numberOfInputComponents);
    std::vector<double> dValuesdEta(numberOfInputComponents);
    std::vector<double> dValuesdZeta(numberOfInputComponents);
    std::vector<data_type> localGradients(numberOfInputComponents*3);

    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

    #pragma omp for collapse(3) nowait
    for (int k=0; k<dims[2]; k++)
    {
      for (int j=0; j<dims[1]; j++)
      {
        for (int i=0; i<dims[0]; i++)
        {
          //  Xi derivatives.
          if ( dims[0] == 1 ) // 2D in this direction
          {
            factor = 1.0;
            for (int ii=0; ii<3; ii++)
            {
              xp[ii] = xm[ii] = 0.0;
            }
            xp[0] = 1.0;
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = minusvalues[inputComponent] = 0;
            }
          }
          else if ( i == 0 )
          {
            factor = 1.0;
            idx = (i+1) + j*dims[0] + k*ijsize;
            idx2 = i + j*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }
          else if ( i == (dims[0]-1) )
          {
            factor = 1.0;
            idx = i + j*dims[0] + k*ijsize;
            idx2 = i-1 + j*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }
          else
          {
            factor = 0.5;
            idx = (i+1) + j*dims[0] + k*ijsize;
            idx2 = (i-1) + j*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }

          xxi = factor * (xp[0] - xm[0]);
          yxi = factor * (xp[1] - xm[1]);
          zxi = factor * (xp[2] - xm[2]);
          for(inputComponent=0;inputComponent<numberOfInputComponents;inputComponent++)
          {
            dValuesdXi[inputComponent] = factor *
              (plusvalues[inputComponent] - minusvalues[inputComponent]);
          }

          //  Eta derivatives.
          if ( dims[1] == 1 ) // 2D in this direction
          {
            factor = 1.0;
            for (int ii=0; ii<3; ii++)
            {
              xp[ii] = xm[ii] = 0.0;
            }
            xp[1] = 1.0;
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = minusvalues[inputComponent] = 0;
            }
          }
          else if ( j == 0 )
          {
            factor = 1.0;
            idx = i + (j+1)*dims[0] + k*ijsize;
            idx2 = i + j*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }
          else if ( j == (dims[1]-1) )
          {
            factor = 1.0;
            idx = i + j*dims[0] + k*ijsize;
            idx2 = i + (j-1)*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }
          else
          {
            factor = 0.5;
            idx = i + (j+1)*dims[0] + k*ijsize;
            idx2 = i + (j-1)*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }

          xeta = factor * (xp[0] - xm[0]);
          yeta = factor * (xp[1] - xm[1]);
          zeta = factor * (xp[2] - xm[2]);
          for(inputComponent=0;inputComponent<numberOfInputComponents;inputComponent++)
          {
            dValuesdEta[inputComponent] = factor *
              (plusvalues[inputComponent] - minusvalues[inputComponent]);
          }

          //  Zeta derivatives.
          if ( dims[2] == 1 ) // 2D in this direction
          {
            factor = 1.0;
            for (int ii=0; ii<3; ii++)
            {
              xp[ii] = xm[ii] = 0.0;
            }
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = minusvalues[inputComponent] = 0;
            }
            xp[2] = 1.0;
          }
          else if ( k == 0 )
          {
            factor = 1.0;
            idx = i + j*dims[0] + (k+1)*ijsize;
            idx2 = i + j*dims[0] + k*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }
          else if ( k == (dims[2]-1) )
          {
            factor = 1.0;
            idx = i + j*dims[0] + k*ijsize;
            idx2 = i + j*dims[0] + (k-1)*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }
          else
          {
            factor = 0.5;
            idx = i + j*dims[0] + (k+1)*ijsize;
            idx2 = i + j*dims[0] + (k-1)*ijsize;
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx, xp);
            GetGridEntityCoordinate(output, fieldAssociation, cell, idx2, xm);
            for(inputComponent=0;inputComponent<numberOfInputComponents;
                inputComponent++)
            {
              plusvalues[inputComponent] = array->GetComponent(idx, inputComponent);
              minusvalues[inputComponent] = array->GetComponent(idx2, inputComponent);
            }
          }

          xzeta = factor * (xp[0] - xm[0]);
          yzeta = factor * (xp[1] - xm[1]);
          zzeta = factor * (xp[2] - xm[2]);
          for(inputComponent=0;inputComponent<numberOfInputComponents;inputComponent++)
          {
            dValuesdZeta[inputComponent] = factor *
              (plusvalues[inputComponent] - minusvalues[inputComponent]);
          }

          // Now calculate the Jacobian.  Grids occasionally have
          // singularities, or points where the Jacobian is infinite (the
          // inverse is zero).  For these cases, we'll set the Jacobian to
          // zero, which will result in a zero derivative.
          //
          aj =  xxi*yeta*zzeta+yxi*zeta*xzeta+zxi*xeta*yzeta
            -zxi*yeta*xzeta-yxi*xeta*zzeta-xxi*zeta*yzeta;
          if (aj != 0.0)
          {
            aj = 1. / aj;
          }

          //  Xi metrics.
          xix  =  aj*(yeta*zzeta-zeta*yzeta);
          xiy  = -aj*(xeta*zzeta-zeta*xzeta);
          xiz  =  aj*(xeta*yzeta-yeta*xzeta);

          //  Eta metrics.
          etax = -aj*(yxi*zzeta-zxi*yzeta);
          etay =  aj*(xxi*zzeta-zxi*xzeta);
          etaz = -aj*(xxi*yzeta-yxi*xzeta);

          //  Zeta metrics.
          zetax=  aj*(yxi*zeta-zxi*yeta);
          zetay= -aj*(xxi*zeta-zxi*xeta);
          zetaz=  aj*(xxi*yeta-yxi*xeta);

          // Finally compute the actual derivatives
          idx = i + j*dims[0] + k*ijsize;
          for(inputComponent=0;inputComponent<numberOfInputComponents;inputComponent++)
          {
            localGradients[inputComponent*3] = static_cast<data_type>(
              xix*dValuesdXi[inputComponent]+etax*dValuesdEta[inputComponent]+
              zetax*dValuesdZeta[inputComponent]);

            localGradients[inputComponent*3+1] = static_cast<data_type>(
              xiy*dValuesdXi[inputComponent]+etay*dValuesdEta[inputComponent]+
              zetay*dValuesdZeta[inputComponent]);

            localGradients[inputComponent*3+2] = static_cast<data_type>(
              dfact*xiz*dValuesdXi[inputComponent]+dfact*etaz*dValuesdEta[inputComponent]+
              dfact*zetaz*dValuesdZeta[inputComponent]);
          }

          if(gradients)
          {
            for(int ii=0;ii<3*numberOfInputComponents;ii++)
            {
              gradients[idx*numberOfInputComponents*3+ii] = localGradients[ii];
            }
          }
          if(vorticity)
          {
            ComputeVorticityFromGradient(&localGradients[0], vorticity+3*idx);
          }
          if(qCriterion)
          {
            ComputeQCriterionFromGradient(&localGradients[0], qCriterion+idx);
          }
          if(L2Criterion)
          {
            ComputeLambda2CriterionFromGradient(&localGradients[0], L2Criterion+idx, 100);
          }
          if(OCriterion)
          {
            ComputeOmegaCriterionFromGradient(&localGradients[0], OCriterion+idx, &aux[idx], epsi);
          }
          if(divergence)
          {
            ComputeDivergenceFromGradient(&localGradients[0], divergence+idx);
          }
        }
      }
    }
    }

    if (OCriterion)
    {
      epsi *= eps0;

      #pragma omp parallel for collapse(3)
      for (int k=0; k<dims[2]; k++)
      {
        for (int j=0; j<dims[1]; j++)
        {
          for (int i=0; i<dims[0]; i++)
          {
            int idx = i + j*dims[0] + k*ijsize;
            OCriterion[idx] /= (aux[idx] + OCriterion[idx] + epsi);
          }
        }
      }
    }
  }

} // end anonymous namespace