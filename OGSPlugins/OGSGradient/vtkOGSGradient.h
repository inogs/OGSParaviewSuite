/*=========================================================================

  Program:   OGS Gradient
  Module:    vtkOGSGradient.h

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

/**
 * @class   vtkOGSGradient
 * @brief   A general filter for gradient estimation.
 *
 *
 * Estimates the gradient of a field in a data set.  The gradient calculation
 * is dependent on the input dataset type.  The created gradient array is
 * of the same type as the array it is calculated from (e.g. point data or cell
 * data) but the data type will be either float or double.  At the boundary
 * the gradient is not central differencing.  The output gradient array has
 * 3*number of components of the input data array.  The ordering for the
 * output gradient tuple will be {du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, dw/dx,
 * dw/dy, dw/dz} for an input array {u, v, w}. There are also the options
 * to additionally compute the vorticity and Q criterion of a vector field.
 * Unstructured grids and polydata can have cells of different dimensions.
 * To handle different use cases in this situation, the user can specify
 * which cells contribute to the computation. The options for this are
 * All (default), Patch and DataSetMax. Patch uses only the highest dimension
 * cells attached to a point. DataSetMax uses the highest cell dimension in
 * the entire data set. For Patch or DataSetMax it is possible that some values
 * will not be computed. The ReplacementValueOption specifies what to use
 * for these values.
*/

#ifndef vtkOGSGradient_h
#define vtkOGSGradient_h

#include "vtkDataSet.h"
#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkDataSetAlgorithm.h"

#include "vtkPVConfig.h" // For PARAVIEW_USE_MPI

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class vtkOGSGradient : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkOGSGradient, vtkDataSetAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent) override;

  /// Options to choose what cells contribute to the gradient calculation
  enum ContributingCellEnum {
    All=0,        //!< All cells
    Patch=1,      //!< Highest dimension cells in the patch of cells contributing to the calculation
    DataSetMax=2  //!< Highest dimension cells in the data set
  };

  /// The replacement value or entities that don't have any gradient computed over them based on
  /// the ContributingCellOption
  enum ReplacementValueEnum {
    Zero=0,        //!< 0
    NaN=1,         //!< NaN
    DataTypeMin=2, //!< The minimum possible value of the input array data type
    DataTypeMax=3  //!< The maximum possible value of the input array data type
  };

  static vtkOGSGradient *New();

  //@{
  /**
   * These are basically a convenience method that calls SetInputArrayToProcess
   * to set the array used as the input scalars.  The fieldAssociation comes
   * from the vtkDataObject::FieldAssociations enum.  The fieldAttributeType
   * comes from the vtkDataSetAttributes::AttributeTypes enum.
   */
  virtual void SetInputScalars(int fieldAssociation, const char *name);
  virtual void SetInputScalars(int fieldAssociation, int fieldAttributeType);
  //@}

  //@{
  /**
   * Get/Set the name of the gradient array to create.  This is only
   * used if ComputeGradient is non-zero. If nullptr (the
   * default) then the output array will be named "Gradients".
   */
  vtkGetStringMacro(ResultArrayName);
  vtkSetStringMacro(ResultArrayName);
  //@}

  //@{
  /**
   * Get/Set the name of the divergence array to create. This is only
   * used if ComputeDivergence is non-zero. If nullptr (the
   * default) then the output array will be named "Divergence".
   */
  vtkGetStringMacro(DivergenceArrayName);
  vtkSetStringMacro(DivergenceArrayName);
  //@}

  //@{
  /**
   * Get/Set the name of the vorticity array to create. This is only
   * used if ComputeVorticity is non-zero. If nullptr (the
   * default) then the output array will be named "Vorticity".
   */
  vtkGetStringMacro(VorticityArrayName);
  vtkSetStringMacro(VorticityArrayName);
  //@}

  //@{
  /**
   * Get/Set the name of the Q criterion array to create. This is only
   * used if ComputeQCriterion is non-zero. If nullptr (the
   * default) then the output array will be named "Q-criterion".
   */
  vtkGetStringMacro(QCriterionArrayName);
  vtkSetStringMacro(QCriterionArrayName);
  //@}

  //@{
  /**
   * Get/Set the name of the Lambda2 criterion array to create. This is only
   * used if ComputeLambda2Criterion is non-zero. If nullptr (the
   * default) then the output array will be named "L2-criterion".
   */
  vtkGetStringMacro(Lambda2CriterionArrayName);
  vtkSetStringMacro(Lambda2CriterionArrayName);
  //@}

  /**
   * Get/Set the name of the Omega criterion array to create. This is only
   * used if ComputeOmegaCriterion is non-zero. If nullptr (the
   * default) then the output array will be named "Omega-criterion".
   */
  vtkGetStringMacro(OmegaCriterionArrayName);
  vtkSetStringMacro(OmegaCriterionArrayName);
  //@}

 //@{
 /**
  * When this flag is on (default is off), the gradient filter will provide a
  * less accurate (but close) algorithm that performs fewer derivative
  * calculations (and is therefore faster).  The error contains some smoothing
  * of the output data and some possible errors on the boundary.  This
  * parameter has no effect when performing the gradient of cell data.
  * This only applies if the input grid is a vtkUnstructuredGrid or a
  * vtkPolyData.
  */
  vtkGetMacro(FasterApproximation, int);
  vtkSetMacro(FasterApproximation, int);
  vtkBooleanMacro(FasterApproximation, int);
 //@}

  //@{
  /**
   * Add the gradient to the output field data.  The name of the array
   * will be ResultArrayName and will be the same type as the input
   * array. The default is on.
   */
  vtkSetMacro(ComputeGradient, int);
  vtkGetMacro(ComputeGradient, int);
  vtkBooleanMacro(ComputeGradient, int);
  //@}

  //@{
  /**
   * Add divergence to the output field data.  The name of the array
   * will be DivergenceArrayName and will be the same type as the input
   * array.  The input array must have 3 components in order to
   * compute this. The default is off.
   */
  vtkSetMacro(ComputeDivergence, int);
  vtkGetMacro(ComputeDivergence, int);
  vtkBooleanMacro(ComputeDivergence, int);
  //@}

  //@{
  /**
   * Add voriticity/curl to the output field data.  The name of the array
   * will be VorticityArrayName and will be the same type as the input
   * array.  The input array must have 3 components in order to
   * compute this. The default is off.
   */
  vtkSetMacro(ComputeVorticity, int);
  vtkGetMacro(ComputeVorticity, int);
  vtkBooleanMacro(ComputeVorticity, int);
  //@}

  //@{
  /**
   * Add Q-criterion to the output field data.  The name of the array
   * will be QCriterionArrayName and will be the same type as the input
   * array.  The input array must have 3 components in order to
   * compute this.  Note that Q-citerion is a balance of the rate
   * of vorticity and the rate of strain. The default is off.
   */
  vtkSetMacro(ComputeQCriterion, int);
  vtkGetMacro(ComputeQCriterion, int);
  vtkBooleanMacro(ComputeQCriterion, int);
  //@}

  /**
   * Add Lambda2-criterion to the output field data.  The name of the array
   * will be Lambda2CriterionArrayName and will be the same type as the input
   * array.  The input array must have 3 components in order to
   * compute this.  Note that Lambda2-criterion the second eigenvalue of
     the strain and rotation tensor squared. The default is off.
   */
  vtkSetMacro(ComputeLambda2Criterion, int);
  vtkGetMacro(ComputeLambda2Criterion, int);
  vtkBooleanMacro(ComputeLambda2Criterion, int);
  //@}

  /**
   * Add Omega-criterion to the output field data.  The name of the array
   * will be OmegaCriterionArrayName and will be the same type as the input
   * array.  The input array must have 3 components in order to
   * compute this. The default is off.
   */
  vtkSetMacro(ComputeOmegaCriterion, int);
  vtkGetMacro(ComputeOmegaCriterion, int);
  vtkBooleanMacro(ComputeOmegaCriterion, int);
  //@}

  /**
   * Small number to avoid division per zero in the Omega criterion.
   */
  vtkSetMacro(epsi, double);
  vtkGetMacro(epsi, double);
  //@}

  //@{
  /**
   * Option to specify what cells to include in the gradient computation.
   * Options are all cells (All, Patch and DataSetMax). The default is All.
   */
  vtkSetClampMacro(ContributingCellOption, int, 0, 2);
  vtkGetMacro(ContributingCellOption, int);
  //@}

  //@{
  /**
   * Option to specify what replacement value or entities that don't have any gradient computed
   * over them based on the ContributingCellOption. Options are (Zero, NaN, DataTypeMin,
   * DataTypeMax). The default is Zero.
   */
  vtkSetClampMacro(ReplacementValueOption, int, 0, 3);
  vtkGetMacro(ReplacementValueOption, int);
  //@}

  #ifdef PARAVIEW_USE_MPI
    // Description:
    // Set the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
  #endif

protected:
  vtkOGSGradient();
  ~vtkOGSGradient() override;

  int RequestUpdateExtent(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  /**
   * Compute the gradients for grids that are not a vtkImageData,
   * vtkRectilinearGrid, or vtkStructuredGrid.
   * Returns non-zero if the operation was successful.
   */
  virtual int ComputeUnstructuredGridGradient(
    vtkDataArray* Array, int fieldAssociation, vtkDataSet* input,
    bool computeVorticity, bool computeQCriterion, bool computeLambda2Criterion,
     bool computeOmegaCriterion, bool computeDivergence, double dfact, vtkDataSet* output);

  /**
   * Compute the gradients for either a vtkImageData, vtkRectilinearGrid or
   * a vtkStructuredGrid.  Computes the gradient using finite differences.
   * Returns non-zero if the operation was successful.
   */
  virtual int ComputeRegularGridGradient(
    vtkDataArray* Array, int fieldAssociation, bool computeVorticity, bool computeQCriterion, 
    bool computeOmegaCriterion, bool computeLambda2Criterion, bool computeDivergence, double dfact, vtkDataSet* output);

  /**
   * Get the proper array type to compute requested derivative quantities for.
   * If the input array is an unsigned type then we switch to a float array
   * for the output but otherwise the output array type is the same as the
   * input array type.
   **/
  int GetOutputArrayType(vtkDataArray* inputArray);

  /**
   * If non-null then it contains the name of the outputted gradient array.
   * By default it is "Gradients".
   */
  char *ResultArrayName;

  /**
   * If non-null then it contains the name of the outputted divergence array.
   * By default it is "Divergence".
   */
  char *DivergenceArrayName;

  /**
   * If non-null then it contains the name of the outputted vorticity array.
   * By default it is "Vorticity".
   */
  char *VorticityArrayName;

  /**
   * If non-null then it contains the name of the outputted Q criterion array.
   * By default it is "Q-criterion".
   */
  char *QCriterionArrayName;

  /**
   * If non-null then it contains the name of the outputted L2 criterion array.
   * By default it is "L2-criterion".
   */
  char *Lambda2CriterionArrayName;

  /**
   * If non-null then it contains the name of the outputted Omega criterion array.
   * By default it is "Omega-criterion".
   */
  char *OmegaCriterionArrayName;

  /**
   * When this flag is on (default is off), the gradient filter will provide a
   * less accurate (but close) algorithm that performs fewer derivative
   * calculations (and is therefore faster).  The error contains some smoothing
   * of the output data and some possible errors on the boundary.  This
   * parameter has no effect when performing the gradient of cell data.
   * This only applies if the input grid is a vtkUnstructuredGrid or a
   * vtkPolyData.
   */
  int FasterApproximation;

  /**
   * Flag to indicate that the gradient of the input vector is to
   * be computed. By default ComputeDivergence is on.
   */
  int ComputeGradient;

  /**
   * Flag to indicate that the divergence of the input vector is to
   * be computed.  The input array to be processed must have
   * 3 components.  By default ComputeDivergence is off.
   */
  int ComputeDivergence;

  /**
   * Flag to indicate that the Q-criterion of the input vector is to
   * be computed.  The input array to be processed must have
   * 3 components.  By default ComputeQCriterion is off.
   */
  int ComputeQCriterion;

  /**
   * Flag to indicate that the Lambda2-criterion of the input vector is to
   * be computed.  The input array to be processed must have
   * 3 components.  By default ComputeLambda2Criterion is off.
   */
  int ComputeLambda2Criterion;

  /**
   * Flag to indicate that the Omega-criterion of the input vector is to
   * be computed.  The input array to be processed must have
   * 3 components.  By default ComputeOmegaCriterion is off.
   */
  int ComputeOmegaCriterion;

  /**
   * Flag to indicate that vorticity/curl of the input vector is to
   * be computed.  The input array to be processed must have
   * 3 components.  By default ComputeVorticity is off.
   */
  int ComputeVorticity;

  /**
   * Option to specify what cells to include in the gradient computation.
   * Options are all cells (All, Patch and DataSet). The default is all.
   */
  int ContributingCellOption;

  /**
   * Option to specify what replacement value or entities that don't have any gradient computed
   * over them based on the ContributingCellOption. Options are (Zero, NaN, DataTypeMin,
   * DataTypeMax). The default is Zero.
   */
  int ReplacementValueOption;

  #ifdef PARAVIEW_USE_MPI
    vtkMultiProcessController* Controller;
  #endif

private:
  vtkOGSGradient(const vtkOGSGradient &) = delete;
  void operator=(const vtkOGSGradient &) = delete;

  int procId, nProcs;
  double epsi;
};

#endif //_vtkGradientFilter_h
