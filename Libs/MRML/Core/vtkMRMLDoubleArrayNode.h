/*=auto=========================================================================

  Portions (c) Copyright 2009 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLCurveAnalysisNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/

#ifndef __vtkMRMLDoubleArrayNode_h
#define __vtkMRMLDoubleArrayNode_h

#include <string>
#include <vector>

#include "vtkMRMLStorableNode.h"
class vtkDoubleArray;
class vtkMRMLStorageNode;

class VTK_MRML_EXPORT vtkMRMLDoubleArrayNode : public vtkMRMLStorableNode
{
public:
  //----------------------------------------------------------------
  /// Constants
  //----------------------------------------------------------------

  /// Interpolation method
  enum {
    INTERP_LINEAR = 0
  };

  //----------------------------------------------------------------
  /// Standard methods for MRML nodes
  //----------------------------------------------------------------

  static vtkMRMLDoubleArrayNode *New();
  vtkTypeMacro(vtkMRMLDoubleArrayNode,vtkMRMLStorableNode);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  virtual vtkMRMLNode* CreateNodeInstance() VTK_OVERRIDE;

  ///
  /// Set node attributes
  virtual void ReadXMLAttributes( const char** atts) VTK_OVERRIDE;

  ///
  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent) VTK_OVERRIDE;

  ///
  /// Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node) VTK_OVERRIDE;

  ///
  /// Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() VTK_OVERRIDE
    {return "DoubleArray";}

  ///
  /// Method to propagate events generated in mrml
  virtual void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, void *callData ) VTK_OVERRIDE;

  //----------------------------------------------------------------
  /// Get and Set Macros
  //----------------------------------------------------------------
  virtual void SetArray(vtkDoubleArray*);
  vtkGetObjectMacro ( Array, vtkDoubleArray );

  //----------------------------------------------------------------
  /// Access methods
  //----------------------------------------------------------------

  ///
  /// Set / change the size of value data (number of points in time dimension)
  void SetSize(unsigned int n);

  ///
  /// Get the size of value data (number of points in time dimension)
  unsigned int GetSize();

  ///
  /// Get Y value by X. If X is between two data points, it interpolates the value
  /// by using the method specified by 'interp'.
  /// NOTE: TO BE IMPLEMENTED.
  double GetYAxisValue(double x, int interp=INTERP_LINEAR);

  /// \bug Fix function GetXYAxisValue(int index, double* x, double* y);
  /// Get X and Y values at the data point specified by 'index'
  /// int GetXYAxisValue(int index, double* x, double* y);

  /// \bug Fix function GetXYAxisValue(int index, double* x, double* y, double*  yerr);
  /// Get X and Y values with error (or standard deviation) of Y,
  /// at the data point specified by 'index'
  /// int GetXYAxisValue(int index, double* x, double* y, double*  yerr);

  ///
  /// Set values at the data point specified by 'index'
  int SetValues(int index, double* values);

  ///
  /// Set value at the data point specified by 'index' at the given 'component'
  int SetValue(int index, int component, double value);

  ///
  /// Set X and Y values at the data point specified by 'index'
  int SetXYValue(int index, double x, double y);

  ///
  /// Set X and Y values with error (or standard deviation) of Y,
  /// at the data point specified by 'index'
  int SetXYValue(int index, double x, double y, double yerr);

  ///
  /// Get values at the data point specified by 'index'
  int GetValues(int index, double* values);

  ///
  /// Get value at the data point specified by 'index' at the given 'component'
  double GetValue(int index, int component, int& success);
  
  ///
  /// Get X and Y values at the data point specified by 'index'
  int GetXYValue(int index, double* x, double* y);

  ///
  /// Get X and Y values with error (or standard deviation) of Y,
  /// at the data point specified by 'index'
  int GetXYValue(int index, double* x, double* y, double* yerr);

  ///
  /// Add values at the end of the array
  int AddValues(double* values);

  ///
  /// Add value at the end of the array at the given 'component'
  int AddValue(int component, double value);

  ///
  /// Add a set of X and Y values at the end of the array
  int AddXYValue(double x, double y);

  ///
  /// Add a set of X and Y values with error (or standard deviation) of Y,
  /// at the end of the array
  int AddXYValue(double x, double y, double yerr);

  ///
  /// Search min and maximum value of X and Y in the array. The result is stored in 'range'.
  /// (range[0]: minimum value, range[1]: maximum value)
  /// if fIncludeError=1 is specified, the range takes account of errors.
  void GetRange(double* rangeX, double* rangeY, int fIncludeError=1);

  ///
  /// Search min and maximum value of X in the array. The result is stored in 'range'.
  /// (range[0]: minimum value, range[1]: maximum value)
  void GetXRange(double* range);

  ///
  /// Search min and maximum value of Y in the array. The result is stored in 'range'.
  /// (range[0]: minimum value, range[1]: maximum value)
  /// if fIncludeError=1 is specified, the range takes account of errors.
  void GetYRange(double* range, int fIncludeError=1);

  // Description:
  //Set labels
  //void SetLabel(std::vector< std::string > labels);
  typedef std::vector< std::string > LabelsVectorType;
  void SetLabels(const LabelsVectorType &labels);

  // Description:
  //Get labels
  const LabelsVectorType & GetLabels() const;

  ///
  /// Create default storage node or NULL if does not have one
  virtual vtkMRMLStorageNode* CreateDefaultStorageNode() VTK_OVERRIDE;

  //----------------------------------------------------------------
  /// Constructor and destroctor
  //----------------------------------------------------------------
 protected:
  vtkMRMLDoubleArrayNode();
  ~vtkMRMLDoubleArrayNode();
  vtkMRMLDoubleArrayNode(const vtkMRMLDoubleArrayNode&);
  void operator=(const vtkMRMLDoubleArrayNode&);


 protected:
  //----------------------------------------------------------------
  /// Data
  //----------------------------------------------------------------

  vtkDoubleArray* Array;

  std::vector< std::string > Unit;
  std::vector< std::string > Labels;


};

#endif
