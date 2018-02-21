/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH)
  All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer

=========================================================================auto=*/

#include "vtkMRMLCoreTestingMacros.h"
#include "vtkMRMLDiffusionTensorVolumeNode.h"

int vtkMRMLDiffusionTensorVolumeNodeTest1(int , char * [] )
{
  vtkNew<vtkMRMLDiffusionTensorVolumeNode> node1;
  EXERCISE_ALL_BASIC_MRML_METHODS(node1.GetPointer());
  return EXIT_SUCCESS;
}