#ifndef VTK_H
#define VTK_H

#ifdef HAVE_VTK
#include <vtkBMPWriter.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellDerivatives.h>
#include <vtkCell.h>
#include <vtkCellType.h>
#include <vtkClipDataSet.h>
#include <vtkContourGrid.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkDoubleArray.h>
#include <vtkExtractVectorComponents.h>
#include <vtkFloatArray.h>
#include <vtkGenericCell.h>
#include <vtkHexahedron.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPolyData.h>
#include <vtkShortArray.h>
#include <vtkStructuredGrid.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkZLibDataCompressor.h>

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

#endif
#endif
