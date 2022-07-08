/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// Qt includes

// qMRML includes
#include "qMRMLUtils.h"
#include "qMRMLLinearTransformSlider.h"

// MRML includes
#include "vtkMRMLTransformNode.h"

// VTK includes
#include <vtkNew.h>
#include <vtkTransform.h>
#include <vtkWeakPointer.h>

//-----------------------------------------------------------------------------
class qMRMLLinearTransformSliderPrivate
{
public:
  qMRMLLinearTransformSliderPrivate();
  qMRMLLinearTransformSlider::TransformType            TypeOfTransform;
  qMRMLLinearTransformSlider::RotationType             TypeOfRotation;
  qMRMLLinearTransformSlider::CoordinateReferenceType  CoordinateReference;
  vtkWeakPointer<vtkMRMLTransformNode>                 MRMLTransformNode;
  double                                               OldPosition;
};

// --------------------------------------------------------------------------
qMRMLLinearTransformSliderPrivate::qMRMLLinearTransformSliderPrivate()
{
  this->TypeOfTransform = qMRMLLinearTransformSlider::TRANSLATION_LR;
  this->CoordinateReference = qMRMLLinearTransformSlider::GLOBAL;
  this->MRMLTransformNode = nullptr;
  this->OldPosition = 0;
}

// --------------------------------------------------------------------------
qMRMLLinearTransformSlider::qMRMLLinearTransformSlider(QWidget* sliderParent)
  : Superclass(sliderParent)
  , d_ptr(new qMRMLLinearTransformSliderPrivate)
{
  this->setQuantity("length");
  this->setUnitAwareProperties(
    qMRMLSliderWidget::Precision | qMRMLSliderWidget::Prefix |
    qMRMLSliderWidget::Scaling | qMRMLSliderWidget::Suffix);
}

// --------------------------------------------------------------------------
qMRMLLinearTransformSlider::~qMRMLLinearTransformSlider() = default;

// --------------------------------------------------------------------------
void qMRMLLinearTransformSlider::setTypeOfTransform(TransformType _typeOfTransform)
{
  Q_D(qMRMLLinearTransformSlider);
  d->TypeOfTransform = _typeOfTransform;
  if (this->isRotation())
    {
    this->setUnitAwareProperties(qMRMLLinearTransformSlider::None);
    this->setSuffix(QString::fromLatin1("\xb0")); // "degree" character
    }
  else
    {
    this->setUnitAwareProperties(
      qMRMLSliderWidget::Precision | qMRMLSliderWidget::Prefix |
      qMRMLSliderWidget::Scaling | qMRMLSliderWidget::Suffix);
    }
  this->onMRMLTransformNodeModified(d->MRMLTransformNode);
}

void qMRMLLinearTransformSlider::setTypeOfRotation(RotationType _typeOfRotation)
{
    Q_D(qMRMLLinearTransformSlider);
    d->TypeOfRotation = _typeOfRotation;
    this->onMRMLTransformNodeModified(d->MRMLTransformNode);
}

// --------------------------------------------------------------------------
qMRMLLinearTransformSlider::TransformType qMRMLLinearTransformSlider::typeOfTransform() const
{
  Q_D(const qMRMLLinearTransformSlider);
  return d->TypeOfTransform;
}

// --------------------------------------------------------------------------
qMRMLLinearTransformSlider::RotationType qMRMLLinearTransformSlider::typeOfRotation() const
{
    Q_D(const qMRMLLinearTransformSlider);
    return d->TypeOfRotation;
}

// --------------------------------------------------------------------------
bool qMRMLLinearTransformSlider::isRotation()const
{
  return (this->typeOfTransform() == ROTATION_LR ||
          this->typeOfTransform() == ROTATION_PA ||
          this->typeOfTransform() == ROTATION_IS);
}

// --------------------------------------------------------------------------
bool qMRMLLinearTransformSlider::isTranslation()const
{
  return (this->typeOfTransform() == TRANSLATION_LR ||
          this->typeOfTransform() == TRANSLATION_PA ||
          this->typeOfTransform() == TRANSLATION_IS);
}

// --------------------------------------------------------------------------
void qMRMLLinearTransformSlider::
setCoordinateReference(CoordinateReferenceType _coordinateReference)
{
  Q_D(qMRMLLinearTransformSlider);
  d->CoordinateReference = _coordinateReference;
  this->onMRMLTransformNodeModified(d->MRMLTransformNode);
}

// --------------------------------------------------------------------------
qMRMLLinearTransformSlider::CoordinateReferenceType qMRMLLinearTransformSlider::coordinateReference() const
{
  Q_D(const qMRMLLinearTransformSlider);
  return d->CoordinateReference;
}

// --------------------------------------------------------------------------
void qMRMLLinearTransformSlider::setMRMLTransformNode(vtkMRMLTransformNode* transformNode)
{
  Q_D(qMRMLLinearTransformSlider);

  if (d->MRMLTransformNode == transformNode) { return; }

  this->qvtkReconnect(d->MRMLTransformNode, transformNode,
    vtkMRMLTransformableNode::TransformModifiedEvent,
    this, SLOT(onMRMLTransformNodeModified(vtkObject*)));

  d->MRMLTransformNode = transformNode;
  this->onMRMLTransformNodeModified(transformNode);
  // If the node is nullptr, any action on the widget is meaningless, this is why
  // the widget is disabled
  this->setEnabled(transformNode != nullptr && transformNode->IsLinear());
}

// --------------------------------------------------------------------------
vtkMRMLTransformNode* qMRMLLinearTransformSlider::mrmlTransformNode()const
{
  Q_D(const qMRMLLinearTransformSlider);
  return d->MRMLTransformNode;
}


void euler_from_matrix(double matrix[3][3], double _NEXT_AXIS[4], double result[3])
{
    /*eturn Euler angles(syxz) from rotation matrix for specified axis sequence.
    : Author :
    `Christoph Gohlke < http://www.lfd.uci.edu/~gohlke/>`_

full library with coplete set of euler triplets(combinations of  s / r x - y - z) at
http ://www.lfd.uci.edu/~gohlke/code/transformations.py.html

Note that many Euler angle triplets can describe one matrix.
"""*/

//epsilon for testing whether a number is close to zero
    double _EPS = 1e-5;

    // axis sequences for Euler angles
    //double _NEXT_AXIS[4] = { 1, 2, 0, 1 };
    int firstaxis = _NEXT_AXIS[3];
    //parity: even(0) if inner axis 'x' is followed by 'y', 'y' is followed by 'z', or 'z' is followed by 'x'.Otherwise odd(1).
    int parity = 1;
    for (int index = 1; index < 3; index++)
    {
        if (_NEXT_AXIS[index] == (_NEXT_AXIS[index + 1] + 1))
            parity = 0;
    }
    if (_NEXT_AXIS[1] == (_NEXT_AXIS[3] + 2))
        parity = 0;
    //repetition : first and last axis are same(1) or different(0).
    int repetition = (_NEXT_AXIS[1] == _NEXT_AXIS[3]);
    int frame = 0;

    int i = firstaxis;
    int j = _NEXT_AXIS[i + parity];
    int k = _NEXT_AXIS[i - parity + 1];

    double M[3][3];
    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
            M[row][col] = matrix[row][col];
    }

    double ax, ay, az, aux;
    if (repetition)
    {
        double sy = sqrt(M[i][j] * M[i][j] + M[i][k] * M[i][k]);
        if (sy > _EPS)
        {
            ax = atan2(M[i][j], M[i][k]);
            ay = atan2(sy, M[i][i]);
            az = atan2(M[j][i], -M[k][i]);
        }
        else
        {
            ax = atan2(-M[j][k], M[j][j]);
            ay = atan2(sy, M[i][i]);
            az = 0.0;
        }
    }
    else
    {
        double cy = sqrt(M[i][i] * M[i][i] + M[j][i] * M[j][i]);
        if (cy > _EPS)
        {
            ax = atan2(M[k][j], M[k][k]);
            ay = atan2(-M[k][i], cy);
            az = atan2(M[j][i], M[i][i]);
        }
        else
        {
            ax = atan2(-M[j][k], M[j][j]);
            ay = atan2(-M[k][i], cy);
            az = 0.0;
        }
    }

    if (parity)
    {
        ax = -ax;
        ay = -ay;
        az = -az;
    }
    if (frame)
    {
        aux = ax;
        ax = az;
        az = aux;
    }

    result[0] = ax;
    result[1] = ay;
    result[2] = az;
    return;
}

void rotationToVtk(double R[3][3], double _NEXT_AXIS[4], double eulerAnglesDeg_xyz[3])
{
    /*
        Concert a rotation matrix into the Mayavi / Vtk rotation paramaters(pitch, roll, yaw)
    */
    double eulerAnglesRad_012[3];
    euler_from_matrix(R, _NEXT_AXIS, eulerAnglesRad_012);

    double eulerAnglesDeg_012[3];
    eulerAnglesDeg_012[0] = eulerAnglesRad_012[0] * 180 / M_PI;
    eulerAnglesDeg_012[1] = eulerAnglesRad_012[1] * 180 / M_PI;
    eulerAnglesDeg_012[2] = eulerAnglesRad_012[2] * 180 / M_PI;

    if ((_NEXT_AXIS[1] == 0) && (_NEXT_AXIS[2] == 1) && (_NEXT_AXIS[3] == 2))
    {
        eulerAnglesDeg_xyz[0] = eulerAnglesDeg_012[0];
        eulerAnglesDeg_xyz[1] = eulerAnglesDeg_012[1];
        eulerAnglesDeg_xyz[2] = eulerAnglesDeg_012[2];
    }
    else
    {
        eulerAnglesDeg_xyz[0] = eulerAnglesDeg_012[1];
        eulerAnglesDeg_xyz[1] = eulerAnglesDeg_012[0];
        eulerAnglesDeg_xyz[2] = eulerAnglesDeg_012[2];
    }

    return;
}


// --------------------------------------------------------------------------
void qMRMLLinearTransformSlider::onMRMLTransformNodeModified(vtkObject* caller)
{
  Q_D(qMRMLLinearTransformSlider);

  vtkMRMLTransformNode* transformNode = vtkMRMLTransformNode::SafeDownCast(caller);
  if (!transformNode)
    {
    return;
    }
  Q_ASSERT(d->MRMLTransformNode == transformNode);

  bool isLinear = transformNode->IsLinear();
  this->setEnabled(isLinear);
  if (!isLinear)
    {
    return;
    }

  vtkNew<vtkTransform> transform;
  if (d->MRMLTransformNode.GetPointer() != nullptr)
    {
    qMRMLUtils::getTransformInCoordinateSystem(d->MRMLTransformNode,
      d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL, transform.GetPointer());
    }

  vtkMatrix4x4 * matrix = transform->GetMatrix();
  Q_ASSERT(matrix);
  if (!matrix) { return; }

  double _value,_ax,_ay,_az = 0.0;
  if (this->typeOfTransform() == TRANSLATION_LR)
    {
    _value = matrix->GetElement(0,3);
    }
  else if (this->typeOfTransform() == TRANSLATION_PA)
    {
    _value = matrix->GetElement(1,3);
    }
  else if (this->typeOfTransform() == TRANSLATION_IS)
    {
    _value = matrix->GetElement(2,3);
    }

  if (this->isTranslation())
    {
    // Slider values only match matrix translation values in case of GLOBAL reference
    if (d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL)
      {
      d->OldPosition = _value;
      // Only attempt to set the current slider value if the current value is reachable
      // (otherwise the value would be truncated to the current slider range).
      // If not reachable then choose closest (truncated) value and block signals to make sure we don't
      // emit signals with the truncated value.
      if (_value<this->minimum())
        {
        bool wasBlocked = this->blockSignals(true);
        this->setValue(this->minimum());
        this->blockSignals(wasBlocked);
        }
      else if (_value>this->maximum())
        {
        bool wasBlocked = this->blockSignals(true);
        this->setValue(this->maximum());
        this->blockSignals(wasBlocked);
        }
      else
        {
        this->setValue(_value);
        }
      }
    else
      {
      d->OldPosition = this->value();
      }
    }
  else if (this->isRotation())
    {
        //d->OldPosition = this->value();
        double matrix3x3[3][3];
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
                matrix3x3[row][col] = matrix->GetElement(row, col);//Transpose since it's a rotation

        double angles_deg_xyz[3];
        if (this->typeOfRotation() == INTRINSIC_ANGLES)
        {
            double configuration[4] = { 1, 2, 0, 1 };
            rotationToVtk(matrix3x3, configuration, angles_deg_xyz);
        }
        else
        {
            double configuration[4] = { 0, 0, 1, 2 };
            rotationToVtk(matrix3x3, configuration, angles_deg_xyz);
        }
        _ax = angles_deg_xyz[0];
        _ay = angles_deg_xyz[1];
        _az = angles_deg_xyz[2];


        // Slider values only match matrix rotation values in case of GLOBAL reference
        if (d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL)
        {
            if (this->typeOfTransform() == ROTATION_LR)
            {
                d->OldPosition = _ax;
                this->setValue(_ax);
            }
            else if (this->typeOfTransform() == ROTATION_PA)
            {
                d->OldPosition = _ay;
                this->setValue(_ay);
            }
            if (this->typeOfTransform() == ROTATION_IS)
            {
                d->OldPosition = _az;
                this->setValue(_az);
            }
        }
        else
        {
            if (this->typeOfTransform() == ROTATION_LR)
            {
                d->OldPosition = _ax;
            }
            else if (this->typeOfTransform() == ROTATION_PA)
            {
                d->OldPosition = _ay;
            }
            if (this->typeOfTransform() == ROTATION_IS)
            {
                d->OldPosition = _az;
            }
        }
    }
}

// --------------------------------------------------------------------------
void qMRMLLinearTransformSlider::applyTransformation(double _sliderPosition)
{
  Q_D(qMRMLLinearTransformSlider);

  if (d->MRMLTransformNode.GetPointer() == nullptr || !d->MRMLTransformNode->IsLinear())
    {
    return;
    }

  vtkNew<vtkTransform> transform;
  qMRMLUtils::getTransformInCoordinateSystem(d->MRMLTransformNode,
    d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL, transform.GetPointer());

  vtkMatrix4x4 * matrix = transform->GetMatrix();
  Q_ASSERT(matrix);
  if (!matrix) { return; }

  bool transformChanged = false;
  const double rotationChangeTolerance = 0.00001;
  const double translationChangeTolerance = 0.00001;

  if (this->typeOfRotation() == EXTRINSIC_ANGLES)
  {
      if (this->typeOfTransform() == ROTATION_LR)
      {
          double angle = _sliderPosition - d->OldPosition;
          transform->RotateX(angle);
          if (fabs(angle) > rotationChangeTolerance)
          {
              transformChanged = true;
          }
      }
      else if (this->typeOfTransform() == ROTATION_PA)
      {
          double angle = _sliderPosition - d->OldPosition;
          transform->RotateY(angle);
          if (fabs(angle) > rotationChangeTolerance)
          {
              transformChanged = true;
          }
      }
      else if (this->typeOfTransform() == ROTATION_IS)
      {
          double angle = _sliderPosition - d->OldPosition;
          transform->RotateZ(angle);
          if (fabs(angle) > rotationChangeTolerance)
          {
              transformChanged = true;
          }
      }
  }
  else if (this->typeOfRotation() == INTRINSIC_ANGLES)
  {
      double matrix3x3[3][3];
      for (int row = 0; row < 3; row++)
          for (int col = 0; col < 3; col++)
              matrix3x3[row][col] = matrix->GetElement(row, col);

      double angles_deg_xyz[3];
      double configuration[4] = { 1, 2, 0, 1 };
      rotationToVtk(matrix3x3, configuration, angles_deg_xyz);

      if (this->typeOfTransform() == ROTATION_LR)
      {
          double angle = _sliderPosition - d->OldPosition;

          transform->RotateY(angles_deg_xyz[1]);
          double rotationAxis[4] = { 1,0,0,0 };
          transform->GetMatrix()->MultiplyPoint(rotationAxis, rotationAxis);
          transform->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);
          rotationAxis[0] = 0;
          rotationAxis[1] = 0;
          rotationAxis[2] = 1;
          transform->GetMatrix()->MultiplyPoint(rotationAxis, rotationAxis);
          transform->RotateWXYZ(angles_deg_xyz[2], rotationAxis[0], rotationAxis[1], rotationAxis[2]);

          if (fabs(angle) > rotationChangeTolerance)
          {
              transformChanged = true;
          }
      }
      else if (this->typeOfTransform() == ROTATION_PA)
      {
          double angle = _sliderPosition - d->OldPosition;

          transform->RotateY(angle);
          double rotationAxis[4] = { 1,0,0,0 };
          transform->GetMatrix()->MultiplyPoint(rotationAxis, rotationAxis);
          transform->RotateWXYZ(angles_deg_xyz[0], rotationAxis[0], rotationAxis[1], rotationAxis[2]);
          rotationAxis[0] = 0;
          rotationAxis[1] = 0;
          rotationAxis[2] = 1;
          transform->GetMatrix()->MultiplyPoint(rotationAxis, rotationAxis);
          transform->RotateWXYZ(angles_deg_xyz[2], rotationAxis[0], rotationAxis[1], rotationAxis[2]);

          if (fabs(angle) > rotationChangeTolerance)
          {
              transformChanged = true;
          }
      }
      else if (this->typeOfTransform() == ROTATION_IS)
      {
          double angle = _sliderPosition - d->OldPosition;

          transform->RotateY(angles_deg_xyz[1]);
          double rotationAxis[4] = { 1,0,0,0 };
          transform->GetMatrix()->MultiplyPoint(rotationAxis, rotationAxis);
          transform->RotateWXYZ(angles_deg_xyz[0], rotationAxis[0], rotationAxis[1], rotationAxis[2]);
          rotationAxis[0] = 0;
          rotationAxis[1] = 0;
          rotationAxis[2] = 1;
          transform->GetMatrix()->MultiplyPoint(rotationAxis, rotationAxis);
          transform->RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2]);

          if (fabs(angle) > rotationChangeTolerance)
          {
              transformChanged = true;
          }
      }
  }
  else if (this->typeOfTransform() == TRANSLATION_LR)
    {
    double vector[] = {0., 0., 0.};
    // Slider values only match matrix translation values in case of GLOBAL reference
    if (d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL)
      {
      vector[0] = _sliderPosition - matrix->GetElement(0,3);
      }
    else
      {
      vector[0] = _sliderPosition - d->OldPosition;
      }
    transform->Translate(vector);
    if (fabs(vector[0])>translationChangeTolerance)
      {
      transformChanged = true;
      }
    }
  else if (this->typeOfTransform() == TRANSLATION_PA)
    {
    double vector[] = {0., 0., 0.};
    // Slider values only match matrix translation values in case of GLOBAL reference
    if (d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL)
      {
      vector[1] = _sliderPosition - matrix->GetElement(1,3);
      }
    else
      {
      vector[1] = _sliderPosition - d->OldPosition;
      }
    transform->Translate(vector);
    if (fabs(vector[1])>translationChangeTolerance)
      {
      transformChanged = true;
      }
    }
  else if (this->typeOfTransform() == TRANSLATION_IS)
    {
    double vector[] = {0., 0., 0.};
    // Slider values only match matrix translation values in case of GLOBAL reference
    if (d->CoordinateReference == qMRMLLinearTransformSlider::GLOBAL)
      {
      vector[2] = _sliderPosition - matrix->GetElement(2,3);
      }
    else
      {
      vector[2] = _sliderPosition - d->OldPosition;
      }
    transform->Translate(vector);
    if (fabs(vector[2])>translationChangeTolerance)
      {
      transformChanged = true;
      }
    }
  d->OldPosition = _sliderPosition;

  if (transformChanged)
    {
  d->MRMLTransformNode->SetMatrixTransformToParent(transform->GetMatrix());
    }
}
