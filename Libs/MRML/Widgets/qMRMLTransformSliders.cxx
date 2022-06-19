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

#include "qMRMLTransformSliders.h"
#include "ui_qMRMLTransformSliders.h"

// Qt includes
#include <QStack>

// qMRML includes
#include <qMRMLUtils.h>

// MRML includes
#include "vtkMRMLTransformNode.h"

// VTK includes
#include <vtkNew.h>
#include <vtkTransform.h>

# define M_PI           3.14159265358979323846  /* pi */
#include <math.h>

//-----------------------------------------------------------------------------
class qMRMLTransformSlidersPrivate: public Ui_qMRMLTransformSliders
{
public:
  qMRMLTransformSlidersPrivate()
    {
    this->TypeOfTransform = -1;
    this->MRMLTransformNode = nullptr;
    }

  int                                    TypeOfTransform;
  vtkMRMLTransformNode*                  MRMLTransformNode;
  QStack<qMRMLLinearTransformSlider*>    ActiveSliders;
};

// --------------------------------------------------------------------------
qMRMLTransformSliders::qMRMLTransformSliders(QWidget* slidersParent)
  : Superclass(slidersParent)
  , d_ptr(new qMRMLTransformSlidersPrivate)
{
  Q_D(qMRMLTransformSliders);

  d->setupUi(this);

  ctkDoubleSpinBox::DecimalsOptions decimalsOptions(ctkDoubleSpinBox::DecimalsByShortcuts | ctkDoubleSpinBox::DecimalsByKey |
    ctkDoubleSpinBox::InsertDecimals | ctkDoubleSpinBox::DecimalsAsMin);
  d->LRSlider->spinBox()->setDecimalsOption(decimalsOptions);
  d->PASlider->spinBox()->setDecimalsOption(decimalsOptions);
  d->ISSlider->spinBox()->setDecimalsOption(decimalsOptions);
  d->LRSlider->setSynchronizeSiblings(ctkSliderWidget::SynchronizeDecimals);
  d->PASlider->setSynchronizeSiblings(ctkSliderWidget::SynchronizeDecimals);
  d->ISSlider->setSynchronizeSiblings(ctkSliderWidget::SynchronizeDecimals);

  this->setCoordinateReference(qMRMLTransformSliders::GLOBAL);
  this->setTypeOfTransform(qMRMLTransformSliders::TRANSLATION);

  this->connect(d->LRSlider, SIGNAL(valueChanged(double)),
                SLOT(onSliderPositionChanged(double)));
  this->connect(d->PASlider, SIGNAL(valueChanged(double)),
                SLOT(onSliderPositionChanged(double)));
  this->connect(d->ISSlider, SIGNAL(valueChanged(double)),
                SLOT(onSliderPositionChanged(double)));

  this->connect(d->MinValueSpinBox, SIGNAL(valueChanged(double)),
                SLOT(onMinimumChanged(double)));
  this->connect(d->MaxValueSpinBox, SIGNAL(valueChanged(double)),
                SLOT(onMaximumChanged(double)));
  // the default values of min and max are set in the .ui file
  this->onMinimumChanged(d->MinValueSpinBox->value());
  this->onMaximumChanged(d->MaxValueSpinBox->value());

  this->connect(d->LRSlider, SIGNAL(decimalsChanged(int)),
                SIGNAL(decimalsChanged(int)));

  // disable as there is not MRML Node associated with the widget
  this->setEnabled(false);
}

// --------------------------------------------------------------------------
qMRMLTransformSliders::~qMRMLTransformSliders() = default;

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setCoordinateReference(CoordinateReferenceType _coordinateReference)
{
  Q_D(qMRMLTransformSliders);

  qMRMLLinearTransformSlider::CoordinateReferenceType ref =
      static_cast<qMRMLLinearTransformSlider::CoordinateReferenceType>(
        _coordinateReference);

  if (this->coordinateReference() != _coordinateReference)
    {
    // reference changed
    if (this->typeOfTransform() == qMRMLTransformSliders::ROTATION
      || (this->typeOfTransform() == qMRMLTransformSliders::TRANSLATION
          && ref == qMRMLLinearTransformSlider::LOCAL) )
      {
      // No one-to-one correspondence between slider and transform matrix values
      bool blocked = false;
      blocked = d->LRSlider->blockSignals(true);
      d->LRSlider->reset();
      d->LRSlider->blockSignals(blocked);
      blocked = d->PASlider->blockSignals(true);
      d->PASlider->reset();
      d->PASlider->blockSignals(blocked);
      blocked = d->ISSlider->blockSignals(true);
      d->ISSlider->reset();
      d->ISSlider->blockSignals(blocked);
      }
    else
      {
      // make sure the current translation values can be set on the slider
      updateRangeFromTransform(d->MRMLTransformNode);
      }
    d->LRSlider->setCoordinateReference(ref);
    d->PASlider->setCoordinateReference(ref);
    d->ISSlider->setCoordinateReference(ref);
    }
}

// --------------------------------------------------------------------------
qMRMLTransformSliders::CoordinateReferenceType qMRMLTransformSliders::coordinateReference() const
{
  Q_D(const qMRMLTransformSliders);

  // Assumes settings of the sliders are all the same
  qMRMLLinearTransformSlider::CoordinateReferenceType ref =
    d->LRSlider->coordinateReference();
  return (ref == qMRMLLinearTransformSlider::GLOBAL) ? GLOBAL : LOCAL;
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setTypeOfTransform(TransformType _typeOfTransform)
{
  Q_D(qMRMLTransformSliders);

  if (d->TypeOfTransform == _typeOfTransform) { return; }
  if (_typeOfTransform == qMRMLTransformSliders::TRANSLATION)
    {
    d->LRSlider->setTypeOfTransform(qMRMLLinearTransformSlider::TRANSLATION_LR);
    d->PASlider->setTypeOfTransform(qMRMLLinearTransformSlider::TRANSLATION_PA);
    d->ISSlider->setTypeOfTransform(qMRMLLinearTransformSlider::TRANSLATION_IS);
    }
  else if (_typeOfTransform == qMRMLTransformSliders::ROTATION)
    {
    d->LRSlider->setTypeOfTransform(qMRMLLinearTransformSlider::ROTATION_LR);
    d->PASlider->setTypeOfTransform(qMRMLLinearTransformSlider::ROTATION_PA);
    d->ISSlider->setTypeOfTransform(qMRMLLinearTransformSlider::ROTATION_IS);

    // Range of Rotation sliders should be fixed to (-180,180)
    //this->setRange(-180.00, 180.00);
    d->LRSlider->setRange(-90.00, 90.00);
    d->PASlider->setRange(-180.00, 180.00);
    d->ISSlider->setRange(-180.00, 180.00);
    }
  d->TypeOfTransform = _typeOfTransform;
}

// --------------------------------------------------------------------------
qMRMLTransformSliders::TransformType qMRMLTransformSliders::typeOfTransform() const
{
  Q_D(const qMRMLTransformSliders);
  return static_cast<qMRMLTransformSliders::TransformType>(d->TypeOfTransform);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setMRMLTransformNode(vtkMRMLNode* node)
{
  this->setMRMLTransformNode(vtkMRMLTransformNode::SafeDownCast(node));
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setMRMLTransformNode(vtkMRMLTransformNode* transformNode)
{
  Q_D(qMRMLTransformSliders);

  if (d->MRMLTransformNode == transformNode)
    {
    // no change
    return;
    }

  this->qvtkReconnect(d->MRMLTransformNode, transformNode,
                      vtkMRMLTransformableNode::TransformModifiedEvent,
                      this, SLOT(onMRMLTransformNodeModified(vtkObject*)));

  bool blocked = d->LRSlider->blockSignals(true);
  d->LRSlider->reset();
  d->LRSlider->blockSignals(blocked);
  blocked = d->PASlider->blockSignals(true);
  d->PASlider->reset();
  d->PASlider->blockSignals(blocked);
  blocked = d->ISSlider->blockSignals(true);
  d->ISSlider->reset();
  d->ISSlider->blockSignals(blocked);

  this->onMRMLTransformNodeModified(transformNode);

  d->LRSlider->setMRMLTransformNode(transformNode);
  d->PASlider->setMRMLTransformNode(transformNode);
  d->ISSlider->setMRMLTransformNode(transformNode);

  // If the node is nullptr, any action on the widget is meaningless, this is why
  // the widget is disabled
  this->setEnabled(transformNode != nullptr && transformNode->IsLinear());
  d->MRMLTransformNode = transformNode;
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::onMRMLTransformNodeModified(vtkObject* caller)
{
  vtkMRMLTransformNode* transformNode = vtkMRMLTransformNode::SafeDownCast(caller);
  if (!transformNode)
    {
    return;
    }
  Q_ASSERT(transformNode);
  bool isLinear = transformNode->IsLinear();
  this->setEnabled(isLinear);
  if (!isLinear)
    {
    return;
    }

  // There is no one-to-one correspondence between matrix values and slider position if transform type is rotation;
  // or transform type is translation and coordinate reference is global. In these cases the slider range must not be updated:
  // it is not necessary (as the slider will be reset to 0 anyway when another slider is moved) and changing the slider range
  // can even cause instability (transform value increasing continuously) when the user drags the slider using the mouse.
  if (this->typeOfTransform() == qMRMLTransformSliders::TRANSLATION && coordinateReference() == LOCAL)
    {
    return;
    }

  this->updateRangeFromTransform(transformNode);

  if (this->typeOfTransform() == qMRMLTransformSliders::ROTATION)
  {
    this->updateAngleValuesFromTransform(transformNode);
  }
}

void euler_from_matrix(double matrix[3][3], double result[3])
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
    double _NEXT_AXIS[4] = { 1, 2, 0, 1 };
    int firstaxis = 1;
    int parity = 1;
    int repetition = 0;
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

void rotationToVtk(double R[3][3], double eulerAnglesDeg_xyz[3])
{
    /*
        Concert a rotation matrix into the Mayavi / Vtk rotation paramaters(pitch, roll, yaw)
    */
    double eulerAnglesRad_yxz[3];
    euler_from_matrix(R, eulerAnglesRad_yxz);

    double eulerAnglesDeg_yxz[3];
    eulerAnglesDeg_yxz[0] = eulerAnglesRad_yxz[0] * 180 / M_PI;
    eulerAnglesDeg_yxz[1] = eulerAnglesRad_yxz[1] * 180 / M_PI;
    eulerAnglesDeg_yxz[2] = eulerAnglesRad_yxz[2] * 180 / M_PI;

    eulerAnglesDeg_xyz[0] = eulerAnglesDeg_yxz[1];
    eulerAnglesDeg_xyz[1] = eulerAnglesDeg_yxz[0];
    eulerAnglesDeg_xyz[2] = eulerAnglesDeg_yxz[2];

    return;
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::updateAngleValuesFromTransform(vtkMRMLTransformNode* transformNode)
{
  Q_D(qMRMLTransformSliders);
  vtkNew<vtkTransform> transform;
  qMRMLUtils::getTransformInCoordinateSystem(transformNode,
      this->coordinateReference() == qMRMLTransformSliders::GLOBAL, transform.GetPointer());

  vtkMatrix4x4 * matrix = transform->GetMatrix();
  Q_ASSERT(matrix);
  if (!matrix) { return; }

  double matrix3x3[3][3];
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 3; col++)
      matrix3x3[col][row] = matrix->GetElement(row,col);//Transpose since it's a rotation

  double angles_deg_xyz[3];
  rotationToVtk(matrix3x3, angles_deg_xyz);

  bool blocked = d->LRSlider->blockSignals(true);
  d->LRSlider->setValue(angles_deg_xyz[0]);
  d->LRSlider->blockSignals(blocked);

  blocked = d->PASlider->blockSignals(true);
  d->PASlider->setValue(angles_deg_xyz[1]);
  d->PASlider->blockSignals(blocked);

  blocked = d->ISSlider->blockSignals(true);
  d->ISSlider->setValue(angles_deg_xyz[2]);
  d->ISSlider->blockSignals(blocked);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::updateRangeFromTransform(vtkMRMLTransformNode* transformNode)
{
    Q_D(qMRMLTransformSliders);
    vtkNew<vtkTransform> transform;
    qMRMLUtils::getTransformInCoordinateSystem(transformNode,
        this->coordinateReference() == qMRMLTransformSliders::GLOBAL, transform.GetPointer());

    vtkMatrix4x4* matrix = transform->GetMatrix();
    Q_ASSERT(matrix);
    if (!matrix) { return; }

    if (this->typeOfTransform() == qMRMLTransformSliders::TRANSLATION)
    {
        QPair<double, double> minmax = this->extractMinMaxTranslationValue(matrix, 0.0);
        if (minmax.first < this->minimum())
        {
            minmax.first = minmax.first - 0.3 * fabs(minmax.first);
            this->setMinimum(minmax.first);
        }
        if (minmax.second > this->maximum())
        {
            minmax.second = minmax.second + 0.3 * fabs(minmax.second);
            this->setMaximum(minmax.second);
        }
    }
}

// --------------------------------------------------------------------------
CTK_GET_CPP(qMRMLTransformSliders, vtkMRMLTransformNode*, mrmlTransformNode, MRMLTransformNode);

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setTitle(const QString& _title)
{
  Q_D(qMRMLTransformSliders);
  d->SlidersGroupBox->setTitle(_title);
}

// --------------------------------------------------------------------------
QString qMRMLTransformSliders::title()const
{
  Q_D(const qMRMLTransformSliders);
  return d->SlidersGroupBox->title();
}

// --------------------------------------------------------------------------
int qMRMLTransformSliders::decimals()const
{
  Q_D(const qMRMLTransformSliders);
  return d->LRSlider->decimals();
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setDecimals(int newDecimals)
{
  Q_D(qMRMLTransformSliders);
  // setting the decimals to LRSlider will propagate to the other widgets.
  d->LRSlider->setDecimals(newDecimals);
}

// --------------------------------------------------------------------------
double qMRMLTransformSliders::minimum()const
{
  Q_D(const qMRMLTransformSliders);
  return d->MinValueSpinBox->value();
}

// --------------------------------------------------------------------------
double qMRMLTransformSliders::maximum()const
{
  Q_D(const qMRMLTransformSliders);
  return d->MaxValueSpinBox->value();
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setMinimum(double min)
{
  Q_D(qMRMLTransformSliders);
  d->MinValueSpinBox->setValue(min);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setMaximum(double max)
{
  Q_D(qMRMLTransformSliders);
  d->MaxValueSpinBox->setValue(max);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setRange(double min, double max)
{
  Q_D(qMRMLTransformSliders);

  // Could be optimized here by blocking signals on spinboxes and manually
  // call the setRange method on the sliders. Does it really worth it ?
  d->MinValueSpinBox->setValue(min);
  d->MaxValueSpinBox->setValue(max);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::onMinimumChanged(double min)
{
  Q_D(qMRMLTransformSliders);

  d->LRSlider->setMinimum(min);
  d->PASlider->setMinimum(min);
  d->ISSlider->setMinimum(min);

  emit this->rangeChanged(min, this->maximum());
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::onMaximumChanged(double max)
{
  Q_D(qMRMLTransformSliders);

  d->LRSlider->setMaximum(max);
  d->PASlider->setMaximum(max);
  d->ISSlider->setMaximum(max);

  emit this->rangeChanged(this->minimum(), max);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setMinMaxVisible(bool visible)
{
  Q_D(qMRMLTransformSliders);
  d->MinMaxWidget->setVisible(visible);
}

// --------------------------------------------------------------------------
bool qMRMLTransformSliders::isMinMaxVisible()const
{
  Q_D(const qMRMLTransformSliders);
  return d->MinMaxWidget->isVisibleTo(
    const_cast<qMRMLTransformSliders*>(this));
}

// --------------------------------------------------------------------------
double qMRMLTransformSliders::singleStep()const
{
  Q_D(const qMRMLTransformSliders);
  // Assumes settings of the sliders are all the same
  return d->PASlider->singleStep();
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setSingleStep(double step)
{
  Q_D(qMRMLTransformSliders);

  d->LRSlider->setSingleStep(step);
  d->PASlider->setSingleStep(step);
  d->ISSlider->setSingleStep(step);
}

// --------------------------------------------------------------------------
QString qMRMLTransformSliders::lrLabel()const
{
  Q_D(const qMRMLTransformSliders);
  return d->LRLabel->text();
}

// --------------------------------------------------------------------------
QString qMRMLTransformSliders::paLabel()const
{
  Q_D(const qMRMLTransformSliders);
  return d->PALabel->text();
}

// --------------------------------------------------------------------------
QString qMRMLTransformSliders::isLabel()const
{
  Q_D(const qMRMLTransformSliders);
  return d->ISLabel->text();
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setLRLabel(const QString& label)
{
  Q_D(qMRMLTransformSliders);
  d->LRLabel->setText(label);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setPALabel(const QString& label)
{
  Q_D(qMRMLTransformSliders);
  d->PALabel->setText(label);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::setISLabel(const QString& label)
{
  Q_D(qMRMLTransformSliders);
  d->ISLabel->setText(label);
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::reset()
{
  Q_D(qMRMLTransformSliders);

  d->LRSlider->reset();
  d->PASlider->reset();
  d->ISSlider->reset();
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::resetUnactiveSliders()
{
  Q_D(qMRMLTransformSliders);

  if (!d->ActiveSliders.contains(d->LRSlider))
    {
    bool blocked = d->LRSlider->blockSignals(true);
    d->LRSlider->reset();
    d->LRSlider->blockSignals(blocked);
    }
  if (!d->ActiveSliders.contains(d->PASlider))
    {
    bool blocked = d->PASlider->blockSignals(true);
    d->PASlider->reset();
    d->PASlider->blockSignals(blocked);
    }
  if (!d->ActiveSliders.contains(d->ISSlider))
    {
    bool blocked = d->ISSlider->blockSignals(true);
    d->ISSlider->reset();
    d->ISSlider->blockSignals(blocked);
    }
}

// --------------------------------------------------------------------------
void qMRMLTransformSliders::onSliderPositionChanged(double position)
{
  Q_D(qMRMLTransformSliders);
  qMRMLLinearTransformSlider* slider =
    qobject_cast<qMRMLLinearTransformSlider*>(this->sender());
  Q_ASSERT(slider);
  d->ActiveSliders.push(slider);
  QWidget* focusWidget = this->focusWidget();

  // If update initiated from spinbox, consider it active, too
  // (when number of decimals are updated then it may change all the sliders
  // one by one, but that should not reset the axis that is currently being changed)
  if (focusWidget)
    {
    if (focusWidget->parent() == d->LRSlider->spinBox())
      {
      d->ActiveSliders.push(d->LRSlider);
      }
    if (focusWidget->parent() == d->PASlider->spinBox())
      {
      d->ActiveSliders.push(d->PASlider);
      }
    if (focusWidget->parent() == d->ISSlider->spinBox())
      {
      d->ActiveSliders.push(d->ISSlider);
      }
    }

  if (this->typeOfTransform() == qMRMLTransformSliders::ROTATION
    || (this->typeOfTransform() == qMRMLTransformSliders::TRANSLATION && coordinateReference() == LOCAL) )
    {
    // When a rotation slider is manipulated, the other rotation sliders are
    // reset to 0. Resetting the other sliders should no fire any event.
    this->resetUnactiveSliders();
    }
  slider->applyTransformation(position);
  emit this->valuesChanged();

  d->ActiveSliders.pop();
}

//-----------------------------------------------------------------------------
QPair<double, double> qMRMLTransformSliders::extractMinMaxTranslationValue(
                                             vtkMatrix4x4 * mat, double pad)
{
  QPair<double, double> minmax;
  if (!mat)
    {
    Q_ASSERT(mat);
    return minmax;
    }
  for (int i=0; i <3; i++)
    {
    minmax.first = qMin(minmax.first, mat->GetElement(i,3));
    minmax.second = qMax(minmax.second, mat->GetElement(i,3));
    }
  double range = minmax.second - minmax.first;
  minmax.first = minmax.first - pad * range;
  minmax.second = minmax.second + pad * range;
  return minmax;
}

