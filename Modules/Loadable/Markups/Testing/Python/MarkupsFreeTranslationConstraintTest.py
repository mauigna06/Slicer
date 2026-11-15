import numpy as np
import vtk

import slicer
from slicer.ScriptedLoadableModule import ScriptedLoadableModuleTest


class MarkupsFreeTranslationConstraintTest(ScriptedLoadableModuleTest):
    """Test that vtkMRMLMarkupsDisplayNode::FreeTranslationConstrainedToHandleAxes restricts the translation of the
    whole markup by middle-button drag, and by shift + left-button drag of a plane, to the enabled translation handle axes.

    The mouse events are sent to the 3D view, so the test goes through the same code path as a user.
    The plane normal is oblique to the view plane, so an unconstrained drag moves the plane along its normal.
    """

    def setUp(self):
        slicer.mrmlScene.Clear(0)
        layoutManager = slicer.app.layoutManager()
        layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUp3DView)
        threeDWidget = layoutManager.threeDWidget(0)
        self.threeDView = threeDWidget.threeDView()

        # Look at the origin along the -Z axis, so the view plane is the XY plane
        self.cameraNode = slicer.modules.cameras.logic().GetViewActiveCameraNode(threeDWidget.mrmlViewNode())
        self.cameraNode.SetPosition(0.0, 0.0, 500.0)
        self.cameraNode.SetFocalPoint(0.0, 0.0, 0.0)
        self.cameraNode.SetViewUp(0.0, 1.0, 0.0)

        # Plane with its normal at 45 degrees from the view direction
        self.planeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsPlaneNode")
        self.planeNode.SetSizeMode(slicer.vtkMRMLMarkupsPlaneNode.SizeModeAbsolute)
        self.planeNode.SetSize(80.0, 80.0)
        planeX = np.array([1.0, 0.0, -1.0]) / np.sqrt(2.0)
        planeY = np.array([0.0, 1.0, 0.0])
        planeZ = np.cross(planeX, planeY)
        self.planeNode.SetAxesWorld(planeX.tolist(), planeY.tolist(), planeZ.tolist())
        self.planeNode.SetOriginWorld(0.0, 0.0, 0.0)

        # Only translation along the X and Y axes of the plane is enabled
        self.displayNode = self.planeNode.GetDisplayNode()
        self.displayNode.SetTranslationHandleComponentVisibility(True, True, False, False)

        self.renderView()

    def tearDown(self):
        slicer.mrmlScene.Clear(0)

    def renderView(self):
        slicer.app.processEvents()
        self.threeDView.forceRender()

    def getDisplayPositionOfWorldPoint(self, point_World):
        renderer = self.threeDView.renderWindow().GetRenderers().GetFirstRenderer()
        renderer.SetWorldPoint(point_World[0], point_World[1], point_World[2], 1.0)
        renderer.WorldToDisplay()
        displayPoint = renderer.GetDisplayPoint()
        return round(displayPoint[0]), round(displayPoint[1])

    def getPlaneOutlinePoint(self):
        # Point at a quarter of an edge of the plane (the middle of the edge is a scale handle)
        cornerPoints_World = vtk.vtkPoints()
        self.planeNode.GetPlaneCornerPointsWorld(cornerPoints_World)
        return 0.75 * np.array(cornerPoints_World.GetPoint(0)) + 0.25 * np.array(cornerPoints_World.GetPoint(1))

    def dragFromWorldPoint(self, point_World, expectedComponentType, button, modifiers=[]):
        """Drag the mouse from the point, after hovering over it, and return the translation of the plane origin
        along the axes of the plane before the drag. The hovered component of the plane must be the expected one.
        """
        planeX, planeY, planeZ = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
        self.planeNode.GetAxesWorld(planeX, planeY, planeZ)
        originBefore = np.array(self.planeNode.GetOriginWorld())

        start = self.getDisplayPositionOfWorldPoint(point_World)
        end = (start[0] + 60, start[1] + 30)
        interactor = self.threeDView.interactorStyle().GetInteractor()
        interactor.SetEventPosition(*start)
        interactor.MouseMoveEvent()
        self.renderView()
        self.assertEqual(self.displayNode.GetActiveComponentType(), expectedComponentType, msg=f"hovered component at display position {start}")
        slicer.util.clickAndDrag(self.threeDView, button=button, start=start, end=end, steps=10, modifiers=modifiers)
        self.renderView()

        translation_World = np.array(self.planeNode.GetOriginWorld()) - originBefore
        return np.array([np.dot(translation_World, axis) for axis in [planeX, planeY, planeZ]])

    def test_MiddleButtonDragIsNotConstrainedByDefault(self):
        self.assertFalse(self.displayNode.GetFreeTranslationConstrainedToHandleAxes())
        translation = self.dragFromWorldPoint(self.planeNode.GetOriginWorld(), slicer.vtkMRMLMarkupsDisplayNode.ComponentControlPoint, "Middle")
        # The drag translates the plane along the view plane, so also along the plane normal
        self.assertGreater(abs(translation[2]), 1.0, msg=f"translation along the plane axes {translation}")

    def test_MiddleButtonDragIsConstrainedToEnabledHandleAxes(self):
        self.displayNode.SetFreeTranslationConstrainedToHandleAxes(True)
        translation = self.dragFromWorldPoint(self.planeNode.GetOriginWorld(), slicer.vtkMRMLMarkupsDisplayNode.ComponentControlPoint, "Middle")
        self.assertGreater(np.linalg.norm(translation[:2]), 1.0, msg=f"translation along the plane axes {translation}")
        self.assertAlmostEqual(translation[2], 0.0, delta=1e-6, msg=f"translation along the plane axes {translation}")

    def test_ShiftLeftButtonDragOfPlaneIsConstrainedToEnabledHandleAxes(self):
        self.displayNode.SetFreeTranslationConstrainedToHandleAxes(True)
        translation = self.dragFromWorldPoint(self.getPlaneOutlinePoint(), slicer.vtkMRMLMarkupsDisplayNode.ComponentPlane, "Left", ["Shift"])
        self.assertGreater(np.linalg.norm(translation[:2]), 1.0, msg=f"translation along the plane axes {translation}")
        self.assertAlmostEqual(translation[2], 0.0, delta=1e-6, msg=f"translation along the plane axes {translation}")

    def test_ViewPlaneComponentKeepsFreeTranslation(self):
        self.displayNode.SetFreeTranslationConstrainedToHandleAxes(True)
        self.displayNode.SetTranslationHandleComponentVisibility(False, False, False, True)
        translation = self.dragFromWorldPoint(self.planeNode.GetOriginWorld(), slicer.vtkMRMLMarkupsDisplayNode.ComponentControlPoint, "Middle")
        self.assertGreater(abs(translation[2]), 1.0, msg=f"translation along the plane axes {translation}")

    def test_FreeTranslationIsDisabledWithoutEnabledComponents(self):
        self.displayNode.SetFreeTranslationConstrainedToHandleAxes(True)
        self.displayNode.SetTranslationHandleComponentVisibility(False, False, False, False)
        focalPointBefore = np.array(self.cameraNode.GetFocalPoint())

        translation = self.dragFromWorldPoint(self.planeNode.GetOriginWorld(), slicer.vtkMRMLMarkupsDisplayNode.ComponentControlPoint, "Middle")
        self.assertTrue(np.allclose(translation, 0.0, atol=1e-9), msg=f"translation along the plane axes {translation}")
        # The markup does not take the event, so the middle-button drag pans the view
        focalPoint = np.array(self.cameraNode.GetFocalPoint())
        self.assertGreater(np.linalg.norm(focalPoint - focalPointBefore), 1.0, msg=f"camera focal point {focalPoint} before {focalPointBefore}")

        translation = self.dragFromWorldPoint(self.getPlaneOutlinePoint(), slicer.vtkMRMLMarkupsDisplayNode.ComponentPlane, "Left", ["Shift"])
        self.assertTrue(np.allclose(translation, 0.0, atol=1e-9), msg=f"translation along the plane axes {translation}")

    def test_FlagIsSavedInScene(self):
        self.displayNode.SetFreeTranslationConstrainedToHandleAxes(True)
        slicer.mrmlScene.SetSaveToXMLString(1)
        slicer.mrmlScene.Commit()
        sceneXML = slicer.mrmlScene.GetSceneXMLString()
        self.assertIn('freeTranslationConstrainedToHandleAxes="true"', sceneXML)

        loadedScene = slicer.vtkMRMLScene()
        loadedScene.SetLoadFromXMLString(1)
        loadedScene.SetSceneXMLString(sceneXML)
        loadedScene.Import()
        loadedDisplayNode = loadedScene.GetNodeByID(self.displayNode.GetID())
        self.assertIsNotNone(loadedDisplayNode)
        self.assertTrue(loadedDisplayNode.GetFreeTranslationConstrainedToHandleAxes())
