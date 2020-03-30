import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import gc
import scipy.optimize
import numpy
import time
import csv
from SlicerDevelopmentToolboxUtils.mixins import ModuleWidgetMixin, ModuleLogicMixin
#
# CryoMC
#

SizeX = 6.0
SizeY = 6.0
SizeZ = 12.0

#SizeX = 10.5
#SizeY = 10.5
#SizeZ = 14.0

class CryoMC(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "CryoMC" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# CryoMCWidget
#

#class CryoMCWidget(ScriptedLoadableModuleWidget):


class CryoMCWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  inputSelector2 = None  # type: object

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputT2image = slicer.qMRMLNodeComboBox()
    self.inputT2image.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputT2image.selectNodeUponCreation = True
    self.inputT2image.addEnabled = False
    self.inputT2image.removeEnabled = False
    self.inputT2image.noneEnabled = False
    self.inputT2image.showHidden = False
    self.inputT2image.showChildNodeTypes = False
    self.inputT2image.setMRMLScene( slicer.mrmlScene )
    self.inputT2image.setToolTip( "Pick the input to the simulation algorithm." )
    parametersFormLayout.addRow("T2 image: ", self.inputT2image)

    self.inputTumor1 = slicer.qMRMLNodeComboBox()
    self.inputTumor1.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.inputTumor1.selectNodeUponCreation = True
    self.inputTumor1.addEnabled = False
    self.inputTumor1.removeEnabled = False
    self.inputTumor1.noneEnabled = False
    self.inputTumor1.showHidden = False
    self.inputTumor1.showChildNodeTypes = False
    self.inputTumor1.setMRMLScene( slicer.mrmlScene )
    self.inputTumor1.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Taget Volume: ", self.inputTumor1)

    self.inputTumor2 = slicer.qMRMLNodeComboBox()
    self.inputTumor2.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.inputTumor2.selectNodeUponCreation = True
    self.inputTumor2.addEnabled = False
    self.inputTumor2.removeEnabled = False
    self.inputTumor2.noneEnabled = False
    self.inputTumor2.showHidden = False
    self.inputTumor2.showChildNodeTypes = False
    self.inputTumor2.setMRMLScene( slicer.mrmlScene )
    self.inputTumor2.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Segmented Urethra: ", self.inputTumor2)


    self.inputTumor3 = slicer.qMRMLNodeComboBox()
    self.inputTumor3.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.inputTumor3.selectNodeUponCreation = True
    self.inputTumor3.addEnabled = False
    self.inputTumor3.removeEnabled = False
    self.inputTumor3.noneEnabled = False
    self.inputTumor3.showHidden = False
    self.inputTumor3.showChildNodeTypes = False
    self.inputTumor3.setMRMLScene( slicer.mrmlScene )
    self.inputTumor3.setToolTip( "Pick the label input to the algorithm." )
    parametersFormLayout.addRow("Obstacle #2: ", self.inputTumor3)

    self.SegmentedUrethra = qt.QCheckBox()
    self.SegmentedUrethra.setChecked(True)
    parametersFormLayout.addRow("Segmented urethra:", self.SegmentedUrethra)

    #
    # threshold value
    #
    self.imageThresholdSliderWidget = ctk.ctkSliderWidget()
    self.imageThresholdSliderWidget.singleStep = 1
    self.imageThresholdSliderWidget.minimum = 0
    self.imageThresholdSliderWidget.maximum = 5
    self.imageThresholdSliderWidget.value = 2
    self.imageThresholdSliderWidget.setToolTip("Set the number of probes")
    parametersFormLayout.addRow("Number of probes", self.imageThresholdSliderWidget)

    self.kdSliderWidget = ctk.ctkSliderWidget()
    self.kdSliderWidget.singleStep = 0.1
    self.kdSliderWidget.minimum = 0
    self.kdSliderWidget.maximum = 10
    self.kdSliderWidget.value = 5
    self.kdSliderWidget.setToolTip("Set value for Kd")
    parametersFormLayout.addRow("Value for Kd", self.kdSliderWidget)

    self.NSliderWidget = ctk.ctkSliderWidget()
    self.NSliderWidget.singleStep = 1
    self.NSliderWidget.minimum = 0
    self.NSliderWidget.maximum = 1000
    self.NSliderWidget.value = 100
    self.NSliderWidget.setToolTip("Set the value of N")
    parametersFormLayout.addRow("Value for N up", self.NSliderWidget)


    self.imageStdev = ctk.ctkSliderWidget()
    self.imageStdev.singleStep = 0.1
    self.imageStdev.minimum = 0
    self.imageStdev.maximum = 5
    self.imageStdev.value = 1.5
    self.imageStdev.setToolTip("Set placement error stdev")
    parametersFormLayout.addRow("placement error", self.imageStdev)


    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)

    outputWindow = ctk.ctkCollapsibleButton()
    outputWindow.text = 'Output'
    self.layout.addWidget(outputWindow)
    self.outputLayout = qt.QFormLayout(outputWindow)

    # Initial output text
    self.outputLabel = qt.QLabel("")
    self.outputLayout.addRow(self.outputLabel)
    self.outputLabel.setText("Results")


    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputT2image.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputT2image.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputT2image.currentNode()

  def onApplyButton(self):
    logic = CryoMCLogic()
    imageThreshold = self.imageThresholdSliderWidget.value
    logic.run(self.inputT2image.currentNode(),self.inputTumor1.currentNode(), self.inputTumor2.currentNode(), self.inputTumor3.currentNode(), self.imageStdev.value, self.imageThresholdSliderWidget.value, self.kdSliderWidget.value, self.NSliderWidget.value, self.outputLabel)


class CryoMCLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """


  def run(self, inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel):
    """
    Run the actual algorithm

    inputVolume2 is the label
    """
    self.outputLabel = outputLabel

    if NoOfProbes == 1:
      self.OneProbes = Cryo1Probes()
      self.OneProbes.TestFunction(inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel)
      return
    elif NoOfProbes == 2:
      self.TwoProbes = Cryo2Probes()
      self.TwoProbes.TestFunction(inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel)
      return
    elif NoOfProbes == 3:
      self.ThreeProbes = Cryo3Probes()
      self.ThreeProbes.TestFunction(inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel)
      return
    elif NoOfProbes == 4:
      self.FourProbes = Cryo4Probes()
      self.FourProbes.TestFunction(inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel)
      return
    else:
      self.outputLabel.setText("Invalid number of Probes")
      return

    return

class Cryo1Probes(ScriptedLoadableModuleLogic):

  def ComputeDice1IceBall(self):

    self.Spheres2 = self.sphere1.GetOutput()

    self.b = self.Spheres2.GetPointData()
    PixelImage2 = numpy.count_nonzero(self.b.GetArray(0))

    self.logicFilter.SetOperationToAnd()
    self.logicFilter.SetInputData(0, self.imData1) #check here.
    self.logicFilter.SetInputData(1, self.Spheres2)
    self.logicFilter.Update()
    self.OutPut = self.logicFilter.GetOutput()

    self.b = self.OutPut.GetPointData()
    intersection = numpy.count_nonzero(self.b.GetArray(0))

    self.b = self.imData1.GetPointData()
    PixelImage1 = numpy.count_nonzero(self.b.GetArray(0))

    return [(2.0 * intersection) / (PixelImage1 + PixelImage2), (1.0*intersection / PixelImage1)]

  def optimizationCryo(self, N, Kd):
    ############
    def Fun(x):
        Metric = self.Probability(x, N)
        temp = Metric[0] / 100.0
        return -temp * temp * Metric[1]
    ############
    x0 = [self.PositionIceBall1[0], self.PositionIceBall1[1], self.PositionIceBall1[2]]
    res = scipy.optimize.minimize(Fun, x0, method='Nelder-Mead', options={'fatol': 0.02, 'xatol': 0.1, 'maxiter': 100, 'maxfev': 100})
    print("==opt output===")
    print(res.x)
    print("=======")
    self.FinalProbePlacement1 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[0],res.x[1],res.x[2],1])
    self.Metric = self.Probability(res.x, 100)
    print(self.Metric)

  def Probability(self, X, No):
    NofInt = int(No)
    PC_vector = numpy.zeros(NofInt)
    DSC_vector = numpy.zeros(NofInt)
    Metric = numpy.zeros(NofInt)
    noiseGaussianX = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianY = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianZ = numpy.random.normal(0.0, self.sdterr/10.0, 1000)

    for i in range(0, NofInt - 1):
      self.sphere1.SetCenter(X[0] + noiseGaussianX[i], X[1] + noiseGaussianY[i], X[2] + noiseGaussianZ[i])
      self.sphere1.Update()
      DICE_temp = self.ComputeDice1IceBall()
      Metric[i] = DICE_temp[0] + self.DiceGain * DICE_temp[1]
      PC_vector[i] = DICE_temp[1]
      DSC_vector[i] = DICE_temp[0]
    count = 0
    for i in range(0, NofInt - 1):
      if PC_vector[i] >= 0.96:
        count = count + 1
    if count == 0:
      count = 1
    return [count, numpy.mean(Metric), numpy.mean(PC_vector), numpy.mean(DSC_vector)]

  def InitDataForPlanning(self,inputVolume2):
    self.sphere1 = vtk.vtkImageEllipsoidSource()
    self.sphere1.SetOutputScalarTypeToShort()
    self.PositionIceBall1 = [0, 0, 0, 1]
    self.PositionIceBall1_temp = self.PositionIceBall1
    self.Spacing = inputVolume2.GetSpacing()
    self.sphere1.SetRadius(SizeX/ self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    size_image = inputVolume2.GetImageData().GetDimensions()
    self.sphere1.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere1.Update()
    self.logicFilter = vtk.vtkImageLogic()
    self.logicFilter2 = vtk.vtkImageLogic()
    self.DiceGain = 5.0
    self.imData1 = inputVolume2.GetImageData()
    self.IjkToRasMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetIJKToRASMatrix(self.IjkToRasMatrix)


  def CreateModels(self, labelMapNode, modelHierarchyNode):

    modelMakerCLI = slicer.modules.modelmaker
    # tf = tempfile.NamedTemporaryFile(prefix='Slicer/Models-', suffix='.mrml')

    modelMakerParameters = {}
    # modelMakerParameters['ColorTable'] = 'vtkMRMLColorTableNodeFileGenericAnatomyColors.txt'
    modelMakerParameters['ModelSceneFile'] = modelHierarchyNode.GetID()
    modelMakerParameters['Name'] = 'Model'
    modelMakerParameters['GenerateAll'] = True
    modelMakerParameters['StartLabel'] = -1
    modelMakerParameters['EndLabel'] = -1
    modelMakerParameters['KkipUnNamed'] = True
    modelMakerParameters['JointSmooth'] = True
    modelMakerParameters['Smooth'] = 15
    modelMakerParameters['FilterType'] = 'Sinc'
    modelMakerParameters['Decimate'] = 0.25
    modelMakerParameters['SplitNormals'] = True
    modelMakerParameters['PointNormals'] = True
    modelMakerParameters['InputVolume'] = labelMapNode.GetID()

    slicer.cli.run(modelMakerCLI, None, modelMakerParameters, True)

  def GetCenterTumor(self, inputVolume2):
    modelHierarchyNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelHierarchyNode")
    slicer.mrmlScene.AddNode(modelHierarchyNode)
    self.CreateModels(inputVolume2, modelHierarchyNode)
    nOfModels = modelHierarchyNode.GetNumberOfChildrenNodes()
    if (nOfModels > 1):
      slicer.util.errorDisplay("More than one segmented ablation volume")
      return
    chnode = modelHierarchyNode.GetNthChildNode(0)
    mnode = chnode.GetAssociatedNode()
    objectPoly = mnode.GetPolyData()

    centerOfmass = vtk.vtkCenterOfMass()
    centerOfmass.SetInputData(objectPoly)
    centerOfmass.SetUseScalarsAsWeights(False)
    centerOfmass.Update()
    self.Center = centerOfmass.GetCenter()

    RasToIjkMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetRASToIJKMatrix(RasToIjkMatrix)
    b_Ijk = objectPoly.GetBounds()

    p_ijk = [self.Center[0], self.Center[1], self.Center[2], 1]
    p_ijk = RasToIjkMatrix.MultiplyDoublePoint(p_ijk)

    self.PositionIceBall1 = p_ijk
    return

  def GetIceMAtrixRAS(self,image,pos):

    Probe1 = pos

    IceMatrix = vtk.vtkMatrix4x4()
    IceMatrix.SetElement(0, 0, 1.0)
    IceMatrix.SetElement(0, 1, 0.0)
    IceMatrix.SetElement(0, 2, 0.0)
    IceMatrix.SetElement(0, 3, Probe1[0])

    IceMatrix.SetElement(1, 0, 0.0)
    IceMatrix.SetElement(1, 1, 1.0)
    IceMatrix.SetElement(1, 2, 0.0)
    IceMatrix.SetElement(1, 3, Probe1[1])

    IceMatrix.SetElement(2, 0, 0.0)
    IceMatrix.SetElement(2, 1, 0.0)
    IceMatrix.SetElement(2, 2, 1.0)
    IceMatrix.SetElement(2, 3, Probe1[2])

    IceMatrix.SetElement(3, 0, 0.0)
    IceMatrix.SetElement(3, 1, 0.0)
    IceMatrix.SetElement(3, 2, 0.0)
    IceMatrix.SetElement(3, 3, 1.0)

    return IceMatrix

  def TestFunction(self, inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel):
    savetime = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMetric = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveDSC = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveHTA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMass = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    self.sdterr = StdevError
    self.DiceGain = Kd

    resultFileName = "Results-N%02d-K%02d-P1-r42" % (N, self.DiceGain)
    resultFilePath = '/Users/pedro/Documents/Test-slicer' + '/' + resultFileName
    resultFile = open(resultFilePath, 'a')

    resultFile.write("Probe; N ; Kd ; Err ; PoS ; PTC ; DSC \n")

    for j in range(0, 13):
      for i in range(0, 5):
        self.sdterr = j
        start_time = time.time()
        self.InitDataForPlanning(inputVolume2)
        self.GetCenterTumor(inputVolume2)
        self.output = outputLabel

        self.optimizationCryo(N, Kd)
        elapsed_time = time.time() - start_time
        print(elapsed_time)
        savetime[i] = elapsed_time
        saveMetric[i] = self.Metric[0]
        saveDSC[i] = self.Metric[3]
        saveHTA[i] = self.Metric[2]
        saveMass[i] = 1
        resultFile.write("01 ; %02d ; %02d ; %02f ; %02d ; %02f ; %02f \n" % (
        N, self.DiceGain, self.sdterr, self.Metric[0], self.Metric[2], self.Metric[3]))


    print("RESULTS TIME:")
    print(numpy.mean(savetime))
    print(numpy.std(savetime))
    print("RESULTS PofS %d:" % N)
    print(numpy.mean(saveMetric))
    print(numpy.std(saveMetric))
    print("RESULTS DSC %d:" % N)
    print(numpy.mean(saveDSC))
    print(numpy.std(saveDSC))
    print(self.DiceGain)
    print(saveDSC)

    self.affectedAreaModelNode = ModuleLogicMixin.createModelNode("AffectedArea")
    print("gfg")

    IceMatrix1 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement1)

    affectedBallArea = vtk.vtkParametricEllipsoid()
    affectedBallArea.SetXRadius(SizeX)
    affectedBallArea.SetYRadius(SizeY)
    affectedBallArea.SetZRadius(SizeZ)
    affectedBallAreaSource = vtk.vtkParametricFunctionSource()
    affectedBallAreaSource.SetParametricFunction(affectedBallArea)
    affectedBallAreaSource.SetScalarModeToV()
    affectedBallAreaSource.Update()
    affectedBallAreaAppend = vtk.vtkAppendPolyData()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix1)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()



    ModuleLogicMixin.createAndObserveDisplayNode(self.affectedAreaModelNode,
                                                 displayNodeClass=slicer.vtkMRMLModelDisplayNode)
    self.affectedAreaModelNode.GetDisplayNode().SetOpacity(0.2)
    self.affectedAreaModelNode.GetDisplayNode().SetColor(1.0, 0.0, 0.0)
    self.affectedAreaModelNode.SetAndObservePolyData(affectedBallAreaAppend.GetOutput())
    ModuleLogicMixin.setNodeVisibility(self.affectedAreaModelNode, True)
    ModuleLogicMixin.setNodeSliceIntersectionVisibility(self.affectedAreaModelNode, True)

    self.output.setText("The probe locations are: \n Probe 1: (%0.1f,%0.1f,%0.1f) \n Probability of success: %0.1f\n"
                        " Average tumor volume covered: %0.1f" % (self.FinalProbePlacement1[0],self.FinalProbePlacement1[1],
                                                                  self.FinalProbePlacement1[2], self.Metric[0], self.Metric[2]))

class Cryo2Probes(ScriptedLoadableModuleLogic):

  def ComputeDice2IceBall(self):

    self.imData2 = self.sphere1.GetOutput()
    self.imData3 = self.sphere2.GetOutput()

    self.logicFilter2.SetOperationToOr()
    self.logicFilter2.SetInputData(0, self.imData2)
    self.logicFilter2.SetInputData(1, self.imData3)
    self.logicFilter2.Update()
    self.Spheres2 = self.logicFilter2.GetOutput()

    self.b = self.Spheres2.GetPointData()
    PixelImage2 = numpy.count_nonzero(self.b.GetArray(0))



    self.logicFilter.SetOperationToAnd()
    self.logicFilter.SetInputData(0, self.imData1) #check here.
    self.logicFilter.SetInputData(1, self.Spheres2)
    self.logicFilter.Update()
    self.OutPut = self.logicFilter.GetOutput()

    self.b = self.OutPut.GetPointData()
    intersection = numpy.count_nonzero(self.b.GetArray(0))

    self.b = self.imData1.GetPointData()
    PixelImage1 = numpy.count_nonzero(self.b.GetArray(0))
    return [(2.0 * intersection) / (PixelImage1 + PixelImage2), (1.0*(intersection) / PixelImage1), PixelImage2]# - PixelImage1]

  def optimizationCryo(self, N, Kd):
    self.DiceGain = Kd
    def Fun(x):
      self.sphere1.SetCenter(x[0], x[1], x[2])
      self.sphere1.Update()
      self.sphere2.SetCenter(x[3], x[4], x[5])
      self.sphere2.Update()
      DICE_temp = self.ComputeDice2IceBall()
      return -DICE_temp[1]

    #x0 = [initialP1_ijk[0]-2, initialP1_ijk[1]-2, initialP1_ijk[2]+2, initialP1_ijk[0]+5, initialP1_ijk[1], initialP1_ijk[2]]
    x0 = [self.PositionIceBall1[0], self.PositionIceBall1[1], self.PositionIceBall1[2], self.PositionIceBall2[0], self.PositionIceBall2[1], self.PositionIceBall2[2]]
    res = scipy.optimize.minimize(Fun, x0, method='Nelder-Mead',options={'fatol': 0.01, 'xatol': 0.02, 'maxiter': 1200, 'maxfev': 10900})#options={'ftol': 0.2, 'gtol': 0.007, 'maxiter': 1200, 'maxfev': 10900})
    #res = scipy.optimize.minimize(Fun, x0, method='Nelder-Mead', options={'fatol': 0.1, 'xatol': 0.7, 'maxiter': 1200, 'maxfev': 10900})

    print(x0)
    print("==opt output===")
    print(res.x)
    print(res.fun)
    print(res.nit)
    print(res.success)
    print("=======")

    self.FinalProbePlacement1 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[0],res.x[1],res.x[2],1])
    self.FinalProbePlacement2 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[3], res.x[4], res.x[5], 1])
    self.Metric = self.Probability(res.x, 1000)
    self.converge = res.success
    print(self.Metric)
    return res.x

  def Probability(self, X, No):
    NofInt = int(No)
    PC_vector = numpy.zeros(NofInt)
    DSC_vector = numpy.zeros(NofInt)
    HTA_vector = numpy.zeros(NofInt)
    Metric = numpy.zeros(NofInt)
    noiseGaussianX = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianY = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianZ = numpy.random.normal(0.0, self.sdterr/10.0, 10000)
    noiseGaussianX2 = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianY2 = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianZ2 = numpy.random.normal(0.0, self.sdterr/10.0, 10000)
    #TODO

    if X[0] <= 0 or X[3] <=0:
      return [0, 0.2, 0, 0, 0]
    for i in range(0, NofInt - 1):
      self.sphere1.SetCenter(X[0] + noiseGaussianX[i], X[1] + noiseGaussianY[i], X[2] + noiseGaussianZ[i])
      self.sphere1.Update()
      self.sphere2.SetCenter(X[3] + noiseGaussianX2[i], X[4] + noiseGaussianY2[i], X[5] + noiseGaussianZ2[i])
      self.sphere2.Update()
      DICE_temp = self.ComputeDice2IceBall()
      Metric[i] = DICE_temp[0] + self.DiceGain * DICE_temp[1] #Changed here: 0 -> 1; 0 eh o DSC ; 1 eh o PTC
      PC_vector[i] = DICE_temp[1]
      DSC_vector[i] = DICE_temp[0]
      HTA_vector[i] = DICE_temp[2]
    count = 0
    for i in range(0, NofInt - 1):
      if PC_vector[i] >= 0.99:
        count = count + 1
    if count == 0:
        count = 10
    return [count, numpy.mean(Metric), numpy.mean(PC_vector), numpy.mean(DSC_vector), numpy.mean(HTA_vector)]

  def InitDataForPlanning(self,inputVolume2):
    self.sphere1 = vtk.vtkImageEllipsoidSource()
    self.sphere2 = vtk.vtkImageEllipsoidSource()
    self.sphere1.SetOutputScalarTypeToShort()
    self.sphere2.SetOutputScalarTypeToShort()
    self.PositionIceBall1 = [0, 0, 0, 1]
    self.PositionIceBall2 = [0, 0, 0, 1]
    self.PositionIceBall1_temp = self.PositionIceBall1
    self.PositionIceBall2_temp = self.PositionIceBall2
    self.Spacing = inputVolume2.GetSpacing()
    self.sphere1.SetRadius(SizeX / self.Spacing[0], SizeY/ self.Spacing[1], SizeZ/ self.Spacing[2])
    self.sphere2.SetRadius(SizeX / self.Spacing[0], SizeY/ self.Spacing[1], SizeZ/ self.Spacing[2])
    size_image = inputVolume2.GetImageData().GetDimensions()
    self.sphere1.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere2.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere1.Update()
    self.sphere2.Update()
    self.logicFilter = vtk.vtkImageLogic()
    self.logicFilter2 = vtk.vtkImageLogic()
    self.DiceGain = 5.0
    self.imData1 = inputVolume2.GetImageData()
    self.IjkToRasMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetIJKToRASMatrix(self.IjkToRasMatrix)


  def CreateModels(self, labelMapNode, modelHierarchyNode):

    modelMakerCLI = slicer.modules.modelmaker
    # tf = tempfile.NamedTemporaryFile(prefix='Slicer/Models-', suffix='.mrml')

    modelMakerParameters = {}
    # modelMakerParameters['ColorTable'] = 'vtkMRMLColorTableNodeFileGenericAnatomyColors.txt'
    modelMakerParameters['ModelSceneFile'] = modelHierarchyNode.GetID()
    modelMakerParameters['Name'] = 'Model'
    modelMakerParameters['GenerateAll'] = True
    modelMakerParameters['StartLabel'] = -1
    modelMakerParameters['EndLabel'] = -1
    modelMakerParameters['KkipUnNamed'] = True
    modelMakerParameters['JointSmooth'] = True
    modelMakerParameters['Smooth'] = 15
    modelMakerParameters['FilterType'] = 'Sinc'
    modelMakerParameters['Decimate'] = 0.25
    modelMakerParameters['SplitNormals'] = True
    modelMakerParameters['PointNormals'] = True
    modelMakerParameters['InputVolume'] = labelMapNode.GetID()

    slicer.cli.run(modelMakerCLI, None, modelMakerParameters, True)

  def GetCenterTumor(self, inputVolume2):
    modelHierarchyNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelHierarchyNode")
    slicer.mrmlScene.AddNode(modelHierarchyNode)
    self.CreateModels(inputVolume2, modelHierarchyNode)
    nOfModels = modelHierarchyNode.GetNumberOfChildrenNodes()
    if (nOfModels > 1):
      slicer.util.errorDisplay("More than one segmented ablation volume")
      return
    chnode = modelHierarchyNode.GetNthChildNode(0)
    mnode = chnode.GetAssociatedNode()
    objectPoly = mnode.GetPolyData()

    centerOfmass = vtk.vtkCenterOfMass()
    centerOfmass.SetInputData(objectPoly)
    centerOfmass.SetUseScalarsAsWeights(False)
    centerOfmass.Update()
    self.Center = centerOfmass.GetCenter()

    RasToIjkMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetRASToIJKMatrix(RasToIjkMatrix)
    b_Ijk = objectPoly.GetBounds()

    # Fisrt option:
    p_Ijka = [0, 0, 0, 1]
    p_Ijk2a = [0, 0, 0, 1]
    p_Ijka[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijka[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijka[2] = self.Center[2]
    p_Ijka[3] = 1
    p_Ijk2a[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2a[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2a[2] = self.Center[2]
    p_Ijk2a[3] = 1
    p_Ijka = RasToIjkMatrix.MultiplyDoublePoint(p_Ijka)
    p_Ijk2a = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2a)
    self.sphere1.SetCenter(p_Ijka[0], p_Ijka[1], p_Ijka[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2a[0], p_Ijk2a[1], p_Ijk2a[2])
    self.sphere2.Update()
    DICE1 = self.ComputeDice2IceBall()

    # Second option:
    p_Ijkb = [0, 0, 0, 1]
    p_Ijk2b = [0, 0, 0, 1]
    p_Ijkb[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijkb[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijkb[2] = self.Center[2]
    p_Ijkb[3] = 1
    p_Ijk2b[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2b[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2b[2] = self.Center[2]
    p_Ijk2b[3] = 1
    p_Ijkb = RasToIjkMatrix.MultiplyDoublePoint(p_Ijkb)
    p_Ijk2b = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2b)
    self.sphere1.SetCenter(p_Ijkb[0], p_Ijkb[1], p_Ijkb[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2b[0], p_Ijk2b[1], p_Ijk2b[2])
    self.sphere2.Update()
    DICE2 = self.ComputeDice2IceBall()

    # Third option:
    p_Ijkc = [0, 0, 0, 1]
    p_Ijk2c = [0, 0, 0, 1]
    p_Ijkc[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijkc[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijkc[2] = self.Center[2]
    p_Ijkc[3] = 1
    p_Ijk2c[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2c[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2c[2] = self.Center[2]
    p_Ijk2c[3] = 1
    p_Ijkc = RasToIjkMatrix.MultiplyDoublePoint(p_Ijkc)
    p_Ijk2c = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2c)
    self.sphere1.SetCenter(p_Ijkc[0], p_Ijkc[1], p_Ijkc[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2c[0], p_Ijk2c[1], p_Ijk2c[2])
    self.sphere2.Update()
    DICE3 = self.ComputeDice2IceBall()

    # Fourth option:
    p_Ijkd = [0, 0, 0, 1]
    p_Ijk2d = [0, 0, 0, 1]
    p_Ijkd[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijkd[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijkd[2] = self.Center[2]
    p_Ijkd[3] = 1
    p_Ijk2d[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2d[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2d[2] = self.Center[2]
    p_Ijk2d[3] = 1
    p_Ijkd = RasToIjkMatrix.MultiplyDoublePoint(p_Ijkd)
    p_Ijk2d = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2d)
    self.sphere1.SetCenter(p_Ijkd[0], p_Ijkd[1], p_Ijkd[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2d[0], p_Ijk2d[1], p_Ijk2d[2])
    self.sphere2.Update()
    DICE4 = self.ComputeDice2IceBall()


    if (DICE1[0] > DICE2[0]) and (DICE1[0] > DICE3[0]) and (DICE1[0] > DICE4[0]):
      self.PositionIceBall1 = [p_Ijka[0], p_Ijka[1], p_Ijka[2], 1]
      self.PositionIceBall2 = [p_Ijk2a[0], p_Ijk2a[1], p_Ijk2a[2], 1]
      return
    elif (DICE2[0] > DICE1[0]) and (DICE2[0] > DICE3[0]) and (DICE2[0] > DICE4[0]):
      self.PositionIceBall1 = [p_Ijkb[0], p_Ijkb[1], p_Ijkb[2], 1]
      self.PositionIceBall2 = [p_Ijk2b[0], p_Ijk2b[1], p_Ijk2b[2], 1]
      return
    elif (DICE3[0] > DICE1[0]) and (DICE3[0] > DICE2[0]) and (DICE3[0] > DICE4[0]):
      self.PositionIceBall1 = [p_Ijkc[0], p_Ijkc[1], p_Ijkc[2], 1]
      self.PositionIceBall2 = [p_Ijk2c[0], p_Ijk2c[1], p_Ijk2c[2], 1]
      return
    elif (DICE4[0] > DICE1[0]) and (DICE4[0] > DICE2[0]) and (DICE4[0] > DICE3[0]):
      self.PositionIceBall1 = [p_Ijkd[0], p_Ijkd[1], p_Ijkd[2], 1]
      self.PositionIceBall2 = [p_Ijk2d[0], p_Ijk2d[1], p_Ijk2d[2], 1]
      return
    else:
      self.PositionIceBall1 = [p_Ijka[0], p_Ijka[1], p_Ijka[2], 1]
      self.PositionIceBall2 = [p_Ijk2a[0], p_Ijk2a[1], p_Ijk2a[2], 1]
      return


  def GetIceMAtrixRAS(self,image,pos):

    Probe1 = pos

    IceMatrix = vtk.vtkMatrix4x4()
    IceMatrix.SetElement(0, 0, 1.0)
    IceMatrix.SetElement(0, 1, 0.0)
    IceMatrix.SetElement(0, 2, 0.0)
    IceMatrix.SetElement(0, 3, Probe1[0])

    IceMatrix.SetElement(1, 0, 0.0)
    IceMatrix.SetElement(1, 1, 1.0)
    IceMatrix.SetElement(1, 2, 0.0)
    IceMatrix.SetElement(1, 3, Probe1[1])

    IceMatrix.SetElement(2, 0, 0.0)
    IceMatrix.SetElement(2, 1, 0.0)
    IceMatrix.SetElement(2, 2, 1.0)
    IceMatrix.SetElement(2, 3, Probe1[2])

    IceMatrix.SetElement(3, 0, 0.0)
    IceMatrix.SetElement(3, 1, 0.0)
    IceMatrix.SetElement(3, 2, 0.0)
    IceMatrix.SetElement(3, 3, 1.0)

    return IceMatrix

  def TestFunction(self, inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel):
    savetime = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMetric = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveDSC = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveHTA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMass = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    self.sdterr = StdevError
    self.DiceGain = 5.0

    cases = ["028","029", "059", "032", "035"]

    for i in range(0, 1):
      path = '/Users/pedro/Dropbox (Partners HealthCare)/DataForCryoPaper/ralp-%s' % (cases[i])
      #path = '/Users/pedro/Dropbox (Partners HealthCare)/DataForCryoPaper/Clinical%s' % (cases[i])
      imageName = 'TumorRed'
      labelName = 'TumorRed-label'
      modelImageFileName = '%s/%s.nrrd' % (path, imageName)
      modelLabelFileName = '%s/%s.nrrd' % (path, labelName)

      (r, modelImageNode) = slicer.util.loadVolume(modelImageFileName, {}, True)
      (r, modelLabelNode) = slicer.util.loadLabelVolume(modelLabelFileName, {}, True)

      print(path)
      inputVolume2 = modelLabelNode
      print(modelLabelNode.GetSpacing())



      resultFileName = "Results-2Probes-2020a-Case0%s" % (cases[i])
      resultFilePath = '/Users/pedro/Projects/MonteCarlo' + '/' + resultFileName
      resultFile = open(resultFilePath, 'a')

      resultFile.write("Probe; Err ; PoF ; PTC \n")


      start_time = time.time()
      self.InitDataForPlanning(inputVolume2)
      self.GetCenterTumor(inputVolume2)
      self.sdterr = 0
      X = self.optimizationCryo(N, Kd)

      for i in range(0, 1):
        self.sdterr = i
        self.output = outputLabel
        result = self.Probability(X, 10000)
        resultFile.write("2; %02d ; %02d ; %02f \n" %(i,result[0],result[2]))

      elapsed_time = time.time() - start_time
      print(elapsed_time)
      print(self.Spacing)






    print(self.FinalProbePlacement1)
    print(self.FinalProbePlacement2)
    print(self.Spacing)
    self.affectedAreaModelNode = ModuleLogicMixin.createModelNode("AffectedArea")




    # Test copy from longquan Begin of affectedBallArea

    IceMatrix1 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement1)
    IceMatrix2 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement2)

    affectedBallArea = vtk.vtkParametricEllipsoid()
    affectedBallArea.SetXRadius(SizeX)
    affectedBallArea.SetYRadius(SizeY)
    affectedBallArea.SetZRadius(SizeZ)
    affectedBallAreaSource = vtk.vtkParametricFunctionSource()
    affectedBallAreaSource.SetParametricFunction(affectedBallArea)
    affectedBallAreaSource.SetScalarModeToV()
    affectedBallAreaSource.Update()
    affectedBallAreaAppend = vtk.vtkAppendPolyData()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix1)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix2)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()


    ModuleLogicMixin.createAndObserveDisplayNode(self.affectedAreaModelNode,
                                                 displayNodeClass=slicer.vtkMRMLModelDisplayNode)
    self.affectedAreaModelNode.GetDisplayNode().SetOpacity(0.2)
    self.affectedAreaModelNode.GetDisplayNode().SetColor(1.0, 0.0, 0.0)
    self.affectedAreaModelNode.SetAndObservePolyData(affectedBallAreaAppend.GetOutput())
    ModuleLogicMixin.setNodeVisibility(self.affectedAreaModelNode, True)
    ModuleLogicMixin.setNodeSliceIntersectionVisibility(self.affectedAreaModelNode, True)

    Distance = numpy.sqrt((IceMatrix1.GetElement(0, 3)-IceMatrix2.GetElement(0, 3))**2.0+(IceMatrix1.GetElement(1, 3)-IceMatrix2.GetElement(1, 3))**2.0)
    self.output.setText("The distance between probes is: %0.1f mm \n The probe locations are: \n Probe 1: (%0.1f,%0.1f,%0.1f) \n Probe 2:"
                        " (%0.1f,%0.1f,%0.1f) \n Probability of success: %0.1f\n Average tumor volume covered: %0.1f" % (Distance, self.FinalProbePlacement1[0],self.FinalProbePlacement1[1],self.FinalProbePlacement1[2],
                                                                                                               self.FinalProbePlacement2[0], self.FinalProbePlacement2[1], self.FinalProbePlacement2[2], self.Metric[0], self.Metric[2]))


class Cryo3Probes(ScriptedLoadableModuleLogic):

#
  def ComputeDice3IceBall(self):

    self.imData2 = self.sphere1.GetOutput()
    self.imData3 = self.sphere2.GetOutput()
    self.imData4 = self.sphere3.GetOutput()

    self.logicFilter2.SetOperationToOr()
    self.logicFilter2.SetInputData(0, self.imData2)
    self.logicFilter2.SetInputData(1, self.imData3)
    self.logicFilter2.Update()
    self.Spheres2 = self.logicFilter2.GetOutput()

    self.logicFilter3.SetOperationToOr()
    self.logicFilter3.SetInputData(0, self.Spheres2)
    self.logicFilter3.SetInputData(1, self.imData4)
    self.logicFilter3.Update()
    self.Spheres3 = self.logicFilter3.GetOutput()

    self.b = self.Spheres3.GetPointData()
    PixelImage2 = numpy.count_nonzero(self.b.GetArray(0))

    self.logicFilter.SetOperationToAnd()
    self.logicFilter.SetInputData(0, self.imData1) #check here.
    self.logicFilter.SetInputData(1, self.Spheres3)
    self.logicFilter.Update()
    self.OutPut = self.logicFilter.GetOutput()

    self.b = self.OutPut.GetPointData()
    intersection = numpy.count_nonzero(self.b.GetArray(0))

    self.b = self.imData1.GetPointData()
    PixelImage1 = numpy.count_nonzero(self.b.GetArray(0))

    return [(2.0 * intersection) / (PixelImage1 + PixelImage2), (1.0*intersection / PixelImage1), PixelImage2]# - PixelImage1]]

  def optimizationCryo(self, N, Kd):
    self.DiceGain = Kd
    def Fun(x):
      self.sphere1.SetCenter(x[0], x[1], x[2])
      self.sphere1.Update()
      self.sphere2.SetCenter(x[3], x[4], x[5])
      self.sphere2.Update()
      self.sphere3.SetCenter(x[6], x[7], x[8])
      self.sphere3.Update()
      DICE_temp = self.ComputeDice3IceBall()
      return -DICE_temp[1]

    x0 = [self.PositionIceBall1[0], self.PositionIceBall1[1], self.PositionIceBall1[2],
          self.PositionIceBall2[0], self.PositionIceBall2[1], self.PositionIceBall2[2],
         self.PositionIceBall3[0], self.PositionIceBall3[1], self.PositionIceBall3[2]]


    res = scipy.optimize.minimize(Fun, x0, method='Nelder-Mead', options={'fatol': 0.001, 'xatol': 0.001, 'maxiter': 1200, 'maxfev': 10900})
    #res = scipy.optimize.minimize(Fun, x0, method='Nelder-Mead',
    #                              options={'fatol': 0.1, 'xatol': 0.7, 'maxiter': 1000, 'maxfev': 19200})

    print("==opt output===")
    print(res.x)
    print(res.nit)
    print(res.success)
    print(res.fun)
    print("=======")

    self.FinalProbePlacement1 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[0],res.x[1],res.x[2],1])
    self.FinalProbePlacement2 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[3], res.x[4], res.x[5], 1])
    self.FinalProbePlacement3 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[6], res.x[7], res.x[8], 1])
    self.Metric = self.Probability(res.x, 1000)
    return res.x

  def Probability(self, X, No):
    NofInt = int(No)
    PC_vector = numpy.zeros(NofInt)
    DSC_vector = numpy.zeros(NofInt)
    HTA_vector = numpy.zeros(NofInt)
    Metric = numpy.zeros(NofInt)
    noiseGaussianX = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianY = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianZ = numpy.random.normal(0.0, self.sdterr/10.0, 10000)
    noiseGaussianX2 = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianY2 = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianZ2 = numpy.random.normal(0.0, self.sdterr/10.0, 10000)
    noiseGaussianX3 = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianY3 = numpy.random.normal(0.0, self.sdterr, 10000)
    noiseGaussianZ3 = numpy.random.normal(0.0, self.sdterr/10.0, 10000)
    if X[0] <= 0 or X[3] <=0 or X[6] <=0:
      return [0, 0.2, 0, 0, 0]
    for i in range(0, NofInt - 1):
      self.sphere1.SetCenter(X[0] + noiseGaussianX[i], X[1] + noiseGaussianY[i], X[2] + noiseGaussianZ[i])
      self.sphere1.Update()
      self.sphere2.SetCenter(X[3] + noiseGaussianX2[i], X[4] + noiseGaussianY2[i], X[5] + noiseGaussianZ2[i])
      self.sphere2.Update()
      self.sphere3.SetCenter(X[6] + noiseGaussianX3[i], X[7] + noiseGaussianY3[i], X[8] + noiseGaussianZ3[i])
      self.sphere3.Update()
      DICE_temp = self.ComputeDice3IceBall()
      Metric[i] = DICE_temp[0] + self.DiceGain * DICE_temp[1]
      PC_vector[i] = DICE_temp[1]
      DSC_vector[i] = DICE_temp[0]
      HTA_vector[i] = DICE_temp[2]
    count = 0
    for i in range(0, NofInt - 1):
      if PC_vector[i] >= 0.90:
        count = count + 1
    if count == 0:
        count = 10
    return [count, numpy.mean(Metric), numpy.mean(PC_vector), numpy.mean(DSC_vector), numpy.mean(HTA_vector)]

  def InitDataForPlanning(self,inputVolume2):
    self.sphere1 = vtk.vtkImageEllipsoidSource()
    self.sphere2 = vtk.vtkImageEllipsoidSource()
    self.sphere3 = vtk.vtkImageEllipsoidSource()
    self.sphere1.SetOutputScalarTypeToShort()
    self.sphere2.SetOutputScalarTypeToShort()
    self.sphere3.SetOutputScalarTypeToShort()
    self.PositionIceBall1 = [0, 0, 0, 1]
    self.PositionIceBall2 = [0, 0, 0, 1]
    self.PositionIceBall3 = [0, 0, 0, 1]
    self.PositionIceBall1_temp = self.PositionIceBall1
    self.PositionIceBall2_temp = self.PositionIceBall2
    self.PositionIceBall3_temp = self.PositionIceBall2
    self.Spacing = inputVolume2.GetSpacing()
    self.sphere1.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    self.sphere2.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    self.sphere3.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    size_image = inputVolume2.GetImageData().GetDimensions()
    self.sphere1.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere2.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere3.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere1.Update()
    self.sphere2.Update()
    self.sphere3.Update()
    self.logicFilter = vtk.vtkImageLogic()
    self.logicFilter2 = vtk.vtkImageLogic()
    self.logicFilter3 = vtk.vtkImageLogic()
    self.DiceGain = 5.0
    self.imData1 = inputVolume2.GetImageData()
    self.IjkToRasMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetIJKToRASMatrix(self.IjkToRasMatrix)


  def CreateModels(self, labelMapNode, modelHierarchyNode):

    modelMakerCLI = slicer.modules.modelmaker
    # tf = tempfile.NamedTemporaryFile(prefix='Slicer/Models-', suffix='.mrml')

    modelMakerParameters = {}
    # modelMakerParameters['ColorTable'] = 'vtkMRMLColorTableNodeFileGenericAnatomyColors.txt'
    modelMakerParameters['ModelSceneFile'] = modelHierarchyNode.GetID()
    modelMakerParameters['Name'] = 'Model'
    modelMakerParameters['GenerateAll'] = True
    modelMakerParameters['StartLabel'] = -1
    modelMakerParameters['EndLabel'] = -1
    modelMakerParameters['KkipUnNamed'] = True
    modelMakerParameters['JointSmooth'] = True
    modelMakerParameters['Smooth'] = 15
    modelMakerParameters['FilterType'] = 'Sinc'
    modelMakerParameters['Decimate'] = 0.25
    modelMakerParameters['SplitNormals'] = True
    modelMakerParameters['PointNormals'] = True
    modelMakerParameters['InputVolume'] = labelMapNode.GetID()

    slicer.cli.run(modelMakerCLI, None, modelMakerParameters, True)

  def GetCenterTumor(self, inputVolume2):
    modelHierarchyNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelHierarchyNode")
    slicer.mrmlScene.AddNode(modelHierarchyNode)
    self.CreateModels(inputVolume2, modelHierarchyNode)
    nOfModels = modelHierarchyNode.GetNumberOfChildrenNodes()
    if (nOfModels > 1):
      slicer.util.errorDisplay("More than one segmented ablation volume")
      return
    chnode = modelHierarchyNode.GetNthChildNode(0)
    mnode = chnode.GetAssociatedNode()
    objectPoly = mnode.GetPolyData()

    centerOfmass = vtk.vtkCenterOfMass()
    centerOfmass.SetInputData(objectPoly)
    centerOfmass.SetUseScalarsAsWeights(False)
    centerOfmass.Update()
    self.Center = centerOfmass.GetCenter()

    RasToIjkMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetRASToIJKMatrix(RasToIjkMatrix)
    b_Ijk = objectPoly.GetBounds()
    # Fisrt option:
    p_Ijka = [0, 0, 0, 1]
    p_Ijk2a = [0, 0, 0, 1]
    p_Ijka[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijka[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijka[2] = self.Center[2]
    p_Ijka[3] = 1
    p_Ijk2a[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2a[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2a[2] = self.Center[2]
    p_Ijk2a[3] = 1
    p_Ijka = RasToIjkMatrix.MultiplyDoublePoint(p_Ijka)
    p_Ijk2a = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2a)
    self.sphere1.SetCenter(p_Ijka[0], p_Ijka[1], p_Ijka[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2a[0], p_Ijk2a[1], p_Ijk2a[2])
    self.sphere2.Update()
    self.sphere3.SetCenter((p_Ijka[0]+p_Ijk2a[0])/2, (p_Ijka[1]+p_Ijk2a[1])/2, (p_Ijka[2]+p_Ijk2a[2])/2)
    self.sphere3.Update()
    DICE1 = self.ComputeDice3IceBall()

    # Second option:
    p_Ijkb = [0, 0, 0, 1]
    p_Ijk2b = [0, 0, 0, 1]
    p_Ijkb[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijkb[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijkb[2] = self.Center[2]
    p_Ijkb[3] = 1
    p_Ijk2b[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2b[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2b[2] = self.Center[2]
    p_Ijk2b[3] = 1
    p_Ijkb = RasToIjkMatrix.MultiplyDoublePoint(p_Ijkb)
    p_Ijk2b = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2b)
    self.sphere1.SetCenter(p_Ijkb[0], p_Ijkb[1], p_Ijkb[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2b[0], p_Ijk2b[1], p_Ijk2b[2])
    self.sphere2.Update()
    self.sphere3.SetCenter((p_Ijkb[0]+p_Ijk2b[0])/2, (p_Ijkb[1]+p_Ijk2b[1])/2, (p_Ijkb[2]+p_Ijk2b[2])/2)
    self.sphere3.Update()
    DICE2 = self.ComputeDice3IceBall()

    # Third option:
    p_Ijkc = [0, 0, 0, 1]
    p_Ijk2c = [0, 0, 0, 1]
    p_Ijkc[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijkc[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijkc[2] = self.Center[2]
    p_Ijkc[3] = 1
    p_Ijk2c[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2c[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2c[2] = self.Center[2]
    p_Ijk2c[3] = 1
    p_Ijkc = RasToIjkMatrix.MultiplyDoublePoint(p_Ijkc)
    p_Ijk2c = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2c)
    self.sphere1.SetCenter(p_Ijkc[0], p_Ijkc[1], p_Ijkc[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2c[0], p_Ijk2c[1], p_Ijk2c[2])
    self.sphere2.Update()
    self.sphere3.SetCenter((p_Ijkc[0]+p_Ijk2c[0])/2, (p_Ijkc[1]+p_Ijk2c[1])/2, (p_Ijkc[2]+p_Ijk2c[2])/2)
    self.sphere3.Update()
    DICE3 = self.ComputeDice3IceBall()

    # Fourth option:
    p_Ijkd = [0, 0, 0, 1]
    p_Ijk2d = [0, 0, 0, 1]
    p_Ijkd[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijkd[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijkd[2] = self.Center[2]
    p_Ijkd[3] = 1
    p_Ijk2d[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2d[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2d[2] = self.Center[2]
    p_Ijk2d[3] = 1
    p_Ijkd = RasToIjkMatrix.MultiplyDoublePoint(p_Ijkd)
    p_Ijk2d = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2d)
    self.sphere1.SetCenter(p_Ijkd[0], p_Ijkd[1], p_Ijkd[2])
    self.sphere1.Update()
    self.sphere2.SetCenter(p_Ijk2d[0], p_Ijk2d[1], p_Ijk2d[2])
    self.sphere2.Update()
    self.sphere3.SetCenter((p_Ijkd[0]+p_Ijk2d[0])/2, (p_Ijkd[1]+p_Ijk2d[1])/2, (p_Ijkd[2]+p_Ijk2d[2])/2)
    self.sphere3.Update()
    DICE4 = self.ComputeDice3IceBall()

    if (DICE1 > DICE2) and (DICE1 > DICE3) and (DICE1 > DICE4):
      self.PositionIceBall1 = [p_Ijka[0], p_Ijka[1], p_Ijka[2], 1]
      self.PositionIceBall2 = [p_Ijk2a[0], p_Ijk2a[1], p_Ijk2a[2], 1]
      self.PositionIceBall3 = [(p_Ijka[0]+p_Ijk2a[0])/2, (p_Ijka[1]+p_Ijk2a[1])/2, (p_Ijka[2]+p_Ijk2a[2])/2, 1]
      return
    elif (DICE2 > DICE1) and (DICE2 > DICE3) and (DICE2 > DICE4):
      self.PositionIceBall1 = [p_Ijkb[0], p_Ijkb[1], p_Ijkb[2], 1]
      self.PositionIceBall2 = [p_Ijk2b[0], p_Ijk2b[1], p_Ijk2b[2], 1]
      self.PositionIceBall3 = [(p_Ijkb[0]+p_Ijk2b[0])/2, (p_Ijkb[1]+p_Ijk2b[1])/2, (p_Ijkb[2]+p_Ijk2b[2])/2, 1]
      return
    elif (DICE3 > DICE1) and (DICE3 > DICE2) and (DICE3 > DICE4):
      self.PositionIceBall1 = [p_Ijkc[0], p_Ijkc[1], p_Ijkc[2], 1]
      self.PositionIceBall2 = [p_Ijk2c[0], p_Ijk2c[1], p_Ijk2c[2], 1]
      self.PositionIceBall3 = [(p_Ijkb[0]+p_Ijk2b[0])/2, (p_Ijkb[1]+p_Ijk2b[1])/2, (p_Ijkb[2]+p_Ijk2b[2])/2, 1]
      return
    elif (DICE4 > DICE1) and (DICE4 > DICE2) and (DICE4 > DICE3):
      self.PositionIceBall1 = [p_Ijkd[0], p_Ijkd[1], p_Ijkd[2], 1]
      self.PositionIceBall2 = [p_Ijk2d[0], p_Ijk2d[1], p_Ijk2d[2], 1]
      self.PositionIceBall3 = [(p_Ijkb[0]+p_Ijk2b[0])/2, (p_Ijkb[1]+p_Ijk2b[1])/2, (p_Ijkb[2]+p_Ijk2b[2])/2, 1]
      return

  def GetIceMAtrixRAS(self,image,pos):

    Probe1 = pos

    IceMatrix = vtk.vtkMatrix4x4()
    IceMatrix.SetElement(0, 0, 1.0)
    IceMatrix.SetElement(0, 1, 0.0)
    IceMatrix.SetElement(0, 2, 0.0)
    IceMatrix.SetElement(0, 3, Probe1[0])

    IceMatrix.SetElement(1, 0, 0.0)
    IceMatrix.SetElement(1, 1, 1.0)
    IceMatrix.SetElement(1, 2, 0.0)
    IceMatrix.SetElement(1, 3, Probe1[1])

    IceMatrix.SetElement(2, 0, 0.0)
    IceMatrix.SetElement(2, 1, 0.0)
    IceMatrix.SetElement(2, 2, 1.0)
    IceMatrix.SetElement(2, 3, Probe1[2])

    IceMatrix.SetElement(3, 0, 0.0)
    IceMatrix.SetElement(3, 1, 0.0)
    IceMatrix.SetElement(3, 2, 0.0)
    IceMatrix.SetElement(3, 3, 1.0)

    return IceMatrix

  def TestFunction(self, inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel):
    savetime = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMetric = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveDSC = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveHTA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMass = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    self.DiceGain = Kd


    cases = ["020", "029", "028", "032", "035", "036", "038", "042"]

    for i in range(0, 1):
      slicer.mrmlScene.Clear(0)
      #path = '/Users/pedro/Dropbox (Partners HealthCare)/DataForCryoPaper/Clinical%s' % (cases[i])
      path = '/Users/pedro/Dropbox (Partners HealthCare)/DataForCryoPaper/ralp-%s' % (cases[i])
      imageName = 'TumorRed'
      labelName = 'TumorRed-label'
      modelImageFileName = '%s/%s.nrrd' % (path, imageName)
      modelLabelFileName = '%s/%s.nrrd' % (path, labelName)

      (r, modelImageNode) = slicer.util.loadVolume(modelImageFileName, {}, True)
      (r, modelLabelNode) = slicer.util.loadLabelVolume(modelLabelFileName, {}, True)

      print(path)
      inputVolume2 = modelLabelNode
      print("==========")
      print(inputVolume2)
      print("==========")

      resultFileName = "Results-3Probes100-Case0%s" % (cases[i])
      resultFilePath = '/Users/pedro/Projects/MonteCarlo' + '/' + resultFileName
      resultFile = open(resultFilePath, 'a')

      resultFile.write("Probe; Err ; PoF ; PTC \n")

      start_time = time.time()
      self.InitDataForPlanning(inputVolume2)
      self.GetCenterTumor(inputVolume2)
      self.sdterr = 0
      X = self.optimizationCryo(N, Kd)

      for i in range(0, 16):
        self.sdterr = i
        self.output = outputLabel
        result = self.Probability(X, 10000)
        resultFile.write("3; %02d ; %02d ; %02f \n" % (i, result[0], result[2]))

      elapsed_time = time.time() - start_time
      print(elapsed_time)


    print(self.FinalProbePlacement1)
    print(self.FinalProbePlacement2)
    print(self.FinalProbePlacement3)
    self.affectedAreaModelNode = ModuleLogicMixin.createModelNode("AffectedArea")

    IceMatrix1 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement1)
    IceMatrix2 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement2)
    IceMatrix3 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement3)

    affectedBallArea = vtk.vtkParametricEllipsoid()
    affectedBallArea.SetXRadius(SizeX)
    affectedBallArea.SetYRadius(SizeY)
    affectedBallArea.SetZRadius(SizeZ)
    affectedBallAreaSource = vtk.vtkParametricFunctionSource()
    affectedBallAreaSource.SetParametricFunction(affectedBallArea)
    affectedBallAreaSource.SetScalarModeToV()
    affectedBallAreaSource.Update()
    affectedBallAreaAppend = vtk.vtkAppendPolyData()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix1)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix2)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix3)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    ModuleLogicMixin.createAndObserveDisplayNode(self.affectedAreaModelNode,
                                                 displayNodeClass=slicer.vtkMRMLModelDisplayNode)
    self.affectedAreaModelNode.GetDisplayNode().SetOpacity(0.2)
    self.affectedAreaModelNode.GetDisplayNode().SetColor(1.0, 0.0, 0.0)
    self.affectedAreaModelNode.SetAndObservePolyData(affectedBallAreaAppend.GetOutput())
    ModuleLogicMixin.setNodeVisibility(self.affectedAreaModelNode, True)
    ModuleLogicMixin.setNodeSliceIntersectionVisibility(self.affectedAreaModelNode, True)

    Distance = numpy.sqrt((IceMatrix1.GetElement(0, 3)-IceMatrix2.GetElement(0, 3))**2.0+(IceMatrix1.GetElement(1, 3)-IceMatrix2.GetElement(1, 3))**2.0)
    self.output.setText("The distance between probes is: %0.1f mm \n The probe locations are: \n Probe 1: (%0.1f,%0.1f,%0.1f) \n Probe 2:"
                        " (%0.1f,%0.1f,%0.1f) \n Probe 3:"
                        " (%0.1f,%0.1f,%0.1f) \n Probability of success: %0.1f\n Average tumor volume covered: %0.1f" % (Distance, self.FinalProbePlacement1[0],self.FinalProbePlacement1[1],self.FinalProbePlacement1[2],
                                                                                                               self.FinalProbePlacement2[0], self.FinalProbePlacement2[1], self.FinalProbePlacement2[2], self.FinalProbePlacement3[0], self.FinalProbePlacement3[1], self.FinalProbePlacement3[2],
                                                                                                                         self.Metric[0], self.Metric[2]))



class Cryo4Probes(ScriptedLoadableModuleLogic):

# RECOMECAR DAQUI.... fazer com 3 probes!!!
  def ComputeDice4IceBall(self):

    self.imData2 = self.sphere1.GetOutput()
    self.imData3 = self.sphere2.GetOutput()
    self.imData4 = self.sphere3.GetOutput()
    self.imData5 = self.sphere4.GetOutput()

    self.logicFilter2.SetOperationToOr()
    self.logicFilter2.SetInputData(0, self.imData2)
    self.logicFilter2.SetInputData(1, self.imData3)
    self.logicFilter2.Update()
    self.Spheres2 = self.logicFilter2.GetOutput()

    self.logicFilter3.SetOperationToOr()
    self.logicFilter3.SetInputData(0, self.imData4)
    self.logicFilter3.SetInputData(1, self.imData5)
    self.logicFilter3.Update()
    self.Spheres2b = self.logicFilter3.GetOutput()

    self.logicFilter4.SetOperationToOr()
    self.logicFilter4.SetInputData(0, self.Spheres2)
    self.logicFilter4.SetInputData(1, self.Spheres2b)
    self.logicFilter4.Update()
    self.Spheres4probes = self.logicFilter3.GetOutput()

    self.logicFilter5.SetOperationToOr()
    self.logicFilter5.SetInputData(0, self.Spheres4probes)
    self.logicFilter5.SetInputData(1, self.imData4)
    self.logicFilter5.Update()
    self.Spheres3 = self.logicFilter5.GetOutput()

    self.b = self.Spheres3.GetPointData()
    PixelImage2 = numpy.count_nonzero(self.b.GetArray(0))

    self.logicFilter.SetOperationToAnd()
    self.logicFilter.SetInputData(0, self.imData1) #check here.
    self.logicFilter.SetInputData(1, self.Spheres3)
    self.logicFilter.Update()
    self.OutPut = self.logicFilter.GetOutput()

    self.b = self.OutPut.GetPointData()
    intersection = numpy.count_nonzero(self.b.GetArray(0))

    self.b = self.imData1.GetPointData()
    PixelImage1 = numpy.count_nonzero(self.b.GetArray(0))

    return [(2.0 * intersection) / (PixelImage1 + PixelImage2), (1.0*intersection / PixelImage1)]

  def optimizationCryo(self, N, Kd):
    self.DiceGain = Kd
    def Fun(x):
        Metric = self.Probability(x, N)
        temp = Metric[0] / float(N)
        return -temp * temp * Metric[1]

    x0 = [self.PositionIceBall1[0], self.PositionIceBall1[1], self.PositionIceBall1[2],
          self.PositionIceBall2[0], self.PositionIceBall2[1], self.PositionIceBall2[2],
         self.PositionIceBall3[0], self.PositionIceBall3[1], self.PositionIceBall3[2],
         self.PositionIceBall4[0], self.PositionIceBall4[1], self.PositionIceBall4[2]]
    res = scipy.optimize.minimize(Fun, x0, method='Nelder-Mead', options={'fatol': 0.02, 'xatol': 0.1, 'maxiter': 300, 'maxfev': 300})
    print("==opt output===")
    print(res.x)
    print(res.nit)
    print(res.success)
    print(res.fun)
    print("=======")

    self.FinalProbePlacement1 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[0],res.x[1],res.x[2],1])
    self.FinalProbePlacement2 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[3], res.x[4], res.x[5], 1])
    self.FinalProbePlacement3 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[6], res.x[7], res.x[8], 1])
    self.FinalProbePlacement4 = self.IjkToRasMatrix.MultiplyDoublePoint([res.x[9], res.x[10], res.x[11], 1])
    self.Metric = self.Probability(res.x, 1000)

  def Probability(self, X, No):
    NofInt = int(No)
    PC_vector = numpy.zeros(NofInt)
    DSC_vector = numpy.zeros(NofInt)
    HTA_vector = numpy.zeros(NofInt)
    Metric = numpy.zeros(NofInt)
    noiseGaussianX = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianY = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianZ = numpy.random.normal(0.0, self.sdterr/10.0, 1000)
    noiseGaussianX2 = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianY2 = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianZ2 = numpy.random.normal(0.0, self.sdterr/10.0, 1000)
    noiseGaussianX3 = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianY3 = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianZ3 = numpy.random.normal(0.0, self.sdterr/10.0, 1000)
    noiseGaussianX4 = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianY4 = numpy.random.normal(0.0, self.sdterr, 1000)
    noiseGaussianZ4 = numpy.random.normal(0.0, self.sdterr/10.0, 1000)
    for i in range(0, NofInt - 1):
      self.sphere1.SetCenter(X[0] + noiseGaussianX[i], X[1] + noiseGaussianY[i], X[2] + noiseGaussianZ[i])
      self.sphere1.Update()
      self.sphere2.SetCenter(X[3] + noiseGaussianX2[i], X[4] + noiseGaussianY2[i], X[5] + noiseGaussianZ2[i])
      self.sphere2.Update()
      self.sphere3.SetCenter(X[6] + noiseGaussianX3[i], X[7] + noiseGaussianY3[i], X[8] + noiseGaussianZ3[i])
      self.sphere3.Update()
      self.sphere4.SetCenter(X[9] + noiseGaussianX4[i], X[10] + noiseGaussianY4[i], X[11] + noiseGaussianZ4[i])
      self.sphere4.Update()
      DICE_temp = self.ComputeDice4IceBall()
      Metric[i] = DICE_temp[0] + self.DiceGain * DICE_temp[1]
      PC_vector[i] = DICE_temp[1]
      DSC_vector[i] = DICE_temp[0]
    count = 0
    for i in range(0, NofInt - 1):
      if PC_vector[i] >= 0.97:
        count = count + 1
    if count == 0:
        count = 10
    return [count, numpy.mean(Metric), numpy.mean(PC_vector), numpy.mean(DSC_vector)]

  def InitDataForPlanning(self,inputVolume2):
    self.sphere1 = vtk.vtkImageEllipsoidSource()
    self.sphere2 = vtk.vtkImageEllipsoidSource()
    self.sphere3 = vtk.vtkImageEllipsoidSource()
    self.sphere4 = vtk.vtkImageEllipsoidSource()
    self.sphere1.SetOutputScalarTypeToShort()
    self.sphere2.SetOutputScalarTypeToShort()
    self.sphere3.SetOutputScalarTypeToShort()
    self.sphere4.SetOutputScalarTypeToShort()
    self.PositionIceBall1 = [0, 0, 0, 1]
    self.PositionIceBall2 = [0, 0, 0, 1]
    self.PositionIceBall3 = [0, 0, 0, 1]
    self.PositionIceBall4 = [0, 0, 0, 1]
    self.PositionIceBall1_temp = self.PositionIceBall1
    self.PositionIceBall2_temp = self.PositionIceBall2
    self.PositionIceBall3_temp = self.PositionIceBall2
    self.PositionIceBall4_temp = self.PositionIceBall2
    self.Spacing = inputVolume2.GetSpacing()
    self.sphere1.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    self.sphere2.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    self.sphere3.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    self.sphere4.SetRadius(SizeX / self.Spacing[0], SizeY / self.Spacing[1], SizeZ / self.Spacing[2])
    size_image = inputVolume2.GetImageData().GetDimensions()
    self.sphere1.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere2.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere3.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere4.SetWholeExtent(0, size_image[0] - 1, 0, size_image[1] - 1, 0, size_image[2] - 1)
    self.sphere1.Update()
    self.sphere2.Update()
    self.sphere3.Update()
    self.sphere4.Update()
    self.logicFilter = vtk.vtkImageLogic()
    self.logicFilter2 = vtk.vtkImageLogic()
    self.logicFilter3 = vtk.vtkImageLogic()
    self.logicFilter4 = vtk.vtkImageLogic()
    self.logicFilter5 = vtk.vtkImageLogic()
    self.DiceGain = 5.0
    self.imData1 = inputVolume2.GetImageData()
    self.IjkToRasMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetIJKToRASMatrix(self.IjkToRasMatrix)


  def CreateModels(self, labelMapNode, modelHierarchyNode):

    modelMakerCLI = slicer.modules.modelmaker
    # tf = tempfile.NamedTemporaryFile(prefix='Slicer/Models-', suffix='.mrml')

    modelMakerParameters = {}
    # modelMakerParameters['ColorTable'] = 'vtkMRMLColorTableNodeFileGenericAnatomyColors.txt'
    modelMakerParameters['ModelSceneFile'] = modelHierarchyNode.GetID()
    modelMakerParameters['Name'] = 'Model'
    modelMakerParameters['GenerateAll'] = True
    modelMakerParameters['StartLabel'] = -1
    modelMakerParameters['EndLabel'] = -1
    modelMakerParameters['KkipUnNamed'] = True
    modelMakerParameters['JointSmooth'] = True
    modelMakerParameters['Smooth'] = 15
    modelMakerParameters['FilterType'] = 'Sinc'
    modelMakerParameters['Decimate'] = 0.25
    modelMakerParameters['SplitNormals'] = True
    modelMakerParameters['PointNormals'] = True
    modelMakerParameters['InputVolume'] = labelMapNode.GetID()

    slicer.cli.run(modelMakerCLI, None, modelMakerParameters, True)

  def GetCenterTumor(self, inputVolume2):
    modelHierarchyNode = slicer.mrmlScene.CreateNodeByClass("vtkMRMLModelHierarchyNode")
    slicer.mrmlScene.AddNode(modelHierarchyNode)
    self.CreateModels(inputVolume2, modelHierarchyNode)
    nOfModels = modelHierarchyNode.GetNumberOfChildrenNodes()
    if (nOfModels > 1):
      slicer.util.errorDisplay("More than one segmented ablation volume")
      return
    chnode = modelHierarchyNode.GetNthChildNode(0)
    mnode = chnode.GetAssociatedNode()
    objectPoly = mnode.GetPolyData()

    centerOfmass = vtk.vtkCenterOfMass()
    centerOfmass.SetInputData(objectPoly)
    centerOfmass.SetUseScalarsAsWeights(False)
    centerOfmass.Update()
    self.Center = centerOfmass.GetCenter()

    RasToIjkMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetRASToIJKMatrix(RasToIjkMatrix)
    b_Ijk = objectPoly.GetBounds()
    # Fisrt option:
    p_Ijka = [0, 0, 0, 1]
    p_Ijk2a = [0, 0, 0, 1]
    p_Ijk3a = [0, 0, 0, 1]
    p_Ijk4a = [0, 0, 0, 1]
    p_Ijka[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijka[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijka[2] = self.Center[2]
    p_Ijka[3] = 1
    p_Ijk2a[0] = b_Ijk[0] + 1.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk2a[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk2a[2] = self.Center[2]
    p_Ijk2a[3] = 1
    p_Ijk3a[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk3a[1] = b_Ijk[2] + 1.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk3a[2] = self.Center[2]
    p_Ijk3a[3] = 1
    p_Ijk4a[0] = b_Ijk[0] + 2.0 * (b_Ijk[1] - b_Ijk[0]) / 3.0
    p_Ijk4a[1] = b_Ijk[2] + 2.0 * (b_Ijk[3] - b_Ijk[2]) / 3.0
    p_Ijk4a[2] = self.Center[2]
    p_Ijk4a[3] = 1

    p_Ijka = RasToIjkMatrix.MultiplyDoublePoint(p_Ijka)
    p_Ijk2a = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk2a)
    p_Ijk3a = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk3a)
    p_Ijk4a = RasToIjkMatrix.MultiplyDoublePoint(p_Ijk4a)

    self.PositionIceBall1 = [p_Ijka[0], p_Ijka[1], p_Ijka[2], 1]
    self.PositionIceBall2 = [p_Ijk2a[0], p_Ijk2a[1], p_Ijk2a[2], 1]
    self.PositionIceBall3 = [p_Ijk3a[0], p_Ijk3a[1], p_Ijk3a[2], 1]
    self.PositionIceBall4 = [p_Ijk4a[0], p_Ijk4a[1], p_Ijk4a[2], 1]

    return

  def GetIceMAtrixRAS(self,image,pos):

    Probe1 = pos

    IceMatrix = vtk.vtkMatrix4x4()
    IceMatrix.SetElement(0, 0, 1.0)
    IceMatrix.SetElement(0, 1, 0.0)
    IceMatrix.SetElement(0, 2, 0.0)
    IceMatrix.SetElement(0, 3, Probe1[0])

    IceMatrix.SetElement(1, 0, 0.0)
    IceMatrix.SetElement(1, 1, 1.0)
    IceMatrix.SetElement(1, 2, 0.0)
    IceMatrix.SetElement(1, 3, Probe1[1])

    IceMatrix.SetElement(2, 0, 0.0)
    IceMatrix.SetElement(2, 1, 0.0)
    IceMatrix.SetElement(2, 2, 1.0)
    IceMatrix.SetElement(2, 3, Probe1[2])

    IceMatrix.SetElement(3, 0, 0.0)
    IceMatrix.SetElement(3, 1, 0.0)
    IceMatrix.SetElement(3, 2, 0.0)
    IceMatrix.SetElement(3, 3, 1.0)

    return IceMatrix

  def TestFunction(self, inputVolume, inputVolume2, inputVolume3, inputVolume4, StdevError, NoOfProbes, Kd, N, outputLabel):
    savetime = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMetric = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveDSC = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveHTA = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    saveMass = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    self.sdterr = 0.5
    print("==========")
    print(inputVolume2)
    print("==========")
    self.DiceGain = Kd

    resultFileName = "Results-N%02d-K%02d-P4-E-1to5" % (N, self.DiceGain)
    resultFilePath = '/Users/pedro/Documents/Test-slicer' + '/' + resultFileName
    resultFile = open(resultFilePath, 'a')

    resultFile.write("Probe rew; N ; Kd ; Err ; PoS ; PTC ; DSC \n")

    for j in range(0, 13):
      for i in range(0,5):
        self.sdterr = j
        start_time = time.time()
        self.InitDataForPlanning(inputVolume2)
        self.GetCenterTumor(inputVolume2)
        self.output = outputLabel
        self.optimizationCryo(N,Kd)
        elapsed_time = time.time() - start_time
        print(elapsed_time)
        saveMetric[i] = self.Metric[0]
        saveDSC[i] = self.Metric[3]
        saveHTA[i] = self.Metric[2]
        saveMass[i] = 1
        resultFile.write("03 ; %02d ; %02d ; %02f ; %02d ; %02f ; %02f ; %02f ; %02f ; %02f ; %02f ; %02f ; %02f \n" % (
        N, self.DiceGain, self.sdterr, self.Metric[0], self.Metric[2], self.Metric[3], self.FinalProbePlacement1[0], self.FinalProbePlacement1[1], self.FinalProbePlacement1[2],self.FinalProbePlacement2[0],self.FinalProbePlacement2[1],self.FinalProbePlacement2[2]))



    start_time = time.time()
    self.InitDataForPlanning(inputVolume2)
    self.GetCenterTumor(inputVolume2)

    self.output = outputLabel

    RasToIjkMatrix = vtk.vtkMatrix4x4()
    inputVolume2.GetRASToIJKMatrix(RasToIjkMatrix)

    self.optimizationCryo(N,Kd)
    elapsed_time = time.time() - start_time
    print(elapsed_time)

    savetime[0] = elapsed_time
    saveMetric[0] = self.Metric[0]




    self.affectedAreaModelNode = ModuleLogicMixin.createModelNode("AffectedArea")

    IceMatrix1 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement1)
    IceMatrix2 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement2)
    IceMatrix3 = self.GetIceMAtrixRAS(inputVolume2, self.FinalProbePlacement3)

    affectedBallArea = vtk.vtkParametricEllipsoid()
    affectedBallArea.SetXRadius(SizeX)
    affectedBallArea.SetYRadius(SizeY)
    affectedBallArea.SetZRadius(SizeZ)
    affectedBallAreaSource = vtk.vtkParametricFunctionSource()
    affectedBallAreaSource.SetParametricFunction(affectedBallArea)
    affectedBallAreaSource.SetScalarModeToV()
    affectedBallAreaSource.Update()
    affectedBallAreaAppend = vtk.vtkAppendPolyData()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix1)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix2)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    transform = vtk.vtkTransform()
    transform.SetMatrix(IceMatrix3)
    tFilter2 = vtk.vtkTransformPolyDataFilter()
    tFilter2.SetTransform(transform)
    tFilter2.SetInputData(affectedBallAreaSource.GetOutput())
    tFilter2.Update()
    affectedBallAreaAppend.AddInputData(tFilter2.GetOutput())
    affectedBallAreaAppend.Update()

    ModuleLogicMixin.createAndObserveDisplayNode(self.affectedAreaModelNode,
                                                 displayNodeClass=slicer.vtkMRMLModelDisplayNode)
    self.affectedAreaModelNode.GetDisplayNode().SetOpacity(0.2)
    self.affectedAreaModelNode.GetDisplayNode().SetColor(1.0, 0.0, 0.0)
    self.affectedAreaModelNode.SetAndObservePolyData(affectedBallAreaAppend.GetOutput())
    ModuleLogicMixin.setNodeVisibility(self.affectedAreaModelNode, True)
    ModuleLogicMixin.setNodeSliceIntersectionVisibility(self.affectedAreaModelNode, True)

    Distance = numpy.sqrt((IceMatrix1.GetElement(0, 3)-IceMatrix2.GetElement(0, 3))**2.0+(IceMatrix1.GetElement(1, 3)-IceMatrix2.GetElement(1, 3))**2.0)
    self.output.setText("The distance between probes is: %0.1f mm \n The probe locations are: \n Probe 1: (%0.1f,%0.1f,%0.1f) \n Probe 2:"
                        " (%0.1f,%0.1f,%0.1f) \n Probe 3:"
                        " (%0.1f,%0.1f,%0.1f) \n Probability of success: %0.1f\n Average tumor volume covered: %0.1f" % (Distance, self.FinalProbePlacement1[0],self.FinalProbePlacement1[1],self.FinalProbePlacement1[2],
                                                                                                               self.FinalProbePlacement2[0], self.FinalProbePlacement2[1], self.FinalProbePlacement2[2], self.FinalProbePlacement3[0], self.FinalProbePlacement3[1], self.FinalProbePlacement3[2],
                                                                                                                         self.Metric[0], self.Metric[2]))

    print(IceMatrix1)
    print(IceMatrix2)

class teste1Test(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_CryoMC1()

  def test_CryoMC1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = CryoMCLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
