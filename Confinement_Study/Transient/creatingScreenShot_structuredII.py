''' how to use
우선 각 폴더에 foam file이 있어야 한다.   -> type 'python3 makingFoamFile.py'
이 파일을 실행하기 전에 오픈폼을 언급해야됨 -> type 'of6'
이 파일 실행                             -> type 'python3 creatingScreenShot.py'
'''

import os
import subprocess
import numpy as np

Voltages = [Volt for Volt in range(20, 111, 30)] # [Volt for Volt in range(21, 41)]
boundaries = ['slip', 'noSlip'] # ['slip', 'Smoluchowski']
Ratios = [2, 4, 8, 16] # [1, 2, 8, 16]
Diameter = 20
addr = '/media/sh/drive3/24_04'
names = ['Diameter'] # ['Diameter', 'Square', 'Square_recombine']

dl = 1/3 # 1/3 um
dr = np.sqrt(4/np.sqrt(3)) * dl
Dn = int(Diameter*np.pi/dr/4)
Dn = int(Dn/2) + 1 # coarse: 기존의 1/2 만큼

for Ratio in Ratios:
  for name in names:
    for Voltage in Voltages:
      for boundary in boundaries:
        print("##################################################################")
        center_x = 0.0
        center_y = 0.0

        if 'Square' in name:
          center_x = 0.5
          center_y = 0.5
        seedHeight = 1/2

        radius = Ratio/2
        meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
        addr_foam = f'{addr}/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/{meshTitle}_{boundary}_V{Voltage}.foam'
        addr_png = f'{addr}/Result/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/{meshTitle}_{boundary}_V{Voltage}.png'
        print(f"{addr}/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}")

        str_head = f"""
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
OpenFOAM = OpenDataFile('{addr_foam}') # GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1439, 802]

# show data in view
OpenFOAMDisplay = Show(OpenFOAM, renderView1)
# trace defaults for the display properties.
OpenFOAMDisplay.Representation = 'Surface'
OpenFOAMDisplay.ColorArrayName = [None, '']
OpenFOAMDisplay.OSPRayScaleArray = 'U'
OpenFOAMDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
OpenFOAMDisplay.SelectOrientationVectors = 'None'
OpenFOAMDisplay.ScaleFactor = 0.19993360042572023
OpenFOAMDisplay.SelectScaleArray = 'None'
OpenFOAMDisplay.GlyphType = 'Arrow'
OpenFOAMDisplay.GlyphTableIndexArray = 'None'
OpenFOAMDisplay.DataAxesGrid = 'GridAxesRepresentation'
OpenFOAMDisplay.PolarAxes = 'PolarAxesRepresentation'
OpenFOAMDisplay.GaussianRadius = 0.09996680021286011
OpenFOAMDisplay.SetScaleArray = ['POINTS', 'p']
OpenFOAMDisplay.ScaleTransferFunction = 'PiecewiseFunction'
OpenFOAMDisplay.OpacityArray = ['POINTS', 'p']
OpenFOAMDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(OpenFOAMDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
OpenFOAMDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# Properties modified on OpenFOAMDisplay
OpenFOAMDisplay.Opacity = 0.2

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=OpenFOAM,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'U']
streamTracer1.MaximumStreamlineLength = 1.9993360042572021

# toggle 3D widget visibility (only when running from the GUI)
# Show3DWidgets(proxy=streamTracer1.SeedType)


# Properties modified on streamTracer1
streamTracer1.SeedType = 'Point Source'

# Properties modified on streamTracer1.SeedType
streamTracer1.SeedType.Center = [{center_x}, {center_y}, {seedHeight}]
streamTracer1.SeedType.NumberOfPoints = {pow(Ratio, 2)*2000}
streamTracer1.SeedType.Radius = {radius}

# show data in view
streamTracer1Display = Show(streamTracer1, renderView1)
# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'
streamTracer1Display.ColorArrayName = [None, '']
streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer1Display.SelectOrientationVectors = 'Normals'
streamTracer1Display.ScaleFactor = 0.19712519049644472
streamTracer1Display.SelectScaleArray = 'AngularVelocity'
streamTracer1Display.GlyphType = 'Arrow'
streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'
streamTracer1Display.GaussianRadius = 0.09856259524822236
streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'

# update the view to ensure updated data information
renderView1.Update()

# current camera placement for renderView1
renderView1.CameraPosition = [0.0006630122661590576, 0.0, 6.292992149887996]
renderView1.CameraFocalPoint = [0.0006630122661590576, 0.0, 0.5]
renderView1.CameraParallelScale = 1.4993366965204098
Hide3DWidgets(proxy=streamTracer1.SeedType)
# save animation

tube1 = Tube(Input=streamTracer1)
tube1.Scalars = ['POINTS', 'AngularVelocity']
tube1.Vectors = ['POINTS', 'Normals']
tube1.Radius = 0.006 # tube1.Radius = 0.006


# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [513, 802]

# show data in view
tube1Display = Show(tube1, renderView1)
# trace defaults for the display properties.
tube1Display.Representation = 'Surface'
tube1Display.ColorArrayName = [None, '']
tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display.SelectOrientationVectors = 'None'
tube1Display.ScaleFactor = -2.0000000000000002e+298
tube1Display.SelectScaleArray = 'None'
tube1Display.GlyphType = 'Arrow'
tube1Display.GlyphTableIndexArray = 'None'
tube1Display.DataAxesGrid = 'GridAxesRepresentation'
tube1Display.PolarAxes = 'PolarAxesRepresentation'
tube1Display.GaussianRadius = -1.0000000000000001e+298
tube1Display.SetScaleArray = [None, '']
tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display.OpacityArray = [None, '']
tube1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(streamTracer1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.0006630122661590576, 0.0, 6.292992149887996]
renderView1.CameraFocalPoint = [0.0006630122661590576, 0.0, 0.5]
renderView1.CameraParallelScale = 2.3439922623964304



animationScene1 = GetAnimationScene()

animationScene1.GoToLast()

SaveScreenshot('{addr_png}', renderView1, ImageResolution=[1439, 802], TransparentBackground=1)

"""
        with open(f'paraFoam.py', 'w') as f:
          f.write(str_head)
        os.system(f"pvpython paraFoam.py") # msh 파일 옮기기
