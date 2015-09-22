#Copyright 2009-2014 Thomas Paviot (tpaviot@gmail.com)
#
#This file is part of pythonOCC.
#
#pythonOCC is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#pythonOCC is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public License
#along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeTorus, BRepPrimAPI_MakeSphere
from OCC.SMESH import SMESH_Gen, SMESH_MeshVSLink

from OCC.NETGENPlugin import *
from OCC.TetgenPlugin import *

from OCC.StdMeshers import (StdMeshers_Arithmetic1D,
                            StdMeshers_QuadranglePreference,
                            StdMeshers_Regular_1D,
                            StdMeshers_Quadrangle_2D)
from OCC.MeshVS import MeshVS_MeshPrsBuilder, MeshVS_Mesh

from OCC.Display.SimpleGui import init_display

# create the shape to mesh
aShape = BRepPrimAPI_MakeTorus(40, 20).Shape()
#aShape = BRepPrimAPI_MakeBox(20, 20, 20).Shape()
#aShape = BRepPrimAPI_MakeSphere(10).Shape()

# create mesh
aMeshGen = SMESH_Gen()

# netgen hypotheses
netgenhyp_1 = NETGENPlugin_SimpleHypothesis_2D(1,0,aMeshGen)
netgenhyp_1.SetMaxElementArea(100)

netgenalgo = NETGENPlugin_NETGEN_2D(2, 0, aMeshGen)


# create hypotheses
tetgenhyp_1 = TetgenPlugin_Hypothesis(3,0,aMeshGen)
tetgenhyp_1.SetMaxVolume(50)

tetgenalgo = TetgenPlugin_Tetgen_3D(4, 0, aMeshGen)
tetgenalgo.setKeepIntermediateFiles(True)

# set shape to mesh
aMesh = aMeshGen.CreateMesh(0, True)
aMesh.ShapeToMesh(aShape)

# assign hyptothesis to mesh
aMesh.AddHypothesis(aShape, 1)
aMesh.AddHypothesis(aShape, 2)

#Compute the mesh
#tetgenalgo.Compute(aMesh, aShape)
aMeshGen.Compute(aMesh, aShape)

aMesh.AddHypothesis(aShape, 3)
aMesh.AddHypothesis(aShape, 4)
tetgenalgo.Compute(aMesh, aShape)

#
# Display the Mesh
display, start_display, add_menu, add_function_to_menu = init_display()
aDS = SMESH_MeshVSLink(aMesh)
aMeshVS = MeshVS_Mesh(True)
DMF = 1  # to wrap!
MeshVS_BP_Mesh = 5  # To wrap!

aPrsBuilder = MeshVS_MeshPrsBuilder(aMeshVS.GetHandle(),
                                    DMF,
                                    aDS.GetHandle(),
                                    0,
                                    MeshVS_BP_Mesh)
aMeshVS.SetDataSource(aDS.GetHandle())
aMeshVS.AddBuilder(aPrsBuilder.GetHandle(), True)

# display both shape and mesh
context = display.Context
context.Display(aMeshVS.GetHandle())
context.Deactivate(aMeshVS.GetHandle())
#display.DisplayShape(aShape)
start_display()
