# Split hexahedra in paraview into separate polytopes
# Allows visualizing single hexahedra, toggling on/off hexahedra.

from paraview.simple import *
from pathlib import Path
import sys
import random

hexmesh_dir = str(Path.home()) + '\Data\OFF\\3D\hexMesh'
default_filename = 'hexMesh.off'

try:
    # Must enter the full path since the directory is currently
    #   wherever this script is stored.
    filename = input("Enter hex mesh .off file:\n" + \
                     "(Directory: " + hexmesh_dir + "." \
                     "  Default filename: " + default_filename + ")")
except EOFError:
    filename = ""

# If filename is blank...
if not filename or filename.isspace():
    # Hard wire input file.
    filename = default_filename
    
filepath = hexmesh_dir + "\\" + filename;
print("\n\nUsing input file: ", filepath)

# read in .off file
with open(filepath) as f:
    header = next(f)
    numPts,numCells,nul = [int(x) for x in next(f).split()] # read first line
    array = []
    line_number = 0
    for line in f: # read rest of lines
        if line_number<numPts:
            array.append([float(x) for x in line.split()])  #read in points
        else:
            array.append([int(x) for x in line.split()])  # read in connectivity
        line_number = line_number + 1

# extract vertex coordinates and cell connectivity from input file
ptsArray = array[0:numPts]
cellsArray = array[numPts+1:numPts+1+numCells]  # '+1' is hard-coded, which represent blank lines in input file

# convert the connectivity to 8 points cube format if it is read in as 4 points square format
if cellsArray[0][0] == 4:
    numCells /= 2
    cellsArray4D = []
    for i in range(numCells):
        cellsArray4D.append(cellsArray[2*i][:] + cellsArray[2*i+1][1:])
        cellsArray4D[i][0] *= 2
        cellsArray4D[i][3], cellsArray4D[i][4] = cellsArray4D[i][4], cellsArray4D[i][3]
        cellsArray4D[i][7], cellsArray4D[i][8] = cellsArray4D[i][8], cellsArray4D[i][7]
    cellsArray = cellsArray4D


def programmable_source(pts, cells, cubeNum):

    numPts = len(pts)
    numCells = len(cells)

    # output script
    output = ''

    output = output + 'numPts = ' + str(numPts) + '\n'
    output = output + 'numCells = ' + str(numCells) + '\n\n'

    output = output + '#Get a vtk.PolyData object for the output \n'
    output = output + 'pdo = self.GetPolyDataOutput()\n'
    output = output + '#This will store the points for the octahedron\n'
    output = output + 'newPts = vtk.vtkPoints()\n\n'

    output = output + '#Set vertices\n'

    for i in range(cells[cubeNum][0]):
        output = output + 'newPts.InsertPoint('+str(i)+','+ \
                     str(pts[cells[cubeNum][i+1]][0])+','+ \
                     str(pts[cells[cubeNum][i+1]][1])+','+ \
                     str(pts[cells[cubeNum][i+1]][2])+') \n'

    output = output + '\n#Add the points to the vtkPolyData object\n'
    output = output + 'pdo.SetPoints(newPts)\n\n'

    output = output + '#Allocate number of cells\n'
    output = output + 'pdo.Allocate(1,1)\n\n'

    output = output + '#Add cells\n'
    output = output + 'quad = vtk.vtkQuad()\n'
    output = output + 'quad.GetPointIds().SetId(0,0)\n'
    output = output + 'quad.GetPointIds().SetId(1,1)\n'
    output = output + 'quad.GetPointIds().SetId(2,3)\n'
    output = output + 'quad.GetPointIds().SetId(3,2)\n'
    output = output + 'pdo.InsertNextCell(quad.GetCellType(),quad.GetPointIds())\n\n'

    output = output + 'quad = vtk.vtkQuad()\n'
    output = output + 'quad.GetPointIds().SetId(0,0)\n'
    output = output + 'quad.GetPointIds().SetId(1,4)\n'
    output = output + 'quad.GetPointIds().SetId(2,6)\n'
    output = output + 'quad.GetPointIds().SetId(3,2)\n'
    output = output + 'pdo.InsertNextCell(quad.GetCellType(),quad.GetPointIds())\n'

    output = output + 'quad = vtk.vtkQuad()\n'
    output = output + 'quad.GetPointIds().SetId(0,0)\n'
    output = output + 'quad.GetPointIds().SetId(1,1)\n'
    output = output + 'quad.GetPointIds().SetId(2,5)\n'
    output = output + 'quad.GetPointIds().SetId(3,4)\n'
    output = output + 'pdo.InsertNextCell(quad.GetCellType(),quad.GetPointIds())\n'

    output = output + 'quad = vtk.vtkQuad()\n'
    output = output + 'quad.GetPointIds().SetId(0,3)\n'
    output = output + 'quad.GetPointIds().SetId(1,7)\n'
    output = output + 'quad.GetPointIds().SetId(2,5)\n'
    output = output + 'quad.GetPointIds().SetId(3,1)\n'
    output = output + 'pdo.InsertNextCell(quad.GetCellType(),quad.GetPointIds())\n'

    output = output + 'quad = vtk.vtkQuad()\n'
    output = output + 'quad.GetPointIds().SetId(0,3)\n'
    output = output + 'quad.GetPointIds().SetId(1,2)\n'
    output = output + 'quad.GetPointIds().SetId(2,6)\n'
    output = output + 'quad.GetPointIds().SetId(3,7)\n'
    output = output + 'pdo.InsertNextCell(quad.GetCellType(),quad.GetPointIds())\n'

    output = output + 'quad = vtk.vtkQuad()\n'
    output = output + 'quad.GetPointIds().SetId(0,7)\n'
    output = output + 'quad.GetPointIds().SetId(1,6)\n'
    output = output + 'quad.GetPointIds().SetId(2,4)\n'
    output = output + 'quad.GetPointIds().SetId(3,5)\n'
    output = output + 'pdo.InsertNextCell(quad.GetCellType(),quad.GetPointIds())\n\n'

    return output

for i in range(numCells):
    output = programmable_source(ptsArray, cellsArray, i)
    obj = ProgrammableSource(Script=output)

    # Change colors and opacity by calling property
    prpty = GetDisplayProperties(obj)
    prpty.Opacity = 0.7
    prpty.DiffuseColor = [random.uniform(0.0,1.0), random.uniform(0.0,1.0), random.uniform(0.0,1.0)]

    # Show
    Show()
    Render()
