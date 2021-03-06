import os
from bellbird.equation import *
from bellbird.operators import *

class BoundaryCondition:
	bcDict = {"DIRICHLET_BOUNDARY_CONDITION":"Dirichlet", "NEUMANN_BOUNDARY_CONDITION":"Neumann"}
	def __init__(self, variableName, boundaryConditionType, boundaryName, value):
		self.variableName 			= variableName
		self.boundaryConditionType 	= boundaryConditionType
		self.boundaryName 			= boundaryName
		self.value 					= value
	@property
	def condition(self):
		return self.bcDict[self.boundaryConditionType]

class InitialCondition:
	def __init__(self, variableName, boundaryConditionType, value):
		self.variableName 			= variableName
		self.boundaryConditionType 	= boundaryConditionType
		self.value 					= value

class Model:
	def __init__(self, name, equationStr, variables, properties, boundaryConditions, meshPath=""):
		self.name 				= name
		self.equationStr 		= equationStr
		self.variables 			= variables
		self.properties 		= properties
		self.boundaryConditions = boundaryConditions
		self.meshPath 			= meshPath

		self.compiled = False
		self.equation = Equation(equationStr)

		self.applyEbFVM()

	def applyEbFVM(self):
		self.discretizedEquation = Equation(self.equationStr)

		self.discretizedEquation.updatePropertyVars(self.properties)
		self.discretizedEquation.updateVariables(self.variables)
		self.discretizedEquation.integrateInSpace()				# a = b -> iiint(a) = iiint(b)
		self.discretizedEquation.applyDivergenceTheorem()		# iiint(div(f)) = iint(f)
		self.discretizedEquation.integrateInTime()				# d/dt(f) = (f-f_old)*(1/Δt)
		self.discretizedEquation.isolateVariables(self.variables) # x+a=y+c -> x-y=-a+c
		self.discretizedEquation.unapplyLinearProperty()

	def compile(self):
		name = self.name.lower()[0] + "".join(self.name.split())[1:]
		self.fileName = "_".join(self.name.lower().split()) + ".py"

		self.text = ""
		def write(ln="", t=0, nl=1):
			self.text += t*'\t' + ln + nl*'\n'

		def writeHeader():
			write("import sys,os")
			write("sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'PyEFVLib'))")
			write("import numpy as np")
			write("import PyEFVLib", nl=3)
		writeHeader()

		def writeProblemData(t):
			write(t=t,   ln=f"def {name}(problemData):")
			write(t=t+1, ln="propertyData 	 = problemData.propertyData")
			write(t=t+1, ln="timeStep 		 = problemData.timeStep")
			write(t=t+1, ln="grid 			 = problemData.grid")
			write(t=t+1, ln="numberOfVertices = grid.numberOfVertices")
			write(t=t+1, ln="dimension 		 = grid.dimension", nl=2)
			write(t=t+1, ln="saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension='xdmf')\n")
		writeProblemData(0)

		def writeFields(t):
			for fieldName in self.variables:
				write(t=t, ln=f"old{fieldName.capitalize()}Field = problemData.initialValues['{fieldName}'].copy()")
				write(t=t, ln=f"{fieldName}Field    = np.repeat(0.0, numberOfVertices)", nl=2)
		writeFields(1)

		def writeMatrix(t):
			write(t=t, ln="# matrix 		= np.zeros((numberOfVertices, numberOfVertices))")
			write(t=t, ln="# independent = np.zeros(numberOfVertices)", nl=2)
		writeMatrix(1)

		def declareProperties(t):
			for propertyName in self.properties:
				write(t=t, ln=f"{propertyName:<10} = propertyData.get(0, '{propertyName}')")
			write("")
		declareProperties(1)

		def assembleMatrix(t):
			write(t=t, ln="def assembleMatrix():")
			write(t=t+1, ln="matrix = np.zeros((numberOfVertices, numberOfVertices))", nl=2)

			def writeMatrixVolumeIntegrals(t):
				for term in self.discretizedEquation.rhs:
					if term.__class__ == VolumetricIntegral and term.arg.__class__ == Multiplication:
						hasVar = False
						for varName in self.variables:
							if varName in map(lambda v:v.name, term.arg.args):
								hasVar = True
								termStr = str(term.arg).replace("Δt", "timeStep")
								coeffStr = "vertex.volume * " + termStr.replace(f"* {varName}", "")

								# In the future implement regions here!!!
								write(t=t, ln="# " + termStr)
								write(t=t, ln="for vertex in grid.vertices:")
								write(t=t+1, ln="matrix[vertex.handle][vertex.handle] += " + coeffStr, nl=2)

						if not hasVar:
							raise Exception(f"Warning: Term {term} is being ignored")
					elif term.__class__ == VolumetricIntegral:
						raise Exception(f"Warning: Term {term} is being ignored")
			writeMatrixVolumeIntegrals(t+1)

			def writeMatrixSurfaceGradIntegrals(t):
				for term in self.discretizedEquation.rhs:
					if term.__class__ == SurfaceIntegral and term.arg.__class__ == Multiplication:
						gradTerms = [arg for arg in term.arg.args if hasSubclass(arg, Gradient)]
						if gradTerms:
							grad = gradTerms[0]

							coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad])
							coeffStr = str(coeff)

							write(t=t+0, ln="# " + str(term.arg))
							write(t=t+0, ln="for element in grid.elements:")
							write(t=t+1, ln="for innerFace in element.innerFaces:")
							write(t=t+2, ln="area = innerFace.area.getCoordinates()[:dimension]")
							write(t=t+2, ln=f"coefficients = {coeffStr} * np.matmul( area.T, innerFace.globalDerivatives )")
							write(t=t+2, ln="backwardVertexHandle, forwardVertexHandle = innerFace.getNeighborVerticesHandles()", nl=2)
							write(t=t+2, ln="for local, vertex in enumerate(element.vertices):")
							write(t=t+3, ln="matrix[backwardVertexHandle][vertex.handle] += coefficients[local]")
							write(t=t+3, ln="matrix[forwardVertexHandle][vertex.handle]  -= coefficients[local]", nl=2)
			writeMatrixSurfaceGradIntegrals(t+1)

			def writeMatrixDirichletBoundaryConditions(t):
				for variableName in self.variables:
					write(t=t+0, ln="# Dirichlet Boundary Conditions")
					write(t=t+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
					write(t=t+1, ln="for vertex in bCondition.boundary.vertices:")
					write(t=t+2, ln="matrix[vertex.handle] = np.zeros(numberOfVertices)")
					write(t=t+2, ln="matrix[vertex.handle][vertex.handle] = 1.0", nl=2)
			writeMatrixDirichletBoundaryConditions(t+1)

			def inverseMatrix(t):
				write(t=t, ln="# Invert Matrix")
				write(t=t, ln="inverseMatrix = np.linalg.inv(matrix)", nl=2)
				write(t=t, ln="return inverseMatrix", nl=2)
			inverseMatrix(t+1)
		assembleMatrix(1)

		def assembleIndependent(t):
			write(t=t, ln="def assembleIndependent():")
			write(t=t+1, ln="independent = np.zeros(numberOfVertices)", nl=2)

			def writeIndependentVolumeIntegral(t):
				for term in self.discretizedEquation.lhs:
					if term.__class__ == VolumetricIntegral and (term.arg.__class__ == Multiplication or hasSubclass(term.arg, Variable)):

						# AINDA PRECISA TRATAR AS VARIÁVEIS DE CAMPO

						termStr = str(term.arg).replace("Δt", "timeStep")
						coeffStr = "vertex.volume * " + termStr

						for field in getTermFields(term.arg):
							if "_old" in field.name:
								fieldName = field.name.replace("_old", "").capitalize()
								coeffStr = coeffStr.replace(f" {field.name}", f" old{fieldName}Field[vertex.handle]")
							else:
								coeffStr = coeffStr.replace(f" {field.name}", f" {field.name}Field[vertex.handle]")

						write(t=t, ln="# " + termStr)
						write(t=t, ln="for vertex in grid.vertices:")
						write(t=t+1, ln="independent[vertex.handle] += " + coeffStr, nl=2)

					elif term.__class__ == VolumetricIntegral:
						raise Exception(f"Warning: Term {term} is being ignored")
			writeIndependentVolumeIntegral(t+1)

			def writeDirichletBoundaryConditions(t):
				for variableName in self.variables:
					write(t=t+0, ln="# Dirichlet Boundary Condition")
					write(t=t+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
					write(t=t+1, ln="for vertex in bCondition.boundary.vertices:")
					write(t=t+2, ln="independent[vertex.handle] = bCondition.getValue(vertex.handle)", nl=2)
			writeDirichletBoundaryConditions(t+1)

			def writeNeumannBoundaryConditions(t):
				for variableName in self.variables:
					write(t=t+0, ln="# Neumann Boundary Condition")
					write(t=t+0, ln=f"for bCondition in problemData.neumannBoundaries['{variableName}']:")
					write(t=t+1, ln="for facet in bCondition.boundary.facets:")
					write(t=t+2, ln="for outerFace in facet.outerFaces:")
					write(t=t+3, ln="independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())", nl=2)
			writeNeumannBoundaryConditions(t+1)

			write(t=t+1, ln="return independent")
		assembleIndependent(1)

		def writeSolverLoop(t):
			# Written very specifically for problems with only one variable
			write(t=t+0, ln="tolerance = problemData.tolerance")
			write(t=t+0, ln="difference = 2*tolerance")
			write(t=t+0, ln="iteration = 0")
			write(t=t+0, ln="currentTime = 0.0")
			write(t=t+0, ln="converged = False", nl=2)
			write(t=t+0, ln="inverseMatrix = assembleMatrix()", nl=2)
			write(t=t+0, ln="while not converged:")
			write(t=t+1, ln="independent = assembleIndependent()")
			write(t=t+1, ln=f"{self.variables[0]}Field = np.matmul(inverseMatrix, independent)", nl=2)
			write(t=t+1, ln=f"difference = max( abs({self.variables[0]}Field - old{self.variables[0].capitalize()}Field) )")
			write(t=t+1, ln=f"old{self.variables[0].capitalize()}Field = {self.variables[0]}Field.copy()")
			write(t=t+1, ln="currentTime += timeStep")
			write(t=t+1, ln="iteration += 1", nl=2)
			write(t=t+1, ln=f"saver.save('{self.variables[0]}', {self.variables[0]}Field, currentTime)", nl=2)
			write(t=t+1, ln="print('{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}'.format(iteration, currentTime, timeStep, difference))")
			write(t=t+1, ln="converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )", nl=2)
			write(t=t+1, ln="if iteration >= problemData.maxNumberOfIterations:")
			write(t=t+2, ln="break", nl=2)
			write(t=t+0, ln="saver.finalize()", nl=2)
		writeSolverLoop(1)

		def writeMainFunction(t):
			def getRegionNames():
				try:
					import sys,os
					sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'PyEFVLib'))
					import PyEFVLib
					return PyEFVLib.read(self.meshPath).gridData.regionsNames
				except:
					return ["Body"]
			regionNames = getRegionNames()
			write(t=t+0, ln="def main():")
			write(t=t+1, ln="problemData = PyEFVLib.ProblemData(")
			write(t=t+1, ln=f"	meshFilePath = '{self.meshPath}',")
			write(t=t+2, ln="outputFilePath = 'results',")
			write(t=t+2, ln="numericalSettings = PyEFVLib.NumericalSettings( timeStep = 0.1, tolerance = 1e-4, maxNumberOfIterations = 300 ),")
			write(t=t+2, ln="propertyData = PyEFVLib.PropertyData({")
			for regionName in regionNames:
				write(t=t+3, ln=f"'{regionName}':")
				write(t=t+3, ln="{")
				for propertyName, propertyValue in self.properties.items():
					write(t=t+4, ln=f"'{propertyName}': {propertyValue},")
				write(t=t+3, ln="},")
			write(t=t+2, ln="}),")
			write(t=t+2, ln="boundaryConditions = PyEFVLib.BoundaryConditions({")
			for variableName in self.variables:
				write(t=t+3, ln=f"'{variableName}': {{")
				for bc in self.boundaryConditions:
					if bc.__class__ == InitialCondition and bc.variableName == variableName:
						write(t=t+4, ln=f"'InitialValue': {bc.value},")
					if bc.__class__ == BoundaryCondition and bc.variableName == variableName:
						write(t=t+4, ln=f"'{bc.boundaryName}': {{ 'condition' : PyEFVLib.{bc.condition}, 'type' : PyEFVLib.Constant, 'value' : {bc.value} }},")
			write(t=t+1, ln="\t\t}\n\t\t}),\n\t)\n")
			write(t=t+1, ln=f"{name}( problemData )", nl=2)
			write(t=t+0, ln="if __name__ == '__main__':")
			write(t=t+1, ln="main()")
		writeMainFunction(0)

		with open(self.fileName, "w") as f:
			f.write(self.text)

		self.compiled = True

	def run(self):
		if not self.compiled:
			self.compile()
		try:
			os.rename(self.fileName, "results.py")
			from results import main
			main()
			os.rename("results.py", self.fileName)
		except:
			Exception("Error while compiling file!")
