from bellbird.equation import *
from bellbird.operators import *

class Model:
	def __init__(self, name, equationStr, variables, properties, boundaryConditions=None, meshPath=""):
		self.name 				= name
		self.equationStr 		= equationStr
		self.variables 			= variables
		self.properties 		= properties
		self.boundaryConditions = boundaryConditions
		self.meshPath 			= meshPath

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
			write(t=t+0, ln="saver.finalize()")
		writeSolverLoop(1)

		self.text += "\n\n\nif __name__ == '__main__':\n\tproblemData = PyEFVLib.ProblemData(\n\t\tmeshFilePath = '{MESHES}/msh/2D/Square.msh',\n\t\toutputFilePath = '{RESULTS}/../../BellBird2/results',\n\t\tnumericalSettings = PyEFVLib.NumericalSettings( timeStep = 0.1, tolerance = 1e-4, maxNumberOfIterations = 300 ),\n\t\tpropertyData = PyEFVLib.PropertyData({\n\t\t    'Body':\n\t\t    {\n\t\t    \t'k': 1,\n\t\t    \t'rho': 1,\n\t\t    \t'cp': 1,\n\t\t    \t'q': 0,\n\t\t    }\n\t\t}),\n\t\tboundaryConditions = PyEFVLib.BoundaryConditions({\n\t\t\t'temperature': {\n\t\t\t\t'InitialValue': 0.0,\n\t\t\t\t'West':\t { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant,'value' : 20.0 },\n\t\t\t\t'East':\t { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant,'value' : 50.0 },\n\t\t\t\t'South': { 'condition' : PyEFVLib.Neumann,   'type' : PyEFVLib.Constant,'value' : 0.0 },\n\t\t\t\t'North': { 'condition' : PyEFVLib.Neumann,   'type' : PyEFVLib.Constant,'value' : 0.0 },\n\t\t\t}\n\t\t}),\n\t)\n\theatTransfer( problemData )\n"

		with open("results.py", "w") as f:
			f.write(self.text)