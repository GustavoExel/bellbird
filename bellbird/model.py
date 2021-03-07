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
	def __init__(self, variableName, value):
		self.variableName 			= variableName
		self.value 					= value

class Model:
	def __init__(self, name, equationsStr, variables, properties, boundaryConditions, meshDimension=2, definitions=[], meshPath=""):
		self.name 				= name
		self.equationsStr 		= equationsStr
		self.variables 			= variables
		self.properties 		= properties
		self.boundaryConditions = boundaryConditions
		self.definitions		= definitions
		self.meshPath 			= meshPath
		self.meshDimension 		= meshDimension

		self.compiled = False

		self.parseEquations()
		self.applyEbFVM()

	def parseEquations(self):
		definedNames = [ definition.split(" = ")[0] for definition in self.definitions ]
		self.definedVars = [ Constant(termStr) if not "[" in definition.split(" = ")[1] else (Vector(termStr) if not "[[" in definition.split(" = ")[1].replace(" ", "") else Matrix(termStr)) for termStr, definition in zip(definedNames, self.definitions)]

		self.equations = []
		for equationStr in self.equationsStr:
			equation = Equation(equationStr)
			equation.updatePropertyVars(self.properties)
			equation.updateVariables(self.variables)
			equation.updateDefinitions(self.definedVars)
			self.equations.append(equation)

		# Arrange equations and variables
		eqDimensions = [self.meshDimension if eq.order else 1 for eq in self.equations]

		self.arranjementDict = {"var": dict(), "eq": dict()}
		self.arranjementDict["var"] = { varName: idx for idx,varName in enumerate(self.variables) }
		self.arranjementDict["eq"]  = { idx: sum(eqDimensions[:idx]) for idx in range(len(self.equations)) }

	def applyEbFVM(self):
		self.discretizedEquations = []
		for equationStr in self.equationsStr:
			discretizedEquation = Equation(equationStr)

			discretizedEquation.updatePropertyVars(self.properties)
			discretizedEquation.updateVariables(self.variables)
			discretizedEquation.updateDefinitions(self.definedVars)
			discretizedEquation.integrateInSpace()				# a = b -> iiint(a) = iiint(b)
			discretizedEquation.applyDivergenceTheorem()		# iiint(div(f)) = iint(f)
			discretizedEquation.integrateInTime()				# d/dt(f) = (f-f_old)*(1/Δt)
			discretizedEquation.isolateVariables(self.variables) # x+a=y+c -> x-y=-a+c
			discretizedEquation.unapplyLinearProperty()

			self.discretizedEquations.append(discretizedEquation)

	def compile(self, fileName="results.py"):
		name = self.name.lower()[0] + "".join(self.name.split())[1:]
		self.fileName = fileName

		self.transient = True in [equation.isTransient() for equation in self.equations]
		self.transientFields = [ var.name for equation in self.equations for term in equation.terms if hasTermWithSubclass(term, TimeDerivative) for var in getTermFields(term) ]
		sizeStr = str(len(self.variables)) + " * numberOfVertices" if len(self.variables) > 1 else "numberOfVertices"

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
			if self.transient:
				write(t=t+1, ln="timeStep 		 = problemData.timeStep")
			write(t=t+1, ln="grid 			 = problemData.grid")
			write(t=t+1, ln="numberOfVertices = grid.numberOfVertices")
			write(t=t+1, ln="dimension 		 = grid.dimension", nl=2)
			write(t=t+1, ln="saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension='xdmf')\n")
		writeProblemData(0)

		def writeFields(t):
			for fieldName in self.variables:
				write(t=t, ln=f"{fieldName}Field    = np.repeat(0.0, numberOfVertices)")
				if fieldName.split("_")[0] in self.transientFields:
					write(t=t, ln=f"old{fieldName.capitalize()}Field = problemData.initialValues['{fieldName}'].copy()", nl=2)
		writeFields(1)

		def writeMatrix(t):
			write(t=t, ln=f"# matrix 		= np.zeros(({sizeStr}, {sizeStr}))")
			write(t=t, ln=f"# independent = np.zeros({sizeStr})", nl=2)
		writeMatrix(1)

		def declareProperties(t):
			for propertyName in self.properties:
				write(t=t, ln=f"{propertyName:<10} = propertyData.get(0, '{propertyName}')")
			write("")
		declareProperties(1)

		def writeDefinitions(t):
			for definition in self.definitions:
				write(t=t, ln=definition)
			write("")
		writeDefinitions(1)

		def writeUtilFuncs(t):
			def writeSymmetricUtilFuncs(t):
				write(t=t+0, ln="def getTransposedVoigtArea(innerFace):")
				write(t=t+1, ln="Sx, Sy, Sz = innerFace.area.getCoordinates()")
				write(t=t+1, ln="return np.array([[Sx,0,Sy],[0,Sy,Sx]]) if dimension==2 else np.array([[Sx,0,0,Sy,0,Sz],[0,Sy,0,Sx,Sz,0],[0,0,Sz,0,Sy,Sx]])", nl=2)
				write(t=t+0, ln="def getVoigtGradientOperator(globalDerivatives):")
				write(t=t+1, ln="if len(globalDerivatives) == 2:")
				write(t=t+2, ln="Nx,Ny = globalDerivatives")
				write(t=t+2, ln="zero=np.zeros(Nx.size)")
				write(t=t+2, ln="return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])", nl=2)
				write(t=t+1, ln="if len(globalDerivatives) == 3:")
				write(t=t+2, ln="Nx,Ny,Nz = globalDerivatives")
				write(t=t+2, ln="zero=np.zeros(Nx.size)")
				write(t=t+2, ln="return np.array([[Nx,zero,zero],[zero,Ny,zero],[zero,zero,Nz],[Ny,Nx,zero],[zero,Nz,Ny],[Nz,zero,Nx]])", nl=2)
			if True in [hasTermWithSubclass(equation.term, SymmetricGradient) for equation in self.equations]:
				writeSymmetricUtilFuncs(t)
		writeUtilFuncs(1)

		def assembleMatrix(t):
			write(t=t, ln="def assembleMatrix():")
			write(t=t+1, ln=f"matrix = np.zeros(({sizeStr}, {sizeStr}))", nl=2)

			def writeMatrixVolumeIntegrals(t):
				for discretizedEquation in self.discretizedEquations:
					for term in discretizedEquation.lhs:
						if term.__class__ == VolumetricIntegral and term.arg.__class__ == Multiplication:
							hasVar = False
							for varName in self.variables:
								if varName in map(lambda v:v.name, term.arg.args):
									hasVar = True
									termStr = str(term.arg).replace("Δt", "timeStep")
									coeffStr = "vertex.volume * " + termStr.replace(f"* {varName}", "")

									# In the future implement regions here!!!
									idxStr = "+" + str(self.arranjementDict["var"][varName]) + "*numberOfVertices" if len(self.variables) > 1 else ""
									write(t=t, ln="# " + termStr)
									write(t=t, ln="for vertex in grid.vertices:")
									write(t=t+1, ln=f"matrix[vertex.handle{idxStr}][vertex.handle{idxStr}] += " + coeffStr, nl=2)

							if not hasVar:
								raise Exception(f"Warning: Term {term} is being ignored")
						elif term.__class__ == VolumetricIntegral:
							raise Exception(f"Warning: Term {term} is being ignored")
			writeMatrixVolumeIntegrals(t+1)

			def writeMatrixSurfaceGradIntegrals(t):
				for discretizedEquation in self.discretizedEquations:
					for term in discretizedEquation.lhs:
						if term.__class__ == SurfaceIntegral and term.arg.__class__ == Multiplication:
							gradTerms = [arg for arg in term.arg.args if hasSubclass(arg, Gradient)]
							if gradTerms and hasSubclass(gradTerms[0].arg, Scalar):
								grad = gradTerms[0]
								var = getTermFields(grad)[0]
								idxStr = "+" + str(self.arranjementDict["var"][var.name]) + "*numberOfVertices" if len(self.variables) > 1 else ""

								coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad])
								coeffStr = str(coeff)

								write(t=t+0, ln="# " + str(term.arg))
								write(t=t+0, ln="for element in grid.elements:")
								write(t=t+1, ln="for innerFace in element.innerFaces:")
								write(t=t+2, ln="area = innerFace.area.getCoordinates()[:dimension]")
								write(t=t+2, ln=f"coefficients = {coeffStr} * np.matmul( area.T, innerFace.globalDerivatives )")
								write(t=t+2, ln="backwardVertexHandle, forwardVertexHandle = innerFace.getNeighborVerticesHandles()", nl=2)
								write(t=t+2, ln="for local, vertex in enumerate(element.vertices):")
								write(t=t+3, ln=f"matrix[backwardVertexHandle{idxStr}][vertex.handle{idxStr}] += coefficients[local]")
								write(t=t+3, ln=f"matrix[forwardVertexHandle{idxStr}][vertex.handle{idxStr}]  -= coefficients[local]", nl=2)
			writeMatrixSurfaceGradIntegrals(t+1)

			def writeMatrixSurfaceSymmetricGradVecIntegrals(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.lhs:
						if term.__class__ == SurfaceIntegral and term.arg.__class__ == Multiplication:
							gradTerms = [arg for arg in term.arg.args if hasSubclass(arg, SymmetricGradient)]
							if gradTerms and hasSubclass(gradTerms[0].arg, Vector):
								grad = gradTerms[0]

								coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad])
								coeffStr = str(coeff)

								arrIdx = self.arranjementDict["eq"][eqIdx]
								idxStr = lambda ref: f"numberOfVertices*{ref}" if (len(self.equations)==1 or arrIdx==0) else f"numberOfVertices*({ref}+{arrIdx})"

								write(t=t, ln=f"# {term.arg}")
								write(t=t, ln="for element in grid.elements:")
								write(t=t+1, ln="for innerFace in element.innerFaces:")
								write(t=t+2, ln="transposedVoigtArea = getTransposedVoigtArea(innerFace)")
								write(t=t+2, ln="voigtGradientOperator = getVoigtGradientOperator(innerFace.globalDerivatives)")
								write(t=t+2, ln=f"coeff = {coeffStr}")
								write(t=t+2, ln="matrixCoefficient = np.einsum('ij,jk,kmn->imn', transposedVoigtArea, coeff, voigtGradientOperator)", nl=2)
								write(t=t+2, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()", nl=2)
								write(t=t+2, ln="for local, vertex in enumerate(element.vertices):")
								write(t=t, ln="			for i in range(dimension):")
								write(t=t, ln="				for j in range(dimension):")
								write(t=t, ln=f"					matrix[backwardsHandle + {idxStr('i')}][vertex.handle + {idxStr('j')}] += matrixCoefficient[i][j][local]")
								write(t=t, ln=f"					matrix[forwardHandle   + {idxStr('i')}][vertex.handle + {idxStr('j')}] -= matrixCoefficient[i][j][local]", nl=2)


			writeMatrixSurfaceSymmetricGradVecIntegrals(t+1)

			def writeMatrixDirichletBoundaryConditions(t):
				write(t=t+0, ln="# Dirichlet Boundary Conditions")
				for variableName in self.variables:
					idxStr = "+" + str(self.arranjementDict["var"][variableName]) + "*numberOfVertices" if len(self.variables) > 1 else ""
					write(t=t+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
					write(t=t+1, ln="for vertex in bCondition.boundary.vertices:")
					write(t=t+2, ln=f"matrix[vertex.handle{idxStr}] = np.zeros({sizeStr})")
					write(t=t+2, ln=f"matrix[vertex.handle{idxStr}][vertex.handle{idxStr}] = 1.0", nl=2)
			writeMatrixDirichletBoundaryConditions(t+1)

			def inverseMatrix(t):
				write(t=t, ln="# Invert Matrix")
				write(t=t, ln="inverseMatrix = np.linalg.inv(matrix)")
				write(t=t, ln="return inverseMatrix", nl=2)
			inverseMatrix(t+1)
		assembleMatrix(1)

		def assembleIndependent(t):
			write(t=t, ln="def assembleIndependent():")
			write(t=t+1, ln=f"independent = np.zeros({sizeStr})", nl=2)

			def writeIndependentVolumeIntegral(t):
				for eqIdx, (equation, discretizedEquation) in enumerate(zip(self.equations, self.discretizedEquations)):
					for term in discretizedEquation.rhs:
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

							if equation.order == 0:
								idxStr = "+" + str(self.arranjementDict["eq"][eqIdx]) + "*numberOfVertices" if len(self.variables) > 1 else ""
								write(t=t, ln="# " + termStr)
								write(t=t, ln="for vertex in grid.vertices:")
								write(t=t+1, ln=f"independent[vertex.handle{idxStr}] += " + coeffStr, nl=2)
							elif equation.order == 1:
								arrIdx = self.arranjementDict["eq"][eqIdx]
								idxStr = "+" + (f"(coord+{arrIdx})" if arrIdx else "coord") + "*numberOfVertices"
								for arg in (term.arg.args if hasSubclass(term.arg, Operator) else term.arg):
									if hasSubclass(arg, Vector):
										coeffStr = coeffStr.replace(str(arg), f"{arg}[coord]")
								write(t=t, ln="# " + termStr)
								write(t=t, ln="for vertex in grid.vertices:")
								write(t=t+1, ln="for coord in range(dimension):")							
								write(t=t+2, ln=f"independent[vertex.handle{idxStr}] += " + coeffStr, nl=2)

						elif term.__class__ == VolumetricIntegral:
							raise Exception(f"Warning: Term {term} is being ignored")
			writeIndependentVolumeIntegral(t+1)

			def writeNeumannBoundaryConditions(t):
				# This signal implies fluxes inwards
				plusVars = [ var.name for equation in self.equations for term in equation.lhs if hasTermWithSubclass(term, Divergence) for var in getTermFields(term) ]
				write(t=t+0, ln="# Neumann Boundary Condition")
				for variableName in self.variables:
					idxStr = "+" + str(self.arranjementDict["var"][variableName]) + "*numberOfVertices" if len(self.variables) > 1 else ""
					signal = "+" if variableName.split("_")[0] in plusVars else "-"
					write(t=t+0, ln=f"for bCondition in problemData.neumannBoundaries['{variableName}']:")
					write(t=t+1, ln="for facet in bCondition.boundary.facets:")
					write(t=t+2, ln="for outerFace in facet.outerFaces:")
					write(t=t+3, ln=f"independent[outerFace.vertex.handle{idxStr}] {signal}= bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())", nl=2)
			writeNeumannBoundaryConditions(t+1)

			def writeDirichletBoundaryConditions(t):
				write(t=t+0, ln="# Dirichlet Boundary Condition")
				for variableName in self.variables:
					idxStr = "+" + str(self.arranjementDict["var"][variableName]) + "*numberOfVertices" if len(self.variables) > 1 else ""
					write(t=t+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
					write(t=t+1, ln="for vertex in bCondition.boundary.vertices:")
					write(t=t+2, ln=f"independent[vertex.handle{idxStr}] = bCondition.getValue(vertex.handle)", nl=2)
			writeDirichletBoundaryConditions(t+1)

			write(t=t+1, ln="return independent", nl=2)
		assembleIndependent(1)

		def writeTransientSolverLoop(t):
			# Written very specifically for problems with only one variable
			write(t=t+0, ln="tolerance = problemData.tolerance")
			write(t=t+0, ln="difference = 2*tolerance")
			write(t=t+0, ln="iteration = 0")
			write(t=t+0, ln="currentTime = 0.0")
			write(t=t+0, ln="converged = False", nl=2)
			write(t=t+0, ln="inverseMatrix = assembleMatrix()", nl=2)
			write(t=t+0, ln="while not converged:")
			write(t=t+1, ln="independent = assembleIndependent()")
			if len(self.variables) == 1:
				write(t=t+1, ln=f"{self.variables[0]}Field = np.matmul(inverseMatrix, independent)", nl=2)
			else:
				write("")
				write(t=t+1, ln=f"results = np.matmul(inverseMatrix, independent)")
				for i, variableName in enumerate(self.variables):
					write(t=t+1, ln=f"{variableName}Field = results[{i}*numberOfVertices:{i+1}*numberOfVertices]")
				write("")
			fieldDifferences = ['max(abs('+varName+'Field - old'+varName.capitalize()+'Field))' for varName in self.variables if varName.split('_')[0] in self.transientFields]
			if len(fieldDifferences) == 1:
				write(t=t+1, ln=f"difference = {fieldDifferences[0]}", nl=2)
			else:
				write(t=t+1, ln=f"difference = max( {', '.join(fieldDifferences)} )", nl=2)
			for variableName in self.variables:
				if variableName.split("_")[0] in self.transientFields:
					write(t=t+1, ln=f"old{variableName.capitalize()}Field = {variableName}Field.copy()")
			write("")
			write(t=t+1, ln="currentTime += timeStep")
			write(t=t+1, ln="iteration += 1", nl=2)
			for variableName in self.variables:
				write(t=t+1, ln=f"saver.save('{variableName}', {variableName}Field, currentTime)")
			write("")
			write(t=t+1, ln="print('{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}'.format(iteration, currentTime, timeStep, difference))")
			write(t=t+1, ln="converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )", nl=2)
			write(t=t+1, ln="if iteration >= problemData.maxNumberOfIterations:")
			write(t=t+2, ln="break", nl=2)
			write(t=t+0, ln="saver.finalize()", nl=2)
		def writePermanentSolverLoop(t):
			# Written very specifically for problems with only one variable
			write(t=t, ln="inverseMatrix = assembleMatrix()")
			write(t=t, ln="independent = assembleIndependent()", nl=2)
			if len(self.variables) == 1:
				write(t=t, ln=f"{self.variables[0]}Field = np.matmul(inverseMatrix, independent)")
			else:
				write(t=t, ln=f"results = np.matmul(inverseMatrix, independent)")
				for i, variableName in enumerate(self.variables):
					write(t=t, ln=f"{variableName}Field = results[{i}*numberOfVertices:{i+1}*numberOfVertices]")
			write(t=t, ln="")
			for variable in self.variables:
				write(t=t, ln=f"saver.save('{variable}', {variable}Field, 0.0)")
			write(t=t, ln="saver.finalize()", nl=2)
		if self.transient:
			writeTransientSolverLoop(1)
		else:
			writePermanentSolverLoop(1)

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
				write(t=t+1, ln="\t\t},")
			write(t=t+2, ln="}),\n\t)\n")
			write(t=t+1, ln=f"{name}( problemData )", nl=2)
			write(t=t+0, ln="if __name__ == '__main__':")
			write(t=t+1, ln="main()")
		writeMainFunction(0)

		with open(self.fileName, "w") as f:
			f.write(self.text)

		self.compiled = True

	def run(self):
		if not (self.compiled and self.fileName=="results.py"):
			self.compile()
		try:
			from results import main
			main()
		except:
			raise Exception("Error while compiling file!")
