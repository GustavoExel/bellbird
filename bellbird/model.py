import os
from bellbird.equation import *
from bellbird.operators import *

class Model:
	def __init__(self, name, equationsStr, variables, properties, boundaryConditions, definitions=[], meshPath="", sparse=False, timeStep=0.1, tolerance=1e-4, maxNumberOfIterations=300):
		self.name 				= name
		self.equationsStr 		= equationsStr
		self.variables 			= variables
		self.properties 		= properties
		self.boundaryConditions = boundaryConditions
		self.definitions		= definitions
		self.meshPath 			= meshPath
		self.sparse				= sparse
		self.timeStep			= timeStep
		self.tolerance			= tolerance
		self.maxNumberOfIterations = maxNumberOfIterations

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

		def _rearrange(terms):
			return str(Equation(f"x={'+'.join(map(str,terms))}").rearrange().term.args[1])

		self.varRadicals = [var.split("_")[0] for var in self.variables]
		self.varRadicals = [var for idx, var in enumerate(self.varRadicals) if not var in self.varRadicals[:idx]]
		varDims = ["1" if var in self.variables else "dimension" for var in self.varRadicals]
		varOffsets = ["0"] + varDims[:-1]

		eqDimensions = ["dimension" if eq.order else "1" for eq in self.equations]
		eqOffsets = ["0"] + [_rearrange(eqDimensions[:idx+1]) for idx,eq_dim in enumerate(eqDimensions[:-1])]

		self.arranjementDict = {"var": dict(), "eq": dict()}
		self.arranjementDict["var"] = {varRadical+("" if varDimension=="1" else "_"+"xyz"[idx]) : _rearrange([varOffset,idx]) for varOffset,varRadical,varDimension in zip(varOffsets,self.varRadicals,varDims) for idx in range({"1":1,"dimension":3}[varDimension])}
		self.arranjementDict["eq"]  = { idx:arr for idx,arr in enumerate(eqOffsets) }

		numberOfVarsStr = _rearrange([varDims.count('1'), f"{varDims.count('dimension')}*dimension"]).replace(" ","")
		self.sizeStr = f"({numberOfVarsStr})*numberOfVertices"

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
		# self.sizeStr = str(len(self.variables)) + " * numberOfVertices" if len(self.variables) > 1 else "numberOfVertices"

		self.text = ""
		def write(ln="", t=0, nl=1):
			self.text += t*'\t' + ln + nl*'\n'

		def add(t, i, j, val, nl=1):
			if self.sparse:
				write(t=t, ln=f"add({i}, {j}, {val})", nl=nl)
			else:
				write(t=t, ln=f"matrix[{i}][{j}] += {val}", nl=nl)

		def set(t, i, nl):
			if self.sparse:
				write(t=t, ln=f"set({i})", nl=nl)
			else:
				write(t=t, ln=f"matrix[{i}] = np.zeros({self.sizeStr})")
				write(t=t, ln=f"matrix[{i}][{i}] = 1.0", nl=nl)

		def writeHeader():
			write("import sys,os")
			write("sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'PyEFVLib'))")
			write("import PyEFVLib")
			write("import numpy as np")
			if self.sparse:
				write("from scipy import sparse")
				write("import scipy.sparse.linalg")
			write("")
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
			for fieldName in self.varRadicals:
				write(t=t, ln=f"{fieldName}Field    = np.repeat(0.0, {'' if fieldName in self.variables else 'dimension*'}numberOfVertices)")
				if fieldName.split("_")[0] in self.transientFields:
					if fieldName in self.variables:
						write(t=t, ln=f"old{fieldName.capitalize()}Field = problemData.initialValues['{fieldName}'].copy()", nl=2)
					else:
						write(t=t, ln=f"old{fieldName.capitalize()}Field = np.concatenate((problemData.initialValues['{fieldName}_x'], problemData.initialValues['{fieldName}_y'], problemData.initialValues['{fieldName}_z'])) if dimension==3 else np.concatenate((problemData.initialValues['{fieldName}_x'], problemData.initialValues['{fieldName}_y']))", nl=2)
		writeFields(1)

		def declareMatrix(t):
			if self.sparse:
				write(t=t, ln="matrixCoords = []")
				write(t=t, ln="matrixVals = []", nl=2)
			else:
				write(t=t, ln=f"matrix 		= np.zeros(({self.sizeStr}, {self.sizeStr}))", nl=2)
			# write(t=t, ln=f"independent = np.zeros({self.sizeStr})", nl=2)
		declareMatrix(1)

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
			def writeSparseAddFunc(t):
				write(t=t+0, ln="def add(i, j, val):")
				write(t=t+1, ln="matrixCoords.append((i, j))")
				write(t=t+1, ln="matrixVals.append(val)", nl=2)
			def writeSparseSetFunc(t):
				write(t=t+0, ln="def set(i):")
				write(t=t+1, ln="_matrixCoords = [coord for coord in matrixCoords if coord[0]!=i] + [(i,i)]")
				write(t=t+1, ln="_matrixVals   = [val for coord,val in zip(matrixCoords,matrixVals) if coord[0]!=i] + [1.0]")
				write(t=t+1, ln="matrixCoords.clear(); matrixCoords.extend(_matrixCoords)")
				write(t=t+1, ln="matrixVals.clear(); matrixVals.extend(_matrixVals)", nl=2)
			if self.sparse:
				writeSparseAddFunc(t)
				writeSparseSetFunc(t)

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

			def writeMatrixScalarFieldVolumeIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.lhs:
						if term.__class__ == VolumetricIntegral:
							for field in getTermFields(term):
								if field.__class__ == ScalarField and (field == term.arg or ( term.arg.__class__ == Multiplication and field in term.arg.args )):
									termStr = str(term.arg)
									coeffStr = "vertex.volume * " + termStr.replace(f"* {field.name}", "")

									arrIdx1 = self.arranjementDict["eq"][eqIdx]
									idxStr1 = f"+({arrIdx1})*numberOfVertices" if arrIdx1!="0" else ""

									arrIdx2 = self.arranjementDict["var"][field.name]
									idxStr2 = f"+({arrIdx2})*numberOfVertices" if arrIdx2!="0" else ""

									# In the future implement regions here!!!
									write(t=t, ln="# " + termStr)
									write(t=t, ln="for vertex in grid.vertices:")

									add(t=t+1, i=f"vertex.handle{idxStr1}", j=f"vertex.handle{idxStr2}", val=coeffStr, nl=2)
			writeMatrixScalarFieldVolumeIntegral(t+1)

			def writeMatrixGradScalarVolumeIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.lhs:
						if term.__class__ == VolumetricIntegral:
							grad, var = None, None							
							if term.arg.__class__ == Gradient:
								grad = term.arg
							elif term.arg.__class__ == Multiplication:
								gradTerms = [arg for arg in term.arg.args if arg.__class__ == Gradient]
								grad = gradTerms[0] if gradTerms else None
							if grad:
								fields = getTermFields(grad)
								var = fields[0] if (fields and fields[0].__class__ == ScalarField) else None

							if grad and var:
								arrIdx1 = self.arranjementDict["eq"][eqIdx]
								idxStr1 = f"+{arrIdx1}" if arrIdx1!="0" else ""

								arrIdx2 = self.arranjementDict["var"][var.name]
								idxStr2 = f"+({arrIdx2})*numberOfVertices" if arrIdx2!="0" else ""

								coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad]) if term.arg.__class__ == Multiplication else "1"
								coeffStr = str(coeff)

								write(t=t+0, ln="# " + str(term.arg))
								write(t=t+0, ln="for element in grid.elements:")
								write(t=t+1, ln="for innerFace in element.innerFaces:")
								write(t=t+2, ln="m = element.vertices.size")
								write(t=t+2, ln="transposedVoigtArea = getTransposedVoigtArea(innerFace)")
								write(t=t+2, ln="shapeFunctions = innerFace.getShapeFunctions()")
								write(t=t+2, ln="identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)]) if dimension == 3 else np.array([shapeFunctions, shapeFunctions, np.zeros(m)])", nl=2)
								write(t=t+2, ln=f"matrixCoefficients = {coeffStr} * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)", nl=2)
								write(t=t+2, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()")
								write(t=t+2, ln="for coord in range(dimension):")
								write(t=t+3, ln="for local, vertex in enumerate(element.vertices):")
								add(t=t+4, i=f"backwardsHandle+numberOfVertices*(coord{idxStr1})", j=f"vertex.handle{idxStr2}", val="matrixCoefficients[coord][local]")
								add(t=t+4, i=f"forwardHandle+numberOfVertices*(coord{idxStr1})", j=f"vertex.handle{idxStr2}", val="-matrixCoefficients[coord][local]", nl=2)
			writeMatrixGradScalarVolumeIntegral(t+1)

			def writeMatrixGradScalarSurfaceIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.lhs:
						if term.__class__ == SurfaceIntegral:
							grad, var = None, None							
							if term.arg.__class__ == Gradient:
								grad = term.arg
							elif term.arg.__class__ == Multiplication:
								gradTerms = [arg for arg in term.arg.args if arg.__class__ == Gradient]
								grad = gradTerms[0] if gradTerms else None
							if grad:
								fields = getTermFields(grad)
								var = fields[0] if (fields and fields[0].__class__ == ScalarField) else None

							if grad and var:
								arrIdx1 = self.arranjementDict["eq"][eqIdx]
								idxStr1 = f"+({arrIdx1})*numberOfVertices" if arrIdx1!="0" else ""

								arrIdx2 = self.arranjementDict["var"][var.name]
								idxStr2 = f"+({arrIdx2})*numberOfVertices" if arrIdx2!="0" else ""

								coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad]) if term.arg.__class__ == Multiplication else "1"
								coeffStr = str(coeff)

								write(t=t+0, ln="# " + str(term.arg))
								write(t=t+0, ln="for element in grid.elements:")
								write(t=t+1, ln="for innerFace in element.innerFaces:")
								write(t=t+2, ln="area = innerFace.area.getCoordinates()[:dimension]")
								write(t=t+2, ln=f"matrixCoefficients = {coeffStr} * np.matmul( area.T, innerFace.globalDerivatives )")
								write(t=t+2, ln="backwardVertexHandle, forwardVertexHandle = innerFace.getNeighborVerticesHandles()")
								write(t=t+2, ln="for local, vertex in enumerate(element.vertices):")
								add(t=t+3, i=f"backwardVertexHandle{idxStr1}", j=f"vertex.handle{idxStr2}", val="matrixCoefficients[local]")
								add(t=t+3, i=f"forwardVertexHandle{idxStr1}", j=f"vertex.handle{idxStr2}", val="-matrixCoefficients[local]", nl=2)
			writeMatrixGradScalarSurfaceIntegral(t+1)

			def writeMatrixSymmetricGradVecSurfaceIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.lhs:
						if term.__class__ == SurfaceIntegral:
							grad, field = None, None
							if term.arg.__class__ == SymmetricGradient:
								grad = term.arg
							elif term.arg.__class__ == Multiplication:
								gradTerms = [arg for arg in term.arg.args if arg.__class__ == SymmetricGradient]
								grad = gradTerms[0] if gradTerms else None

							if grad:
								fields = getTermFields(grad)
								vecField = fields[0] if (fields and fields[0].__class__ == VectorField) else None
								vecField = vecField if (grad.arg == vecField or (grad.arg.__class__ == Multiplication and vecField in grad.arg.args)) else None
							
							if grad and vecField:
								coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad])
								coeffStr = str(coeff)

								arrIdx1 = self.arranjementDict["eq"][eqIdx]
								idxStr1 = "" if arrIdx1=="0" else f"+{arrIdx1}"

								arrIdx2 = self.arranjementDict["var"][vecField.name+"_x"]
								idxStr2 = "" if arrIdx2=="0" else f"+{arrIdx2}"

								write(t=t, ln=f"# {term.arg}")
								write(t=t, ln="for element in grid.elements:")
								write(t=t+1, ln="for innerFace in element.innerFaces:")
								write(t=t+2, ln="transposedVoigtArea = getTransposedVoigtArea(innerFace)")
								write(t=t+2, ln="voigtGradientOperator = getVoigtGradientOperator(innerFace.globalDerivatives)")
								write(t=t+2, ln=f"coeff = {coeffStr}")
								write(t=t+2, ln="matrixCoefficients = np.einsum('ij,jk,kmn->imn', transposedVoigtArea, coeff, voigtGradientOperator)", nl=2)
								write(t=t+2, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()")
								write(t=t+2, ln="for local, vertex in enumerate(element.vertices):")
								write(t=t+3, ln="for i in range(dimension):")
								write(t=t+4, ln="for j in range(dimension):")								
								add(t=t+5, i=f"backwardsHandle + numberOfVertices*(i{idxStr1})", j=f"vertex.handle + numberOfVertices*(j{idxStr2})", val=f" matrixCoefficients[i][j][local]")
								add(t=t+5, i=f"forwardHandle   + numberOfVertices*(i{idxStr1})", j=f"vertex.handle + numberOfVertices*(j{idxStr2})", val=f"-matrixCoefficients[i][j][local]", nl=2)
			writeMatrixSymmetricGradVecSurfaceIntegral(t+1)

			def writeMatrixVecFieldSurfaceIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.lhs:
						if term.__class__ == SurfaceIntegral:
							for field in getTermFields(term):
								if field.__class__ == VectorField and (field == term.arg or ( term.arg.__class__ == Multiplication and field in term.arg.args )):
									coeff = Multiplication(*[arg for arg in term.arg.args if arg != field]) if (term.arg.__class__ == Multiplication and len(term.arg.args) > 1) else "1"
									coeffStr = str(coeff)

									arrIdx1 = self.arranjementDict["eq"][eqIdx]
									idxStr1 = f"+numberOfVertices*({arrIdx1})" if arrIdx1 != "0" else ""
									idxStr2 = self.arranjementDict["var"][field.name+"_x"]

									write(t=t+0, ln="# " + str(term.arg))
									write(t=t+0, ln="for element in grid.elements:")
									write(t=t+1, ln="for innerFace in element.innerFaces:")
									write(t=t+2, ln="area = innerFace.area.getCoordinates()[:dimension]")
									write(t=t+2, ln="shapeFunctions = innerFace.getShapeFunctions()", nl=2)
									write(t=t+2, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()")
									write(t=t+2, ln="for coord in range(dimension):")
									write(t=t+3, ln="for local, vertex in enumerate(element.vertices):")
									add(t=t+4, i=f"backwardsHandle{idxStr1}", j=f"vertex.handle + numberOfVertices*(coord+{idxStr2})", val=f"{coeffStr} * shapeFunctions[local] * area[coord]")
									add(t=t+4, i=f"forwardHandle{idxStr1}", j=f"vertex.handle + numberOfVertices*(coord+{idxStr2})", val=f"-{coeffStr} * shapeFunctions[local] * area[coord]", nl=2)
			writeMatrixVecFieldSurfaceIntegral(t+1)

			def writeMatrixDirichletBoundaryConditions(t):
				write(t=t+0, ln="# Dirichlet Boundary Conditions")
				for variableName in self.variables:
					idxStr = f"+({self.arranjementDict['var'][variableName]})*numberOfVertices" if len(self.variables) > 1 else ""
					zVar = (variableName[-2:] == "_z")
					a = (1 if zVar else 0)
					if zVar: write(t=t, ln="if dimension == 3:")
					write(t=t+a+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
					write(t=t+a+1, ln="for vertex in bCondition.boundary.vertices:")
					set(t=t+a+2, i=f"vertex.handle{idxStr}", nl=2)
			writeMatrixDirichletBoundaryConditions(t+1)

			def inverseMatrix(t):
				write(t=t, ln="# Inverse Matrix")
				if self.sparse:
					write(t=t, ln=f"matrix = sparse.csc_matrix( (matrixVals, zip(*matrixCoords)), shape=({self.sizeStr}, {self.sizeStr}) )")
					write(t=t, ln="inverseMatrix = sparse.linalg.inv( matrix )", nl=2)
				else:
					write(t=t, ln="inverseMatrix = np.linalg.inv(matrix)")
				write(t=t, ln="return inverseMatrix", nl=2)

			inverseMatrix(t+1)
		assembleMatrix(1)

		def assembleIndependent(t):
			write(t=t, ln="def assembleIndependent():")
			write(t=t+1, ln=f"independent = np.zeros({self.sizeStr})", nl=2)

			def writeIndependentVolumeIntegral(t):
				# Tem muita coisa pra revisar aqui
				for eqIdx, (equation, discretizedEquation) in enumerate(zip(self.equations, self.discretizedEquations)):
					for term in discretizedEquation.rhs:
						if term.__class__ == VolumetricIntegral and (term.arg.__class__ == Multiplication or hasSubclass(term.arg, Variable)):
							termStr = str(term.arg)
							coeffStr = "vertex.volume * " + termStr

							for field in getTermFields(term):
								if "_old" in field.name:
									fieldName = field.name.replace("_old", "").capitalize()
									coeffStr = coeffStr.replace(f" {field.name}", f" old{fieldName}Field[vertex.handle]")
								else:
									coeffStr = coeffStr.replace(f" {field.name}", f" {field.name}Field[vertex.handle]")

							if equation.order == 0:
								idxStr = f"+({self.arranjementDict['eq'][eqIdx]})*numberOfVertices" if len(self.variables) > 1 else ""
								write(t=t, ln="# " + termStr)
								write(t=t, ln="for vertex in grid.vertices:")
								write(t=t+1, ln=f"independent[vertex.handle{idxStr}] += " + coeffStr, nl=2)
							elif equation.order == 1:
								arrIdx = self.arranjementDict["eq"][eqIdx]
								idxStr = "+" + (f"(coord+{arrIdx})" if arrIdx!="0" else "coord") + "*numberOfVertices"
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

			def writeIndependentVecConstSurfaceIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.rhs:
						if term.__class__ == SurfaceIntegral and not getTermFields(term):
							coeffStr = str(term.arg)

							arrIdx = self.arranjementDict["eq"][eqIdx]
							idxStr = f"+numberOfVertices*({arrIdx})" if arrIdx != "0" else ""

							write(t=t+0, ln=f"# {coeffStr}")
							write(t=t+0, ln="for element in grid.elements:")
							write(t=t+1, ln="for innerFace in element.innerFaces:")
							write(t=t+2, ln="area = innerFace.area.getCoordinates()[:dimension]")
							write(t=t+2, ln=f"coefficient = {coeffStr}")
							write(t=t+2, ln="coefficient = np.matmul(area.T, coefficient)", nl=2)
							write(t=t+2, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()")
							write(t=t+2, ln=f"independent[backwardsHandle{idxStr}] += coefficient")
							write(t=t+2, ln=f"independent[forwardHandle{idxStr}]   -= coefficient", nl=2)
			writeIndependentVecConstSurfaceIntegral(t+1)

			def writeIndependentVecFieldSurfaceIntegral(t):
				for eqIdx, discretizedEquation in enumerate(self.discretizedEquations):
					for term in discretizedEquation.rhs:
						if term.__class__ == SurfaceIntegral:
							for field in getTermFields(term):
								if field.__class__ == VectorField and (field == term.arg or ( term.arg.__class__ == Multiplication and field in term.arg.args )):
									coeff = Multiplication(*[arg for arg in term.arg.args if arg != field]) if (term.arg.__class__ == Multiplication and len(term.arg.args) > 1) else "1"
									coeffStr = str(coeff)

									fieldName = f"old{field.name.replace('_old','').capitalize()}Field"

									arrIdx1 = self.arranjementDict["eq"][eqIdx]
									idxStr1 = f"+numberOfVertices*({arrIdx1})" if arrIdx1 != "0" else ""

									arrIdx2 = self.arranjementDict["var"][field.name.replace("_old","_x")]
									idxStr2 = f"+{arrIdx2}" if arrIdx2 != "0" else ""

									write(t=t+0, ln=f"# {term.arg}")
									write(t=t+0, ln="for element in grid.elements:")
									write(t=t+1, ln="for innerFace in element.innerFaces:")
									write(t=t+2, ln="area = innerFace.area.getCoordinates()[:dimension]")
									write(t=t+2, ln="shapeFunctions = innerFace.getShapeFunctions()", nl=2)
									write(t=t+2, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()")
									write(t=t+2, ln="for coord in range(dimension):")
									write(t=t+3, ln="for local, vertex in enumerate(element.vertices):")
									write(t=t+4, ln=f"independent[backwardsHandle{idxStr1}] += {coeffStr} * shapeFunctions[local] * area[coord] * {fieldName}[vertex.handle + numberOfVertices*(coord{idxStr2})]")
									write(t=t+4, ln=f"independent[forwardHandle{idxStr1}] -= {coeffStr} * shapeFunctions[local] * area[coord] * {fieldName}[vertex.handle + numberOfVertices*(coord{idxStr2})]", nl=2)


			writeIndependentVecFieldSurfaceIntegral(t+1)

			def writeNeumannBoundaryConditions(t):
				# This signal implies fluxes inwards
				plusVars = [ var.name for equation in self.equations for term in equation.lhs if hasTermWithSubclass(term, Divergence) for var in getTermFields(term) ]
				write(t=t+0, ln="# Neumann Boundary Condition")
				for variableName in self.variables:
					idxStr = f"+({self.arranjementDict['var'][variableName]})*numberOfVertices" if len(self.variables) > 1 else ""
					signal = "+" if variableName.split("_")[0] in plusVars else "-"
					zVar = (variableName[-2:] == "_z")
					a = (1 if zVar else 0)
					if zVar: write(t=t, ln="if dimension == 3:")
					write(t=t+a+0, ln=f"for bCondition in problemData.neumannBoundaries['{variableName}']:")
					write(t=t+a+1, ln="for facet in bCondition.boundary.facets:")
					write(t=t+a+2, ln="for outerFace in facet.outerFaces:")
					write(t=t+a+3, ln=f"independent[outerFace.vertex.handle{idxStr}] {signal}= bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())", nl=2)
			writeNeumannBoundaryConditions(t+1)

			def writeDirichletBoundaryConditions(t):
				write(t=t+0, ln="# Dirichlet Boundary Condition")
				for variableName in self.variables:
					idxStr = f"+({self.arranjementDict['var'][variableName]})*numberOfVertices" if len(self.variables) > 1 else ""
					zVar = (variableName[-2:] == "_z")
					a = (1 if zVar else 0)
					if zVar: write(t=t, ln="if dimension == 3:")
					write(t=t+a+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
					write(t=t+a+1, ln="for vertex in bCondition.boundary.vertices:")
					write(t=t+a+2, ln=f"independent[vertex.handle{idxStr}] = bCondition.getValue(vertex.handle)", nl=2)
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
				if self.sparse:
					write(t=t+1, ln=f"{self.variables[0]}Field = inverseMatrix.dot(independent)", nl=2)
				else:
					write(t=t+1, ln=f"{self.variables[0]}Field = np.matmul(inverseMatrix, independent)", nl=2)
			else:
				write("")
				if self.sparse:
					write(t=t+1, ln=f"results = inverseMatrix.dot(independent)")
				else:
					write(t=t+1, ln=f"results = np.matmul(inverseMatrix, independent)")
				# for i, variableName in enumerate(self.variables):
				# 	idxStr = self.arranjementDict["var"][variableName]
				# 	zVar = (variableName[-2:] == "_z")
				# 	a = (1 if zVar else 0)
				# 	if zVar: write(t=t+1, ln="if dimension == 3:")
				# 	write(t=t+a+1, ln=f"{variableName}Field = results[({idxStr})*numberOfVertices:({idxStr}+1)*numberOfVertices]")
				for variableName in self.varRadicals:
					if variableName in self.variables:
						idxStr = self.arranjementDict["var"][variableName]
						write(t=t+1, ln=f"{variableName}Field = results[({idxStr})*numberOfVertices:({idxStr}+1)*numberOfVertices]")
					else:
						idxStr = self.arranjementDict["var"][variableName+"_x"]
						write(t=t+1, ln=f"{variableName}Field = results[({idxStr})*numberOfVertices:({idxStr}+dimension)*numberOfVertices]")
				write("")
			fieldDifferences = ['max(abs('+varName+'Field - old'+varName.capitalize()+'Field))' for varName in self.transientFields]
			if len(fieldDifferences) == 1:
				write(t=t+1, ln=f"difference = {fieldDifferences[0]}", nl=2)
			else:
				write(t=t+1, ln=f"difference = max( {', '.join(fieldDifferences)} )", nl=2)
			for variableName in self.transientFields:
				write(t=t+1, ln=f"old{variableName.capitalize()}Field = {variableName}Field.copy()")
			write("")
			write(t=t+1, ln="currentTime += timeStep")
			write(t=t+1, ln="iteration += 1", nl=2)
			for variableName in self.variables:
				if "_" in variableName:
					zVar = (variableName[-2:] == "_z")
					a = (1 if zVar else 0)
					if zVar: write(t=t+1, ln="if dimension == 3:")
					i = "xyz".index(variableName.split("_")[1])
					write(t=t+a+1, ln=f"saver.save('{variableName}', {variableName.split('_')[0]}Field[{i}*numberOfVertices:{i+1}*numberOfVertices], currentTime)")
				else:
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
				if self.sparse:
					write(t=t, ln=f"{self.variables[0]}Field = inverseMatrix.dot(independent)")
				else:
					write(t=t, ln=f"{self.variables[0]}Field = np.matmul(inverseMatrix, independent)")
			else:
				if self.sparse:
					write(t=t, ln=f"results = inverseMatrix.dot(independent)")
				else:
					write(t=t, ln=f"results = np.matmul(inverseMatrix, independent)")
				for i, variableName in enumerate(self.variables):
					zVar = (variableName[-2:] == "_z")
					a = (1 if zVar else 0)
					if zVar: write(t=t, ln="if dimension == 3:")
					write(t=t+a, ln=f"{variableName}Field = results[{i}*numberOfVertices:{i+1}*numberOfVertices]")
			write(t=t, ln="")
			for variable in self.variables:
				zVar = (variableName[-2:] == "_z")
				a = (1 if zVar else 0)
				if zVar: write(t=t, ln="if dimension == 3:")
				write(t=t+a, ln=f"saver.save('{variable}', {variable}Field, 0.0)")
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
			write(t=t+2, ln=f"meshFilePath = '{self.meshPath}',")
			write(t=t+2, ln="outputFilePath = 'results',")
			write(t=t+2, ln=f"numericalSettings = PyEFVLib.NumericalSettings( timeStep = {self.timeStep}, tolerance = {self.tolerance}, maxNumberOfIterations = {self.maxNumberOfIterations} ),")
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
				if [bc for bc in self.boundaryConditions if bc.variableName==variableName]:
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

		self.text = self.text.replace("Δt", "timeStep")
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