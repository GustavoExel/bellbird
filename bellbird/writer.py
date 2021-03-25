from bellbird.equation import *
from bellbird.operators import *

class Writer:
	def __init__(self, model):
		self.model = model
		self.compiled = False

	def write(self, ln="", t=0, nl=1):
		self.text += t*'\t' + ln + nl*'\n'

	def compile(self, fileName):
		self.name = self.model.name.lower()[0] + "".join(self.model.name.split())[1:]
		self.fileName = fileName

		self.transient = True in [equation.isTransient() for equation in self.model.equations]
		self.transientFields = [ var.name for equation in self.model.equations for term in equation.terms if hasTermWithSubclass(term, TimeDerivative) for var in getTermFields(term) ]
		self.notImplementedTerms = [term for eq in self.model.discretizedEquations for term in eq.terms]

		self.text = ""

		self.writeHeader(t=0)
		self.writeProblemData(0)
		self.writeFields(1)
		self.declareMatrix(1)
		# self.declareProperties(1)
		# self.writeDefinitions(1)
		self.writeUtilFuncs(1)

		self.assembleMatrix(1)
		self.assembleIndependent(1)

		if self.transient:
			self.writeTransientSolverLoop(1)
		else:
			self.writePermanentSolverLoop(1)

		self.writeMainFunction(0)


		self.text = self.text.replace("Δt", "timeStep").replace("'","\"")
		with open(self.fileName, "w") as f:
			f.write(self.text)

		self.compiled = True

	# ------------------------------------------------------------------------------------------------------

	def assembleMatrix(self, t):
		def add(t, i, j, val, nl=1):
			if self.model.sparse:
				self.write(t=t, ln=f"add({i}, {j}, {val})", nl=nl)
			else:
				self.write(t=t, ln=f"matrix[{i}][{j}] += {val}", nl=nl)

		def set(t, i, nl):
			if self.model.sparse:
				self.write(t=t, ln=f"set({i})", nl=nl)
			else:
				self.write(t=t, ln=f"matrix[{i}] = np.zeros({self.model.matrixSize})")
				self.write(t=t, ln=f"matrix[{i}][{i}] = 1.0", nl=nl)

		def writeMatrixSurfaceIntegral(t, term, eqIdx, field, defs, val, vecEq, vecVar, ax, outerFace):
			self.notImplementedTerms.remove(term)

			arrIdx1 = self.model.arranjementDict["eq"][eqIdx]
			arrIdx2 = self.model.arranjementDict["var"][field.name+('_x' if vecVar else '')]
			idxStr1 = ('coord+' if vecEq and ax==0 else ('i+' if vecEq and ax==3 else ''))
			idxStr2 = 'coord+' if vecEq and ax==1 else ('j+' if vecEq and ax==3 else '')

			vT = (2 if ax==3 else 1) if vecEq else 0

			self.write(t=t+0, ln="# " + str(term.arg))
			self.write(t=t+0, ln=f"for face in element.{'faces' if outerFace else 'innerFaces'}:")
			self.write(t=t+1, ln="area = face.area.getCoordinates()[:dimension]")

			for _def in defs:
				self.write(t=t+1, ln=_def)

			if not outerFace:
				self.write(t=t+1, ln="backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()")
			
			if vecEq and ax==3:
				self.write(t=t+1, ln="for i in range(dimension):")
				self.write(t=t+2, ln="for j in range(dimension):")
			elif vecEq:
				self.write(t=t+1, ln="for coord in range(dimension):")


			self.write(t=t+1+vT, ln="for local, vertex in enumerate(element.vertices):")

			if not outerFace:
				add(t=t+2+vT, i=f"backwardsHandle+({idxStr1}{arrIdx1})*numberOfVertices", j=f"vertex.handle+({idxStr2}{arrIdx2})*numberOfVertices", val=f"{val}")
				add(t=t+2+vT, i=f"forwardHandle+({idxStr1}{arrIdx1})*numberOfVertices", j=f"vertex.handle+({idxStr2}{arrIdx2})*numberOfVertices", val=f"-{val}", nl=2)
			else:
				self.write(t=t+2+vT, ln="if type(face) == PyEFVLib.InnerFace:")
				self.write(t=t+3+vT, ln="backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()", nl=2)
				add(t=t+3+vT, i=f"backwardsHandle+({idxStr1}{arrIdx1})*numberOfVertices", j=f"vertex.handle+({idxStr2}{arrIdx2})*numberOfVertices", val=f"{val}")
				add(t=t+3+vT, i=f"forwardHandle+({idxStr1}{arrIdx1})*numberOfVertices", j=f"vertex.handle+({idxStr2}{arrIdx2})*numberOfVertices", val=f"-{val}", nl=2)
				self.write(t=t+2+vT, ln="elif type(face) == PyEFVLib.OuterFace:")
				add(t=t+3+vT, i=f"face.vertex.handle+({idxStr1}{arrIdx1})*numberOfVertices", j=f"vertex.handle+({idxStr2}{arrIdx2})*numberOfVertices", val=f"{val}", nl=2)

		self.write(t=t, ln="def assembleMatrix():")
		self.write(t=t+1, ln="for region in grid.regions:")
		self.declareProperties(t=t+2, region=True)
		self.writeDefinitions(t=t+2)
		self.write(t=t+2, ln="for element in region.elements:")

		def writeMatrixScalarFieldVolumeIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
				for term in discretizedEquation.lhs:
					if term.__class__ == VolumetricIntegral:
						for field in getTermFields(term):
							if field.__class__ == ScalarField and (field == term.arg or ( term.arg.__class__ == Multiplication and field in term.arg.args )):
								self.notImplementedTerms.remove(term)

								termStr = str(term.arg)
								# coeffStr = "vertex.volume * " + termStr.replace(f"* {field.name}", "")
								coeffStr = "vertex.getSubElementVolume(element) * " + termStr.replace(f"* {field.name}", "")

								arrIdx1 = self.model.arranjementDict["eq"][eqIdx]
								idxStr1 = f"+({arrIdx1})*numberOfVertices" if arrIdx1!="0" else ""

								arrIdx2 = self.model.arranjementDict["var"][field.name]
								idxStr2 = f"+({arrIdx2})*numberOfVertices" if arrIdx2!="0" else ""

								self.write(t=t, ln="# " + termStr)
								self.write(t=t, ln="for vertex in element.vertices:")

								add(t=t+1, i=f"vertex.handle{idxStr1}", j=f"vertex.handle{idxStr2}", val=coeffStr, nl=2)
		writeMatrixScalarFieldVolumeIntegral(t+3)

		def writeMatrixGradScalarVolumeIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
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
							coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad]) if term.arg.__class__ == Multiplication else "1"
							coeffStr = str(coeff)

							writeMatrixSurfaceIntegral(
								t = t, term = term, eqIdx = eqIdx, field = var,
								defs = [
									"m = element.vertices.size",
									"transposedVoigtArea = getTransposedVoigtArea(face)",
									"shapeFunctions = face.getShapeFunctions()",
									"identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)]) if dimension == 3 else np.array([shapeFunctions, shapeFunctions, np.zeros(m)])",
									f"matrixCoefficients = {coeffStr} * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)",
								],
								val = "matrixCoefficients[coord][local]",
								vecEq = True, vecVar = False, ax = 0, outerFace = False,
							)
		writeMatrixGradScalarVolumeIntegral(t+3)

		def writeMatrixGradScalarSurfaceIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
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
							coeff = Multiplication(*[arg for arg in term.arg.args if arg != grad]) if term.arg.__class__ == Multiplication else "1"
							coeffStr = str(coeff)

							writeMatrixSurfaceIntegral(
								t = t, term = term, eqIdx = eqIdx, field = var,
								defs = [f"matrixCoefficients = {coeffStr} * np.matmul( area.T, face.globalDerivatives )"],
								val = "matrixCoefficients[local]",
								vecEq = False, vecVar = False, ax = None, outerFace = False,
							)
		writeMatrixGradScalarSurfaceIntegral(t+3)

		def writeMatrixSymmetricGradVecSurfaceIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
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

							writeMatrixSurfaceIntegral(
								t = t, term = term, eqIdx = eqIdx, field = vecField,
								defs = [
									"transposedVoigtArea = getTransposedVoigtArea(face)",
									"voigtGradientOperator = getVoigtGradientOperator(face.globalDerivatives)",
									f"coeff = {coeffStr}",
									"matrixCoefficients = np.einsum('ij,jk,kmn->imn', transposedVoigtArea, coeff, voigtGradientOperator)",
								],
								val = "matrixCoefficients[i][j][local]",
								vecEq = True, vecVar = True, ax = 3, outerFace = False,
							)
		writeMatrixSymmetricGradVecSurfaceIntegral(t+3)

		def writeMatrixVecFieldSurfaceIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
				for term in discretizedEquation.lhs:
					if term.__class__ == SurfaceIntegral:
						for field in getTermFields(term):
							if field.__class__ == VectorField and (field == term.arg or ( term.arg.__class__ == Multiplication and field in term.arg.args )):
								coeff = Multiplication(*[arg for arg in term.arg.args if arg != field]) if (term.arg.__class__ == Multiplication and len(term.arg.args) > 1) else "1"
								coeffStr = str(coeff)

								writeMatrixSurfaceIntegral(
									t = t, term = term, eqIdx = eqIdx, field = field,
									defs = ["shapeFunctions = face.getShapeFunctions()"],
									val = f"{coeffStr} * shapeFunctions[local] * area[coord]",
									vecEq = True, vecVar = True, ax = 1, outerFace = True,
								)
		writeMatrixVecFieldSurfaceIntegral(t+3)

		def writeMatrixDirichletBoundaryConditions(t):
			self.write(t=t+0, ln="# Dirichlet Boundary Conditions")
			for variableName in self.model.variables:
				idxStr = f"+({self.model.arranjementDict['var'][variableName]})*numberOfVertices" if len(self.model.variables) > 1 else ""
				zVar = (variableName[-2:] == "_z")
				a = (1 if zVar else 0)
				if zVar: self.write(t=t, ln="if dimension == 3:")
				self.write(t=t+a+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
				self.write(t=t+a+1, ln="for vertex in bCondition.boundary.vertices:")
				set(t=t+a+2, i=f"vertex.handle{idxStr}", nl=2)
		writeMatrixDirichletBoundaryConditions(t+1)

		def inverseMatrix(t):
			self.write(t=t, ln="# Inverse Matrix")
			if self.model.sparse:
				self.write(t=t, ln=f"matrix = sparse.csc_matrix( (matrixVals, zip(*matrixCoords)), shape=({self.model.matrixSize}, {self.model.matrixSize}) )")
				self.write(t=t, ln="inverseMatrix = sparse.linalg.inv( matrix )", nl=2)
			else:
				self.write(t=t, ln="inverseMatrix = np.linalg.inv(matrix)")
			self.write(t=t, ln="return inverseMatrix", nl=2)

		inverseMatrix(t+1)

	def assembleIndependent(self, t):
		self.write(t=t, ln="def assembleIndependent():")
		self.write(t=t+1, ln=f"independent = np.zeros({self.model.matrixSize})", nl=2)

		self.write(t=t+1, ln="for region in grid.regions:")
		self.declareProperties(t=t+2, region=True)
		self.writeDefinitions(t=t+2)
		self.write(t=t+2, ln="for element in region.elements:")

		def writeIndependentVolumeIntegral(t):
			for eqIdx, (equation, discretizedEquation) in enumerate(zip(self.model.equations, self.model.discretizedEquations)):
				for term in discretizedEquation.rhs:
					if term.__class__ == VolumetricIntegral and (term.arg.__class__ == Multiplication or hasSubclass(term.arg, Variable)):
						self.notImplementedTerms.remove(term)

						termStr = str(term.arg)
						coeffStr = "vertex.getSubElementVolume(element) * " + termStr

						for field in getTermFields(term):
							if "_old" in field.name:
								fieldName = field.name.replace("_old", "").capitalize()
								coeffStr = coeffStr.replace(f" {field.name}", f" old{fieldName}Field[vertex.handle]")
							else:
								coeffStr = coeffStr.replace(f" {field.name}", f" {field.name}Field[vertex.handle]")

						if equation.order == 0:
							idxStr = f"+({self.model.arranjementDict['eq'][eqIdx]})*numberOfVertices" if len(self.model.variables) > 1 else ""
							self.write(t=t+0, ln="# " + termStr)
							self.write(t=t+0, ln="for vertex in element.vertices:")
							self.write(t=t+1, ln=f"independent[vertex.handle{idxStr}] += " + coeffStr, nl=2)
						elif equation.order == 1:
							arrIdx = self.model.arranjementDict["eq"][eqIdx]
							idxStr = "+" + (f"(coord+{arrIdx})" if arrIdx!="0" else "coord") + "*numberOfVertices"
							for arg in (term.arg.args if hasSubclass(term.arg, Operator) else term.arg):
								if hasSubclass(arg, Vector):
									coeffStr = coeffStr.replace(f" {arg}", f" {arg}[coord]")
							self.write(t=t+0, ln="# " + termStr)
							self.write(t=t+0, ln="for vertex in element.vertices:")
							self.write(t=t+1, ln="for coord in range(dimension):")							
							self.write(t=t+2, ln=f"independent[vertex.handle{idxStr}] += " + coeffStr, nl=2)

					elif term.__class__ == VolumetricIntegral:
						raise Exception(f"Warning: Term {term} is being ignored")
		writeIndependentVolumeIntegral(t+3)

		def writeIndependentVecConstSurfaceIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
				for term in discretizedEquation.rhs:
					if term.__class__ == SurfaceIntegral and not getTermFields(term):
						self.notImplementedTerms.remove(term)

						coeffStr = str(term.arg)

						arrIdx = self.model.arranjementDict["eq"][eqIdx]
						idxStr = f"+({arrIdx})*numberOfVertices" if arrIdx != "0" else ""

						self.write(t=t+0, ln=f"# {coeffStr}")
						self.write(t=t+0, ln="for innerFace in element.innerFaces:")
						self.write(t=t+1, ln="area = innerFace.area.getCoordinates()[:dimension]")
						self.write(t=t+1, ln=f"coefficient = {coeffStr}")
						self.write(t=t+1, ln="coefficient = np.matmul(area.T, coefficient)", nl=2)
						self.write(t=t+1, ln="backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()")
						self.write(t=t+1, ln=f"independent[backwardsHandle{idxStr}] += coefficient")
						self.write(t=t+1, ln=f"independent[forwardHandle{idxStr}]   -= coefficient", nl=2)
		writeIndependentVecConstSurfaceIntegral(t+3)

		def writeIndependentVecFieldSurfaceIntegral(t):
			for eqIdx, discretizedEquation in enumerate(self.model.discretizedEquations):
				for term in discretizedEquation.rhs:
					if term.__class__ == SurfaceIntegral:
						for field in getTermFields(term):
							if field.__class__ == VectorField and (field == term.arg or ( term.arg.__class__ == Multiplication and field in term.arg.args )):
								self.notImplementedTerms.remove(term)

								coeff = Multiplication(*[arg for arg in term.arg.args if arg != field]) if (term.arg.__class__ == Multiplication and len(term.arg.args) > 1) else "1"
								coeffStr = str(coeff)

								fieldName = f"old{field.name.replace('_old','').capitalize()}Field"

								arrIdx1 = self.model.arranjementDict["eq"][eqIdx]
								idxStr1 = f"+({arrIdx1})*numberOfVertices" if arrIdx1 != "0" else ""

								arrIdx2 = self.model.arranjementDict["var"][field.name.replace("_old","_x")]
								idxStr2 = f"+{arrIdx2}" if arrIdx2 != "0" else ""

								self.write(t=t+0, ln=f"# {term.arg}")
								self.write(t=t+0, ln="for face in element.faces:")
								self.write(t=t+1, ln="area = face.area.getCoordinates()[:dimension]")
								self.write(t=t+1, ln="shapeFunctions = face.getShapeFunctions()", nl=2)
								self.write(t=t+1, ln="for coord in range(dimension):")
								self.write(t=t+2, ln="for local, vertex in enumerate(element.vertices):")
								self.write(t=t+3, ln="if type(face) == PyEFVLib.InnerFace:")
								self.write(t=t+4, ln="backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()")
								self.write(t=t+4, ln=f"independent[backwardsHandle{idxStr1}] += {coeffStr} * shapeFunctions[local] * area[coord] * {fieldName}[vertex.handle + (coord{idxStr2})*numberOfVertices]")
								self.write(t=t+4, ln=f"independent[forwardHandle{idxStr1}] -= {coeffStr} * shapeFunctions[local] * area[coord] * {fieldName}[vertex.handle + (coord{idxStr2})*numberOfVertices]", nl=2)
								self.write(t=t+3, ln="elif type(face) == PyEFVLib.OuterFace:")
								self.write(t=t+4, ln=f"independent[face.vertex.handle{idxStr1}] += {coeffStr} * shapeFunctions[local] * area[coord] * {fieldName}[vertex.handle + (coord{idxStr2})*numberOfVertices]", nl=2)
		writeIndependentVecFieldSurfaceIntegral(t+3)

		def writeNeumannBoundaryConditions(t):
			# This signal implies fluxes inwards
			plusVars = [ var.name for equation in self.model.equations for term in equation.lhs if hasTermWithSubclass(term, Divergence) for var in getTermFields(term) ]
			self.write(t=t+0, ln="# Neumann Boundary Condition")
			for variableName in self.model.variables:
				idxStr = f"+({self.model.arranjementDict['var'][variableName]})*numberOfVertices" if len(self.model.variables) > 1 else ""
				signal = "+" if variableName.split("_")[0] in plusVars else "-"
				zVar = (variableName[-2:] == "_z")
				a = (1 if zVar else 0)
				if zVar: self.write(t=t, ln="if dimension == 3:")
				self.write(t=t+a+0, ln=f"for bCondition in problemData.neumannBoundaries['{variableName}']:")
				self.write(t=t+a+1, ln="for facet in bCondition.boundary.facets:")
				self.write(t=t+a+2, ln="for outerFace in facet.outerFaces:")
				self.write(t=t+a+3, ln=f"independent[outerFace.vertex.handle{idxStr}] {signal}= bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())", nl=2)
		writeNeumannBoundaryConditions(t+1)

		def writeDirichletBoundaryConditions(t):
			self.write(t=t+0, ln="# Dirichlet Boundary Condition")
			for variableName in self.model.variables:
				idxStr = f"+({self.model.arranjementDict['var'][variableName]})*numberOfVertices" if len(self.model.variables) > 1 else ""
				zVar = (variableName[-2:] == "_z")
				a = (1 if zVar else 0)
				if zVar: self.write(t=t, ln="if dimension == 3:")
				self.write(t=t+a+0, ln=f"for bCondition in problemData.dirichletBoundaries['{variableName}']:")
				self.write(t=t+a+1, ln="for vertex in bCondition.boundary.vertices:")
				self.write(t=t+a+2, ln=f"independent[vertex.handle{idxStr}] = bCondition.getValue(vertex.handle)", nl=2)
		writeDirichletBoundaryConditions(t+1)

		self.write(t=t+1, ln="return independent", nl=2)
	
	# ------------------------------------------------------------------------------------------------------

	def writeTransientSolverLoop(self, t):
		# Written very specifically for problems with only one variable
		self.write(t=t+0, ln="tolerance = problemData.tolerance")
		self.write(t=t+0, ln="difference = 2*tolerance")
		self.write(t=t+0, ln="iteration = 0")
		self.write(t=t+0, ln="currentTime = 0.0")
		self.write(t=t+0, ln="converged = False", nl=2)
		self.write(t=t+0, ln="inverseMatrix = assembleMatrix()", nl=2)
		self.write(t=t+0, ln="while not converged:")
		self.write(t=t+1, ln="independent = assembleIndependent()")
		if len(self.model.variables) == 1:
			if self.model.sparse:
				self.write(t=t+1, ln=f"{self.model.variables[0]}Field = inverseMatrix.dot(independent)", nl=2)
			else:
				self.write(t=t+1, ln=f"{self.model.variables[0]}Field = np.matmul(inverseMatrix, independent)", nl=2)
		else:
			self.write("")
			if self.model.sparse:
				self.write(t=t+1, ln=f"results = inverseMatrix.dot(independent)")
			else:
				self.write(t=t+1, ln=f"results = np.matmul(inverseMatrix, independent)")
			# for i, variableName in enumerate(self.model.variables):
			# 	idxStr = self.model.arranjementDict["var"][variableName]
			# 	zVar = (variableName[-2:] == "_z")
			# 	a = (1 if zVar else 0)
			# 	if zVar: self.write(t=t+1, ln="if dimension == 3:")
			# 	self.write(t=t+a+1, ln=f"{variableName}Field = results[({idxStr})*numberOfVertices:({idxStr}+1)*numberOfVertices]")
			for variableName in self.model.varRadicals:
				if variableName in self.model.variables:
					idxStr = self.model.arranjementDict["var"][variableName]
					self.write(t=t+1, ln=f"{variableName}Field = results[({idxStr})*numberOfVertices:({idxStr}+1)*numberOfVertices]")
				else:
					idxStr = self.model.arranjementDict["var"][variableName+"_x"]
					self.write(t=t+1, ln=f"{variableName}Field = results[({idxStr})*numberOfVertices:({idxStr}+dimension)*numberOfVertices]")
			self.write("")
		fieldDifferences = ['max(abs('+varName+'Field - old'+varName.capitalize()+'Field))' for varName in self.transientFields]
		if len(fieldDifferences) == 1:
			self.write(t=t+1, ln=f"difference = {fieldDifferences[0]}", nl=2)
		else:
			self.write(t=t+1, ln=f"difference = max( {', '.join(fieldDifferences)} )", nl=2)
		for variableName in self.transientFields:
			self.write(t=t+1, ln=f"old{variableName.capitalize()}Field = {variableName}Field.copy()")
		self.write("")
		self.write(t=t+1, ln="currentTime += timeStep")
		self.write(t=t+1, ln="iteration += 1", nl=2)
		for variableName in self.model.variables:
			if "_" in variableName:
				zVar = (variableName[-2:] == "_z")
				a = (1 if zVar else 0)
				if zVar: self.write(t=t+1, ln="if dimension == 3:")
				i = "xyz".index(variableName.split("_")[1])
				self.write(t=t+a+1, ln=f"saver.save('{variableName}', {variableName.split('_')[0]}Field[{i}*numberOfVertices:{i+1}*numberOfVertices], currentTime)")
			else:
				self.write(t=t+1, ln=f"saver.save('{variableName}', {variableName}Field, currentTime)")
		self.write("")
		self.write(t=t+1, ln="print('{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}'.format(iteration, currentTime, timeStep, difference))")
		self.write(t=t+1, ln="converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )", nl=2)
		self.write(t=t+1, ln="if iteration >= problemData.maxNumberOfIterations:")
		self.write(t=t+2, ln="break", nl=2)
		self.write(t=t+0, ln="saver.finalize()", nl=2)
	
	def writePermanentSolverLoop(self, t):
		# Written very specifically for problems with only one variable
		self.write(t=t, ln="inverseMatrix = assembleMatrix()")
		self.write(t=t, ln="independent = assembleIndependent()", nl=2)
		if len(self.model.variables) == 1:
			if self.model.sparse:
				self.write(t=t, ln=f"{self.model.variables[0]}Field = inverseMatrix.dot(independent)")
			else:
				self.write(t=t, ln=f"{self.model.variables[0]}Field = np.matmul(inverseMatrix, independent)")
		else:
			if self.model.sparse:
				self.write(t=t, ln=f"results = inverseMatrix.dot(independent)")
			else:
				self.write(t=t, ln=f"results = np.matmul(inverseMatrix, independent)")
			for i, variableName in enumerate(self.model.variables):
				zVar = (variableName[-2:] == "_z")
				a = (1 if zVar else 0)
				if zVar: self.write(t=t, ln="if dimension == 3:")
				self.write(t=t+a, ln=f"{variableName}Field = results[{i}*numberOfVertices:{i+1}*numberOfVertices]")
		self.write(t=t, ln="")
		for variableName in self.model.variables:
			zVar = (variableName[-2:] == "_z")
			a = (1 if zVar else 0)
			if zVar: self.write(t=t, ln="if dimension == 3:")
			self.write(t=t+a, ln=f"saver.save('{variableName}', {variableName}Field, 0.0)")
		self.write(t=t, ln="saver.finalize()", nl=2)

	def writeMainFunction(self, t):
		# def getRegionNames():
		# 	try:
		# 		import sys,os
		# 		sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'PyEFVLib'))
		# 		import PyEFVLib
		# 		return PyEFVLib.read(self.model.meshPath).gridData.regionsNames
		# 	except:
		# 		return ["Body"]
		# regionNames = getRegionNames()
		self.write(t=t+0, ln="def main():")
		self.write(t=t+1, ln="problemData = PyEFVLib.ProblemData(")
		self.write(t=t+2, ln=f"meshFilePath = '{self.model.meshPath}',")
		self.write(t=t+2, ln="outputFilePath = 'results',")
		self.write(t=t+2, ln=f"numericalSettings = PyEFVLib.NumericalSettings( timeStep = {self.model.timeStep}, tolerance = {self.model.tolerance}, maxNumberOfIterations = {self.model.maxNumberOfIterations} ),")
		self.write(t=t+2, ln="propertyData = PyEFVLib.PropertyData({")
		for regionName in self.model.properties:
			self.write(t=t+3, ln=f"'{regionName}':")
			self.write(t=t+3, ln="{")
			for propertyName, propertyValue in self.model.properties[regionName].items():
				self.write(t=t+4, ln=f"'{propertyName}': {propertyValue},")
			self.write(t=t+3, ln="},")
		self.write(t=t+2, ln="}),")
		self.write(t=t+2, ln="boundaryConditions = PyEFVLib.BoundaryConditions({")
		for variableName in self.model.variables:
			if [bc for bc in self.model.boundaryConditions if bc.variableName==variableName]:
				self.write(t=t+3, ln=f"'{variableName}': {{")
				for bc in self.model.boundaryConditions:
					if bc.__class__.__name__ == "InitialCondition" and bc.variableName == variableName:
						self.write(t=t+4, ln=f"'InitialValue': {bc.value},")
					if bc.__class__.__name__ == "BoundaryCondition" and bc.variableName == variableName:
						self.write(t=t+4, ln=f"'{bc.boundaryName}': {{ 'condition' : PyEFVLib.{bc.condition}, 'type' : PyEFVLib.Constant, 'value' : {bc.value} }},")
				self.write(t=t+1, ln="\t\t},")
		self.write(t=t+2, ln="}),\n\t)\n")
		self.write(t=t+1, ln=f"{self.name}( problemData )", nl=2)
		self.write(t=t+0, ln="if __name__ == '__main__':")
		self.write(t=t+1, ln="main()")

	# ------------------------------------------------------------------------------------------------------

	def writeHeader(self, t):
		self.write(t=t, ln="import sys,os")
		self.write(t=t, ln="sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'PyEFVLib'))")
		self.write(t=t, ln="import PyEFVLib")
		self.write(t=t, ln="import numpy as np")
		if self.model.sparse:
			self.write(t=t, ln="from scipy import sparse")
			self.write(t=t, ln="import scipy.sparse.linalg")
		self.write("")

	def writeProblemData(self, t):
		self.write(t=t,   ln=f"def {self.name}(problemData):")
		self.write(t=t+1, ln="propertyData 	 = problemData.propertyData")
		if self.transient:
			self.write(t=t+1, ln="timeStep 		 = problemData.timeStep")
		self.write(t=t+1, ln="grid 			 = problemData.grid")
		self.write(t=t+1, ln="numberOfVertices = grid.numberOfVertices")
		self.write(t=t+1, ln="dimension 		 = grid.dimension", nl=2)
		self.write(t=t+1, ln="saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension='xdmf')\n")

	def writeFields(self, t):
		for fieldName in self.model.varRadicals:
			self.write(t=t, ln=f"{fieldName}Field    = np.repeat(0.0, {'' if fieldName in self.model.variables else 'dimension*'}numberOfVertices)")
			if fieldName.split("_")[0] in self.transientFields:
				if fieldName in self.model.variables:
					self.write(t=t, ln=f"old{fieldName.capitalize()}Field = problemData.initialValues['{fieldName}'].copy()", nl=2)
				else:
					self.write(t=t, ln=f"old{fieldName.capitalize()}Field = np.concatenate((problemData.initialValues['{fieldName}_x'], problemData.initialValues['{fieldName}_y'], problemData.initialValues['{fieldName}_z'])) if dimension==3 else np.concatenate((problemData.initialValues['{fieldName}_x'], problemData.initialValues['{fieldName}_y']))", nl=2)

	def declareMatrix(self, t):
		if self.model.sparse:
			self.write(t=t, ln="matrixCoords = []")
			self.write(t=t, ln="matrixVals = []", nl=2)
		else:
			self.write(t=t, ln=f"matrix 		= np.zeros(({self.model.matrixSize}, {self.model.matrixSize}))", nl=2)
		# self.write(t=t, ln=f"independent = np.zeros({self.model.matrixSize})", nl=2)

	def declareProperties(self, t, region=False):
		for propertyName in self.model.propertyNames:
			self.write(t=t, ln=f"{propertyName:<10} = propertyData.get({'region.handle' if region else '0'}, '{propertyName}')")
		self.write("")

	def writeDefinitions(self, t):
		for definition in self.model.definitions:
			self.write(t=t, ln=definition)
		self.write("")

	def writeUtilFuncs(self, t):
		def writeSparseAddFunc(t):
			self.write(t=t+0, ln="def add(i, j, val):")
			self.write(t=t+1, ln="matrixCoords.append((i, j))")
			self.write(t=t+1, ln="matrixVals.append(val)", nl=2)
		def writeSparseSetFunc(t):
			self.write(t=t+0, ln="def set(i):")
			self.write(t=t+1, ln="_matrixCoords = [coord for coord in matrixCoords if coord[0]!=i] + [(i,i)]")
			self.write(t=t+1, ln="_matrixVals   = [val for coord,val in zip(matrixCoords,matrixVals) if coord[0]!=i] + [1.0]")
			self.write(t=t+1, ln="matrixCoords.clear(); matrixCoords.extend(_matrixCoords)")
			self.write(t=t+1, ln="matrixVals.clear(); matrixVals.extend(_matrixVals)", nl=2)
		if self.model.sparse:
			writeSparseAddFunc(t)
			writeSparseSetFunc(t)

		def writeSymmetricUtilFuncs(t):
			self.write(t=t+0, ln="def getTransposedVoigtArea(innerFace):")
			self.write(t=t+1, ln="Sx, Sy, Sz = innerFace.area.getCoordinates()")
			self.write(t=t+1, ln="return np.array([[Sx,0,Sy],[0,Sy,Sx]]) if dimension==2 else np.array([[Sx,0,0,Sy,0,Sz],[0,Sy,0,Sx,Sz,0],[0,0,Sz,0,Sy,Sx]])", nl=2)
			self.write(t=t+0, ln="def getVoigtGradientOperator(globalDerivatives):")
			self.write(t=t+1, ln="if len(globalDerivatives) == 2:")
			self.write(t=t+2, ln="Nx,Ny = globalDerivatives")
			self.write(t=t+2, ln="zero=np.zeros(Nx.size)")
			self.write(t=t+2, ln="return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])", nl=2)
			self.write(t=t+1, ln="if len(globalDerivatives) == 3:")
			self.write(t=t+2, ln="Nx,Ny,Nz = globalDerivatives")
			self.write(t=t+2, ln="zero=np.zeros(Nx.size)")
			self.write(t=t+2, ln="return np.array([[Nx,zero,zero],[zero,Ny,zero],[zero,zero,Nz],[Ny,Nx,zero],[zero,Nz,Ny],[Nz,zero,Nx]])", nl=2)
		if True in [hasTermWithSubclass(equation.term, SymmetricGradient) for equation in self.model.equations]:
			writeSymmetricUtilFuncs(t)
