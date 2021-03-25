import os
from bellbird.equation import *
from bellbird.operators import *
from bellbird.writer import Writer

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
		self.writer = Writer(self)
		self.compiled = False

		self.parseEquations()
		self.applyEbFVM()

	def parseEquations(self):
		self.propertyNames = list(self.properties.values())[0].keys()
		definedNames = [ definition.split(" = ")[0] for definition in self.definitions ]
		self.definedVars = [ Constant(termStr) if not "[" in definition.split(" = ")[1] else (Vector(termStr) if not "[[" in definition.split(" = ")[1].replace(" ", "") else Matrix(termStr)) for termStr, definition in zip(definedNames, self.definitions)]

		self.equations = []
		for equationStr in self.equationsStr:
			equation = Equation(equationStr)
			equation.updatePropertyVars(self.propertyNames)
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
		self.matrixSize = f"({numberOfVarsStr})*numberOfVertices"

	def applyEbFVM(self):
		self.discretizedEquations = []
		for equationStr in self.equationsStr:
			discretizedEquation = Equation(equationStr)

			discretizedEquation.updatePropertyVars(self.propertyNames)
			discretizedEquation.updateVariables(self.variables)
			discretizedEquation.updateDefinitions(self.definedVars)
			discretizedEquation.integrateInSpace()				# a = b -> iiint(a) = iiint(b)
			discretizedEquation.applyDivergenceTheorem()		# iiint(div(f)) = iint(f)
			discretizedEquation.integrateInTime()				# d/dt(f) = (f-f_old)*(1/Δt)
			discretizedEquation.isolateVariables(self.variables) # x+a=y+c -> x-y=-a+c
			discretizedEquation.unapplyLinearProperty()

			self.discretizedEquations.append(discretizedEquation)

	def compile(self, filename="results.py"):
		self.writer.compile(filename)
		if self.writer.notImplementedTerms:
			print(f"Warning: it was not possible to implement the term(s) {', '.join(map(str,self.writer.notImplementedTerms))}, implement it manually in results.py ")

	def run(self):
		if not (self.writer.compiled and self.writer.fileName=="results.py"):
			self.compile()
		if not self.writer.notImplementedTerms:
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