from bellbird.equation import Equation
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
		self.discretizedEquation.integrateInSpace()				# a = b -> iiint(a) = iiint(b)
		self.discretizedEquation.applyDivergenceTheorem()		# iiint(div(f)) = iint(f)
		self.discretizedEquation.integrateInTime()				# d/dt(f) = (f-f_old)*(1/Δt)
		self.discretizedEquation.isolateVariables(self.variables) # x+a=y+c -> x-y=-a+c

	def compile(self):
		