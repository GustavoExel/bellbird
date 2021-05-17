from .variables import *
from .operators import *
from .equation import *
from .model import *
from .writer import *

Dirichlet = "DIRICHLET_BOUNDARY_CONDITION"
Neumann   = "NEUMANN_BOUNDARY_CONDITION"


class BoundaryCondition:
	bcDict = {"DIRICHLET_BOUNDARY_CONDITION":"Dirichlet", "NEUMANN_BOUNDARY_CONDITION":"Neumann"}
	conditionType = "BOUNDARY_CONDITION"
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