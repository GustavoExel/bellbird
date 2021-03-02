class Variable:
	symbol = "var"
	def __init__(self, name=None):
		self.name = name

	def __repr__(self):
		# return f"{self.symbol}({self.name})"
		return self.name

	def setVarName(self, name):
		self.name = name

class Scalar(Variable):
	symbol = "scalar"

class Vector(Variable):
	symbol = "vec"

class Constant(Variable):
	symbol = "const"


varTypes = {
	"var":		Variable,
	"scalar": 	Scalar,
	"vector":	Vector,
	"const":	Constant,
}

timeDifferential 		= Scalar("Δt")
overTimeDifferential 	= Scalar("(1/Δt)")
minusOne				= Constant("(-1)")
zero					= Constant("0")

variablesDict = {
	"timeDifferential"		: timeDifferential,
	"overTimeDifferential"	: overTimeDifferential,
	"minusOne"				: minusOne,
	"zero"					: zero,
	"0"						: zero,
}