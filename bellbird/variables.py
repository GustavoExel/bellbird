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

class Matrix(Variable):
	symbol = "mat"

class Tensor(Variable):
	symbol = "tens"

class Constant(Variable):
	symbol = "const"

class Field(Variable):
	symbol = "field"

class ScalarField(Scalar, Field):
	symbol = "scalarF"

class VectorField(Vector, Field):
	symbol = "vecF"

varTypes = {
	"var"	:	Variable,
	"scalar": 	Scalar,
	"vec"	:	Vector,
	"mat"	:	Matrix,
	"tens"	:	Tensor,
	"const"	:	Constant,
	"field"	:	Field,
	"scalarF":	ScalarField,
	"vecF"	:	VectorField,
}

timeDifferential 		= Scalar("Δt")
overTimeDifferential 	= Scalar("(1/Δt)")
minusOne				= Constant("(-1)")
zero					= Constant("0")
one 					= Constant("1")
zeroVec					= Vector("0")


variablesDict = {
	"timeDifferential"		: timeDifferential,
	"dt"					: timeDifferential,
	"Δt"					: timeDifferential,
	"overTimeDifferential"	: overTimeDifferential,
	"(1/Δt)"				: overTimeDifferential,
	"(1/dt)"				: overTimeDifferential,
	"minusOne"				: minusOne,
	"(-1)"					: minusOne,
	"zero"					: zero,
	"0"						: zero,
	"one"					: one,
	"1"						: one,
	"zeroVec"				: zeroVec,
}