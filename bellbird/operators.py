# --------- Operators ---------
class Operator:
	def __init__(self, *args):
		self.args = list(args)

	def setArgs(self, args):
		self.args = args

	def setArg(self, arg):
		self.args = [ arg ]

class Function(Operator):
	def __repr__(self):
		return f"{self.func}({self.args[0]})"

	@property
	def arg(self):
		return self.args[0]
	

class BinaryOperator(Operator):
	def __repr__(self):
		return f" {self.symbol} ".join(map(str, self.args))
		# return f"{self.func}({', '.join(map(str, self.args))})"


# --------- Functions --------- 
class Divergence(Function):
	func = "div"
	linear = True

class Gradient(Function):
	func = "grad"
	linear = True

class SymmetricGradient(Function):
	# (1/2) * ( grad(u) + grad(u)^T )
	func = "grad_s"
	linear = True

class TimeDerivative(Function):
	func = "d/dt"
	linear = True

class TimeIntegral(Function):
	func = "intT"
	linear = True

class VolumetricIntegral(Function):
	func = "iiint"
	linear = True

class SurfaceIntegral(Function):
	func = "iint"
	linear = True

class VolumetricSummatory(Function):
	func = "summ3D"
	linear = True

class SurfaceSummatory(Function):
	func = "summ2D"
	linear = True

# --------- BinaryOperator --------- 
class Equals(BinaryOperator):
	symbol = "="
	func   = "equals"
	commutative = True

class Sum(BinaryOperator):
	symbol = "+"
	func   = "sum"
	commutative = True

class Subtraction(BinaryOperator):
	symbol = "-"
	func   = "subtraction"
	commutative = False

class Multiplication(BinaryOperator):
	symbol = "*"
	func   = "multiplication"
	commutative = True

class Division(BinaryOperator):
	symbol = "/"
	func   = "division"
	commutative = False

class DotProduct(BinaryOperator):
	symbol = "·"
	func   = "dotProduct"
	commutative = True

class CrossProduct(BinaryOperator):
	symbol = "×"
	func   = "crossProduct"
	commutative = False


# --------- Dicts ---------

functions = {
	"div"	 : Divergence,
	"grad"	 : Gradient,
	"grad_s" : SymmetricGradient,
	"d/dt"	 : TimeDerivative,
	"intT"	 : TimeIntegral,
	"iiint"	 : VolumetricIntegral,
	"iint"	 : SurfaceIntegral,
	"summ3D" : VolumetricSummatory,
	"summ2D" : SurfaceSummatory,
}

binaryOperators = {
	"*"		: Multiplication,
	"/"		: Division,
	"·"		: DotProduct,
	"×"		: CrossProduct,
	"+"		: Sum,
	"-"		: Subtraction,
	"="		: Equals,	
}

linearOperators = [Sum, Subtraction, Equals]

other = ["(", ")"]

symbols = list(functions.keys()) + list(binaryOperators.keys()) + other