from bellbird.operators import *
from bellbird.variables import *

class Equation:
	def __init__(self, equationStr):
		self.equationStr = equationStr
		self.term = parseEquationStr(self.equationStr)

	def __repr__(self):
		return str(self.term)

	@property
	def lhs(self):
		if self.term.args[0].__class__ == Sum:
			return self.term.args[0].args
		else:
			return [ self.term.args[0] ]

	@property
	def rhs(self):
		if self.term.args[1].__class__ == Sum:
			return self.term.args[1].args
		else:
			return [ self.term.args[1] ]

	@property
	def terms(self):
		return self.lhs + self.rhs

	def applyLinearProperty(self):
		applyLinearProperty(self.term)

	def applyDistributiveProperty(self):
		applyDistributiveProperty(self.term)

	def clusterOperations(self):
		clusterOperations(self.term)

	def unapplyLinearProperty(self):
		unapplyLinearProperty(self.term)
		self.clusterOperations()
		self.applyDistributiveProperty()

	def rewrite(self):
		self.applyDistributiveProperty()	# (a+b)*c = a*c + b*c
		self.clusterOperations()			# times(a, times(b,c)) = times(a, b, c)

		self.applyLinearProperty()			# f(a*(x+y)) = a*f(x) + a*f(y)

		self.applyDistributiveProperty()	# (a+b)*c = a*c + b*c
		self.clusterOperations()			# times(a, times(b,c)) = times(a, b, c)

		self.applyLinearProperty()			# f(a*(x+y)) = a*f(x) + a*f(y)

	def updatePropertyVars(self, properties):
		propertyVars = properties.keys()
		updatePropertyVars(self.term, propertyVars)

	def updateVariables(self, variableNames):
		updateVariables(self.term, variableNames)

	def updateDefinitions(self, definedVars):
		updateDefinitions(self.term, definedVars)

	def integrateInTime(self):
		integrateInTime(self.term)
		self.rewrite()

	def integrateInSpace(self):
		self.term.setArgs([ VolumetricIntegral(self.term.args[0]), VolumetricIntegral(self.term.args[1]) ])
		self.rewrite()

	def applyDivergenceTheorem(self):
		applyDivergenceTheorem(self.term)

	def isolateVariables(self, variables):
		lhsVars = [ term for term in self.lhs if 	 containsVariables(term, variables) ]
		lhsInds = [ term for term in self.lhs if not containsVariables(term, variables) ]
		rhsVars = [ term for term in self.rhs if 	 containsVariables(term, variables) ]
		rhsInds = [ term for term in self.rhs if not containsVariables(term, variables) ]

		self.term = Equals( 
			Sum( *lhsVars, Multiplication(minusOne, *rhsVars) ) if rhsVars else Sum(*lhsVars),
			Sum( *rhsInds, Multiplication(minusOne, *lhsInds) ) if lhsInds else Sum(*rhsInds)
		)

		self.rewrite()

	_order = None
	@property
	def order(self):
		if not self._order:
			self._order = self.getOrder()
		return self._order
	def getOrder(self):
		return getOrder(self.term)

	def isTransient(self):
		return isTransient(self.term)

def hasSubclass(obj, subclass):
	return subclass in obj.__class__.__mro__

def parseEquationStr(eqString):
	# Separate symbols
	for symbol in symbols:
		eqString = eqString.replace(symbol, f" {symbol} ")
	while "  " in eqString:
		eqString = eqString.replace("  ", " ")
	eqString = eqString.replace("d / dt", "d/dt")
	eqString = eqString.replace("grad _s", "grad_s")
	eqComponentsStr = eqString.split()

	# Separate parenthesis
	while "(" in eqComponentsStr:
		i0 = eqComponentsStr.index("(")
		i1 = [(eqComponentsStr[:i].count("(")-eqComponentsStr[:i].count(")") == 1) and (c == ")") for i,c in enumerate(eqComponentsStr)].index(1)
		eqComponentsStr = eqComponentsStr[:i0] + [ eqComponentsStr[i0+1:i1] ] + eqComponentsStr[i1+1:]
	eqComponents = []

	# Create objects
	for idx, compStr in enumerate(eqComponentsStr):
		# Parenthesis
		if type(compStr) == list:
			arg = parseEquationStr( " ".join(compStr) )
			eqComponents.append(arg)

		# Function
		elif compStr in functions:
			comp = functions[compStr]()
			eqComponents.append(comp)

		# Binary Operator
		elif compStr in binaryOperators:
			comp = binaryOperators[compStr]()
			eqComponents.append(comp)

		# Listed Variables
		elif compStr in variablesDict:
			comp = variablesDict[compStr]
			eqComponents.append(comp)

		# Variable
		elif compStr in varTypes:
			name = eqComponentsStr[idx+1][0]
			comp = varTypes[compStr](name)
			eqComponents.append(comp)
			eqComponentsStr.pop(idx+1)
		else:
			comp = Variable(compStr)
			eqComponents.append(comp)

	eqFunctions = [comp for idx, comp in enumerate(eqComponents) if hasSubclass(comp, Function) and not comp.args]
	eqOperators = [comp for idx, comp in enumerate(eqComponents) if hasSubclass(comp, BinaryOperator) and not comp.args]

	eqFunctions.sort( key=lambda func:list(functions.values()).index(func.__class__) )
	eqOperators.sort( key=lambda op:list(binaryOperators.values()).index(op.__class__) )

	for function in eqFunctions:
		idx = eqComponents.index(function)
		arg = eqComponents[idx+1]
		function.setArg(arg)

		eqComponents.pop( idx+1 )

	for operator in eqOperators:
		idx = eqComponents.index(operator)
		args = [ eqComponents[idx-1], eqComponents[idx+1] ]
		operator.setArgs(args)

		eqComponents.pop( idx+1 )
		eqComponents.pop( idx-1 )

	return eqComponents[0]

def updatePropertyVars(term, propertyVars):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Variable):
			if arg.name in propertyVars:
				term.args[idx] = Constant(arg.name)

		if hasSubclass(arg, Operator):
			updatePropertyVars(arg, propertyVars)

def updateVariables(term, variableNames):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Variable):
			if arg.name in variableNames:
				term.args[idx] = ScalarField(arg.name)
			elif arg.name in [varName.replace("_x","") for varName in variableNames if "_x" in varName]:
				term.args[idx] = VectorField(arg.name)

		if hasSubclass(arg, Operator):
			updateVariables(arg, variableNames)	

def updateDefinitions(term, definedVars):
	definedNames = [ var.name for var in definedVars ]
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Variable):
			if arg.name in definedNames:
				term.args[idx] = definedVars[ definedNames.index(arg.name) ]

		if hasSubclass(arg, Operator):
			updateDefinitions(arg, definedVars)	

def applyLinearProperty(term, flag=1):
	for idx, arg in enumerate(term.args):
		function = arg.__class__
		if hasSubclass(arg, Function) and arg.linear:
			if hasSubclass(arg.arg, Multiplication):
				constants = [ mArg for mArg in arg.arg.args if hasSubclass(mArg, Constant) ]
				variables = [ mArg for mArg in arg.arg.args if not hasSubclass(mArg, Constant) ]
				varArg    = Multiplication(*variables) if len(variables) > 1 else ( variables[0] if len(variables) == 1 else None )
				
				if constants and varArg:
					term.args[idx] = Multiplication( *constants, function(varArg) )
				elif len(constants) > 1:
					term.args[idx] = Multiplication( *constants )
				elif len(constants) == 1:
					term.args[idx] = constants[0]

			if hasSubclass(arg.arg, Variable) and arg.arg in [zero, zeroVec]:
				term.args[idx] = arg.arg

			for operator in linearOperators:
				if hasSubclass(arg.arg, operator):
					args = arg.arg.args
					term.args[idx] = operator( *[function(arg) for arg in args] )
		if hasSubclass(arg, Operator):
			applyLinearProperty(arg, flag+1)
	if flag > 0:
		applyLinearProperty(term, flag-1)

def unapplyLinearProperty(term):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Multiplication):
			for mArg in arg.args:
				if hasSubclass(mArg, Function) and mArg.linear:
					function = mArg.__class__

					variables = [mmArg for mmArg in arg.args if not hasSubclass(mmArg, Constant) and mmArg != mArg]
					constants = [mmArg for mmArg in arg.args if hasSubclass(mmArg, Constant)]

					if constants and variables:
						term.args[idx] = Multiplication(*variables, function(Multiplication(*constants, *mArg.args)))
					elif constants and not variables:
						term.args[idx] = function(Multiplication(*constants, *mArg.args))
					break

		if hasSubclass(arg, Operator):
			unapplyLinearProperty(arg)

def applyDivergenceTheorem(term):
	for idx, arg in enumerate(term.args):
		if arg.__class__ == VolumetricIntegral and arg.arg.__class__ == Divergence:
			divArg = arg.arg.arg
			term.args[idx] = SurfaceIntegral(divArg)

		elif hasSubclass(arg, BinaryOperator):
			applyDivergenceTheorem(arg)

def integrateInTime(term):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, TimeDerivative) and hasSubclass(arg.arg, Variable):
			var = arg.arg
			varOld = var.__class__(f"{var.name}_old")
			term.args[idx] = Multiplication( Subtraction(var, varOld), overTimeDifferential )

		elif hasSubclass(arg, Operator):
			integrateInTime(arg)

def applyDistributiveProperty(term, flag=1):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Multiplication):
			for sArg in arg.args:
				if hasSubclass(sArg, Sum) or hasSubclass(sArg, Subtraction):
					plus_func = sArg.__class__
					args = [ mArg for mArg in arg.args if not mArg == sArg ]

					term.args[idx] = plus_func(*[Multiplication(*args, ssArg) for ssArg in sArg.args])

					# only one sArg
					break

		if hasSubclass(arg, Operator):
			applyDistributiveProperty(arg, flag+1)

	if flag > 0:
		applyDistributiveProperty(term, flag-1)

def clusterOperations(term, flag=5):
	for idx, arg in enumerate(term.args):
		# Multiplication(A, Multiplication(B, C)) -> Multiplication(A, B, C)
		if hasSubclass(arg, Multiplication):
			args = []
			for mArg in arg.args:
				if hasSubclass(mArg, Multiplication):
					args += mArg.args
				else:
					args.append(mArg)
			arg.args = args

		# Sum(A, Sum(B, C)) -> Sum(A, B, C)
		if hasSubclass(arg, Sum):
			args = []
			for sArg in arg.args:
				if hasSubclass(sArg, Sum):
					args += sArg.args
				else:
					args.append(sArg)
			arg.args = args

		# Subtraction(A, B) -> Sum(A, Multiplication(minusOne, B))
		if hasSubclass(arg, Subtraction):
			term.args[idx] = Sum(arg.args[0], Multiplication(minusOne, arg.args[1]))

		# Multiplication(minusOne, minusOne, A, B, ...) -> Multiplication(A, B, ...)
		if hasSubclass(arg, Multiplication):
			if arg.args.count(minusOne) >= 2:
				arg.args.remove(minusOne)
				arg.args.remove(minusOne)

		# Sum(0, A, B, ...) -> Sum(A, B, ...)
		if hasSubclass(arg, Sum):
			if zero in arg.args:
				arg.args.remove(zero)
			elif zeroVec in arg.args:
				arg.args.remove(zeroVec)

		# Multiplication(1, A, B, ...) -> Multiplication(A, B, ...)
		if hasSubclass(arg, Multiplication):
			if one in arg.args:
				arg.args.remove(one)

		# Multiplication(0, A, B, ...) -> 0
		if hasSubclass(arg, Multiplication):
			if zero in arg.args:
				term.args[idx] = zero

		# Sum() -> 0
		if hasSubclass(arg, Sum):
			if len(arg.args) == 0:
				term.args[idx] = zero

		# Multiplication() -> 1
		if hasSubclass(arg, Multiplication):
			if len(arg.args) == 0:
				term.args[idx] = one

		# binaryOperator(A) -> A
		if hasSubclass(arg, BinaryOperator):
			if len(arg.args) == 1:
				term.args[idx] = arg.args[0]

		# d/dX( const(A) ) -> 0
		if arg.__class__ in [TimeDerivative, Gradient]:
			if hasSubclass(arg.arg, Constant):
				term.args[idx] = zero

		# vec(0) -> zeroVec
		if hasSubclass(arg, Vector):
			if arg.name == "0":
				term.args[idx] = zeroVec

		if hasSubclass(arg, Operator):
			clusterOperations(arg, flag+1)

	if flag > 0:
		clusterOperations(term, flag-1)

def containsVariables(term, variables):
	if hasSubclass(term, Variable):
		return term.name in variables + [varName.replace("_x", "") for varName in variables if "_x" in varName]

	elif hasSubclass(term, Operator):
		for arg in term.args:
			if containsVariables(arg, variables):
				return True
		return False

def getTermFields(term):
	if hasSubclass(term, Field):
		return [term]
	
	fields = []
	if hasSubclass(term, Operator):
		for arg in term.args:
			fields += getTermFields(arg)

	return fields

def getOrder(term):
	if hasSubclass(term, Scalar):
		return 0
	elif hasSubclass(term, Vector):
		return 1
	elif hasSubclass(term, Matrix):
		return 2
	elif hasSubclass(term, Tensor):
		return 3
	elif hasSubclass(term, Function):
		if term.__class__ in [TimeDerivative, TimeIntegral, VolumetricIntegral, SurfaceIntegral, VolumetricSummatory, SurfaceSummatory]:
			return getOrder(term.arg)
		elif term.__class__ == Divergence:
			return getOrder(term.arg) - 1
		elif term.__class__ in [Gradient, SymmetricGradient]:
			return getOrder(term.arg) + 1

	elif hasSubclass(term, BinaryOperator):
		if term.__class__ in [Sum, Subtraction, Equals]:
			argsOrders = [ getOrder(arg) for arg in term.args ]
			order = argsOrders[0]
			if [argOrder for argOrder in argsOrders if argOrder != order ]:
				raise Exception(f"Invalid operation {term}")
			return order

		elif term.__class__ == Multiplication:
			argsOrders = [ getOrder(arg) for arg in term.args ]
			order = argsOrders[0]

			for argOrder in argsOrders[1:]:
				if order == 0:
					order = argOrder
				elif order == 1 and argOrder == 1:
					order = 0
				elif order == 2 and argOrder == 1:
					order = 1
				# else:
				# 	order = order
			return order

		elif term.__class__ == DotProduct:
			return sum([ getOrder(arg) for arg in term.args ]) - 2

		elif term.__class__ == CrossProduct:
			# May need revision
			return getOrder(term.args[0])

	# print("-----------------")
	# print(term)
	# exit()

def isTransient(term):
	if hasSubclass(term, TimeDerivative):
		return True

	elif hasSubclass(term, Operator):
		for arg in term.args:
			if isTransient(arg):
				return True

	return False

def hasTermWithSubclass(term, subclass):
	if hasSubclass(term, subclass):
		return True
	elif hasSubclass(term, Operator):
		for arg in term.args:
			if hasTermWithSubclass(arg, subclass):
				return True
	return False