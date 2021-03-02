from bellbird.operators import *
from bellbird.variables import *

class Equation:
	def __init__(self, equationStr):
		self.equationStr = equationStr

		self.term = parseEquationStr(self.equationStr)[0]

	def __repr__(self):
		return str(self.term)

	def applyLinearProperty(self):
		applyLinearProperty(self.term)

	def applyDistributiveProperty(self):
		applyDistributiveProperty(self.term)

	def clusterOperations(self):
		clusterOperations(self.term)

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

	def integrateInTime(self):
		integrateInTime(self.term)
		self.rewrite()

	def integrateInSpace(self):
		self.term.setArgs([ VolumetricIntegral(self.term.args[0]), VolumetricIntegral(self.term.args[1]) ])
		self.rewrite()

	def applyDivergenceTheorem(self):
		applyDivergenceTheorem(self.term)

	def isolateVariables(self, variables):
		rhsTerms = self.term.args[0].args if self.term.args[0].__class__ == Sum else [ self.term.args[0] ]
		lhsTerms = self.term.args[1].args if self.term.args[1].__class__ == Sum else [ self.term.args[1] ]

		rhsVars = [ term for term in rhsTerms if 	 containsVariables(term, variables) ]
		rhsInds = [ term for term in rhsTerms if not containsVariables(term, variables) ]
		lhsVars = [ term for term in lhsTerms if 	 containsVariables(term, variables) ]
		lhsInds = [ term for term in lhsTerms if not containsVariables(term, variables) ]

		self.term = Equals( 
			Sum( *rhsVars, Multiplication(minusOne, *lhsVars) ) if lhsVars else Sum(*rhsVars),
			Sum( *lhsInds, Multiplication(minusOne, *rhsInds) ) if rhsInds else Sum(*lhsInds)
		)

		self.rewrite()

def hasSubclass(obj, subclass):
	return subclass in obj.__class__.__mro__

def parseEquationStr(eqString):
	# Separate symbols
	for symbol in symbols:
		eqString = eqString.replace(symbol, f" {symbol} ")
	while "  " in eqString:
		eqString = eqString.replace("  ", " ")
	eqString = eqString.replace("d / dt", "d/dt")
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
			# print(arg[0])
			eqComponents.append(arg[0])

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

	return eqComponents

def updatePropertyVars(term, propertyVars):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Variable):
			if arg.name in propertyVars:
				term.args[idx] = Constant(arg.name)

		if hasSubclass(arg, Operator):
			updatePropertyVars(arg, propertyVars)

def applyLinearProperty(term, flag=1):
	for idx, arg in enumerate(term.args):
		function = arg.__class__
		if hasSubclass(arg, Function) and arg.linear:
			if hasSubclass(arg.args[0], Multiplication):
				constants = [ mArg for mArg in arg.args[0].args if hasSubclass(mArg, Constant) ]
				variables = [ mArg for mArg in arg.args[0].args if not hasSubclass(mArg, Constant) ]
				varArg    = Multiplication(*variables) if len(variables) > 1 else ( variables[0] if len(variables) == 1 else None )
				
				if constants and varArg:
					term.args[idx] = Multiplication( *constants, function(varArg) )
				elif len(constants) > 1:
					term.args[idx] = Multiplication( *constants )
				elif len(constants) == 1:
					term.args[idx] = constants[0]

			if hasSubclass(arg.args[0], Constant) and arg.args[0].name == "0":
				term.args[idx] = zero

			for operator in linearOperators:
				if hasSubclass(arg.args[0], operator):
					args = arg.args[0].args
					# term.args[idx] = operator( function(args[0]), function(args[1]) )
					term.args[idx] = operator( *[function(arg) for arg in args] )
		elif hasSubclass(arg, Operator):
			applyLinearProperty(arg, flag+1)
	if flag > 0:
		applyLinearProperty(term, flag-1)

def applyDivergenceTheorem(term):
	for idx, arg in enumerate(term.args):
		if arg.__class__ == VolumetricIntegral and arg.args[0].__class__ == Divergence:
			divArg = arg.args[0].args[0]
			term.args[idx] = SurfaceIntegral(divArg)

		elif hasSubclass(arg, BinaryOperator):
			applyDivergenceTheorem(arg)

def integrateInTime(term):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, TimeDerivative) and hasSubclass(arg.args[0], Variable):
			var = arg.args[0]
			varOld = var.__class__(f"{var.name}_old")
			term.args[idx] = Multiplication( Subtraction(var, varOld), overTimeDifferential )

		elif hasSubclass(arg, Operator):
			integrateInTime(arg)

def applyDistributiveProperty(term, flag=1):
	for idx, arg in enumerate(term.args):
		if hasSubclass(arg, Multiplication):
			# mult_func( a, b, ..., sum(c,d,e), sum(f,g,h) )
			# sArg = [ mArg for mArg in arg.args if mArg is Sum ][0]
			# X = [mArg for mArg in arg.args if mArg is not sArg]
			# sum(mult_func(*X, f), mult_func(*X, g), mult_func(*X, h))
			for sArg in arg.args:
				if hasSubclass(sArg, Sum) or hasSubclass(sArg, Subtraction):
					plus_func = sArg.__class__
					args = [ mArg for mArg in arg.args if not mArg == sArg ]

					term.args[idx] = plus_func(*[Multiplication(*args, ssArg) for ssArg in sArg.args])

					# only one sArg
					break




		# if hasSubclass(arg, Multiplication) or hasSubclass(arg, Division):
		# 	if hasSubclass(arg, Multiplication) and ( hasSubclass(arg.args[1], Sum) or hasSubclass(arg.args[1], Subtraction) ):
		# 		arg.args = arg.args[::-1]
		# 	if hasSubclass(arg.args[0], Sum) or hasSubclass(arg.args[0], Subtraction):
		# 		mult_func = arg.__class__
		# 		plus_func = arg.args[0].__class__
		# 		a,b = arg.args[0].args
		# 		c = arg.args[1]

		# 		term.args[idx] = plus_func( mult_func(a,c), mult_func(b,c) )

		if hasSubclass(arg, Operator):
			applyDistributiveProperty(arg, flag+1)

	if flag > 0:
		applyDistributiveProperty(term, flag-1)

def clusterOperations(term, flag=4):
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

		# binaryOperator(A) -> A
		if hasSubclass(arg, BinaryOperator):
			if len(arg.args) == 1:
				term.args[idx] = arg.args[0]

		if hasSubclass(arg, Operator):
			clusterOperations(arg, flag+1)

	if flag > 0:
		clusterOperations(term, flag-1)

def containsVariables(term, variables):
	if hasSubclass(term, Variable):
		return term.name in variables
	
	elif hasSubclass(term, Operator):
		for arg in term.args:
			if containsVariables(arg, variables):
				return True
		return False
