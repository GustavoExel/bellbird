import bellbird

# Heat Transfer
model = bellbird.Model(
	name = "Heat Transfer",
	equationsStr = ["rho * cp * d/dt(temperature) = const(k) * div( grad(temperature) ) + q"],
	# equationStr = "const(k) * div( grad(temperature) ) + q = rho * cp * d/dt(temperature)",
	variables   = ["temperature"],
	properties = {
		"k" : 1.0,
		"q" : 0.0,
		"rho" : 1.0,
		"cp" : 1.0,
	},
	boundaryConditions = [
		bellbird.InitialCondition("temperature", 0.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "West", 100.0),
		bellbird.BoundaryCondition("temperature", bellbird.Dirichlet, "East", 100.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "North", 0.0),
	],
	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",
	meshDimension = 2,
)

# print(model.discretizedEquations[0])
model.compile()
# model.run()