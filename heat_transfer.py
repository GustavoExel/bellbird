import bellbird

# Heat Transfer
model = bellbird.Model(
	name = "Heat Transfer",
	equationsStr = ["k * div( grad(T) ) + q = rho * cp * d/dt(T)"],
	variables   = ["T"],
	properties = {
		"k" : 1.0,
		"q" : 0.0,
		"rho" : 1.0,
		"cp" : 1.0,
	},
	boundaryConditions = [
		bellbird.InitialCondition("T", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "West", 100.0),
		bellbird.BoundaryCondition("T", bellbird.Dirichlet, "East", 100.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "North", 0.0),
	],
	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",
	timeStep = 10,
	tolerance = 1e-4,
	maxNumberOfIterations = 300,
	sparse = True,
)

# print(model.discretizedEquations[0])
# model.compile()
model.run()