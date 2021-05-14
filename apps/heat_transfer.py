import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import bellbird

# Heat Transfer
model = bellbird.Model(
	name = "Heat Transfer",
	equationsStr = ["k * div( grad(T) ) + q = rho * cp * d/dt(T)"],
	variables   = ["T"],
	properties = {
		"Outer":{
			"k" : 1.0,
			"q" : 0.0,
			"rho" : 1.0,
			"cp" : 1.0,
		},
		"Inner":{
			"k" : 10000.0,
			"q" : 0.0,
			"rho" : 1.0,
			"cp" : 1.0,
		},
	},
	boundaryConditions = [
		bellbird.InitialCondition("T", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Dirichlet, "West", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Dirichlet, "East", 100.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "North", 0.0),
	],
	meshPath = "../../PyEFVLib/meshes/msh/2D/vug 15x15.msh",
	timeStep = 0.05,
	tolerance = 1e-4,
	maxNumberOfIterations = 300,
	sparse = False,
)

model.run()