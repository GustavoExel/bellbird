import bellbird

# Stress Equilibrium
model = bellbird.Model(
	name = "Stress Equilibrium",
	equationsStr = ["div(Ce * grad_s(vec(u))) + rho * g = vec(0)"],
	variables   = ["u_x", "u_y"],
	definitions = [
		"g = [0.0, -9.81]",
		"lame = 2*G*nu/(1-2*nu)",
		"Ce = [[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]",
	],
	properties = {
    	"rho": 1800.0,
        "nu": 0.4,
        "G": 6.0e+06,
	},
	boundaryConditions = [
		bellbird.InitialCondition("u_x", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Dirichlet, "West", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Dirichlet, "East", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Neumann, "North", 0.0),
		bellbird.InitialCondition("u_y", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "West", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "East", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Dirichlet, "South", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "North", -1e4),
	],
	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",
	meshDimension = 2,
)

# print(model.discretizedEquations[0])
model.compile()
model.run()