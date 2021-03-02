import bellbird

name = "Heat Transfer"
# equationStr = "S * d/dt(p) + a * div(d/dt(u)) + div(v) = 0"
equationStr = "rho * cp * d/dt(T) = const(k) * div( grad(T) ) + q"

variables = ["T", "p", "u", "v"]

properties = {
	"k" : 1.0,
	"q" : 0.0,
	"rho" : 1000.0,
	"cp" : 120.0,
}

# boundaryConditions = [
# 	bellbird.BoundaryCondition("T", bellbird.Dirichlet, "Body", 0.0),
# 	bellbird.BoundaryCondition("T", bellbird.Dirichlet, "West", 0.0),
# 	bellbird.BoundaryCondition("T", bellbird.Dirichlet, "East", 100.0),
# 	bellbird.BoundaryCondition("T", bellbird.Neumann, "South", 0.0),
# 	bellbird.BoundaryCondition("T", bellbird.Neumann, "North", 0.0),
# ]

meshPath = "./mesh.msh"

model = bellbird.Model(
	name = name,
	equationStr = equationStr,
	variables   = variables,
	properties = properties,
	# boundaryConditions = boundaryConditions,
	meshPath = meshPath,
)

print(model.discretizedEquation)