import bellbird

# Stress Equilibrium + Heat Transfer
model = bellbird.Model(
	name = "Stress Equilibrium",
	equationsStr = [
		"div(Ce * grad_s(vec(u))) + rho * g = vec(0)",
	],
	variables   = ["u_x", "u_y", "u_z"],
	definitions = [
		"g = np.array([0.0, -9.81, 0.0])[:dimension]",
		"lame = 2*G*nu/(1-2*nu)",
		"Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])",
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
		# bellbird.BoundaryCondition("u_x", bellbird.Neumann, "Top", 0.0),				# 3D
		# bellbird.BoundaryCondition("u_x", bellbird.Neumann, "Bottom", 0.0),			# 3D
		bellbird.InitialCondition("u_y", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "West", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "East", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Dirichlet, "South", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "North", 1e4),
		# bellbird.BoundaryCondition("u_y", bellbird.Neumann, "Top", 0.0),				# 3D
		# bellbird.BoundaryCondition("u_y", bellbird.Neumann, "Bottom", 0.0),			# 3D
		# bellbird.InitialCondition("u_z", 0.0),										# 3D
		# bellbird.BoundaryCondition("u_z", bellbird.Neumann, "West", 0.0),				# 3D
		# bellbird.BoundaryCondition("u_z", bellbird.Neumann, "East", 0.0),				# 3D
		# bellbird.BoundaryCondition("u_z", bellbird.Neumann, "South", 0.0),			# 3D
		# bellbird.BoundaryCondition("u_z", bellbird.Neumann, "North", 0.0),			# 3D
		# bellbird.BoundaryCondition("u_z", bellbird.Dirichlet, "Top", 0.0),			# 3D
		# bellbird.BoundaryCondition("u_z", bellbird.Dirichlet, "Bottom", 0.0),			# 3D
		bellbird.InitialCondition("temperature", 0.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "West", 100.0),
		bellbird.BoundaryCondition("temperature", bellbird.Dirichlet, "East", 100.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "North", 0.0),
		# bellbird.BoundaryCondition("temperature", bellbird.Neumann, "Top", 0.0),		# 3D
		# bellbird.BoundaryCondition("temperature", bellbird.Neumann, "Bottom", 0.0),	# 3D
	],
	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",									# 2D
	# meshPath = "../PyEFVLib/meshes/msh/3D/Hexas.msh",									# 3D
	sparse = True,
)

model.compile()
model.run()