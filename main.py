import bellbird

# Stress Equilibrium
model = bellbird.Model(
	name = "Stress Equilibrium and Heat Transfer",
	equationsStr = [
		"div(Ce * grad_s(vec(u))) + rho * g = vec(0)",
		"rho * cp * d/dt(temperature) = const(k) * div( grad(temperature) ) + q",
	],
	variables   = ["u_x", "u_y", "u_z", "temperature"],
	definitions = [
		"g = np.array([0.0, -9.81, 0.0])[:dimension]",
		"lame = 2*G*nu/(1-2*nu)",
		"Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])",
	],
	properties = {
    	"rho": 1800.0,
        "nu": 0.4,
        "G": 6.0e+06,
		"k" : 1.0,
		"q" : 0.0,
		"rho" : 1.0,
		"cp" : 1.0,
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
)

model.compile()
model.run()

"""
Plano:
	- Implementar o HeatTransfer com bastantes assumptions e sem muita preocupação com a
		idependência do programa, como calcular as condições de contorno, ...
		+ Não se preocupar com termos diferentes (como função de forma..., campos na integral de superfície)
		+ Não se preocupar com termos especiais como matriz constitutiva
		+ Não se preocupar com o arranjo da matriz, ou com a dimensão dos termos (como grad, div, matrizes)

	- Implementar o Stress Equilibrium
		+ Começar a se preocupar com calcular as condições de contorno
		+ Se preocupar um pouco mais com o arranjo da matriz (mas não taanto)
		+ Se preocupar com as dimensões dos termos (vetores, grads, divs, e matrizes)
		+ Se preocupar com campos na iint (que precisa de funções de forma)
		+ Não se preocupar com várias equações
		+ Não se preocupar com arranjo na matriz

	- Implementar a geomecânica
		+ Se preocupar com mais de uma equação

"""

"""
--> PENDÊNCIAS
	--> Sparse
	--> Implementar funções do clusterTerms como:
		--> (1+1) = 2
		--> A+A+B = 2*A+B
		--> n*A+k*A+B = (n+k)*A + B
"""