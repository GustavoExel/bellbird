import bellbird

# Stress Equilibrium
model = bellbird.Model(
	name = "Stress Equilibrium and Heat Transfer",
	equationsStr = [
		"div(Ce * grad_s(vec(u))) + rho * g = vec(0)",
		"rho * cp * d/dt(temperature) = const(k) * div( grad(temperature) ) + q",
	],
	variables   = ["u_x", "u_y", "temperature"],
	definitions = [
		"g = [0.0, -9.81]",
		"lame = 2*G*nu/(1-2*nu)",
		"Ce = [[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]",
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
		bellbird.InitialCondition("u_y", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "West", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "East", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Dirichlet, "South", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "North", -1e4),
		bellbird.InitialCondition("temperature", 0.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "West", 100.0),
		bellbird.BoundaryCondition("temperature", bellbird.Dirichlet, "East", 100.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "North", 0.0),
	],
	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",
	meshDimension = 2,
)

# model.compile()
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
	--> Resolver o problema do arranjo e do sinal do fluxo (Neumann)
	--> Declarar old só pras que precisam
	--> Largar de forecer a meshDimension e fazer tudo no programa, mais automático
"""