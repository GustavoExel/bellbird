import bellbird

name = "Heat Transfer"
# equationStr = "S * d/dt(p) + a * div(d/dt(u)) + div(v) = 0"
equationStr = "rho * cp * d/dt(temperature) = const(k) * div( grad(temperature) ) + q"

variables = ["temperature"]

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

# print(model.discretizedEquation)
model.compile()


"""
Plano:
	- Implementar o HeatTransfer com bastantes assumptions e sem muita preocupação com a
		idependência do programa, como calcular as condições de contorno, ...
		+ Não se preocupar com o arranjo da matriz, ou com a dimensão dos termos (como grad, div, matrizes)
		+ Não se preocupar com termos diferentes (como função de forma..., campos na integral de superfície)
		+ Não se preocupar com termos especiais como matriz constitutiva

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