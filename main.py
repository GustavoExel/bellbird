import bellbird

# Stress Equilibrium
model = bellbird.Model(
	name = "Stress Equilibrium",
	equationStr = "div(Ce * grad_s(vec(u))) + rho * g = vec(0)",
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
		bellbird.InitialCondition("u_x",  bellbird.Dirichlet, 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Dirichlet, "West", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Dirichlet, "East", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("u_x", bellbird.Neumann, "North", 0.0),
		bellbird.InitialCondition("u_y",  bellbird.Dirichlet, 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "West", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "East", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Dirichlet, "South", 0.0),
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "North", -1e6),
	],
	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",
)

print(model.discretizedEquation)
# print(model.equation.getDimension())
model.compile()
# model.run()

# # Heat Transfer
# model = bellbird.Model(
# 	name = "Heat Transfer",
# 	equationStr = "rho * cp * d/dt(temperature) = const(k) * div( grad(temperature) ) + q",
# 	variables   = ["temperature"],
# 	properties = {
# 		"k" : 1.0,
# 		"q" : 0.0,
# 		"rho" : 1.0,
# 		"cp" : 1.0,
# 	},
# 	boundaryConditions = [
# 		bellbird.InitialCondition("temperature",  bellbird.Dirichlet, 0.0),
# 		bellbird.BoundaryCondition("temperature", bellbird.Dirichlet, "West", 0.0),
# 		bellbird.BoundaryCondition("temperature", bellbird.Dirichlet, "East", 100.0),
# 		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "South", 0.0),
# 		bellbird.BoundaryCondition("temperature", bellbird.Neumann, "North", 0.0),
# 	],
# 	meshPath = "../PyEFVLib/meshes/msh/2D/Square.msh",
# )

# # print(model.discretizedEquation)
# # model.compile()
# # model.run()

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
!!!!!!!!!!!!!!!!!!!
--> PENDÊNCIAS
!!!!!!!!!!!!!!!!!!!
	--> Falta calcular automaticamente as condições de contorno
	--> falta implementar todos os esquemas pro stress equilibrium
	--> AINDA TEM MUUUITO PROBLEMA COM O ARRANJO DAS VARIÁVEIS NA MATRIZ [MUUUUUUUUUUUUUUUITO]
		--> Principalmente porque não podemos tratar as componentes dos vetores muito independentemente
		--> E também compoenentes de vetores como iiint(rho * g)

	--> Gostaria de fazer uma anotação pertinente sobre condições de contorno:
		Quando se discretiza uma equação teremos:
			div(k*X) -> iiint(div(k*X)) -> iint(k*X) --> summ2D(k*X * s)_i
			-> summ2D(k*X * s)_IN + summ2D(k*X * s)_Neumann + summ2D(k*X * s)_Dirichlet
			A(X) + div(k*X) = B
			iiint(A) + iint(k*X)_IN + iint(k*X)_N + iint(k*X)_D = iiint(B)
				iint_N passa para o rhs pois na cond de neumann o fluxo é dado e n depende de X
				iint_D desaparece pois se a condição for de dirichlet a eq toda muda
			iiint(A) + iint(k*X)_IN = iiint(B) - iint(k*X)_N
			iiint(A) + iint(k*X)_IN = iiint(B) - dot<flux, area>

			por isso no geral fazemos o somatório apenas com as innerFaces
			a de dirichlet é muito simples...



"""

"""
DISCRETIZAÇÃO DO STRESS EQUILIBRIUM

01) div(Ce * grad_s(u)) + rho * g = 0
02) iiint( div(Ce * grad_s(u)) ) + iiint( rho * g ) = 0
03) iint( Ce * grad_s(u) ) = - iiint( rho * g )


"""