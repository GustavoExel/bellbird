import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import bellbird

# Geomechanics
model = bellbird.Model(
	name = "Geomechanics",
	equationsStr = [
		"div( Ce*grad_s(u) ) - grad(alpha*p) + rho * g = vec(0)",
		"S * d/dt(p) + alpha*d/dt(div(u)) + div( (k*mu) * ( rho*g - grad(p) ) ) = 0",
	],
	variables   = ["u_x", "u_y", "u_z", "p"],
	definitions = [
		"g = np.array([0.0, 0.0, 0.0])[:dimension]",
		"lame = 2*G*nu/(1-2*nu)",
		"Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])",

		"rho = phi * rhof + (1-phi) * rhos",
		"K = 2*G*(1 + nu) / 3*(1-2*nu)",
		"cb = 1 / K",
		"alpha = 1 - cs / cb",
		"S = (phi * cf + (alpha-phi) * cs)",
	],
	properties = {
		"Body":{
			"nu": 2.0e-1,
			"G": 6.0e+9,
			"cs": 0.0,
			"phi": 1.9e-1,
			"k": 1.9e-15,
			"cf": 3.0303e-10,
			"mu": 1/1.0e-03,
			"rhos": 2.7e+3,
			"rhof": 1.0e+3,
		},
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
		bellbird.BoundaryCondition("u_y", bellbird.Neumann, "North", 1e5),
		bellbird.InitialCondition("p", 0.0),
		bellbird.BoundaryCondition("p", bellbird.Neumann, "West", 0.0),
		bellbird.BoundaryCondition("p", bellbird.Neumann, "East", 0.0),
		bellbird.BoundaryCondition("p", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("p", bellbird.Neumann, "North", 0.0),
	],
	sparse = False,
	meshPath = "../../PyEFVLib/meshes/msh/2D/10x10.msh",
	maxNumberOfIterations = 200,
	timeStep = 10,
	tolerance = 1e-4,
)

model.compile()
model.run()


# Urgente
"""
	--> Ver se o self.variables é compatível com as condições de contorno
	--> Quem sabe eliminar a necessidade de escrever vec(0)
	--> Precisamos melhorar o __str__ por causa que é o que vai para o script
		- a * (b + c) pode virar a * b + c
		- a * (b / c) * d pode virar a * b / c * d

	--> Revisar o writeIndependentVolumeIntegral

	--> Adicionar um help
		- Avisar para usar o np.array por causa do (-1) * 
	se n colocar nada tolerancia = 0
"""

# Não Urgente
"""
	--> Tinhamos pensado em espelhar todas as integrais da matriz pro independente, mas não é uma boa fazer isso sem olhar com cuidado para cada uma
		> Mesmo assim é bom pensar que se quisermos aplicar uma derivada temporal em tudo teremos que passar tudo da matrix pro independente
	--> Quem sabe colocar comentários para explicar a indexação
	--> No stress equilibrium a gente não arrumou o u_xField, e talvez n precisaria, mas é bom

"""


# Dúvida
"""
	--> Verificar as condições de contorno (se vai outerFace ou não, cond de Neumann, fluxo...)
	--> Se não der certo o esquema da Boundary Condition alertar o usuário que ele pode ter que revisar as boundary conditions


# Novo
	--> Implementar timeStep variável
	--> Declarar matrix dentro da função de matrix
	--> Usar gravidade como propriedade (isso é no uso)
	--> Condição de Neumann errada (independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) ao invés de outerFace.vertex.handle)

"""