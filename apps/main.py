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
	--> Ver se o self.variables ?? compat??vel com as condi????es de contorno
	--> Quem sabe eliminar a necessidade de escrever vec(0)
	--> Precisamos melhorar o __str__ por causa que ?? o que vai para o script
		- a * (b + c) pode virar a * b + c
		- a * (b / c) * d pode virar a * b / c * d

	--> Revisar o writeIndependentVolumeIntegral

	--> Adicionar um help
		- Avisar para usar o np.array por causa do (-1) * 
	se n colocar nada tolerancia = 0
"""

# N??o Urgente
"""
	--> Tinhamos pensado em espelhar todas as integrais da matriz pro independente, mas n??o ?? uma boa fazer isso sem olhar com cuidado para cada uma
		> Mesmo assim ?? bom pensar que se quisermos aplicar uma derivada temporal em tudo teremos que passar tudo da matrix pro independente
	--> Quem sabe colocar coment??rios para explicar a indexa????o
	--> No stress equilibrium a gente n??o arrumou o u_xField, e talvez n precisaria, mas ?? bom

"""


# D??vida
"""
	--> Verificar as condi????es de contorno (se vai outerFace ou n??o, cond de Neumann, fluxo...)
	--> Se n??o der certo o esquema da Boundary Condition alertar o usu??rio que ele pode ter que revisar as boundary conditions


# Novo
	--> Implementar timeStep vari??vel
	--> Declarar matrix dentro da fun????o de matrix
	--> Usar gravidade como propriedade (isso ?? no uso)
	--> Condi????o de Neumann errada (independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) ao inv??s de outerFace.vertex.handle)

"""