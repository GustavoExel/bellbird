import bellbird

# Geomechanics
model = bellbird.Model(
	name = "Geomechanics",
	equationsStr = [
		"div( Ce*grad_s(vec(u)) ) - grad(alpha*p) + rho * g = vec(0)",
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
	sparse = True,
	meshPath = "../PyEFVLib/meshes/msh/2D/10x10.msh",
	maxNumberOfIterations = 70,
	timeStep = 10,
	tolerance = 1e-4,
)

model.compile()
model.run()

"""
--> PENDÊNCIAS
	--> Tá demorando muuito, precisamos revisar os esquemas com as flags

	--> Avisar se o programa não der conta de implementar algum termo

	--> Adicionar um help
		- Avisar para usar o np.array por causa do (-1) * 

	--> Precisamos melhorar o __str__ por causa que é o que vai para o script
		- a * (b + c) pode virar a * b + c
		- a * (b / c) * d pode virar a * b / c * d

	--> Fica mais organizado se fizermos uma classe separada do tipo writer
	--> Espelhar tudo o que foi implementado na matrix pro independente
	--> Revisar a rigorosidade dos filtros no writer, ELES TÃO UM LIXO
	--> Padronizar a indexação
		--> coeff antes ou depois do number of vertices?
		--> Padronizar quando esconde ou não o zero
		--> Quem sabe colocar comentários para explicar a indexação
	--> Urgentíssimo verificar as condições de contorno e ver o esquema da outer face

	--> Revisar o writeIndependentVolumeIntegral

	--> No stress equilibrium a gente não arrumou o u_xField, e talvez n precisaria, mas é bom

	--> Ver dimension==3 no saver sem a pressão

	--> Ver direitinho o que ainda não tá certo

	--> Ver a geomecânica

	--> Ver CERTINHO a boundary condition (tipo se tem outerFace ou não)

	--> Se não der certo o esquema da Boundary Condition alertar o usuário que ele pode ter que revisar as boundary conditions
"""