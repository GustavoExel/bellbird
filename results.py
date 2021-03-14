import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'PyEFVLib'))
import PyEFVLib
import numpy as np

def geomechanics(problemData):
	propertyData 	 = problemData.propertyData
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension='xdmf')

	uField    = np.repeat(0.0, dimension*numberOfVertices)
	matrix 		= np.zeros(((dimension)*numberOfVertices, (dimension)*numberOfVertices))

	nu         = propertyData.get(0, 'nu')
	G          = propertyData.get(0, 'G')
	cs         = propertyData.get(0, 'cs')
	phi        = propertyData.get(0, 'phi')
	k          = propertyData.get(0, 'k')
	cf         = propertyData.get(0, 'cf')
	mu         = propertyData.get(0, 'mu')
	rhos       = propertyData.get(0, 'rhos')
	rhof       = propertyData.get(0, 'rhof')
	p          = propertyData.get(0, 'p')

	g = np.array([0.0, 0.0, 0.0])[:dimension]
	lame = 2*G*nu/(1-2*nu)
	Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
	rho = phi * rhof + (1-phi) * rhos
	K = 2*G*(1 + nu) / 3*(1-2*nu)
	cb = 1 / K
	alpha = 1 - cs / cb
	S = (phi * cf + (alpha-phi) * cs)

	def getTransposedVoigtArea(innerFace):
		Sx, Sy, Sz = innerFace.area.getCoordinates()
		return np.array([[Sx,0,Sy],[0,Sy,Sx]]) if dimension==2 else np.array([[Sx,0,0,Sy,0,Sz],[0,Sy,0,Sx,Sz,0],[0,0,Sz,0,Sy,Sx]])

	def getVoigtGradientOperator(globalDerivatives):
		if len(globalDerivatives) == 2:
			Nx,Ny = globalDerivatives
			zero=np.zeros(Nx.size)
			return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])

		if len(globalDerivatives) == 3:
			Nx,Ny,Nz = globalDerivatives
			zero=np.zeros(Nx.size)
			return np.array([[Nx,zero,zero],[zero,Ny,zero],[zero,zero,Nz],[Ny,Nx,zero],[zero,Nz,Ny],[Nz,zero,Nx]])

	def assembleMatrix():
		# Ce * grad_s(u)
		for element in grid.elements:
			for face in element.innerFaces:
				area = face.area.getCoordinates()[:dimension]
				transposedVoigtArea = getTransposedVoigtArea(face)
				voigtGradientOperator = getVoigtGradientOperator(face.globalDerivatives)
				coeff = Ce
				matrixCoefficients = np.einsum('ij,jk,kmn->imn', transposedVoigtArea, coeff, voigtGradientOperator)
				backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
				for i in range(dimension):
					for j in range(dimension):
						for local, vertex in enumerate(element.vertices):
							matrix[backwardsHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += matrixCoefficients[i][j][local]
							matrix[forwardHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += -matrixCoefficients[i][j][local]

		# Dirichlet Boundary Conditions
		for bCondition in problemData.dirichletBoundaries['u_x']:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle+(0)*numberOfVertices] = np.zeros((dimension)*numberOfVertices)
				matrix[vertex.handle+(0)*numberOfVertices][vertex.handle+(0)*numberOfVertices] = 1.0

		for bCondition in problemData.dirichletBoundaries['u_y']:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle+(1)*numberOfVertices] = np.zeros((dimension)*numberOfVertices)
				matrix[vertex.handle+(1)*numberOfVertices][vertex.handle+(1)*numberOfVertices] = 1.0

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries['u_z']:
				for vertex in bCondition.boundary.vertices:
					matrix[vertex.handle+(2)*numberOfVertices] = np.zeros((dimension)*numberOfVertices)
					matrix[vertex.handle+(2)*numberOfVertices][vertex.handle+(2)*numberOfVertices] = 1.0

		# Inverse Matrix
		inverseMatrix = np.linalg.inv(matrix)
		return inverseMatrix

	def assembleIndependent():
		independent = np.zeros((dimension)*numberOfVertices)

		# rho * g
		for vertex in grid.vertices:
			for coord in range(dimension):
				independent[vertex.handle+coord*numberOfVertices] += vertex.volume * rho * g[coord]

		# Neumann Boundary Condition
		for bCondition in problemData.neumannBoundaries['u_x']:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle+(0)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		for bCondition in problemData.neumannBoundaries['u_y']:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle+(1)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		if dimension == 3:
			for bCondition in problemData.neumannBoundaries['u_z']:
				for facet in bCondition.boundary.facets:
					for outerFace in facet.outerFaces:
						independent[outerFace.vertex.handle+(2)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in problemData.dirichletBoundaries['u_x']:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle+(0)*numberOfVertices] = bCondition.getValue(vertex.handle)

		for bCondition in problemData.dirichletBoundaries['u_y']:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle+(1)*numberOfVertices] = bCondition.getValue(vertex.handle)

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries['u_z']:
				for vertex in bCondition.boundary.vertices:
					independent[vertex.handle+(2)*numberOfVertices] = bCondition.getValue(vertex.handle)

		return independent

	inverseMatrix = assembleMatrix()
	independent = assembleIndependent()

	results = np.matmul(inverseMatrix, independent)
	u_xField = results[0*numberOfVertices:1*numberOfVertices]
	u_yField = results[1*numberOfVertices:2*numberOfVertices]
	if dimension == 3:
		u_zField = results[2*numberOfVertices:3*numberOfVertices]
	
	saver.save('u_x', u_xField, 0.0)
	saver.save('u_y', u_yField, 0.0)
	if dimension == 3:
		saver.save('u_z', u_zField, 0.0)
	saver.finalize()

def main():
	problemData = PyEFVLib.ProblemData(
		meshFilePath = '../PyEFVLib/meshes/msh/2D/10x10.msh',
		outputFilePath = 'results',
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 10, tolerance = 0.0001, maxNumberOfIterations = 70 ),
		propertyData = PyEFVLib.PropertyData({
			'Body':
			{
				'nu': 0.2,
				'G': 6000000000.0,
				'cs': 0.0,
				'phi': 0.19,
				'k': 1.9e-15,
				'cf': 3.0303e-10,
				'mu': 1000.0,
				'rhos': 2700.0,
				'rhof': 1000.0,
				'p': 52000.0,
			},
		}),
		boundaryConditions = PyEFVLib.BoundaryConditions({
			'u_x': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
			'u_y': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 100000.0 },
			},
		}),
	)

	geomechanics( problemData )

if __name__ == '__main__':
	main()
