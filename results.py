import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, 'PyEFVLib'))
import PyEFVLib
import numpy as np
from scipy import sparse
import scipy.sparse.linalg

def geomechanics(problemData):
	propertyData 	 = problemData.propertyData
	timeStep 		 = problemData.timeStep
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension='xdmf')

	uField    = np.repeat(0.0, dimension*numberOfVertices)
	oldUField = np.concatenate((problemData.initialValues['u_x'], problemData.initialValues['u_y'], problemData.initialValues['u_z'])) if dimension==3 else np.concatenate((problemData.initialValues['u_x'], problemData.initialValues['u_y']))

	pField    = np.repeat(0.0, numberOfVertices)
	oldPField = problemData.initialValues['p'].copy()

	matrixCoords = []
	matrixVals = []

	nu         = propertyData.get(0, 'nu')
	G          = propertyData.get(0, 'G')
	cs         = propertyData.get(0, 'cs')
	phi        = propertyData.get(0, 'phi')
	k          = propertyData.get(0, 'k')
	cf         = propertyData.get(0, 'cf')
	mu         = propertyData.get(0, 'mu')
	rhos       = propertyData.get(0, 'rhos')
	rhof       = propertyData.get(0, 'rhof')

	g = np.array([0.0, 0.0, 0.0])[:dimension]
	lame = 2*G*nu/(1-2*nu)
	Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
	rho = phi * rhof + (1-phi) * rhos
	K = 2*G*(1 + nu) / 3*(1-2*nu)
	cb = 1 / K
	alpha = 1 - cs / cb
	S = (phi * cf + (alpha-phi) * cs)

	def add(i, j, val):
		matrixCoords.append((i, j))
		matrixVals.append(val)

	def set(i):
		_matrixCoords = [coord for coord in matrixCoords if coord[0]!=i] + [(i,i)]
		_matrixVals   = [val for coord,val in zip(matrixCoords,matrixVals) if coord[0]!=i] + [1.0]
		matrixCoords.clear(); matrixCoords.extend(_matrixCoords)
		matrixVals.clear(); matrixVals.extend(_matrixVals)

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
		# S * (1/timeStep) * p
		for vertex in grid.vertices:
			add(vertex.handle+(dimension)*numberOfVertices, vertex.handle+(dimension)*numberOfVertices, vertex.volume * S * (1/timeStep) )

		# (-1) * alpha * grad(p)
		for element in grid.elements:
			for innerFace in element.innerFaces:
				m = element.vertices.size
				transposedVoigtArea = getTransposedVoigtArea(innerFace)
				shapeFunctions = innerFace.getShapeFunctions()
				identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)]) if dimension == 3 else np.array([shapeFunctions, shapeFunctions, np.zeros(m)])

				matrixCoefficients = (-1) * alpha * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)

				backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()
				for coord in range(dimension):
					for local, vertex in enumerate(element.vertices):
						add(backwardsHandle+numberOfVertices*(coord), vertex.handle+(dimension)*numberOfVertices, matrixCoefficients[coord][local])
						add(forwardHandle+numberOfVertices*(coord), vertex.handle+(dimension)*numberOfVertices, -matrixCoefficients[coord][local])

		# (-1) * k * mu * grad(p)
		for element in grid.elements:
			for innerFace in element.innerFaces:
				area = innerFace.area.getCoordinates()[:dimension]
				matrixCoefficients = (-1) * k * mu * np.matmul( area.T, innerFace.globalDerivatives )
				backwardVertexHandle, forwardVertexHandle = innerFace.getNeighborVerticesHandles()
				for local, vertex in enumerate(element.vertices):
					add(backwardVertexHandle+(dimension)*numberOfVertices, vertex.handle+(dimension)*numberOfVertices, matrixCoefficients[local])
					add(forwardVertexHandle+(dimension)*numberOfVertices, vertex.handle+(dimension)*numberOfVertices, -matrixCoefficients[local])

		# Ce * grad_s(u)
		for element in grid.elements:
			for innerFace in element.innerFaces:
				transposedVoigtArea = getTransposedVoigtArea(innerFace)
				voigtGradientOperator = getVoigtGradientOperator(innerFace.globalDerivatives)
				coeff = Ce
				matrixCoefficients = np.einsum('ij,jk,kmn->imn', transposedVoigtArea, coeff, voigtGradientOperator)

				backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()
				for local, vertex in enumerate(element.vertices):
					for i in range(dimension):
						for j in range(dimension):
							add(backwardsHandle + numberOfVertices*(i), vertex.handle + numberOfVertices*(j),  matrixCoefficients[i][j][local])
							add(forwardHandle   + numberOfVertices*(i), vertex.handle + numberOfVertices*(j), -matrixCoefficients[i][j][local])

		# alpha * (1/timeStep) * u
		for element in grid.elements:
			for face in element.faces:
				area = face.area.getCoordinates()[:dimension]
				shapeFunctions = face.getShapeFunctions()

				for coord in range(dimension):
					for local, vertex in enumerate(element.vertices):
						if type(face) == PyEFVLib.InnerFace:
							backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()

							add(backwardsHandle+numberOfVertices*(dimension), vertex.handle + numberOfVertices*(coord+0), alpha * (1/timeStep) * shapeFunctions[local] * area[coord])
							add(forwardHandle+numberOfVertices*(dimension), vertex.handle + numberOfVertices*(coord+0), -alpha * (1/timeStep) * shapeFunctions[local] * area[coord])

						elif type(face) == PyEFVLib.OuterFace:
							add(face.vertex.handle+numberOfVertices*(dimension), vertex.handle + numberOfVertices*(coord+0), alpha * (1/timeStep) * shapeFunctions[local] * area[coord])

		# Dirichlet Boundary Conditions
		for bCondition in problemData.dirichletBoundaries['u_x']:
			for vertex in bCondition.boundary.vertices:
				set(vertex.handle+(0)*numberOfVertices)

		for bCondition in problemData.dirichletBoundaries['u_y']:
			for vertex in bCondition.boundary.vertices:
				set(vertex.handle+(1)*numberOfVertices)

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries['u_z']:
				for vertex in bCondition.boundary.vertices:
					set(vertex.handle+(2)*numberOfVertices)

		for bCondition in problemData.dirichletBoundaries['p']:
			for vertex in bCondition.boundary.vertices:
				set(vertex.handle+(dimension)*numberOfVertices)

		# Inverse Matrix
		matrix = sparse.csc_matrix( (matrixVals, zip(*matrixCoords)), shape=((1+dimension)*numberOfVertices, (1+dimension)*numberOfVertices) )
		inverseMatrix = sparse.linalg.inv( matrix )

		return inverseMatrix

	def assembleIndependent():
		independent = np.zeros((1+dimension)*numberOfVertices)

		# rho * g
		for vertex in grid.vertices:
			for coord in range(dimension):
				independent[vertex.handle+coord*numberOfVertices] += vertex.volume * rho * g[coord]

		# S * (1/timeStep) * p_old
		for vertex in grid.vertices:
			independent[vertex.handle+(dimension)*numberOfVertices] += vertex.volume * S * (1/timeStep) * oldPField[vertex.handle]

		# (-1) * k * mu * rho * g
		for element in grid.elements:
			for innerFace in element.innerFaces:
				area = innerFace.area.getCoordinates()[:dimension]
				coefficient = (-1) * k * mu * rho * g
				coefficient = np.matmul(area.T, coefficient)

				backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()
				independent[backwardsHandle+numberOfVertices*(dimension)] += coefficient
				independent[forwardHandle+numberOfVertices*(dimension)]   -= coefficient

		# alpha * (1/timeStep) * u_old
		for element in grid.elements:
			for face in element.faces:
				area = face.area.getCoordinates()[:dimension]
				shapeFunctions = face.getShapeFunctions()

				for coord in range(dimension):
					for local, vertex in enumerate(element.vertices):
						if type(face) == PyEFVLib.InnerFace:
							backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
							independent[backwardsHandle+numberOfVertices*(dimension)] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * oldUField[vertex.handle + numberOfVertices*(coord)]
							independent[forwardHandle+numberOfVertices*(dimension)] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * oldUField[vertex.handle + numberOfVertices*(coord)]

						elif type(face) == PyEFVLib.OuterFace:
							independent[face.vertex.handle+numberOfVertices*(dimension)] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * oldUField[vertex.handle + numberOfVertices*(coord)]
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

		for bCondition in problemData.neumannBoundaries['p']:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle+(dimension)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

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

		for bCondition in problemData.dirichletBoundaries['p']:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle+(dimension)*numberOfVertices] = bCondition.getValue(vertex.handle)

		return independent

	tolerance = problemData.tolerance
	difference = 2*tolerance
	iteration = 0
	currentTime = 0.0
	converged = False

	inverseMatrix = assembleMatrix()

	while not converged:
		independent = assembleIndependent()

		results = inverseMatrix.dot(independent)
		uField = results[(0)*numberOfVertices:(0+dimension)*numberOfVertices]
		pField = results[(dimension)*numberOfVertices:(dimension+1)*numberOfVertices]

		difference = max( max(abs(pField - oldPField)), max(abs(uField - oldUField)) )

		oldPField = pField.copy()
		oldUField = uField.copy()

		currentTime += timeStep
		iteration += 1

		saver.save('u_x', uField[0*numberOfVertices:1*numberOfVertices], currentTime)
		saver.save('u_y', uField[1*numberOfVertices:2*numberOfVertices], currentTime)
		if dimension == 3:
			saver.save('u_z', uField[2*numberOfVertices:3*numberOfVertices], currentTime)
		saver.save('p', pField, currentTime)

		print('{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}'.format(iteration, currentTime, timeStep, difference))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )

		if iteration >= problemData.maxNumberOfIterations:
			break

	saver.finalize()

def main():
	problemData = PyEFVLib.ProblemData(
		meshFilePath = '../PyEFVLib/meshes/msh/2D/10x10.msh',
		outputFilePath = 'results',
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 10, tolerance = 0.001, maxNumberOfIterations = 70 ),
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
			'p': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
		}),
	)

	geomechanics( problemData )

if __name__ == '__main__':
	main()
