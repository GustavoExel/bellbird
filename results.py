import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, "PyEFVLib"))
import PyEFVLib
import numpy as np

def heatTransfer(problemData):
	propertyData 	 = problemData.propertyData
	timeStep 		 = problemData.timeStep
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension="xdmf")

	TField    = np.repeat(0.0, numberOfVertices)
	oldTField = problemData.initialValues["T"].copy()

	matrix 		= np.zeros(((1)*numberOfVertices, (1)*numberOfVertices))

	def assembleMatrix():
		for region in grid.regions:
			k          = propertyData.get(region.handle, "k")
			q          = propertyData.get(region.handle, "q")
			rho        = propertyData.get(region.handle, "rho")
			cp         = propertyData.get(region.handle, "cp")


			for element in region.elements:
				# (-1) * rho * cp * (1/timeStep) * T
				for vertex in element.vertices:
					matrix[vertex.handle][vertex.handle] += vertex.getSubElementVolume(element) * (-1) * rho * cp * (1/timeStep) 

				# k * grad(T)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					matrixCoefficients = k * np.matmul( area.T, face.globalDerivatives )
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for local, vertex in enumerate(element.vertices):
						matrix[backwardsHandle+(0)*numberOfVertices][vertex.handle+(0)*numberOfVertices] += matrixCoefficients[local]
						matrix[forwardHandle+(0)*numberOfVertices][vertex.handle+(0)*numberOfVertices] += -matrixCoefficients[local]

		# Dirichlet Boundary Conditions
		for bCondition in problemData.dirichletBoundaries["T"]:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle] = np.zeros((1)*numberOfVertices)
				matrix[vertex.handle][vertex.handle] = 1.0

		# Inverse Matrix
		inverseMatrix = np.linalg.inv(matrix)
		return inverseMatrix

	def assembleIndependent():
		independent = np.zeros((1)*numberOfVertices)

		for region in grid.regions:
			k          = propertyData.get(region.handle, "k")
			q          = propertyData.get(region.handle, "q")
			rho        = propertyData.get(region.handle, "rho")
			cp         = propertyData.get(region.handle, "cp")


			for element in region.elements:
				# rho * cp * (-1) * (1/timeStep) * T_old
				for vertex in element.vertices:
					independent[vertex.handle] += vertex.getSubElementVolume(element) * rho * cp * (-1) * (1/timeStep) * oldTField[vertex.handle]

				# (-1) * q
				for vertex in element.vertices:
					independent[vertex.handle] += vertex.getSubElementVolume(element) * (-1) * q

		# Neumann Boundary Condition
		for bCondition in problemData.neumannBoundaries["T"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in problemData.dirichletBoundaries["T"]:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle] = bCondition.getValue(vertex.handle)

		return independent

	tolerance = problemData.tolerance
	difference = 2*tolerance
	iteration = 0
	currentTime = 0.0
	converged = False

	inverseMatrix = assembleMatrix()

	while not converged:
		independent = assembleIndependent()
		TField = np.matmul(inverseMatrix, independent)

		difference = max(abs(TField - oldTField))

		oldTField = TField.copy()

		currentTime += timeStep
		iteration += 1

		saver.save("T", TField, currentTime)

		print("{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}".format(iteration, currentTime, timeStep, difference))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )

		if iteration >= problemData.maxNumberOfIterations:
			break

	saver.finalize()

def main():
	problemData = PyEFVLib.ProblemData(
		meshFilePath = "../PyEFVLib/meshes/msh/2D/vug 15x15.msh",
		outputFilePath = "results",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 0.05, tolerance = 0.0001, maxNumberOfIterations = 300 ),
		propertyData = PyEFVLib.PropertyData({
			"Outer":
			{
				"k": 1.0,
				"q": 0.0,
				"rho": 1.0,
				"cp": 1.0,
			},
			"Inner":
			{
				"k": 10000.0,
				"q": 0.0,
				"rho": 1.0,
				"cp": 1.0,
			},
		}),
		boundaryConditions = PyEFVLib.BoundaryConditions({
			"T": {
				"InitialValue": 0.0,
				"West": { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant, "value" : 0.0 },
				"East": { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant, "value" : 100.0 },
				"South": { "condition" : PyEFVLib.Neumann, "type" : PyEFVLib.Constant, "value" : 0.0 },
				"North": { "condition" : PyEFVLib.Neumann, "type" : PyEFVLib.Constant, "value" : 0.0 },
			},
		}),
	)

	heatTransfer( problemData )

if __name__ == "__main__":
	main()
