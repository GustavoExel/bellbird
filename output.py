import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, "PyEFVLib"))
import numpy as np
import PyEFVLib

def heatTransfer(problemData):
	propertyData 	 = problemData.propertyData
	timeStep 		 = problemData.timeStep
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension="xdmf")

	oldTemperatureField = problemData.initialValues["temperature"].copy()
	temperatureField    = np.repeat(0.0, numberOfVertices)

	matrix 		= np.zeros((numberOfVertices, numberOfVertices))
	independent = np.zeros(numberOfVertices)

	k 	= propertyData.get(0, "k")
	q 	= propertyData.get(0, "q")
	rho = propertyData.get(0, "rho")
	cp 	= propertyData.get(0, "cp")

	def assembleMatrix():
		for region in grid.regions:
			k = propertyData.get(region.handle, "k")
			for element in region.elements:
				for innerFace in element.innerFaces:
					area = innerFace.area.getCoordinates()[:dimension]
					coefficients = -1.0 * k * np.matmul( area.T, innerFace.globalDerivatives )
					backwardVertexHandle, forwardVertexHandle = innerFace.getNeighborVerticesHandles()

					for local, vertex in enumerate(element.vertices):
						matrix[backwardVertexHandle][vertex.handle] += coefficients[local]
						matrix[forwardVertexHandle][vertex.handle]  -= coefficients[local]

		for region in grid.regions:
			rho = propertyData.get(region.handle, "rho")
			cp = propertyData.get(region.handle, "cp")
			coefficients = rho * cp / timeStep

			for element in region.elements:
				for local, vertex in enumerate(element.vertices):
					matrix[vertex.handle][vertex.handle] += element.subelementVolumes[local] * coefficients

		for bCondition in problemData.dirichletBoundaries["temperature"]:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle] = np.zeros(numberOfVertices)
				matrix[vertex.handle][vertex.handle] = 1.0

		# Invert Matrix
		inverseMatrix = np.linalg.inv(matrix)

		return inverseMatrix


	def assembleVector():
		independent = np.zeros(grid.numberOfVertices)

		# Generation Term
		for region in grid.regions:
			q = propertyData.get(region.handle, "q")
			for element in region.elements:
				for local, vertex in enumerate(element.vertices):
					independent[vertex.handle] += element.subelementVolumes[local] * q

		# Transient Term
		for region in grid.regions:
			rho = propertyData.get(region.handle, "rho")
			cp = propertyData.get(region.handle, "cp")
			coefficients = rho * cp / timeStep

			for element in region.elements:
				for local, vertex in enumerate(element.vertices):
					independent[vertex.handle] += element.subelementVolumes[local] * coefficients * oldTemperatureField[vertex.handle]

		# Neumann Boundary Condition
		for bCondition in problemData.neumannBoundaries["temperature"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in problemData.dirichletBoundaries["temperature"]:
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
		independent = assembleVector()
		temperatureField = np.matmul(inverseMatrix, independent)
		
		difference = max( abs(temperatureField - oldTemperatureField) )
		temperatureField = oldTemperatureField.copy()
		currentTime += timeStep
		iteration += 1

		saver.save("temperature", temperatureField, currentTime)

		print("{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}".format(iteration, currentTime, timeStep, difference))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )

		if iteration >= problemData.maxNumberOfIterations:
			break

	saver.finalize()

if __name__ == '__main__':
	problemData = PyEFVLib.ProblemData(
		meshFilePath = "{MESHES}/msh/2D/Square.msh",
		outputFilePath = "{RESULTS}/../../BellBird2/results",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 0.1, tolerance = 1e-4, maxNumberOfIterations = 300 ),
		propertyData = PyEFVLib.PropertyData({
		    "Body":
		    {
		    	"k": 1,
		    	"rho": 1,
		    	"cp": 1,
		    	"q": 0,
		    }
		}),
		boundaryConditions = PyEFVLib.BoundaryConditions({
			"temperature": {
				"InitialValue": 0.0,
				"West":	 { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant,"value" : 20.0 },
				"East":	 { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant,"value" : 50.0 },
				"South": { "condition" : PyEFVLib.Neumann,   "type" : PyEFVLib.Constant,"value" : 0.0 },
				"North": { "condition" : PyEFVLib.Neumann,   "type" : PyEFVLib.Constant,"value" : 0.0 },
			}
		}),
	)


	heatTransfer( problemData )