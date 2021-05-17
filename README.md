# bellbird

Library that parses differential equations (PDEs) strings and automates the creation of scripts for PyEFVLib using the element-based finite volume method. 

## Dependencies & Installation

- [Python 3](https://www.python.org/downloads/) (3.8.2);
- [meshio](https://pypi.org/project/meshio/) (4.0.15);
- [numpy](https://numpy.org/) (1.17.4);
- [pandas](https://pandas.pydata.org/)(1.1.3);
- [scipy](https://www.scipy.org/) (1.5.3);
- [PyEFVLib](https://pypi.org/project/PyEFVLib/) (1.0.5);


## Usage

```python
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import bellbird

# Heat Transfer
model = bellbird.Model(
	name = "Heat Transfer",
	equationsStr = ["k * div( grad(T) ) + q = rho * cp * d/dt(T)"],
	variables   = ["T"],
	properties = {
		"Body":{
			"k" : 10000.0,
			"q" : 0.0,
			"rho" : 1.0,
			"cp" : 1.0,
		},
	},
	boundaryConditions = [
		bellbird.InitialCondition("T", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Dirichlet, "West", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Dirichlet, "East", 100.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "South", 0.0),
		bellbird.BoundaryCondition("T", bellbird.Neumann, "North", 0.0),
	],
	meshPath = "mesh.msh",
	timeStep = 0.05,
	tolerance = 1e-4,
	maxNumberOfIterations = 300,
	sparse = False,
)

model.run()
```