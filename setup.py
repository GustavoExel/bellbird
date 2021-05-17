from setuptools import setup

with open("README.md", "r") as f:
	long_description = f.read()

setup(
	name="bellbird",
	version="0.0.1",
	description="Library that parses differential equations (PDEs) strings and automates the creation of scripts for PyEFVLib using the element-based finite volume method. ",
	py_modules=["bellbird/__init__", "bellbird/variables", "bellbird/operators", "bellbird/equation", "bellbird/model", "bellbird/writer"],
	classifiers=[
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent"
	],
	long_description=long_description,
	long_description_content_type="text/markdown",
	install_requires=[
		"numpy ~= 1.19.2",
		"PyEFVLib",
	],
	url="https://github.com/GustavoExel/bellbird",
	author="Gustavo Exel",
	author_email="gustavoexelgpe@gmail.com",

)