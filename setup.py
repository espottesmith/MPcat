from setuptools import setup, find_packages

"""
Setup script for mpcat.
"""

setup(
    name='mpcat',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/espottesmith/MPcat',
    author='Evan Walter Clark Spotte-Smith',
    author_email='espottesmith@gmail.com',
    python_requires=">=3.8",
    description="""An interface between the Materials Project software suite and the Schrodinger
Python API, designed to allow for high-throughput execution of Jaguar 
and AutoTS calculations for molecular thermodynamics and kinetics.

Let the cat out of the box!""",
    keywords=["schrodinger", "jaguar", "autots", "chemistry", "kinetics", "materials",
              "project", "electronic", "structure", "analysis", "automation",
              "transition", "state"]
)
