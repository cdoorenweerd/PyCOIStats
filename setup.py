import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pycoistats",
    version="1.4.2",
    author="Camiel Doorenweerd",
    author_email="c.doorenweerd@gmail.com",
    description="A python package to calculate statistics on haploid DNA sequence data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cdoorenweerd/PyCOIStats",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: GPL-3.0 License",
    ),
)