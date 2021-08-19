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
    scripts=['aln_filter','aln_hapcounter','aln_pdistancer','aln_renamer','aln_splitspecies','aln_summary','basefunctions'],
    packages=setuptools.find_packages(include=['pycoistats', 'pycoistats.*']),
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: GPL-3.0 License",
    ),
)