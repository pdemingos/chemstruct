import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(    
    name="chemstruct",
    version="2019.06.2",
    author="Pedro Guerra Demingos",
    author_email="pdemingos@gmail.com",
    description="Package for chemical structure analysis of atomic files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pdemingos/chemstruct",
    packages=["chemstruct"],
    package_dir={"chemstruct": "chemstruct"},
    package_data={"chemstruct": ["*"]},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)

