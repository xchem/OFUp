import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ofup-tools",
    version="0.0.1",
    author="Rachael Skyner",
    author_email="rachael.skyner@diamond.ac.uk",
    description="Open Follow-Up tools used by XChem",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xchem/OFUp/tree/master/tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE",
        "Operating System :: OS Independent",
    ],
)