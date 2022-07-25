from setuptools import setup, find_packages

setup(
    name = "CellDART",
    version = "0.1.0",
    description = "Cell type inference by domain adaptation of single-cell and spatial transcriptomic data",
    url = "https://github.com/mexchy1000/CellDART.git",
    author = "Hongyoon Choi, Sungwoo Bae",
    packages=find_packages(include=['CellDART', 'CellDART.*']),
    install_requires = ["tensorflow~=2.9.0","tensorflow-gpu~=2.9.0", 
                        "pandas~=1.4.0","numpy~=1.20.0",
                        "scanpy","leidenalg","python-igraph",
                        "jupyter","ply","pytest"]
)