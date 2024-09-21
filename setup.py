from setuptools import setup, find_packages

setup(
    name = "CellDART",
    version = "0.1.2",
    description = "Cell type inference by domain adaptation of single-cell and spatial transcriptomic data",
    url = "https://github.com/mexchy1000/CellDART.git",
    author = "Hongyoon Choi, Sungwoo Bae",
    packages=find_packages(include=['CellDART', 'CellDART.*']),
    install_requires = ["tensorflow~=2.9","tensorflow-gpu~=2.9", 
                        "pandas~=1.4","numpy~=1.20",
                        "scanpy","leidenalg","igraph",
                        "jupyter","ply","pytest"]
)
