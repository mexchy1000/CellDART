from setuptools import setup

setup(
    name = "CellDART",
    version = "0.0.1",
    description = "Cell type inference by domain adaptation of single-cell and spatial transcriptomic data",
    url = "https://github.com/mexchy1000/CellDART.git",
    author = "Hongyoon Choi, Sungwoo Bae",
    install_requires = ["scanpy==1.5.1","pandas==1.3.5","numpy==1.21.6",
                        "h5py==2.10.0", "jupyter",
                        "keras==2.3.1", "tensorflow==1.14.0", "tensorflow-gpu==1.14.0"]
)