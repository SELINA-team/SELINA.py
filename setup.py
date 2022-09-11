import setuptools

setuptools.setup(
    name="selina",
    version="0.1",
    author="Pengfei Ren",
    author_email="pfren@tongji.edu.cn",
    description="A annotation tool based on large-scale reference data",
    url="https://github.com/pfren1998/SELINA",
    packages=setuptools.find_packages(),
    package_dir={'selina': 'selina'},
    package_data={
        'selina':
        ['predict/downstream.R', 'preprocess/r/*', 'preprocess/data/*']
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'torch>=1.10.0', 'datatable>=0.11.1', 'pandas>=1.3.4', 'numpy>=1.21.2',
        'tqdm>=4.62.3', 'h5py>=3.4.0', 'tables>=3.6.1', 'scipy>=1.6.3',
        'imbalanced-learn>=0.8.1'
    ],
    scripts=['bin/selina'])
