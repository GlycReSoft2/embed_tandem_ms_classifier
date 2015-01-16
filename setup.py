from setuptools import setup, find_packages

setup(
    name="GlycReSoft",
    version="1.0.2",
    packages=find_packages(),
    install_requires=[
        "scikit-learn >= 0.14.1",
        "pandas >= 0.14.0",
        "pyyaml >= 3.11",
        "pyteomics >= 2.5",
        "sqlitedict >= 1.1.0",
        "numexpr  >= 2.1",
        "xray >= 0.3.2"
    ],
    zip_safe=False,
    include_package_data=True,
    package_data={
        'glycresoft_ms2_classification': ["*.csv", "*.xml", "*.json", "data/*.csv"],
        'glycresoft_ms2_classification.structure': ["structure/data/*.csv", "structure/data/*.json"]
    },
    entry_points={
        'console_scripts': [
            "glycresoft-ms2 = glycresoft_ms2_classification.entry_point:main",
        ],
        'setuptools.installation': [
            "eggsecutable = glycresoft_ms2_classification.entry_point:main"
        ]
    },
    namespace_packages=["glycresoft_ms2_classification"]
)
