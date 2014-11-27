from setuptools import setup, find_packages

setup(
    name="GlycReSoft",
    version="1.0.a.dev",
    packages=find_packages("glycresoft_ms2_classification"),
    install_requires=[
        "sklearn >= 0.14.1",
        "pandas >= 0.14.0",
        "yaml >= 3.11"
    ],
    include_package_data=True,
    package_data={
        '': ["**.csv", "*.xml", "*.json"],
        'structure': ["data/*.csv"]
    },
    entry_points={
        'console_scripts': [
            "glycresoft-ms2 = glycresoft_ms2_classification.entry_point:main"
        ]
        #,
        # 'setuptools.installation': [
        #     "eggsecutable = glycresoft_ms2_classification.entry_point:main"
        # ]
    }

)
