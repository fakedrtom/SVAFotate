from setuptools import setup
import re

with open("svafotate/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)


setup(
    name="svafotate",
    version=version,
    description="Annotaton tools for structural variant annotations",
    long_description=open("README.md").read(),
    author="",
    author_email="",
    url="https://github.com/mikecormier/SVAFotate",
    packages=["svafotate"],
    package_data={"": ["LICENSE", "README.md"]},
    package_dir={"svafotate": "svafotate"},
    include_package_data=True,
    python_requires=">=3",
    license="MIT",
    zip_safe=False,

    entry_points={
        "console_scripts": [
            "svafotate = svafotate.__main__:main"
        ]
    },

    classifiers=[
        "Programming Language :: Python :: 3"
        "License :: MIT License",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"]
)

