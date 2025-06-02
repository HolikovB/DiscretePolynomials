
import setuptools

setuptools.setup(
    name="discpoly",
    version="0.1.0",
    description="Polynomial‐function tools on finite abelian groups (e.g. Z/4×Z/4 → R/Z)",
    author="Borys Holikov",
    author_email="boris2107g@gmail.com",
    url="https://github.com/yourusername/polyfunc",  # update as appropriate
    packages=setuptools.find_packages(),  # will include the `polyfunc/` folder
    python_requires=">=3.7",

    # No install‐time requirements for core functionality:
    install_requires=[
        # (leave empty if there are really no runtime dependencies)
    ],

    # If you want to allow `python setup.py test`, you can use tests_require:
    tests_require=[
        "pytest>=7.0.0"
    ],

    # Optional metadata:
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)