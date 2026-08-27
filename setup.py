from setuptools import setup, find_packages

setup(
    name="frp-flexure",           # this is the PyPI/distribution name
    version="1.0.0",
    packages=find_packages(),     # will find ["frp_flexure"]
    install_requires=["numpy"],
    extras_require={"examples": ["matplotlib"], "tests": ["pytest"]},
    python_requires=">=3.8",
)
