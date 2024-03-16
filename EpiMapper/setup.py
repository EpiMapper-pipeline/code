# -*- coding: utf-8 -*-


"""setup.py: setuptools control."""


from setuptools import setup

with open("README.txt", "rb") as f:
    long_descr = f.read().decode("utf-8")


setup(
    name = "epimapper",
    packages = [
        "epimapper",
        "epimapper.function_scripts",
        "epimapper.function_scripts.dar_scripts",
        "epimapper.function_scripts.dar_scripts.scripts_high",
        "epimapper.seacr_tool"],
    package_data={'epimapper': ['seacr_tool/*']},
    entry_points = {
        "console_scripts": ['epimapper = epimapper.epimapper:main']
        },
    version = 1.0,
    description = "Python Package for analysing epigenetic data such as ChIP-seq, CUT&RUN, ATAC-seq and CUT&Tag",
    long_description = long_descr,
    author = "Jenny Sofie Dragland",
    author_email = "jennysodr@gmail.com",
    )
