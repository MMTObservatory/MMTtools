import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
        name="MMTtools",
        version="0.1",
        author="Chun Ly",
        author_email="astro.chun@gmail.com",
        description="MMTO tools",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/astrochun/MMTtools",
        packages=setuptools.find_packages(),
        install_requires=[
            'astropy',
            'photutils',
            'ccdproc',
            'pymysql',
        ],
        classifiers=[
            "Programming Language :: Python",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
)
