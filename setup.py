from setuptools import setup, find_packages

setup(
    name='parseq_analyze',
    version='0.1.1',
    packages=find_packages(),
    install_requires=[
        
        'numpy',
        'pandas',
        'seaborn',
        'matplotlib',
        'biopython',
        'edlib'
    ],
    python_requires='>=3.8',
    author='Moustafa Houmani',
    author_email='mhoumani112@gmail.com',
    description='parseq is a python package for the analysis of parSEQ sequencing data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/m-houm/parseq_analyze',
    classifiers=[
    'Development Status :: 4 - Beta',  # Adjust as needed
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3 :: Only',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Framework :: Jupyter',
    'Environment :: Console',
    
],
)
