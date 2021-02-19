import flanker

from setuptools import setup, find_packages


setup(
    name = 'flanker',
    version = flanker.__version__,
    description = ' Gene-flank analysis tool ',
    url = 'https://github.com/wtmatlock/flanker',
    license = '',
    package_dir={'':'flanker'},
    py_modules=['cluster','salami'],
    python_requires='>=3.6',
    install_requires=['pandas', 'biopython', 'pytest'],
    entry_points = {'console_scripts':['flanker=flanker.flanker:main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English'])
