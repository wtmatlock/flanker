import flanker

from setuptools import setup


setup(
    name='flanker',
    version=flanker.__version__,
    description=' Gene-flank analysis tool ',
    url='https://github.com/wtmatlock/flanker',
    license='LICENSE',
    python_requires='>=3.6',
    packages=['flanker'],
    install_requires=['pandas', 'biopython', 'networkx', 'jellyfish'],
    entry_points={'console_scripts':['flanker=flanker.flanker:main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English'])
