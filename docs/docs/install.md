# Installation

## Conda

This is the recommended method and by far the easiest. It will be available soon.

## During development

```
conda create -n flanker python=3 abricate biopython mash  # needs bioconda channel
pip install --editable git+https://github.com/wtmatlock/flanker
pip install --editable /path/to/flanker/package  # where 'package' contains setup.py
```

### Run tests

```
pytest
```

## From GitHub

Ensure you have all dependencies installed:

**Python dependencies:**

* multiprocessing
* sys
* argparse
* pandas
* numpy
* subprocess
* biopython
* pathlib
* time
* logging
* tempfile
* os
* collections
* glob
* networkx

**External software:**

* abricate
* mash

Then simply clone the repository:

```
  git clone https://github.com/wtmatlock/flanker
```

and check everything is working:

```
  python flanker.py --help
```
