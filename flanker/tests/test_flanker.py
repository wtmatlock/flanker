import os
import subprocess

import pytest

from flanker import flanker


cwd = os.getcwd()
data_dir = 'flanker/tests/data'  # Test from project root

def run(cmd, cwd=cwd, check=True):  # Helper for CLI testing
    return subprocess.run(cmd,
                          cwd=cwd,
                          shell=True,
                          check=check,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

def test_cli_blactxm():
	run('flanker -i RHB02-C01_4.fasta -g blaCTX-M-55', cwd=data_dir)
	run('rm -f RHB02-C01_4_blaCTX-M-55_1_1000_both_flank.fasta RHB02-C01_4.fasta_resfinder', cwd=data_dir)  # Tidy up

# def test_default():
# 	run('flanker -i test.fasta -g blaKPC-2 -w 100', cwd=data_dir)
# 	run('rm -f *flank.fasta *_resfinder', cwd=data_dir)

# def test_mm():
# 	run('flanker -i test.fasta -g blaKPC-2 blaKPC-3 -m mm -w 100', cwd=data_dir)
# 	run('rm -f *flank.fasta *_resfinder', cwd=data_dir)

# def test_circ_inc():
# 	run('flanker -i test.fasta -g blaKPC-2 -w 100 -circ -inc', cwd=data_dir)
# 	run('rm -f *flank.fasta *_resfinder', cwd=data_dir)

# def test_cluster():
# 	run('flanker -i test.fasta -g blaKPC-2 -w 100 -cl', cwd=data_dir)
# 	run('rm -f *flank.fasta *_resfinder out*', cwd=data_dir)

# def test_salami():
# 	run('flanker -i test.fasta -g blaKPC-2 -w 100 -m sm', cwd=data_dir)
# 	run('rm -f *flank.fasta *_resfinder', cwd=data_dir)