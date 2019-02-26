# TiltRefiner
Finds a plane that fits the defocus values for all particles and removes particles that fall outside of some distance threshold.

Requires CtfRefine to be run or per-particle-defocus.

# Instructions
This requires pyem to read starfiles. Install that from:
https://github.com/asarnow/pyem/wiki/Install-pyem

Usage is hopefully straight forward:

/path/to/python TiltRefiner.py StarFile.star

see TiltRefiner.py --help for more details

Use the --test flag to see the first 10 micrographs to get a sense for what 
  threshold value you need to use
