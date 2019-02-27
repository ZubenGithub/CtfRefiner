# CtfRefiner

Relion 3 CTF Refine compares projections from a refinement model to each particle
and by applying noise and a CTF can estimate individual defocus values.

If a particle does not resemble the refinement model (i.e., it is not really 
the same type of particle) then it may need to have a much larger defocus value
applied to cause the projection/particle to match.

By finding a plane that fits all particle defocus values for each micrograph 
it is possible to identify particles that are outside some distance threshold 
and remove them (since they have a defocus that is very different from all their
siblings).

# Instructions
This requires pyem to read starfiles. Install that from:
https://github.com/asarnow/pyem/wiki/Install-pyem

Usage is hopefully straight forward:

/path/to/python CtfRefiner.py --star_file StarFile.star

see CtfRefiner.py --help for more details

Use the --test flag to see the first 10 micrographs to get a sense for what 
  threshold value you need to use
