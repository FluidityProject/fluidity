# Python scripts will use Fluidity binaries associated with the present
# directory unless specified otherwise.  When running tests which are
# remote from the branch, the user will need to point FLUIDITYPATH back
# to this branch.
FLUIDITYPATH ?= $(CURDIR)/../../

export PYTHONPATH := $(PYTHONPATH):$(FLUIDITYPATH)/python:$(FLUIDITYPATH)/tools:$(FLUIDITYPATH)/tests/darcy_impes_common
