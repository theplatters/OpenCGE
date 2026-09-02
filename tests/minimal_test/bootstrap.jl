# bootstrap.jl — shared bootstrap for standalone diagnostics in tests/minimal_test.
#
# Activates the repository project (located two directories above this file),
# imports the real BeyondHulten package, and points the working directory at
# the repository root so that read_data(filename) (which joins pwd()/"data"/filename)
# resolves to data/I-O_DE2019_formatiert.csv.
#
# Common dependencies used by the diagnostics are also imported here so each
# script need not repeat them.

using Pkg

const REPO_ROOT = let
    # This file lives in <repo>/tests/minimal_test
    dir = @__DIR__
    dirname(dirname(dir))   # climb two levels to the repository root
end

Pkg.activate(REPO_ROOT)
using BeyondHulten

# Standard-library / package imports used by the diagnostics.
using Printf
using Statistics
using LinearAlgebra
using DataFrames
using NonlinearSolve: NonlinearSolve

# Ensure read_data (which uses pwd()) finds data/I-O_DE2019_formatiert.csv
cd(REPO_ROOT)
