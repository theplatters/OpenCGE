# src/closures/financing/financing.jl — financing closures (relocation placeholder)
#
# Phase 2, Part 1: only the baseline reference exists so the relocated kernel
# stays loadable. Concrete F1/F2/F3 hooks land here in Part 2 (ported from
# the frozen cbase2/src/financing.jl; never edited there).

"""Baseline reference: no programme, household untaxed, no additive demand."""
struct NoFinancing <: AbstractFinancing end
