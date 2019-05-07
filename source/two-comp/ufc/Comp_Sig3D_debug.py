#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
# Copyright (C) 2017 Van-Dang Nguyen (vdnguyen@dtth.se)

# This file is part of FEniCS
#
# FEniCS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

Ve = FiniteElement("CG", tetrahedron, 1)
TH = MixedElement([Ve,Ve])
V_DG = FiniteElement("DG", tetrahedron, 0)

phase = Coefficient(V_DG)

u = Coefficient(TH)
M = u[0]*(1-phase)*dx + u[2]*phase*dx

