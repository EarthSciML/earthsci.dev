# This file was generated, do not modify it. # hide
discretization = MOLFiniteDifference([lon => 50], t, approx_order=2)
prob = discretize(get_mtk(model), discretization)
sol = solve(prob, TRBDF2(), saveat=3600.0)