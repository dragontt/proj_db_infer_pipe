#!/usr/bin/python
import numpy

# load data
dir_expr = 'resources/fly_netowrk_zeke_singles_net_full/'
fn_expr = dir_expr + 'data.singles.expr'
fn_expr_add_pseudocount = dir_expr + 'data.singles.add_pseudocount.expr'

expr = numpy.loadtxt(fn_expr)

print numpy.max(expr)
print numpy.min(expr)

# compute pseudocount
expr_indices = numpy.transpose(numpy.nonzero(expr))
expr_nonzero = numpy.zeros(len(expr_indices))
for i in range(expr_indices.shape[0]):
	expr_nonzero[i] = expr[expr_indices[i,0], expr_indices[i,1]]
pseudocount = numpy.percentile(expr_nonzero, 0.5)
print pseudocount
expr += pseudocount

# write expr with pseudocount added
# numpy.savetxt(fn_expr_add_pseudocount, expr, fmt='%.5e', delimiter=' ', newline='\n')

