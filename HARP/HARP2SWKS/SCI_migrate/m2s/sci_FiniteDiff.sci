function [tree] = sci_FiniteDiff(tree)
// Copyright INRIA (Generated by M2SCI)
// Conversion function for Matlab FiniteDiff()
// Input: tree = Matlab funcall tree
// Output: tree = Scilab equivalent for tree
tree.lhs(1).dims=list(-1,1)
tree.lhs(1).type=Type(Double,Unknown)
endfunction
