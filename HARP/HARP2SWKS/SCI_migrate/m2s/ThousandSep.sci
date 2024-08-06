function [out] = ThousandSep(in)

// Output variables initialisation (not found in input variables)
out=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//THOUSANDSEP adds thousands Separators to a 1x1 array.
// Example:
// ThousandSep(1234567)
// !! L.5: Matlab function import not yet converted, original calling sequence used.
import("java.text.*")
// ! L.6: mtlb(DecimalFormat) can be replaced by DecimalFormat() or DecimalFormat whether DecimalFormat is an M-file or not.
v = mtlb(DecimalFormat);
out = char(v.format(in));

endfunction
