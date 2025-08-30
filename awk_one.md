
## AWK function for absolute value:
awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($9) < 500) print $0}'
https://stackoverflow.com/questions/11184915/absolute-value-in-awk-doesnt-work
