# Plot the estimated functional coefficients with the corresponding credible interval(s).

Plot the estimated functional coefficients with the corresponding
credible interval(s).

## Usage

``` r
# S3 method for class 'refundBayes'
plot(x = NULL, ..., prob = 0.95, include = "both")
```

## Arguments

- x:

  A fitted object returned by
  [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md).

- ...:

  Other parameters

- prob:

  Coverage probability for the credible interval(s). Defaults to 0.95.

- include:

  Type of interval to include. `"pointwise"` produces pointwise credible
  intervals; `"CMA"` produces the CMA credible band; `"both"` produces
  both. Defaults to `"both"`.

## Value

A list of `ggplot` objects, one for each functional coefficient.
