# Generate the summary table for the Bayesian model

Generate the summary table for the Bayesian model

## Usage

``` r
# S3 method for class 'refundBayes'
summary(object = NULL, ..., prob = 0.95)
```

## Arguments

- object:

  A fitted object returned by
  [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md).

- ...:

  Other parameters

- prob:

  Coverage probability for the reported confidence intervals. Defaults
  to 0.95.

## Value

A list of two objects, the first is the summary table for the estimated
scalar coefficients, the second is the plots for the estimated
functional coefficients.
