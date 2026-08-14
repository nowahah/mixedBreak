---
title: "mixedBreak: Mixed-effects Models with Unknown Breakpoints"
---

Modelling segmented relationships with linear mixed-effects model with breakpoint(s) estimates.


## Installation (R >= 4.3.1)

You can use `pak` package to install the latest stable version of `mixedBreak`:

```r
pak::pak("nowahah/mixedBreak@v0.0.0") # TODO - change when stable version out
```


## Fitting a Segmented Linear Mixed-effects Model

The functionalities of the package will be exemplified on the following dataset:

```r
library(mixedBreak)
data(SDIpsilo, package = "lmbreak") # TODO - update before release
str(SDIpsilo)
```

where the experience of 15 individuals after drug intake is monitored over time:

![Subjective Drug Intensity in 15 healthy individuals after drug intake (t=0min)](https://github.com/nowahah/mixedBreak/tree/main/inst/figures/gg-SDI-all.png)

The `mixed.break` function can be used to model the experience of the patients using a segmented mixed model. This model not only estimates the regression coefficients in every segment, but also the breakpoint values themselves, for every individual.

```r
break.fit <- mixed.break(
  score ~ 0 + time + (0 + time|ID),
  pattern = "101",
  data = SDIpsilo
)
```

The syntax was designed to remain close to a `lme4::lmer` call, to which one adds a `pattern` argument, used for specifying whether a non-null slope has to be estimated in the corresponding segment:

<!--- 
- `pattern=11` fits 1 breakpoints and assumes every slope is to be estimated (unconstrained) ;
- `pattern=111` fits 2 breakpoints (unconstrained) ;
- `pattern=1111` fits 3 breakpoints (unconstrained) ; 
--->
- `pattern=10` fits 1 breakpoint with a constant segment in the second phase ;
- `pattern=101` fits 2 breakpoints with a constant segment in the second phase ;
- `pattern=1010` fits 3 breakpoints with a constant segment in phases 2 & 4.

Available patterns either one of `c('10','101', '1010')`. Unconstrained model with 1 to 3 breakpoints are currently under development.



```r
plot(break.fit, fit.avg = TRUE)
```

![Individual fits (orange lines) and average fit (red line) of the segmented mixed model - `pattern 101` ](https://github.com/nowahah/mixedBreak/tree/main/inst/figures/gg-SDI-101-avgfit.png)



## Model formulation and estimation

The estimation procedures is mostly based on the work on 2014 Muggeo's article (see references). The formulation of the model, in its simplest case of `pattern = 11`, is the following:

![equation](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D%20y_%7Bij%7D=%5Cbeta_%7B1i%7Dt_%7Bij%7D&plus;%5Cvarepsilon_%7Bij%7D)

<!--- very ugly but works --->

It is essentially an iterative procedure over approximated linearized mixed-effects models.

## Limitations

For the sake of the model's estimation, one of the assumption for the breakpoints is they are ordered, _i.e._ $ \psi_{1} < $

Switched breakpoints

- Message / Warnings explained for this situation

With real-life data, it comes unsurprisingly that a single pattern cannot always fit every individual for a given phenomena. 
<!--- In the event that many individuals does not fit to the presumed pattern, they may be modeled separately with a different pattern. --->

- The model typically assumes that breakpoints are ordered to make an interpretation of its output. However, there are no mathematical constraints ensuring this order. Therefore, it could happen that some breakpoints of some individuals are switched ( _i.e._ )
- No intercept is allowed in the model for now ;
- No other covariates than the segmented variables is allowed in the model ;
- The model showed limited statistical properties (bias, undercoverage) in a simulation study carried on for a single breakpoint and across realistic Data-Generating Mechanisms.


## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.


## References

1. Segmented mixed models with random changepoints: a maximum likelihood approach with application to
treatment for depression study. Statistical Modelling. 2014 Aug ;14(4):293–313. doi:10.1177/1471082X13504721
2. segmented: Regression Models with Break-Points / Change-Points Estimation (with Possibly Random
Effects) [Internet]. 2003 [cited 2026 Mar 24]. Available from: https://cran.r-project.org/package=segmented
doi:10.32614/CRAN.package.segmented


## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)