<!--- DO NOT EDIT:  AUTOMATICALLY GENERATED from README.Rmd -->

# SOAs

Creates stratum orthogonal arrays (also known as strong orthogonal
arrays).

  - **Author**: Ulrike Groemping, BHT Berlin.
  - **Contributor**: Rob Carnell.

|                                                                                                                                                                Actions                                                                                                                                                                 |                                                                      Coverage                                                                      |                                            Website                                             |
| :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------: |
| [![R-CMD-check](https://github.com/bertcarnell/SOAs/actions/workflows/r_cmd_check.yml/badge.svg)](https://github.com/bertcarnell/SOAs/actions/workflows/r_cmd_check.yml)[![pkgdown](https://github.com/bertcarnell/SOAs/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/bertcarnell/SOAs/actions/workflows/pkgdown.yaml) | [![Codecov test coverage](https://codecov.io/gh/bertcarnell/SOAs/branch/main/graph/badge.svg)](https://codecov.io/gh/bertcarnell/SOAs?branch=main) | [![](https://img.shields.io/badge/pkgdown-SOAs-blue.svg)](https://bertcarnell.github.io/SOAs/) |

## Installation

`SOAs` has not yet been released on [CRAN](https://CRAN.R-project.org).

You can install the development version of `SOAs` from here with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("bertcarnell/SOAs")
```

### System Dependencies

For the `arrangements` package:

  - Ubuntu: `apt-get update; apt-get install libgmp3-dev -y`
  - MacOS: `brew install gmp`

## Details

<p>

This package constructs arrays in <i>s^el</i> levels from orthogonal
arrays in s levels. These are all based on equations of the type<br />
<i>D = s^{el-1} A\_1 + … + s A\_{el-1} + A\_{el}</i>, or <br /> for
<i>s^2</i> levels, <i>D = s A + B</i> and <br /> for <i>s^3</i> levels,
<i>D = s^2 A + s B + C</i>. <br /> The constructions differ in how they
obtain the ingredient matrices, and what properties can be guaranteed
for the resulting D. Where a construction function guarantees orthogonal
columns for all matrices D it produces, its name starts with a OSOA,
otherwise with SOA.<br />

</p>

<p>

If optimization is requested (default TRUE), space filling properties of
D are improved using a level permutation algorithm by Weng (2014). This
algorithm is applied for improving the <code>phi\_p</code> criterion,
which is often a reasonable surrogate for increasing the minimum
distance.

</p>

<p>

Groemping (2021) describes the constructions by He and Tang (2013,
function <code>SOAs</code>), Liu and Liu (2015, function
<code>OSOAs\_LiuLiu</code>), He, Cheng and Tang (2018, function
<code>SOAs2plus\_regular</code>), Zhou and Tang (2019), Shi and Tang
(2020, function <code>SOAs\_8level</code>) and Li, Liu and Yang (2021)
in unified notation. The constructions by Zhou and Tang (2019) and Li et
al. (2021) are very close to each other and are both implemented in the
three functions <code>OSOAs</code>, <code>OSOAs\_hadamard</code> and
<code>OSOAs\_regular</code>.

</p>

<p>

Besides the construction functions, properties of the resulting array D
can be checked using the aforementioned function <code>phi\_p</code> as
well as check functions <code>ocheck</code>, <code>ocheck3</code> for
orthogonality and <code>soacheck2D</code>, <code>soacheck3D</code> for
(O)SOA stratification properties.

</p>

<p>

There is one further construction, maximin distance level expansion
(<code>XiaoXuMDLE</code>, <code>MDLEs</code>), that does not yield
stratum (aka strong) orthogonal arrays and is available for comparison
only (Xiao and Xu 2018).

</p>

## References

<p>

Groemping, U. (2021 - soon). A unified implementation of stratum (aka
strong) orthogonal arrays. Report 2021/01, Reports in Mathematics,
Physics and Chemistry, Berliner Hochschule für Technik.
<a href="http://www1.beuth-hochschule.de/FB_II/reports/Report-2021-001.pdf">http://www1.beuth-hochschule.de/FB\_II/reports/Report-2021-001.pdf</a>.

</p>

<p>

He, Y., Cheng, C.S. and Tang, B. (2018). Strong orthogonal arrays of
strength two plus. <em>The Annals of Statistics</em> <b>46</b>, 457-468.

</p>

<p>

He, Y. and Tang, B. (2013). Strong orthogonal arrays and associated
Latin hypercubes for computer experiments. <em>Biometrika</em>
<b>100</b>, 254-260.

</p>

<p>

Li, W., Liu, M.-Q. and Yang, J.-F. (2021, in press). Construction of
column-orthogonal strong orthogonal arrays. <em>Statistical Papers</em>
<a href="https://doi.org/10.1007/s00362-021-01249-w">https://doi.org/10.1007/s00362-021-01249-w</a>.

</p>

<p>

Liu, H. and Liu, M.-Q. (2015). Column-orthogonal strong orthogonal
arrays and sliced strong orthogonal arrays. <em>Statistica Sinica</em>
<strong>25</strong>, 1713-1734.

</p>

<p>

Shi, L. and Tang, B. (2020). Construction results for strong orthogonal
arrays of strength three. <em>Bernoulli</em> <strong>26</strong>,
418-431.
<a href="https://doi.org/10.3150/19-BEJ1130">https://doi.org/10.3150/19-BEJ1130</a>

</p>

<p>

Weng, J. (2014). Maximin Strong Orthognal Arrays. <em>Master’s
thesis</em> at Simon Fraser University under supervision of Boxin Tang
and Jiguo Cao.
<a href="https://summit.sfu.ca/item/14433">https://summit.sfu.ca/item/14433</a>

</p>

<p>

Xiao, Q. and Xu, H. (2018). Construction of Maximin Distance Designs via
Level Permutation and Expansion. <em>Statistica Sinica</em> <b>28</b>,
1395-1414.

</p>

<p>

Zhou, Y.D. and Tang, B. (2019). Column-orthogonal strong orthogonal
arrays of strength two plus and three minus. <em>Biometrika</em>
<strong>106</strong>, 997-1004.
<a href="https://doi.org/10.1093/biomet/asz043">https://doi.org/10.1093/biomet/asz043</a>

</p>
