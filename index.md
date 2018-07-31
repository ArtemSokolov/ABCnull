# ABCmonster

---

We present companion materials for our paper titled _Characterizing properties of ABC-transporter substrates using a clean-slate genetic background_. All code and data used in the paper is encapsulated inside an R package. The package can be installed by running the following commands within R:

``` r
if( !require(devtools) ) install.packages( "devtools" )
devtools::install_github( "ArtemSokolov/ABCmonster" )
```

The package provides several functions that simplify reproducibility of the analyses described in the paper. These functions are outlined in the [Reference](reference/index.html) section. All other scripts are described across three vignettes:

  * [Overview of raw data](articles/rawdata.html), where we describe the training data and present several unsupervised analyses;
  * [Cross-validation](articles/sv.html), where we evaluate several prediction methods in a cross-validation setting;
  * and [Prospective validation](articles/prosp.html), where we describe novel drug sensitivity data that we collected on the ABC-16 strain and compare predictions made by the model against this data.
