# Read a GMT gene set file

Reads a tab-delimited GMT file (as used by MSigDB/GSEA: one gene set per
line, columns are id, description, gene1, gene2, ...) into a named list
suitable for
[`eem_search`](https://atusiniida.github.io/eemr/reference/eem_search.md).

## Usage

``` r
read_gmt(file)
```

## Arguments

- file:

  path to a GMT file.

## Value

a named list of character vectors (gene set id -\> gene ids).

## Examples

``` r
gmt_file <- system.file("extdata", "test.gmt", package = "eemr")
geneSets <- read_gmt(gmt_file)
length(geneSets)
#> [1] 100
geneSets[[1]]
#>   [1] "HNRPK"    "XRCC6"    "INTS1"    "INTS2"    "MED23"    "MED22"   
#>   [7] "ZNF638"   "INTS3"    "HDAC7A"   "MED21"    "ISG20"    "INTS9"   
#>  [13] "EPC1"     "INTS8"    "INTS5"    "INTS4"    "MED26"    "INTS7"   
#>  [19] "FOXF1"    "INTS6"    "FOXF2"    "U2AF1"    "GIT2"     "PITX2"   
#>  [25] "IKBKAP"   "MCRS1"    "POLI"     "POLH"     "PTBP1"    "MED12"   
#>  [31] "MED11"    "MED14"    "MED13"    "TAF6L"    "PPARGC1A" "HTATIP"  
#>  [37] "PPARGC1B" "SFRS2"    "MED16"    "MED17"    "MED18"    "LRCH4"   
#>  [43] "SUDS3"    "CPSF3L"   "MRE11A"   "FIGLA"    "CHEK2"    "COIL"    
#>  [49] "NKX2-1"   "NUP54"    "WDHD1"    "POLQ"     "NKX2-5"   "YEATS4"  
#>  [55] "SUB1"     "RTCD1"    "PHB"      "TP53BP1"  "SMAD4"    "NXF3"    
#>  [61] "SMAD2"    "C19orf2"  "PMF1"     "HNRNPA1"  "JAZF1"    "SRrp35"  
#>  [67] "CHAF1A"   "CHAF1B"   "NUAK2"    "ELF4"     "EZH2"     "NAP1L1"  
#>  [73] "BNIP3"    "NAP1L3"   "NAP1L2"   "NAP1L4"   "RPS19BP1" "ZNF350"  
#>  [79] "ASH2L"    "ZNF148"   "TOP2B"    "CCNO"     "TOP2A"    "ZNF281"  
#>  [85] "TFIP11"   "PRPF31"   "SYMPK"    "RBBP5"    "SNAPC4"   "PTF1A"   
#>  [91] "GTF2H4"   "TP53"     "APTX"     "NDUFA13"  "TTF2"     "GTF2H1"  
#>  [97] "MED6"     "TAF11"    "TAF10"    "MED7"     "MED4"     "TAF13"   
#> [103] "TAF12"    "MED8"     "RUVBL2"   "RBM39"    "RUVBL1"   "CLOCK"   
#> [109] "MED1"     "POLR2E"   "POLR2L"   "POLR2J"   "SHFM1"    "POLA1"   
#> [115] "PML"      "NUFIP1"   "ZBTB16"   "POLR2D"   "POLR2B"   "POLR2A"  
#> [121] "GATAD2A"  "ACTL6B"   "ACTL6A"   "POLR3G"   "POLR3H"   "POLR3K"  
#> [127] "HCLS1"    "TRIM27"   "PHF12"    "APPL1"    "POLR3C"   "APPL2"   
#> [133] "NARG1"    "PHF21A"   "MPHOSPH1" "OGG1"     "RBM14"    "CIR"     
#> [139] "FOXE3"    "C1orf124" "NCBP1"    "XPO1"     "LMO4"     "TBP"     
#> [145] "CBX1"     "SETX"     "TOP1"     "CDKN2A"   "SMARCD2"  "OIP5"    
#> [151] "SMARCD3"  "SMARCD1"  "EED"      "CIB1"     "MLL2"     "MYO6"    
#> [157] "RELA"     "FMR1"     "HDAC10"   "TOPBP1"   "ARID1A"   "HDAC11"  
#> [163] "ARID1B"   "ARNTL"    "HES6"     "SP140"    "ELL2"     "SUZ12"   
#> [169] "RPAIN"    "HIF1A"    "SMARCE1"  "SMARCA5"  "EDF1"     "SLU7"    
#> [175] "SMARCA2"  "SMARCA4"  "TRRAP"    "HNRNPL"   "SMARCB1"  "TAF9"    
#> [181] "ASF1A"    "TAF2"     "TAF1"     "TAF4"     "TAF6"     "EPAS1"   
#> [187] "TAF5"     "TAF8"     "SPTBN4"   "TAF7"     "SAP18"    "INTS12"  
#> [193] "WRN"      "INTS10"   "HDAC5"    "ATF6"     "HDAC4"    "HDAC3"   
#> [199] "HDAC2"    "PPIG"     "GCM1"     "HDAC1"    "CDKN2AIP" "SMARCC1" 
#> [205] "GTF2F1"   "SMARCC2"  "DYRK1A"   "MYO16"    "WDR4"     "TDG"     
#> [211] "HDAC9"    "HDAC8"    "HDAC6"    "RSF1"     "ARID4A"   "NR6A1"   
#> [217] "MKI67IP"  "TIMM50"   "MEN1"     "SAP30"    "MCM7"     "DKC1"    
#> [223] "NARG1L"   "PAX8"     "ELP4"     "ELP3"     "ELP2"     "SP100"   
#> [229] "CDK8"     "CDK9"     "UBE2I"    "RB1"      "IFI16"    "TOPORS"  
#> [235] "NCR1"     "HMGA1"    "CD3EAP"   "CCDC55"   "HIPK2"    "NCOA6"   
#> [241] "CPSF6"    "MDM2"     "CPSF3"    "KPNA2"    "CPSF1"    "FUSIP1"  
#> [247] "NUP98"    "ING2"     "GLIS2"    "NFYC"     "NFYB"     "NFYA"    
#> [253] "IVNS1ABP" "RPA1"     "ERCC5"    "NPAS2"    "ERCC6"    "NUDT21"  
#> [259] "GTF3C6"   "GTF3C5"   "SFRS2IP"  "ERCC3"    "GTF3C2"   "GEMIN6"  
#> [265] "GTF3C1"   "GEMIN7"   "GTF3C4"   "ERCC2"    "GTF3C3"   "GEMIN5"  
#> [271] "ACTB"     "ATXN1"    "PPP1R9B"  "MPG"      "ATXN3"    "MED30"   
#> [277] "THRAP3"   "PRM1"     "PRM2"    
```
