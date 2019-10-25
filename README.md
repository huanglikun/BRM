# <a name="top"></a>Block Regression Mapping (BRM)
**Block Regression Mapping (BRM)** is a statistical method for QTL mapping based on bulked segregant analysis by deep sequencing. The core function is programmed by R language.

## Content
* [Introduction](#intro)
* [Getting started](#getstart)
* [Input data format](#inputdata)
  * [Example](#inputdataexample)
* [Input configuration files](#inputconf)
  * [Chromosome length file](#inputconfchrlen)
    * [Example](#inputconfchrlenexample)
  * [Block regression configuration file](#inputconfblk)
    * [Example](#inputconfblkexample)
  * [Population information configuration file](#inputconfpop)
    * [Example](#inputconfpopexample)
* [Output files explanation](#output)
  * [BRM step 1 output: the statistics](#outputstep1)
    * [Example](#outputstep1example)
  * [BRM step 2 output: the threshold](#outputstep2)
    * [Example](#outputstep2example)
  * [BRM step 3 output: the peaks and the confidence intervals](#outputstep3)
    * [Example](#outputstep3example)

## <a name="intro"></a>Introduction
Bulked segregant analysis by deep sequencing (BSA-seq) has been widely used for QTL mapping in recent 
years. A number of different statistical methods for BSA-seq have been proposed. However, determination 
of significance threshold, the key point for QTL identification, remains to be a problem that has not been 
well solved due to the difficulty of multiple test correction. In addition, estimation of the confidence interval is also a problem to be solved.
Here, we propose a new statistical method for BSA-seq, named Block Regression Mapping (BRM). BRM is 
robust to sequencing noise and is applicable to the case of low sequencing depth. Significance threshold 
can be reasonably determined by taking multiple test correction into account. Meanwhile, the confidence 
interval of QTL position can also be estimated.

 [back to top](#top)

## <a name="getstart"></a>Getting started
Get all scripts and the examples from github and have a try:
```bash
git clone https://github.com/huanglikun/BRM.git
cd BRM
mkdir result
Rscript BRMstep1.loess.R configureExample/loess_conf.txt configureExample/chr_length.txt dataExample/yeast_markers_dp50.bsa result/AF1.xls AF1
Rscript BRMstep1.loess.R configureExample/loess_conf.txt configureExample/chr_length.txt dataExample/yeast_markers_dp50.bsa result/AF2.xls AF2
Rscript BRMstep1.loess.R configureExample/loess_conf.txt configureExample/chr_length.txt dataExample/yeast_markers_dp50.bsa result/AAF.xls AAF
Rscript BRMstep1.loess.R configureExample/loess_conf.txt configureExample/chr_length.txt dataExample/yeast_markers_dp50.bsa result/AFD.xls AFD
Rscript BRMstep2.threshold.R configureExample/threshold_conf.txt result/AAF.xls result/AF1.xls result/AF2.xls result/threshold.xls
Rscript BRMstep3.peak_and_CI.R result/AFD.xls result/threshold.xls result/allpeaks.xls
```

 [back to top](#top)
 
## <a name="inputdata"></a>Input data format
We use a **tab** separated file named bsa format as input data. There are six columns in this file:

| Chromosome Name | Position | a | b | c | d |
| --- | --- | --- | --- | --- | --- |

a, b, c, d stand for the allele counts of one parent for each pool.

|  | Parent1 allele | Parent2 allele |
| :---: | :---: | :---: |
| High pool | a | b |
| Low pool | c | d |

  ### <a name="inputdataexample"></a>Example
  **dataExample/yeast_markers_dp50.bsa**

      I	33070	30	19	19	26
      I	33147	37	19	22	27
      I	33152	32	20	17	32
      I	33200	21	27	27	31
      I	33293	25	19	35	35
      ......

 [back to top](#top)

## <a name="inputconf"></a>Input configuration files
* ### <a name="inputconfchrlen"></a>Chromosome length file
    It's a **tab** separated file with two columns:
    
    | Chromosome Name | Chromosome Length |
    | --- | --- |
    
    ### <a name="inputconfchrlenexample"></a>Example
    **configureExample/chr_length.txt**
       
        I    230218
        II    813184
        III    316620
        ......

 [back to top](#top)

* ### <a name="inputconfblk"></a>Block regression configuration file
    It's a **key-value** file contains block regression parameters. The separator is "=". And the space between key and value will be ignored.
    There are five parameters needed to be set:
    
  | Key | Value type | Description | Recommend value |
  | --- | --- | --- | --- |
  | UNIT | Integer | The unit of block size (bp).  | default: 1000 |
  | DEG | Integer | The degree of the polynomials to be used in Local Polynomial Regression Fitting. | 2 |
  | BLK | Integer or float | The block size is equal to BLK * UNIT. | e.g.: 0.2 for yeast, 20 for rice |
  | MIN | Integer | Min total depth in one valid block. | 1 |
  | MINVALID | Integer | Min valid block number in one chromosome. | 10 |
    
  ### <a name="inputconfblkexample"></a>Example
   **configureExample/loess_conf.txt**

        UNIT     = 1000
        DEG      = 2
        BLK      = 0.2
        MIN      = 1  # min depth in block
        MINVALID = 10 # min valid blocks in one chromosome (needed to be at least 10)

 [back to top](#top)

*  ### <a name="inputconfpop"></a>Population information configuration file
    It's a **key-value** file contains experimental design information. The separator is "=". And the space between key and value will be ignored.
    There are four parameters needed to be set:
    
   | Key | Value type | Description | Recommend value |
   | --- | --- | --- | --- |
   | n1 | Integer | The individuals in pool1 (high pool).  | it depends on experiment design |
   | n2 | Integer | The individuals in pool2 (low pool).  | it depends on experiment design |
   | t | 0,1 | Population type.  | For DH or RI etc., t=0; ![equation](https://latex.codecogs.com/gif.latex?F_2) or ![equation](https://latex.codecogs.com/gif.latex?F_3) etc., t=1 |
   | ua | float | The ![equation](https://latex.codecogs.com/gif.latex?u_%7b%5calpha%2f2%7d) value.  | see the table below |
   
   The values of ![equation](https://latex.codecogs.com/gif.latex?u_%7b%5calpha%2f2%7d) of some species 
   
      <table>
      <tr>
      <th rowspan=2>Species</th>
    <th rowspan=2>n<sup>a</sup></th>
      <th colspan=2>Genome size<sup>b</sup></th>
      <th rowspan=2>Ratio</br>(kb/cM)</th>
      <th colspan=4><img src="https://latex.codecogs.com/gif.latex?u_%7b%5calpha%2f2%7d"></img><sup>  c</sup></th>
      <th rowspan=2>Ref.</th>
    </tr>
    <tr>
      <th>cM</th>
      <th>Mb</th>
    <th>H/DH/F<sub>2</sub></th>
      <th>F<sub>3</sub></th>
      <th>F<sub>4</sub></th>
      <th>RIL</th>
    </tr>
    <tr>
     <td><i>Arabidopsis</i></td>
      <td>5</td>
      <td>600</td>
      <td>119</td>
      <td>199</td>
      <td>3.41</td>
      <td>3.50</td>
      <td>3.54</td>
      <td>3.57</td>
      <td>[4]</td>
    </tr>
    <tr>
      <td>Cucumber</td>
      <td>7</td>
      <td>1390</td>
      <td>192</td>
      <td>138</td>
      <td>3.62</td>
      <td>3.72</td>
      <td>3.75</td>
      <td>3.78</td>
      <td>[9]</td>
    </tr>
    <tr>
      <td>Maize</td>
      <td>10</td>
      <td>2060</td>
      <td>2106</td>
      <td>1023</td>
      <td>3.72</td>
      <td>3.81</td>
      <td>3.85</td>
      <td>3.87</td>
      <td>[3]</td>
    </tr>
    <tr>
      <td>Rapeseed</td>
      <td>18</td>
      <td>2520</td>
      <td>855</td>
      <td>339</td>
      <td>3.78</td>
      <td>3.87</td>
      <td>3.90</td>
      <td>3.93</td>
      <td>[5]</td>
    </tr>
    <tr>
      <td>Rice</td>
      <td>12</td>
      <td>1530</td>
      <td>382</td>
      <td>250</td>
      <td>3.65</td>
      <td>3.74</td>
      <td>3.78</td>
      <td>3.80</td>
      <td>[6]</td>
    </tr>
    <tr>
      <td>Tobacco</td>
      <td>12</td>
      <td>3270</td>
      <td>3613</td>
      <td>1105</td>
      <td>3.84</td>
      <td>3.92</td>
      <td>3.96</td>
      <td>3.98</td>
      <td>[1]</td>
    </tr>
    <tr>
      <td>Tomato</td>
      <td>12</td>
      <td>1470</td>
      <td>807</td>
      <td>549</td>
      <td>3.65</td>
      <td>3.73</td>
      <td>3.77</td>
      <td>3.80</td>
      <td>[7]</td>
    </tr>
    <tr>
      <td>Wheat</td>
      <td>21</td>
      <td>3140</td>
      <td>14547</td>
      <td>4633</td>
      <td>3.83</td>
      <td>3.92</td>
      <td>3.95</td>
      <td>3.98</td>
      <td>[8]</td>
    </tr>
    <tr>
      <td>Yeast</td>
      <td>16</td>
      <td>4900</td>
      <td>12</td>
      <td>2.5 </td>
      <td>3.93</td>
      <td>4.02</td>
      <td>4.05</td>
      <td>4.08</td>
      <td>[2]</td>
    </tr>
    </table>

   Note: **a.** n, number of chromosomes in haploid. **b.** The genetic map length (cM) of each species was 
from the references listed except for that of yeast, which was calculated by us using the data from [2]. The genome size (Mb) of each species was all from NCBI. **c.** Corresponding to the overall 
significance level of 0.05.

**References**

[1] Bindler,G. *et al.* (2011) A high density genetic map of tobacco (*Nicotiana tabacum* L.) obtained 
from large scale microsatellite marker development. *Theor Appl Genet*, **123**, 219.

[2] Bloom,J.S. *et al.* (2015) Genetic interactions contribute less than additive effects to quantitative 
trait variation in yeast. *Nat Commun*, **6**, 8712.

[3] Civardi,L. *et al.* (1994) The relationship between genetic and physical distances in the cloned 
a1-sh2 interval of the Zea mays L. genome. *Proc Natl Acad Sci U S A*, **91**, 8268–8272.

[4] Garcia-Hernandez,M. *et al.* (2002) TAIR: a resource for integrated Arabidopsis data. *Funct 
Integr Genomics*, **2**, 239–253.

[5] Raman,H. *et al.* (2014) SNP markers-based map construction and genome-wide linkage analysis 
in Brassica napus. *Plant Biotechnology Journal*, **12**, 851–860.

[6] International Rice Genome Sequencing Project (2005) The map-based sequence of the rice 
genome. *Nature*, **436**, 793.

[7] Shirasawa,K. *et al.* (2010) SNP Discovery and linkage map construction in cultivated tomato. 
*DNA Research*, **17**, 381–391.

[8] Yang,Q. *et al.* (2018) High-density genetic map construction and mapping of the homologous 
transformation sterility gene (hts) in wheat using GBS markers. *BMC Plant Biology*, **18**, 301.

[9] Zhou,Q. *et al.* (2015) A sequencing-based linkage map of cucumber. *Molecular Plant*, **8**, 
961–963.
   
   ### <a name="inputconfpopexample"></a>Example
   **configureExample/threshold_conf.txt**
    
        n1 = 300
        n2 = 300
        t  = 0        # For DH or RI etc., t=0; F2 or F3 etc., t=1
        ua = 4.08    # For rice, F2:3.65; F3:3.74; F4:3.78 .
        
 [back to top](#top)
    
## <a name="output"></a>Output files explanation
* ### <a name="outputstep1"></a>BRM step 1 output: the statistics
    In BRM step 1 (BRMstep1.loess.R), four statistics will be calculated.
   * Pool1 (high pool) allele frequency file (AF1.xls)
   * Pool2 (low pool) allele frequency file (AF2.xls)
   * The **A**verage **A**llele **F**requency between Pool1 and Pool2 (AAF.xls)
   * The **A**llele **F**requency **D**ifference between Pool1 and Pool2 (AFD.xls)
         
    The first two files include seven columns:
    
    | Chromosome Name | Position | Block Average | Fitted Average | Std. Error | No. of Markers | Total depth in one block |
    | --- | --- | --- | --- | --- | --- | --- |
    
    The last two files include six columns:
    
    | Chromosome Name | Position | Block Average | Fitted Average | Std. Error | No. of Markers
    | --- | --- | --- | --- | --- | --- |
    
    ### <a name="outputstep1example"></a>Example
    **result/AF1.xls**
    
      ......
      I	178200	0.5625	0.544412954114093	0.00467274528845156	1	48
      I	178400	NA	0.544531420628367	0.00468496414921583	0	0
      I	178600	0.580310880829015	0.544650496651613	0.00469761880431398	4	193
      I	178800	0.618181818181818	0.544770181325821	0.00471071099925492	1	55
      I	179000	0.490066225165563	0.544890473792982	0.00472424237513985	3	151
      I	179200	0.5	0.545011373195085	0.00473821446883456	1	32
      ......
    
    **result/AF2.xls**
    
      ......
      I	178200	0.446428571428571	0.470160746538764	0.00523407229630387	1	56
      I	178400	NA	0.470451065105337	0.00522434327573345	0	0
      I	178600	0.474576271186441	0.470742430395541	0.00521502608500492	4	177
      I	178800	0.510204081632653	0.471034730131143	0.00520612879591068	1	49
      I	179000	0.491017964071856	0.471327852033912	0.00519766006364739	3	167
      I	179200	0.431372549019608	0.471621683825613	0.00518962917866429	1	51
      ......
    
    **result/AAF.xls**
    
      ......
      I	178200	0.504464285714286	0.508484165438124	0.0032480224632241	1
      I	178400	NA	0.508654123407129	0.00325653910141993	0
      I	178600	0.527443576007728	0.508824797046294	0.00326535807973694	4
      I	178800	0.564192949907236	0.508996187470268	0.00327448059757064	1
      I	179000	0.49054209461871	0.509168295793705	0.00328390778191357	3
      I	179200	0.465686274509804	0.509341123131258	0.00329364068748569	1
      ......
    
    **result/AFD.xls**
    
      ......
      I	178200	0.116071428571429	0.071887997632869	0.00656012459522455	1
      I	178400	NA	0.0717842339388587	0.00657723114717712	0
      I	178600	0.105734609642575	0.0716802521073883	0.00659495069191366	4
      I	178800	0.107977736549165	0.0715760481149835	0.00661328570843137	1
      I	179000	-0.000951738906293353	0.0714716179381699	0.00663223852880404	3
      I	179200	0.0686274509803921	0.0713669575534733	0.00665181133840192	1
      ......
    
 [back to top](#top)
    
* ### <a name="outputstep2"></a>BRM step 2 output: the threshold
    The output file of BRM step 2 (BRMstep2.threshold.R) include eight columns:
    
    | Chromosome Name | Position | Allele Frequency | Fitted Average | Sample threshold | Variance of sample | Theoretical threshold | Variance if this is a QTL |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    
    Among these columns, the most important is "Theoretical threshold", which will be used in the next step. 
    
    ### <a name="outputstep2example"></a>Example
    **result/threshold.xls**
    
      .....
      I	178200	0.504464285714286	0.508484165438124	0.166541321687839	0.00166618679291212	0.166565302509256	0.00165712369486579
      I	178400	NA	0.508654123407129	0.166540351206206	0.00166616737432036	0.166565302509256	0.00165714604341137
      I	178600	0.527443576007728	0.508824797046294	0.166539357262139	0.00166614748638061	0.166565302509256	0.00165716775923202
      I	178800	0.564192949907236	0.508996187470268	0.166538339604352	0.00166612712407333	0.166565302509256	0.00165718881335159
      I	179000	0.49054209461871	0.509168295793705	0.166537297980008	0.00166610628234826	0.166565302509256	0.00165720917764551
      I	179200	0.465686274509804	0.509341123131258	0.166536232134715	0.00166608495612431	0.166565302509256	0.00165722882484733
      ......
      
 [back to top](#top)

* ### <a name="outputstep3"></a>BRM step 3 output: the peaks and the confidence intervals
    The script (BRMstep3.peak_and_CI.R) identify all possible peaks and calculate the confidence interval assuming that is a QTL peak. After getting the output, we can filter out most of the peaks whose value below threshold by using Excel etc. first. Then, we can pick up those reasonable peak lines.
    The peaks file include six columns:
    
    | Chromosome Name | Peak Position | Peak Value | Peak Direction | Confidence Interval Start Position | Confidence Interval End Position |
    | --- | --- | --- | --- | --- | --- |
    
    ### <a name="outputstep3example"></a>Example
    The theoretical threshold in example is 0.166565302509256.
    
    So,
    
    Before:
    
      #Chromosome	Peak Position	Peak Value	Peak Direction	Interval Start	Interval End
      I	33000	0.183543648153291	+	33000	207200
      I	207200	0.0527821549800999	-	33000	207200
      II	10200	-0.0198254858036586	+	10200	8e+05
      II	220070.969156435	-0.0849990247083773	-	10200	8e+05
      II	555244.927316666	0.175754210593449	+	10200	8e+05
      II	8e+05	-0.0322540917987707	-	10200	8e+05
      III	3000	-0.0450448433411402	+	3000	303000
      III	196836.714809522	-0.0984815101521141	-	3000	303000
      III	303000	-0.0772517404558697	+	3000	303000
      IV	19600	0.155465732313207	+	19600	1522000
      IV	102753.718264883	0.0324352807289751	-	19600	1522000
      IV	111814.525991933	0.0362548543519746	+	19600	1522000
      IV	162523.296298926	0.0120050028667934	-	19600	1522000
      ......
      
    After:
    
      #Chromosome	Peak Position	Peak Value	Peak Direction	Interval Start	Interval End
      I	33000	0.183543648153291	+	33000	207200
      II	555244.927316666	0.175754210593449	+	10200	8e+05
      ......
    
 [back to top](#top)
