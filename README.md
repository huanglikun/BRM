# <a name="top"></a>Block Regression Mapping (BRM)
**Block Regression Mapping (BRM)** is a statistical method for QTL mapping based on bulked segregant analysis by deep sequencing. The core function is programmed by R language. For the detailed description of the method, please see the original article "BRM: A statistical method for QTL mapping based on bulked segregant analysis by deep sequencing" published in Bioinformatics.

Please cite: Huang L, Tang W, Bu S, et al. BRM: A statistical method for QTL mapping based on bulked segregant analysis by deep sequencing. Bioinformatics, 2019.  https://doi.org/10.1093/bioinformatics/btz861

## Content
* [Introduction](#intro)
* [Getting started](#getstart)
* [Input data file](#inputdata)
  * [Example](#inputdataexample)
* [Input configuration files](#inputconf)
  * [Chromosome length file](#inputconfchrlen)
    * [Example](#inputconfchrlenexample)
  * [Block regression mapping configuration file](#inputconfblk)
    * [Example](#inputconfblkexample)
* [About <i>u</i><sub>α/2</sub>](#aboutua)
  * [The <i>u</i><sub>α/2</sub> values in various populations](#uafk)
    * [The <i>u</i><sub>α/2</sub> values of some species in various populations](#uafktable)
    * [A tool for calculating <i>u</i><sub>α/2</sub> in an F<sub>k</sub> population](#uafktool)
  * [The <i>u</i><sub>α/2</sub> values in random-mating progeny populations](#uark)
    * [The <i>u</i><sub>α/2</sub> values of yeast and maize in random-mating progeny populations](#uarktable)
    * [A tool for calculating <i>u</i><sub>α/2</sub> in an R<sub>k</sub> population](#uarktool)
* [Output files](#output)
  * [Result 1 file](#out1)
    * [Example](#out1example)
  * [Result 2 file](#out2)
    * [Example](#out2example)
* [Q&A](#qa)
  * [1. If I have the VCF file generated by Freebayes/GATK, how can I convert the VCF format into BSA format?](#vcf2bsa)
  * [2. Why the output shows the data size of one chromosome is 0?](#chrchk)
  * [3. What's the proper block size?](#blksize)

## <a name="intro"></a>Introduction

BRM is a method of BSA-seq for mapping QTLs or major genes. It can apply to different populations including recombinant inbred lines (RIL), doubled haploid (DH), haploid (H), F<sub>2</sub> , F<sub>3</sub>, and so on. 

BRM finds out candidate QTL (or gene) peaks in three main steps. The first step is to divide the genome into many small blocks of equal size and calculate the average **a**llele **f**requency (AF) of each block in each pool, the average **a**llele **f**requency in the **p**opulation (AFP) of each block as well as the **a**llele **f**requency **d**ifference (AFD) between two pools in each block. The second step is to figure out the AFD threshold of the 5% overall significance level at every genomic position. The third step is to identify possible QTL positions (significant AFD peaks) and calculate the 95% confidence interval of each QTL. The BRM scripts output the above results in two files. The first file contains the results of step one and step two, and the second file contains the results of step three.

 [back to top](#top)

## <a name="getstart"></a>Getting started
Download all scripts and examples from GitHub, then you can have a try with the example data:

```bash
git clone https://github.com/huanglikun/BRM.git
cd BRM
# Usage: 
# Rscript BRM.R <Block regression mapping configuration file> <Chromosome length file> <Input data in bsa format>
# Design A
Rscript BRM.R configureExample/designA/BRM_conf.txt configureExample/chr_length.tsv dataExample/designA/yeast_markers_dp10.bsa
# Design B (the experiment which has a high selected pool and a random pool)
Rscript BRM.R configureExample/designBH/BRM_conf.txt configureExample/chr_length.tsv dataExample/designBH/yeast_markers_dp10.bsa
# Design B (the experiment which has a low selected pool and a random pool)
Rscript BRM.R configureExample/designBL/BRM_conf.txt configureExample/chr_length.tsv dataExample/designBL/yeast_markers_dp10.bsa
```
 [back to top](#top)
 
## <a name="inputdata"></a>Input data file
The input data file is a **tab-separated values** file named with bsa format. It contains six columns:

| Column 1 | Column 2 | Column 3 | Column 4 | Column 5 | Column 6 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| Chromosome code | Marker position (bp) | a | b | c | d |

Note: a, b, c and d stand for the counts of marker allele of PARENT 1 in pool 1 (high pool in Design A or selected high/low pool in Design B), PARENT 2 in pool 1, PARENT 1 in pool 2 (low pool in Design A or random pool in Design B) and PARENT 2 in pool 2, respectively.

![Pools illustration](https://github.com/huanglikun/BRM/blob/master/img/pools.png)

  * ### <a name="inputdataexample"></a>Example
    **dataExample/designBL/yeast_markers_dp10.bsa**

        I	33070	6	6	6	7
        I	33147	7	7	2	5
        I	33152	4	5	2	5
        I	33200	3	9	4	3
        I	33293	7	7	7	2  
        ......

[back to top](#top)

## <a name="inputconf"></a>Input configuration files
* ### <a name="inputconfchrlen"></a>Chromosome length file
    It's a **tab-separated values** file with two columns:
    
    | Column 1 | Column 2 |
    | --- | --- |
    | Chromosome code | Chromosome length (bp) |
    
    * ### <a name="inputconfchrlenexample"></a>Example
        **configureExample/chr_length.tsv**
          
            I    230218
            II    813184
            III    316620
            ......

Note: The chromosome code in data file should corresponding to the chromosome code in chromosome length file.

 [back to top](#top)

* ### <a name="inputconfblk"></a>Block regression mapping configuration file
    It's a **key-value** file containing parameters required for block regression mapping (see the following table). The separator is "=". The space between key and value will be ignored.
    
  | Key | Value type | Description | Recommend value |
  | --- | --- | --- | --- |
  | Design | A, BH, BL | A: Design A; <br>BH: Design B with selected high pool; <br>BL: Design B with selected low pool | Default: A |
  | n1 | Integer | Number of individuals in pool1 (high pool) | Depending on experiment design |
  | n2 | Integer | Number of individuals in pool2 (low pool) | Depending on experiment design |
  | t | 0, 1 | Population type | t = 0 for DH, RI or H; t = 1 for F<sub>2</sub>, F<sub>3</sub>, etc. |
  | ua | float | <i>u</i><sub>α/2</sub> value | See the table below, or can be calculated by **cal_ua_fk.R** or **cal_ua_rk.R** in the "tools" directory (see explanation below). |
  | UNIT | Integer | Unit block size (bp)  | Default: 1000 |
  | DEG | Integer | The degree of the polynomials used in LOESS (local polynomial regression fitting) | 2 |
  | BLK | Integer or float | Block size = BLK * UNIT | e.g.: yeast, BLK = 0.2 ; rice, BLK = 20 |
  | MIN | Integer | Minimum total depth in a valid block | 10 |
  | MINVALID | Integer | Minimum valid block number in a chromosome | 10 |
  | Result1_File | File path | To define the result 1 output path (optional) | Default: result/result1.xls |
  | Result2_File | File path | To define the result 2 output path (optional) | Default: result/result2.xls |
    
  * ### <a name="inputconfblkexample"></a>Example
      **configureExample/designBL/BRM_conf.txt**

          # version 0.3

          # Experiment
          Design   = BL       # available: A, BH, BL. BH: Design B with HIGH selected pool. BL: Design B with LOW selected pool.
          n1       = 300      # pool 1 size. The high pool for design A or the selected pool for design B.
          n2       = 300      # pool 2 size. The low pool for design A or the random pool for design B.
          t        = 0        # For DH or RI etc., t=0; F2 or F3 etc., t=1
          ua       = 4.08     # For rice, F2:3.65; F3:3.74; F4:3.78.

          # Block regression
          UNIT     = 1000     # block unit(bp)
          DEG      = 2        # The degree of the polynomials to be used in Local Polynomial Regression Fitting.
          BLK      = 0.2      # block size = BLK * UNIT
          MIN      = 10       # min total depth in block
          MINVALID = 10       # min valid blocks in one chromosome (needed to be at least 10)

          # output file determination (optional)
          # Result1_File = result/result1.xls	# including allele frequency and threshold
          # Result2_File = result/result2.xls	# including peak information and confidence interval
          Result1_File = result_random_L/result1.xls	# including allele frequency and threshold
          Result2_File = result_random_L/result2.xls	# including peak information and confidence interval

 [back to top](#top)

## <a name="aboutua"></a>About <i>u</i><sub>α/2</sub>
### <a name="uafk"></a>The <i>u</i><sub>α/2</sub> values in various populations

  * ### <a name="uafktable"></a>The <i>u</i><sub>α/2</sub> values of some species in various populations
   
      <table style="text-align: center"><thead><tr><th rowspan=2 style="text-align: center">Species</th><th rowspan=2 style="text-align: center">n<sup>a</sup></th><th colspan=2 style="text-align: center">Genome size<sup>b</sup></th><th rowspan=2 style="text-align: center">Ratio</br>(kb/cM)</th><th colspan=4 style="text-align: center"><i>u</i><sub>α/2</sub><sup>c</sup></th><th rowspan=2 style="text-align: center">Ref.</th></tr>
    <tr><th style="text-align: center">cM</th><th style="text-align: center">Mb</th><th style="text-align: center">H/DH/<i>F</i><sub>2</sub></th><th style="text-align: center"><i>F</i><sub>3</sub></th><th style="text-align: center"><i>F</i><sub>4</sub></th><th style="text-align: center">RIL</th></tr></thead>
    <tbody><tr><td><i>Arabidopsis</i></td><td>5</td><td>600</td><td>119</td><td>199</td><td>3.41</td><td>3.50</td><td>3.54</td><td>3.57</td><td>[1]</td></tr>
    <tr><td>Cucumber</td><td>7</td><td>1390</td><td>192</td><td>138</td><td>3.62</td><td>3.72</td><td>3.75</td><td>3.78</td><td>[2]</td></tr>
    <tr><td>Maize</td><td>10</td><td>2060</td><td>2106</td><td>1023</td><td>3.72</td><td>3.81</td><td>3.85</td><td>3.87</td><td>[3]</td></tr>
    <tr><td>Rapeseed</td><td>18</td><td>2520</td><td>855</td><td>339</td><td>3.78</td><td>3.87</td><td>3.90</td><td>3.93</td><td>[4]</td></tr>
    <tr><td>Rice</td><td>12</td><td>1530</td><td>382</td><td>250</td><td>3.65</td><td>3.74</td><td>3.78</td><td>3.80</td><td>[5]</td></tr>
    <tr><td>Tobacco</td><td>12</td><td>3270</td><td>3613</td><td>1105</td><td>3.84</td><td>3.92</td><td>3.96</td><td>3.98</td><td>[6]</td></tr>
    <tr><td>Tomato</td><td>12</td><td>1470</td><td>807</td><td>549</td><td>3.65</td><td>3.73</td><td>3.77</td><td>3.80</td><td>[7]</td></tr>
    <tr><td>Wheat</td><td>21</td><td>3140</td><td>14547</td><td>4633</td><td>3.83</td><td>3.92</td><td>3.95</td><td>3.98</td><td>[8]</td></tr>
    <tr><td>Yeast</td><td>16</td><td>4900</td><td>12</td><td>2.5</td><td>3.93</td><td></td><td></td><td></td><td>[9]</td></tr></tbody></table>

    Note: **a.** n, number of chromosomes in haploid. **b.** The genetic map length (cM) of each species was from the references listed except for that of yeast, which was calculated by us using the data from [9]. The genome size (Mb) of each species was all from NCBI. **c.** Corresponding to the overall significance level of 0.05.

    **References**

    [1] Garcia-Hernandez,M. *et al.* (2002) TAIR: a resource for integrated Arabidopsis data. *Funct Integr Genomics*, **2**, 239–253.

    [2] Zhou,Q. *et al.* (2015) A sequencing-based linkage map of cucumber. *Molecular Plant*, **8**, 961–963.
    
    [3] Civardi,L. *et al.* (1994) The relationship between genetic and physical distances in the cloned a1-sh2 interval of the Zea mays L. genome. *Proc Natl Acad Sci U S A*, **91**, 8268–8272.

    [4] Raman,H. *et al.* (2014) SNP markers-based map construction and genome-wide linkage analysis in Brassica napus. *Plant Biotechnology Journal*, **12**, 851–860.

    [5] International Rice Genome Sequencing Project (2005) The map-based sequence of the rice genome. *Nature*, **436**, 793.

    [6] Bindler,G. *et al.* (2011) A high density genetic map of tobacco (*Nicotiana tabacum* L.) obtained from large scale microsatellite marker development. *Theor Appl Genet*, **123**, 219.

    [7] Shirasawa,K. *et al.* (2010) SNP Discovery and linkage map construction in cultivated tomato. *DNA Research*, **17**, 381–391.

    [8] Yang,Q. *et al.* (2018) High-density genetic map construction and mapping of the homologous transformation sterility gene (hts) in wheat using GBS markers. *BMC Plant Biology*, **18**, 301.

    [9] Bloom,J.S. *et al.* (2015) Genetic interactions contribute less than additive effects to quantitative trait variation in yeast. *Nat Commun*, **6**, 8712.

  * ### <a name="uafktool"></a>A tool for calculating <i>u</i><sub>α/2</sub> in an F<sub>k</sub> population

      A tool is provided in **tools/cal_ua_fk.R** as below for calculating <i>u</i><sub>α/2</sub> in an F<sub>k</sub> population given the values of necessary parameters.
    
      * #### Quick start

        ```bash
        # Usage
        # Rscript tools/cal_ua_fk.R <population and species information file>
        Rscript tools/cal_ua_fk.R configureExample/fk_rice/ua_fk_conf.txt
        ```

      * **Example**

        **configureExample/fk_rice/ua_fk_conf.txt**

        If we know the approximately total genetic distance of this species,

            #
            precision = 1e-13
            k         = 2            # k=0(RIL) k=1(H/DH) k=2(F2) k=3(F3) k=4(F4) k>4(F5,F6...)
            #
            chrNum         = 12      # chromosome number
            geneticLength  = 1530    # Genome size (cM)
            #genomeSize    = 382     # Genome size (Mb)
            #ratio         = 249.7   # Ratio (kb/cM)

        Otherwise,

            #
            precision = 1e-13
            k         = 2            # k=0(RIL) k=1(H/DH) k=2(F2) k=3(F3) k=4(F4) k>4(F5,F6...)
            #
            chrNum         = 12      # chromosome number
            #geneticLength = 1530    # Genome size (cM)
            genomeSize     = 382     # Genome size (Mb)
            ratio          = 249.7   # Ratio (kb/cM)

        It is a **key-value** file contains required parameters. The separator is "=". The space between key and value will be ignored. There are at least four parameters needed to be set:

        | Key | Value type | Description | Recommend value | 
        | --- | --- | --- | --- |
        | precision | Float | Precision for numerical calculation. | Default: 1e-13 |
        | k | Integer | Generation. <br>k=0 means RIL<br>k=1 means H/DH<br>k=2 means F<sub>2</sub><br>k=3 means F<sub>3</sub><br>k=4 means F<sub>4</sub><br>k>4 means F<sub>5</sub>, F<sub>6</sub> ..., and is considered as RIL. | Depending on experimental design |
        | chrNum | Integer | Chromosome number. | Depending on species |
        | geneticLength | Integer | Approximately total genetic distance of the species, the unit is cM. Comment this by "#" to use parameters below instead. | Depending on species |
        | genomeSize | Float | Genome size of the species, the unit is Mb. (optional) | Depending on species |
        | ratio | Float | Physical distance/Genetic distance ratio (kb/cM). (optional) | Depending on species |

      * #### Result

        The <i>u</i><sub>α/2</sub> will be displayed as ua on the screen:

            Read global parameters from configureExample/fk_rice/ua_fk_conf.txt .
            Parameter precision=1e-13
            Parameter k=2
            Parameter chrNum=12
            Parameter geneticLength=1530

            ua is  3.654641  

[back to top](#top)

### <a name="uark"></a>The <i>u</i><sub>α/2</sub> values in random-mating progeny populations

  For the case that the progeny population of a cross between two pure-line parents is not generated by selfing but by random mating, the <i>u</i><sub>α/2</sub> value will be larger. The following table shows the <i>u</i><sub>α/2</sub> values of yeast and maize in different random-mating progeny populations, where H<sub>1</sub> (or H) is the gamete (haploid) generated by F<sub>1</sub>, and R<sub>1</sub> (or F<sub>2</sub>) is the sporophyte generated by combination of F<sub>1</sub> gametes; H<sub>2</sub> is the gamete generated by R<sub>1</sub> (F<sub>2</sub>), and R<sub>2</sub> is the sporophyte generated by combination of R<sub>1</sub> (F<sub>2</sub>) gametes; the others can be deduced likewise. 

  * ### <a name="uarktable"></a>The <i>u</i><sub>α/2</sub> values of yeast and maize in random-mating progeny populations

      | Species | H<sub>1</sub>(H)/R<sub>1</sub>(F<sub>2</sub>) | H<sub>2</sub>/R<sub>2</sub> | H<sub>3</sub>/R<sub>3</sub> | H<sub>4</sub>/R<sub>4</sub> | H<sub>5</sub>/R<sub>5</sub> | H<sub>6</sub>/R<sub>6</sub> | H<sub>7</sub>/R<sub>7</sub> | H<sub>8</sub>/R<sub>8</sub> | H<sub>9</sub>/R<sub>9</sub> |
      | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
      | Yeast | 3.93 | 4.02 | 4.08 | 4.13 | 4.17 | 4.21 | 4.24 | 4.26 | 4.29 |
      | Maize | 3.72 | 3.81 | 3.88 | 3.93 | 3.97 | 4.01 | 4.04 | 4.07 | 4.09 |

  * ### <a name="uarktool"></a>A tool for calculating <i>u</i><sub>α/2</sub> in an R<sub>k</sub> population

      A tool is provided in **tools/cal_ua_rk.R** as below for calculating <i>u</i><sub>α/2</sub> in a random-mating progeny population given the values of necessary parameters.

      * #### Quick start

        ```bash
        # Usage
        # Rscript tools/cal_ua_rk.R <population and species information file>
        Rscript tools/cal_ua_rk.R configureExample/rk_yeast/ua_rk_conf.txt
        ```

      * #### Example

        **configureExample/rk_yeast/ua_rk_conf.txt**

        If we know the approximately total genetic distance of this species,

            #
            precision = 1e-13
            k         = 6
            #
            chrNum         = 16      # chromosome number
            geneticLength  = 4900    # Genome size (cM)
            #genomeSize    = 12      # Genome size (Mb)
            #ratio         = 2.5     # Ratio (kb/cM)

        Otherwise,

            #
            precision = 1e-13
            k         = 6
            #
            chrNum         = 16      # chromosome number
            #geneticLength = 4900    # Genome size (cM)
            genomeSize     = 12      # Genome size (Mb)
            ratio          = 2.5     # Ratio (kb/cM)

        It is a **key-value** file contains required parameters. The separator is "=". The space between key and value will be ignored. There are at least four parameters needed to be set:

        | Key | Value type | Description | Recommend value | 
        | --- | --- | --- | --- |
        | precision | Float | Precision for numerical calculation. | Default: 1e-13 |
        | k | Integer | Generation, e.g. k=6 means R<sub>6</sub> for diploid or H<sub>6</sub> for haploid. | Depending on experimental design |
        | chrNum | Integer | Chromosome number. | Depending on species |
        | geneticLength | Integer | Approximately total genetic distance of the species, the unit is cM. Comment this by "#" to use parameters below instead. | Depending on species |
        | genomeSize | Float | Genome size of the species, the unit is Mb. (optional) | Depending on species |
        | ratio | Float | Physical distance/Genetic distance ratio (kb/cM). (optional) | Depending on species |

      * #### Result

        The <i>u</i><sub>α/2</sub> will be displayed as ua on the screen:

            Read global parameters from configureExample/rk_yeast/ua_rk_conf.txt .
            Parameter precision=1e-13
            Parameter k=6
            Parameter chrNum=16
            Parameter geneticLength=4900

            recombination rate is  0.02545
            Genetic distance is  2.5472  cM
            ua is  4.207883  


[back to top](#top)

## <a name="output"></a>Output files
### <a name="out1"></a>Result 1 file

  This file contains one line to show the theoretical threshold (assuming AF = 0.5 in the population) and a 11 columns table following. The AF is defined by the allele from Parent 1.

  The data table described as below:

  | Column | Heading | Description |
  | :---: | --- | --- |
  | 1 | Chr. | Chromosome |
  | 2 | Pos. | Block position (bp) |
  | 3 | AF1-Observed | Observed value of block average **a**llele **f**requency in pool 1 (AF1) |
  | 4 | AF1-Expected | Expected value of block average **a**llele **f**requency in pool 1 (AF1) |
  | 5 | AF2-Observed | Observed value of block average **a**llele **f**requency in pool 2 (AF2) | 
  | 6 | AF2-Expected | Expected value of block average **a**llele **f**requency in pool 2 (AF2) |
  | 7 | AFD-Observed | Observed value of **a**llele **f**requency **d**ifferency (AF1 - AF2) |
  | 8 | AFD-Expected | Expected value of **a**llele **f**requency **d**ifferency (AF1 - AF2) |
  | 9 | AFP-Observed | Observed value of block average **a**llele **f**requency in the **p**opulation (AFP) |
  | 10 | AFP-Expected |Expected value of Observed value of block average **a**llele **f**requency in the **p**opulation (AFP) |
  | 11 | Sample threshold | Threshold estimated based on the expected AFP |

  * ### <a name="out1example"></a>Example

    **result_random_L/result1.xls**

        #Theoretical threshold is ±	0.1666
        #Chr.	Pos.	AF1-Observed	AF1-Expected	AF2-Observed	AF2-Expected	AFD-Observed	AFD-Expected	AFP-Observed	AFP-Expected	Sample threshold
        ......
        I	178300	NA	0.4655	0.6667	0.4993	NA	0.0409	0.6667	0.4993	0.1666
        I	178500	NA	0.4658	NA	0.4993	NA	0.0404	NA	0.4993	0.1666
        I	178700	0.5	0.4661	0.4583	0.4993	-0.0417	0.04	0.4583	0.4993	0.1666
        I	178900	0.4545	0.4664	NA	0.4992	NA	0.0396	NA	0.4992	0.1666
        I	179100	0.4444	0.4667	0.5	0.4992	0.0556	0.0392	0.5	0.4992	0.1666
        I	179300	0.6	0.467	0.2727	0.4992	-0.3273	0.0387	0.2727	0.4992	0.1666
        ......

 [back to top](#top)
    
### <a name="out2"></a>Result 2 file
  This file shows the information of candidate QTLs (significant AFD peaks). It contains six columns. By far, some peaks are needed to filter out manually to get the final candidate QTL list.

  | Column | Heading | Description |
  | --- | --- | --- |
  | 1 | Chr. | Chromosome code | 
  | 2 | Pos. | Position of block center | 
  | 3 | Val. | Peak value of AFD |
  | 4 | Peak Dir. | Peak direction: +, upward; -,downward |
  | 5 | Start | Start point of confidence interval |
  | 6 | End | End point of confidence interval|
  
  * ### <a name="out2example"></a>Example
    **result_random_L/result2.xls**
  
        #Chr.	Pos.	Val.	Peak Dir.	Start	End
        IV	599459	-0.1914	-	534520	680013
        VIII	222625	-0.2848	-	153041	287275
        XI	643900	-0.3316	-	618381	643900
        XII	655812	0.2102	+	550089	733741
        XIII	9900	-0.2488	-	9900	108758
        XIV	377027	0.273	+	321716	450359
      
 [back to top](#top)

## <a name="qa"></a>Q&A
1. <a name="vcf2bsa"></a>**If I have the VCF file generated by Freebayes/GATK, how can I convert the VCF format into BSA format?**
  
    A perl script is provided for transforming the VCF format into BSA format.

      #### Quick start

      ```bash
      # Usage:
      # perl tools/vcf2bsa.pl <samples information file> <markers.vcf> <output file>
      perl tools/vcf2bsa.pl configureExample/vcf2bsa/vcf2bsa_conf.txt dataExample/vcf2bsa/markers.freebayes.vcf result/markers.bsa
      ```
      
      #### Example
      **configureExample/vcf2bsa/vcf2bsa_conf.txt**
        
        # Samples in 2x2 Table
        # Bulk/Pool 1: the first sample
        Table2x2.pool1 = low-pool-RG
        # Bulk/Pool 2: the second sample
        Table2x2.pool2 = random-pool-RG
        # Parent 1
        Table2x2.parent1 = P1-RG
        # Parent 2
        Table2x2.parent2 = P2-RG

      It's a **key-value** file contains samples relationship. The separator is "=". And the space between key and value will be ignored. There are four parameters needed to be set:

      | Key | Value type | Description |
      | --- | --- | --- |
      | Table2x2.pool1 | String | Design A: high selected pool RG tag in VCF file. <br>Design B: selected pool RG tag in VCF file. |
      | Table2x2.pool2 | String | Design A: low selected pool RG tag in VCF file. <br>Design B: random pool RG tag in VCF file. |
      | Table2x2.parent1 | String | High phenotype value parent sample RG tag in VCF file. |
      | Table2x2.parent2 | String | Low phenotype value parent sample RG tag in VCF file. |

      RG tag is the “Read Group” setting at the reads mapping step. For example, the -R parameter value in BWA.

      If one parent is missing, the script will consider the genotype different from the known parent sample as the other parent genotype. If both parents are missing, the script will consider the reference genotype as the known parent genotype.


 [back to top](#top)

2. <a name="chrchk"></a>**Why the output shows the data size of one chromosome is 0?**

    If the chromosome code in the data file and the chromosome code in the chromosome length file didn't match, such report will be shown. The chromosome code match is **case sensitive**.

 [back to top](#top)

3. <a name="blksize"></a>**What's the proper block size?**

    We recommend about 0.1 cM as the block size. For example, as the yeast is about 2.5 kb/cM, we choose 0.2 kb as the block size. While the rice is around 250 kb/cM, we choose 20 kb as the block size.

 [back to top](#top)
