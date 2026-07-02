# phylogeoflow: Citations

## Framework

- **nf-core**
  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020;38(3):276-278. doi: 10.1038/s41587-020-0439-x.

- **Nextflow**
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319. doi: 10.1038/nbt.3820.

## Part 1 — Molecular data

### Data sources

- **BOLD (Barcode of Life Data System)**
  > Ratnasingham S, Hebert PDN. BOLD: The Barcode of Life Data System. Mol Ecol Notes. 2007;7(3):355-364. doi: 10.1111/j.1471-8286.2007.01678.x.
- **GenBank**
  > Benson DA, Cavanaugh M, Clark K, et al. GenBank. Nucleic Acids Res. 2013;41(Database issue):D36-42. doi: 10.1093/nar/gks1195.
- **GBIF**
  > GBIF.org. GBIF Occurrence Download. (Cite the dataset DOI written to `gbif_doi.txt` for each run.)

### Retrieval & cleaning

- **BOLDconnectR** — BOLD Systems Central. https://github.com/boldsystems-central/BOLDconnectR
- **rentrez**
  > Winter DJ. rentrez: an R package for the NCBI eUtils API. R Journal. 2017;9(2):520-526.
- **Entrez Direct (EDirect)**
  > Kans J. Entrez Direct: E-utilities on the Unix Command Line. NCBI. https://www.ncbi.nlm.nih.gov/books/NBK179288/
- **rgbif**
  > Chamberlain S, Boettiger C. R Python, and Ruby clients for GBIF species occurrence data. PeerJ Preprints. 2017. doi: 10.7287/peerj.preprints.3304v1.
- **CoordinateCleaner**
  > Zizka A, Silvestro D, Andermann T, et al. CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases. Methods Ecol Evol. 2019;10(5):744-751. doi: 10.1111/2041-210X.13152.

### Taxonomic classification

- **RDPClassifier**
  > Wang Q, Garrity GM, Tiedje JM, Cole JR. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 2007;73(16):5261-5267. doi: 10.1128/AEM.00062-07.
- **COI eukaryote training set**
  > Porter TM, Gibson JF, Shokralla S, et al. Rapid and accurate taxonomic classification of insect (class Insecta) COI DNA barcode sequences using a naive Bayesian classifier. Mol Ecol Resour. 2014;14(5):929-942. doi: 10.1111/1755-0998.12240.
  > Porter TM, Hajibabaei M. Automated high throughput animal CO1 metabarcode classification. Sci Rep. 2018;8:4226. doi: 10.1038/s41598-018-22505-4.

### Multiple sequence alignment & evaluation

- **MAFFT**
  > Katoh K, Standley DM. MAFFT Multiple Sequence Alignment Software Version 7. Mol Biol Evol. 2013;30(4):772-780. doi: 10.1093/molbev/mst010.
- **PASTA**
  > Mirarab S, Nguyen N, Guo S, et al. PASTA: Ultra-Large Multiple Sequence Alignment. J Comput Biol. 2015;22(5):377-386. doi: 10.1089/cmb.2014.0156.
- **MUSCLE**
  > Edgar RC. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 2004;32(5):1792-1797. doi: 10.1093/nar/gkh340.
- **T-Coffee (CORE / TCS evaluation)**
  > Notredame C, Higgins DG, Heringa J. T-Coffee: A novel method for fast and accurate multiple sequence alignment. J Mol Biol. 2000;302(1):205-217. doi: 10.1006/jmbi.2000.4042.
  > Chang JM, Di Tommaso P, Notredame C. TCS: a new multiple sequence alignment reliability measure. Mol Biol Evol. 2014;31(6):1625-1637. doi: 10.1093/molbev/msu117.
- **SeaView (visualisation)**
  > Gouy M, Guindon S, Gascuel O. SeaView Version 4. Mol Biol Evol. 2010;27(2):221-224. doi: 10.1093/molbev/msp259.
- **BMGE (trimming)**
  > Criscuolo A, Gribaldo S. BMGE (Block Mapping and Gathering with Entropy). BMC Evol Biol. 2010;10:210. doi: 10.1186/1471-2148-10-210.

### Phylogenetic inference & placement

- **ModelFinder**
  > Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. ModelFinder. Nat Methods. 2017;14(6):587-589. doi: 10.1038/nmeth.4285.
- **ModelTest-NG**
  > Darriba D, Posada D, Kozlov AM, et al. ModelTest-NG. Mol Biol Evol. 2020;37(1):291-294. doi: 10.1093/molbev/msz189.
- **IQ-TREE**
  > Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ. IQ-TREE. Mol Biol Evol. 2015;32(1):268-274. doi: 10.1093/molbev/msu300.
- **UFBoot2**
  > Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. UFBoot2. Mol Biol Evol. 2018;35(2):518-522. doi: 10.1093/molbev/msx281.
- **RAxML-NG**
  > Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. RAxML-NG. Bioinformatics. 2019;35(21):4453-4455. doi: 10.1093/bioinformatics/btz305.
- **RAxML**
  > Stamatakis A. RAxML version 8. Bioinformatics. 2014;30(9):1312-1313. doi: 10.1093/bioinformatics/btu033.
- **FastTree**
  > Price MN, Dehal PS, Arkin AP. FastTree 2. PLoS ONE. 2010;5(3):e9490. doi: 10.1371/journal.pone.0009490.
- **EPA-NG (phylogenetic placement)**
  > Barbera P, Kozlov AM, Czech L, et al. EPA-ng. Syst Biol. 2019;68(2):365-369. doi: 10.1093/sysbio/syy054.
- **Bio.Phylo**
  > Talevich E, Invergo BM, Cock PJA, Chapman BA. Bio.Phylo. BMC Bioinformatics. 2012;13:209. doi: 10.1186/1471-2105-13-209.

### Species delimitation

- **mPTP**
  > Kapli P, Lutteropp S, Zhang J, et al. Multi-rate Poisson tree processes for single-locus species delimitation. Bioinformatics. 2017;33(11):1630-1638. doi: 10.1093/bioinformatics/btx025.

### Population structure & phylogeography

- **pegas**
  > Paradis E. pegas: an R package for population genetics with an integrated-modular approach. Bioinformatics. 2010;26(3):419-420. doi: 10.1093/bioinformatics/btp696.
- **ape**
  > Paradis E, Schliep K. ape 5.0. Bioinformatics. 2019;35(3):526-528. doi: 10.1093/bioinformatics/bty633.
- **adegenet**
  > Jombart T. adegenet. Bioinformatics. 2008;24(11):1403-1405. doi: 10.1093/bioinformatics/btn129.
- **diveRsity**
  > Keenan K, McGinnity P, Cross TF, Crozier WW, Prodohl PA. diveRsity. Methods Ecol Evol. 2013;4(8):782-788. doi: 10.1111/2041-210X.12067.
- **FinePop**
  > Kitada S, Nakamichi R, Kishino H. The empirical Bayes estimators of fine-scale population structure in high gene flow species. Mol Ecol Resour. 2017;17(6):1210-1222. doi: 10.1111/1755-0998.12663.
- **Jost's D / FST**
  > Jost L. GST and its relatives do not measure differentiation. Mol Ecol. 2008;17(18):4015-4026. doi: 10.1111/j.1365-294X.2008.03887.x.
  > Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. Evolution. 1984;38(6):1358-1370. doi: 10.2307/2408641.
- **AMOVA**
  > Excoffier L, Smouse PE, Quattro JM. Analysis of molecular variance. Genetics. 1992;131(2):479-491.
- **GENELAND**
  > Guillot G, Mortier F, Estoup A. GENELAND. Mol Ecol Notes. 2005;5(3):712-715. doi: 10.1111/j.1471-8286.2005.01031.x.
- **SAMOVA**
  > Dupanloup I, Schneider S, Excoffier L. A simulated annealing approach to define the genetic structure of populations. Mol Ecol. 2002;11(12):2571-2581. doi: 10.1046/j.1365-294X.2002.01650.x.
- **Neutrality tests**
  > Tajima F. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics. 1989;123(3):585-595.
  > Fu YX. Statistical tests of neutrality of mutations against population growth, hitchhiking and background selection. Genetics. 1997;147(2):915-925.
- **Mantel tests (ade4 / vegan)**
  > Dray S, Dufour AB. The ade4 package. J Stat Softw. 2007;22(4):1-20. doi: 10.18637/jss.v022.i04.
  > Legendre P, Fortin MJ, Borcard D. Should the Mantel test be used in spatial analysis? Methods Ecol Evol. 2015;6(11):1239-1247. doi: 10.1111/2041-210X.12425.
- **BASTA (structured coalescent) / BEAST2**
  > De Maio N, Wu CH, O'Reilly KM, Wilson D. New routes to phylogeography: a Bayesian structured coalescent approximation. PLoS Genet. 2015;11(8):e1005421. doi: 10.1371/journal.pgen.1005421.
  > Bouckaert R, Vaughan TG, Barido-Sottani J, et al. BEAST 2.5. PLoS Comput Biol. 2019;15(4):e1006650. doi: 10.1371/journal.pcbi.1006650.
- **Geocoding helper**
  > Gratton P, Marta S, Bocksberger G, et al. A world of sequences: can we use georeferenced nucleotide databases for a robust automated phylogeography? J Biogeogr. 2017;44(2):475-486. doi: 10.1111/jbi.12786.
- **PGDSpider (format conversion)**
  > Lischer HEL, Excoffier L. PGDSpider. Bioinformatics. 2012;28(2):298-299. doi: 10.1093/bioinformatics/btr642.

## Part 2 — Environmental integration

### Environmental data & access

- **WorldClim 2.1**
  > Fick SE, Hijmans RJ. WorldClim 2. Int J Climatol. 2017;37(12):4302-4315. doi: 10.1002/joc.5086.
- **CHELSA**
  > Karger DN, Conrad O, Bohner J, et al. Climatologies at high resolution for the earth's land surface areas. Sci Data. 2017;4:170122. doi: 10.1038/sdata.2017.122.
- **geodata / terra** — Hijmans RJ, et al. geodata & terra R packages. https://cran.r-project.org/package=geodata
- **MODIS NDVI (MODISTools)**
  > Tuck SL, Phillips HRP, Hintzen RE, et al. MODISTools. Ecol Evol. 2014;4(24):4658-4668. doi: 10.1002/ece3.1273.
- **NASA earthaccess / earthdatalogin** — https://github.com/nsidc/earthaccess ; https://boettiger-lab.github.io/earthdatalogin/
- **TerraClimate**
  > Abatzoglou JT, Dobrowski SZ, Parks SA, Hegewisch KC. TerraClimate. Sci Data. 2018;5:170191. doi: 10.1038/sdata.2017.191.
- **JRC Global Surface Water**
  > Pekel JF, Cottam A, Gorelick N, Belward AS. High-resolution mapping of global surface water and its long-term changes. Nature. 2016;540:418-422. doi: 10.1038/nature20584.
- **WorldPop**
  > Tatem AJ. WorldPop, open data for spatial demography. Sci Data. 2017;4:170004. doi: 10.1038/sdata.2017.4.

### Niche / distribution modelling

- **MaxEnt**
  > Phillips SJ, Anderson RP, Schapire RE. Maximum entropy modeling of species geographic distributions. Ecol Modell. 2006;190(3-4):231-259. doi: 10.1016/j.ecolmodel.2005.03.026.
- **ENMeval**
  > Kass JM, Muscarella R, Galante PJ, et al. ENMeval 2.0. Methods Ecol Evol. 2021;12(9):1602-1608. doi: 10.1111/2041-210X.13628.

### Landscape genetics

- **MMRR**
  > Wang IL. Examining the full effects of landscape heterogeneity on spatial genetic variation: a multiple matrix regression approach. Evolution. 2013;67(12):3403-3411. doi: 10.1111/evo.12134.
- **GDM**
  > Ferrier S, Manion G, Elith J, Richardson K. Using generalized dissimilarity modelling. Divers Distrib. 2007;13(3):252-264. doi: 10.1111/j.1472-4642.2007.00341.x.
  > Fitzpatrick MC, Keller SR. Ecological genomics meets community-level modelling of biodiversity. Ecol Lett. 2015;18(1):1-16. doi: 10.1111/ele.12376.
- **RDA (genotype-environment association)**
  > Forester BR, Lasky JR, Wagner HH, Urban DL. Comparing methods for detecting multilocus adaptation with multivariate genotype-environment associations. Mol Ecol. 2018;27(9):2215-2233. doi: 10.1111/mec.14584.
- **vegan** — Oksanen J, et al. vegan: Community Ecology Package. https://cran.r-project.org/package=vegan

## Software packaging & reproducibility

- **Anaconda** — https://anaconda.com
- **Bioconda**
  > Gruning B, Dale R, Sjodin A, et al. Bioconda. Nat Methods. 2018;15(7):475-476. doi: 10.1038/s41592-018-0046-7.
- **BioContainers**
  > da Veiga Leprevost F, et al. BioContainers. Bioinformatics. 2017;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192.
- **Docker** — Merkel D. Docker. Linux Journal. 2014;2014(239):2.
- **Singularity**
  > Kurtzer GM, Sochat V, Bauer MW. Singularity. PLoS ONE. 2017;12(5):e0177459. doi: 10.1371/journal.pone.0177459.

## Protocols & references from the original meta-analysis

- **COI barcode preparation**
  > Prosser SWJ, deWaard JR, Miller SE, Hebert PDN. DNA barcodes from century-old type specimens. Mol Ecol Resour. 2016;16(2):487-497. doi: 10.1111/1755-0998.12474.
  > Wilson JJ. DNA barcodes for insects. Methods Mol Biol. 2012;858:17-46. doi: 10.1007/978-1-61779-591-6_3.
