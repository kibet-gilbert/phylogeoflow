# **Three Key Databases for Phylogeography: BOLD, GenBank, and GBIF**  

For phylogeography—the study of evolutionary history and biogeographic patterns within species across geographic space—these three databases serve complementary but distinct roles:   

---

## **1. BOLD (Barcode of Life Data System)**  
**The Standardized Genetic Reference**  

BOLD is an online workbench and database that supports the assembly and use of DNA barcode data, functioning as a searchable database of Barcode Index Numbers (BINs)—sequence clusters that closely approximate species.   

**For phylogeography:**  
- BOLD focuses on the mitochondrial gene cytochrome c oxidase I (COI), which serves as the core of a global bioidentification system for animals  
- BOLD requires standardized metadata including species name, voucher data (catalogue number and institution storing), specimen identifier, COI sequences of at least 500 bp, and trace files  
- BOLD data are curated with associated metadata about taxonomy, collection localities, and geographic coordinates, which is essential for phylogeographic analyses  
- **Strength**: High-quality, standardized sequences linked to physical specimens and precise location data  
- **Limitation**: Primarily animals and single-gene (COI) focus  

---

## **2. GenBank (NCBI)**  
**The Comprehensive Sequence Repository**  

GenBank is an annotated collection of all publicly available nucleotide sequences and their protein translations, produced by the National Center for Biotechnology Information as part of an international collaboration with the European Molecular Biology Laboratory and the DNA Data Bank of Japan.  

**For phylogeography:**  
- GenBank contains publicly available nucleotide sequences for 420,000 formally described species, with most submissions made by individual researchers and sequencing centers  
- Supports multiple molecular markers: mitochondrial genes (COI, cytochrome *b*, 16S rRNA), nuclear genes (ITS), and larger genomic regions—ideal for multilocus phylogeographic studies  
- GenBank serves as a central repository enabling comparability of molecular data across research studies, despite variations in data size and complexity  
- **Strength**: Massive scale, diverse markers, integrates with literature via PubMed  
- **Limitation**: Limited georeferencing, uneven geographic representation, and data compatibility issues constrain its utility for phylogeographic studies  

---

## **3. GBIF (Global Biodiversity Information Facility)**  
**The Species Occurrence and Distribution Layer**  

GBIF is an international network and data infrastructure providing global data that document the occurrence of species, currently integrating datasets documenting over 1.6 billion species occurrences from specimen-related data of natural history museums, observations from citizen science networks, and automated environmental surveys.   

**For phylogeography:**  
- Provides georeferenced occurrence data crucial for mapping phylogeographic patterns and understanding species distributions across space and time   
- GBIF provides access to about 15 million iBOL occurrence records (98% of BOLD's public database), plus 500,000 occurrence records from DNA metabarcoding studies   
- GBIF accumulates several hundred million records and serves as the basis for large-scale analyses of macroecological and biogeographic patterns and to document environmental changes over time  
- **Strength**: Massive geographic coverage linking genetic and occurrence data at global scale  
- **Limitation**: Data suffer from spatial and taxonomic biases, data quality issues, and errors that require careful cleaning before use in quantitative analyses   

---

## **Practical Phylogeography Workflow**   

BOLD and GBIF work together by maintaining mapping layers linking DNA sequences with Linnaean taxonomy through Barcode Index Numbers (BINs) and occurrence data taxonomically annotated at the operational taxonomic unit (OTU) level, allowing indexing of species occurrence data.   

**Integrated approach:**   
1. **BOLD** → retrieve validated COI sequences with precise specimen location data   
2. **GenBank** → supplement with additional molecular markers for deeper phylogenetic resolution   
3. **GBIF** → contextualize genetic data with broader species distribution patterns and environmental variables   

Combined use of these databases enables researchers to conduct phylogeographic, macrogenetic, and conservation analyses by integrating mitochondrial sequences with associated metadata about taxonomy, collection localities, and geographic coordinates.

---

## **Other Specialized Databases for Phylogeographic Work**

### **UNITE - Fungal & Eukaryotic ITS Sequences**

UNITE (https://unite.ut.ee) is a web-based database and sequence management environment for molecular identification of eukaryotes, targeting the nuclear ribosomal internal transcribed spacer (ITS) region and offering nearly 10 million such sequences for reference, clustered into ~2.4 million species hypotheses (SHs).

- Especially valuable for fungi, which are under-represented in BOLD
- UNITE integrates with the taxonomic backbone of GBIF and regularly exchanges data with major fungal sequence databases

### **Dryad - Data Repository for Phylogenetic Studies**

Dryad is a major data repository that integrates large-scale genetic datasets in genetics with other large, interoperable databases such as GBIF, WorldClim, the International Nucleotide Sequence Database Collaboration (INSDC), and BOLD, facilitating the emerging field of macrogenetics.

- Stores phylogeographic datasets, genetic alignments, and associated metadata
- Researchers deposit spatial genetic variation data (e.g., mitochondrial DNA sequences) and analytical scripts for reproducibility

### **TreeBASE - Phylogenetic Trees & Matrices**

TreeBASE is a database of phylogenetic trees and sequence alignment matrices submitted voluntarily by researchers. As of July 2018, it contains 15,223 matrices and 20,246 trees covering ~117,231 distinct taxa across diverse organisms (34% plants, 33% fungi, 28% animals, with remaining bacteria, archaea, and viruses).

- Useful for comparative phylogeographic studies and reconciling topologies

### **MorphoBank - Morphological Matrices & Images**

Founded in 2001, MorphoBank is a collaborative Open Source and largely Open Access resource for the creation, preservation, and sharing of phylogenetic matrices, used by researchers in comparative biology, paleontology, and anthropology to construct phylogenetic trees with integrated morphological data and images.

- Complements genetic data with phenotypic traits for integrated phylogeographic analyses

---

## **INSDC - The International Nucleotide Sequence Database Collaboration**

The International Nucleotide Sequence Database Collaboration (INSDC) includes macrogenetic databases containing summaries of genetic variation from publications and genetic data repositories, plus associated contextual metadata like geographic coordinates.

- INSDC partners: GenBank (NCBI), EMBL (Europe), and DDBJ (Japan) exchange data daily
- Broader umbrella for nucleotide sequence data

---

## **Integrated Workflow Summary**

For comprehensive phylogeographic analysis, you might use:

1. **BOLD** or **UNITE** → for standardized barcodes (COI for animals, ITS for fungi)
2. **GenBank/INSDC** → for additional molecular markers (mitochondrial, nuclear genes)
3. **GBIF** → for species occurrence/distribution data
4. **Dryad** → for archiving your own datasets and supplementary phylogeographic data
5. **TreeBASE** → for comparing your phylogenetic trees with published analyses
6. **ZooKeys** → for publishing your findings with integrated genetic and taxonomic data

This integrated approach aligns with the emerging field of macrogenetics, which emphasizes the integration of large-scale genetic datasets with ecological data, collections science, biogeography, and phylogeography.
