# FH-interpret: Automating the Clinical Interpretation of LDLR, APOB, and PCSK9 Variants

[![Preprint](https://img.shields.io/badge/medRxiv-10.64898%2F2026.01.10.26343831-red)](https://doi.org/10.64898/2026.01.10.26343831)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**FH-interpret** is a web-based platform designed to automate the description of evidence for the clinical interpretation of variants in *LDLR*, *APOB*, and *PCSK9* used in FH (Familial Hypercholesterolemia) diagnosis.

<h3>This Github Repository includes databases used by the app (/data) and the Shiny application script itself.</h3>
Latest update: 14-07-2026

<h1>Documentation</h1>
FH-Interpret (v1.5) provides data for classifying variants in the LDLR, APOB, and PCSK9 genes. For the selected variant, the app generates a draft interpretation text. The app describes the available evidence but does not make suggestions for the final classification. The following data is included:

<h2>The Genome Aggregation Database, <a href='https://gnomad.broadinstitute.org/' target='_blank'>gnomAD</a>:</h2>
A database of publicly available aggregate data from large exome and genome sequencing projects from around the world. The interpretation text presents the total allele frequency for the selected variant in the approximately 1.6 million sequenced alleles currently available in gnomAD (v4.1.0). There is also a hyperlink to allele frequencies stratified by ethnicity for the selected variant (click on 'gnomAD').

<h2>UK Biobank:</h2>
UK Biobank is a large cohort of individuals from the British population. The app contains UK Biobank data for variants in the coding regions of the genes LDLR, APOB, and PCSK9 and their association with statin-adjusted plasma LDL-C in approximately 470,000 UK Biobank participants. The number of alleles and heterozygotes (but not homozygotes) as well as the mean LDL-C concentration for each genotype are described in the interpretation text. The effects on LDL-C are also presented graphically (mean +/-95% confidence interval) in the app itself.  Participants on statins have had their measured LDL-C multiplied by 1.43 to correct for the statin-induced reduction in LDL-C of an average of 30%. Statistical testing for differences in LDL-C between genotypes is via linear regression adjusted for sex, age, and ethnicity (principal genetic components 1-10). These analyses were conducted under UK Biobank Application Number: 104807. While rigorous quality control has been applied, variant calling and genotype-phenotype associations are probabilistic processes; users should assume that inherent errors persist.

<h2>Global Lipids Genetics Consortium (GLGC) exome chip GWAS meta-analysis from 2017 (PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/29083408/' target='_blank'>29083408</a>):</h2>
Per-allele effects on LDL-C in approximately 300,000 individuals from 73 different cohorts. The UK Biobank was not included in the meta-analysis, so these data can be considered independent. Only a small group of variants in LDLR , APOB and PCSK9 were measured in this GWAS. For the variants that were not measured, the statement regarding these data is omitted in the app's interpretation text.

<h2>AlphaMissense:</h2>
AlphaMissense is an AI tool developed by Google Deepmind that predicts the harmfulness of missense variants (PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/37733863/' target='_blank'>37733863</a>). AlphaMissense is based on AlphaFold's ability to predict protein structures, as well as evidence for how rarely variation in a specific amino acid occurs across humans and primates; changes in amino acids that show little or no variation are more likely to be harmful. AlphaMissense harmfulness is expressed as a value from 0 (benign) to 1 (pathogenic), with recommended cutoff values: <0.34: benign, 0.34-0.56: uncertain; >0.56: likely pathogenic. The closer to 0 or 1, respectively, the more certain the prediction.

<h2>REVEL:</h2>
<a href='https://sites.google.com/site/revelgenomics/about/' target='_blank'>Rare Exome Variant Ensemble Learner (REVEL v1.3) </a> scores were obtained from the UCSC genome track on genome build 38. Four files were obtained, one for each possible nucleotide substitution (A,C,G,T) at all positions across the exomes restricted to LDLR. Each variant combination of position and nucleotide substitution was merged to our variant database. Pathogenic classification by REVEL scores were flagged using the cutoff of 0.644 suggested by <a href='https://pubmed.ncbi.nlm.nih.gov/40084623/' target='_blank'>Bergquist et al.</a>

<h2>SpliceAI:</h2>
<a href='https://spliceailookup.broadinstitute.org/' target='_blank'> SpliceAI (v1.3)</a> scores were acquired as separate single-nucleotide variant (SNV) and indel files from Illumina’s BaseSpace Sequence Hub. For this application, unmasked (raw) pre-computed scores were used. SpliceAI provides predictive scores for four distinct splicing consequences: acceptor gain, acceptor loss, donor gain, and donor loss. A predicted impact on splicing was flagged and described in the narrative summary text if the maximum delta score was greater than or equal to 0.2, following the thresholds outlined by <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/' target='_blank' >Walker et al.</a>

<h2>Human Gene Mutation Database (<a href="https://www.hgmd.cf.ac.uk/ac/index.php" target="_blank">HGMD</a>):</h2>
Is a curated database of genetic variants in disease-associated genes. For each variant, there is information about which articles the variant has appeared in, which phenotype it associates with, and the overall classification of the variant. The app includes data on all missense, nonsense and splice variants in LDLR, APOB, and PCSK9 that have been strongly associated with familial hypercholesterolemia (Variant Class 'DM', red box in HGMD). Latest extract: December 2025. For the selected variant, the PMIDs of articles referencing the variant are included as evidence in the interpretation text. For variants with more than 4 articles, 4 are mentioned, followed by 'and x other references'.

<h2><a href='https://www.ncbi.nlm.nih.gov/clinvar/' target='_blank' >ClinVar</a>:</h2>
A publicly available database of genetic variants and their association with diseases. All missense, nonsense, and splice variants in LDLR and PCSK9 that have been associated with familial hypercholesterolemia phenotypes have been extracted (latest extract: December 2025) and included as evidence in the app. Variant classification is divided into benign, likely benign, VUS, likely pathogenic, pathogenic, or 'unclear'. The app links to the ClinVar page for the selected variant in the interpretation response (click on ClinVar). ClinVar data for APOB is not included due to inconsistent naming of phenotype associations for variants in this gene (it is not possible to create a rule that reliably captures FH).

<h2>Evaluation of amino acid position:</h2>
To streamline the evaluation of the PM5 criterion (novel missense change at an amino acid residue where a different missense change is determined to be pathogenic), the application cross-references ClinVar and the Human Gene Mutation Database (HGMD). It automatically identifies if alternative substitutions at the exact same residue as the selected variant are classified as 'Pathogenic/Likely Pathogenic' or 'Disease-Causing Mutation' (DM), respectively. Additionally, the application automates PM1 (variant located in a mutational hot spot or critical functional domain without benign variation). Specifically, two high-risk categories defined by <a href='https://pubmed.ncbi.nlm.nih.gov/34906454/' target='_blank'>expert consensus</a> are flagged: 1) LDLR exon 4 Missense Variants (flagged using GRCh38 genomic coordinates (chr19:11,105,220–11,105,600, corresponding to amino acid positions 105–232); and 2) Critical cysteine residues: Missense variants disrupting highly conserved disulfide bonds by substituting a cysteine residue within LDLR exons 2–8 or exon 14.

<h2>Effect on LDL uptake and LDLR abundance on the cell surface:</h2>
In a study from 2025, the functional consequence of 17,000 different missense variants (nearly all theoretically possible) in LDLR were investigated (PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/41166440/' target='_blank'>41166440</a>). The variants were generated one at a time using CRISPR/CAS9 modified HeLa cells. Each cell line with a specific mutation was then examined for cellular uptake of LDL-C (two different assays, one without and one with the presence of VLDL) and for abundance of the LDLR receptor on the cell surface. The function was quantified for each assay on a scale from 0 (no function) to 1 (normal, wild-type function). For each of the three assays, a score <0.5 is interpreted as impaired function. The full dataset can be downloaded from <a href='https://www.science.org/doi/10.1126/science.ady7186#supplementary-materials' target='_blank'>the supplementary files of Tabet al</a>. The app's interpretation text describes the overall result of these assays for the selected variant. Please note that the interpretation text only relates to data from this new study. Thus, older functional data that is not described in the app may exist for a given variant (such data may be identified through manual lookup in HGMD or ClinVar, which can be done using the links in the interpretation text).

<h2>Terms of use:</h2>
Use of this app and its content is at the investigator's sole risk. It is the user’s responsibility to implement necessary security precautions (e.g., virus scanning) for all downloaded content. By using this data, investigators certify compliance with all applicable local, national, and institutional regulations governing human genetics research. We welcome feedback regarding any significant errors identified during your analysis.

<h2>Contact information:</h2>
<h4>FH-Interpret (v1.5) is developed and maintained by MD, PhD <a href='https://scholar.google.com/citations?user=0zcd41YAAAAJ&hl=en&oi=ao/' target='_blank'>Helene Gellert-Kristensen</a>, cand. scient. <a href='https://scholar.google.com/citations?user=4H5xhzgAAAAJ&hl=en&oi=ao/' target='_blank'>Tim Møller Eyrich</a>, and MD, associate professor, PhD <a href='https://scholar.google.com/citations?user=mtgbiKoAAAAJ&hl=en&oi=ao/' target='_blank'>Stefan Stender</a>, all from the Department of Clinical Biochemistry at Rigshospitalet, Denmark. Contact: <a href='mailto:stefan.stender@regionh.dk'>stefan.stender@regionh.dk</a></h4>






