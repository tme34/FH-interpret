<h3>This Github Repository includes databases used by the app (/data) and the Shiny application script itself.</h3>
Latest update: 08-01-2026

<h1>Documentation</h1>
The app contains data for use in classifying variants in the LDLR , APOB , and PCSK9 genes . For the selected variant, the app generates a draft interpretation text. 
The app describes the available evidence but does not make suggestions for the final classification. The following data is included:

<h2>The Genome Aggregation Database, gnomAD:</h2>
A database of publicly available aggregate data from large exome and genome sequencing projects from around the world. The interpretation text presents the total allele frequency for the selected variant in the approximately 1.6 million sequenced alleles currently available in gnomAD (v4.1.0). There is also a hyperlink to allele frequencies stratified by ethnicity for the selected variant (click on 'gnomAD').

<h2>UK Biobank:</h2>
UK Biobank is a large cohort of individuals from the British population. The app contains UK Biobank data for variants in the coding regions of the genes LDLR , APOB , and PCSK9 and their association with statin-adjusted plasma LDL-C in approximately 470,000 UK Biobank participants. The number of non-carriers, heterozygotes, and homozygotes, as well as the mean LDL-C concentration and standard deviation for each genotype are described in the interpretation text. The effects on LDL-C are also presented graphically (mean +/-95% confidence interval) in the app itself.  Participants on statins have had their measured LDL-C multiplied by 1.42 to correct for the statin-induced reduction in LDL-C of an average of 30%. Statistical testing for differences in LDL-C between genotypes is via linear regression adjusted for sex, age, and ethnicity (principal genetic components 1-10). These analyses were conducted under UK Biobank Application Number: 104807.  

<h2>Global Lipids Genetics Consortium (GLGC) exome chip GWAS meta-analysis from 2017 (PMID: 29083408):</h2>
Per-allele effects on LDL-C in approximately 300,000 individuals from 73 different cohorts. UK Biobank was not included in the meta-analysis, so these data can be considered independent. Only a small group of variants in LDLR , APOB and PCSK9 were measured in this GWAS. For the variants that were not measured, the statement regarding these data is omitted in the app's interpretation text.

<h2>AlphaMissense:</h2>
AlphaMissense is an AI tool developed by Google Deepmind that predicts the harmfulness of missense variants (PMID: 37733863). AlphaMissense is based on AlphaFold's ability to predict protein structures, as well as evidence for how rarely variation in a specific amino acid occurs across humans and primates; changes in amino acids that show little or no variation are more likely to be harmful. AlphaMissense harmfulness is expressed as a value from 0 (benign) to 1 (pathogenic), with recommended cutoff values: <0.34: benign, 0.34-0.56: uncertain; >0.56: likely pathogenic. The closer to 0 or 1, respectively, the more certain the prediction.

<h2>Human Gene Mutation Database (HGMD, https://www.hgmd.cf.ac.uk/ac/index.php):</h2>
Is a curated database of genetic variants in disease-associated genes. For each variant, there is information about the scientific articles in which the variant has appeared, which phenotype it is associated with, and the overall classification. The app includes data from all missense, nonsense and splice variants in LDLR , APOB , and PCSK9 that have been strongly associated with familial hypercholesterolemia (Variant Class 'DM', red box in HGMD). Latest extract: December 2025. For the selected variant, any referring articles are included as evidence in the interpretation text, as PMID. For variants with more than 4 referring articles, 4 are mentioned, followed by 'and x other references'.

<h2>ClinVar, https://www.ncbi.nlm.nih.gov/clinvar/:</h2>
A publicly available database of genetic variants and their association with diseases. All missense, nonsense, and splice variants in LDLR and PCSK9 that have been associated with familial hypercholesterolemia phenotypes have been extracted (latest extract: December 2025) and included as evidence in the app. Variant classification is divided into benign, likely benign, VUS, likely pathogenic, pathogenic, or 'unclear'. The app links to the ClinVar page for the selected variant in the interpretation response (click on ClinVar). ClinVar data for APOB is not included due to inconsistent naming of phenotype associations for variants in this gene (it is not possible to create a rule that reliably captures FH).

<h2>Other missense variants at the same amino acid position:</h2>
If another missense variant at the same amino acid position as the patient's variant is known to cause FH, this counts as moderate evidence of pathogenicity (ACMG criterion PM5). The app reports if other missense variants at the same position have been strongly associated with FH in HGMD (variant class 'DM') or in ClinVar (presumably pathogenic or pathogenic).

<h2>Effect on LDL uptake and LDLR transport to the cell surface:</h2>
In a study from 2025, 17,000 different (almost all theoretically possible) missense variants in LDLR were functionally investigated (PMID: 41166440). The variants were generated one at a time via CRISPR/CAS9 in cell lines. Each cell line with a specific mutation was examined for cellular uptake of LDL-C (two different assays, one without and one with the presence of VLDL) and for transport of the LDLR receptor to the cell surface. The function was quantified for each assay on a scale from 0 (no function) to 1 (normal, wild-type function). For each of the three assays, a score <0.5 is interpreted as impaired function. The app's interpretation text describes the overall result of these assays for the selected variant. Please note that the interpretation text only relates to data from this new study. Thus, older functional data may exist for a given variant that is not described in the app (such data may be identified via manual lookup in HGMD or ClinVar interpretation texts).

The app was developed by physician, PhD Helene Gellert-Kristensen, cand. scient. Tim Møller Eyrich, and senior physician, associate professor, PhD Stefan Stender, all Department of Clinical Biochemistry, Rigshospitalet. Contact: stefan.stender@regionh.dk

