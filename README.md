[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blueviolet.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Generic badge](https://img.shields.io/badge/Version-6.0.4-informational.svg)](https://shields.io/)
[![Build and Publish](https://github.com/sirius-ms/sirius/actions/workflows/distribute.yaml/badge.svg?branch=release-4-pre)](https://github.com/sirius-ms/sirius/actions/workflows/distribute.yaml)

*<span style="color: #808080;">Our methods are offered to the scientific community as freely available resources. (Re-)distribution of the
methods, in whole or in part, for commercial purposes is prohibited.
CSI:FingerID and CANOPUS web services hosted by the [Böcker group](https://bio.informatik.uni-jena.de/) are for academic research and education use only.
Please review the [terms of service](https://bio.informatik.uni-jena.de/terms-of-service-fsu-csi) of the academic version for details.
For non-academic users, the [Bright Giant GmbH](https://bright-giant.com) provides licenses and all related services.
We ask that users of our tools cite the corresponding papers in any resulting publications.</span>*

COMET (Combinatorial Mass Encoding decoding Tool) is a java-based software framework for the analysis of LC-MS/MS data 
obtained from an affinity selection-mass spectrometry (AS-MS) experiment where self-encoded libraries were screened. 
The use case focuses mainly on combinatorial libraries where all library compounds represent distinct combinations of predefined building blocks. 
Nevertheless, COMET can also be used in conjunction with natural product libraries as screening libraries.

Currently, the [SIRIUS](https://github.com/sirius-ms/sirius) platform is used in which ASMS-related features are integrated,
e.g. excluding features (or MS/MS spectra) unrelated to compounds of the screened combinatorial library and ranking of candidate structures using EPIMETHEUS.

Main developers of COMET are the [Böcker group](https://bio.informatik.uni-jena.de/) and the [Bright Giant GmbH](https://bright-giant.com)

### Documentation of SIRIUS
- [Online Documentation](https://v6.docs.sirius-ms.io/)
- [Video tutorials](https://www.youtube.com/channel/UCIbW_ZFSADRUQ-T5nmgU4VA/featured)
- [Bookchapter on using SIRIUS 4](https://doi.org/10.1007/978-1-0716-0239-3_11) ([Preprint](https://bio.informatik.uni-jena.de/wp/wp-content/uploads/2020/12/SIRIUS4_book_chapter_preprint-2.pdf)) -- does not cover the new LC-MS/MS processing option
- [Demo data](data/demo.zip)
- [Logos for publications and presentations](https://bio.informatik.uni-jena.de/software/sirius/sirius-logos/)

### Installation and Dependencies
COMET is available for Windows (64bit), MacOS (64bit), and Linux (64bit) and can be installed via the fowllowing links:
<!--begin download-->
- for Windows (64bit): [msi](https://github.com/sirius-ms/comet/releases/download/v6.0.4-SNAPSHOT/sirius-6.0.4-SNAPSHOT-win64.msi) / [zip](https://github.com/sirius-ms/comet/releases/download/v6.0.4-SNAPSHOT/sirius-6.0.4-SNAPSHOT-win64.zip)
- for Mac (64bit): [pkg](https://github.com/sirius-ms/comet/releases/download/v6.0.4-SNAPSHOT/sirius-6.0.4-SNAPSHOT-osx64.pkg) / [zip](https://github.com/sirius-ms/comet/releases/download/v6.0.4-SNAPSHOT/sirius-6.0.4-SNAPSHOT-osx64.zip)
- for Linux (64bit): [zip](https://github.com/sirius-ms/comet/releases/download/v6.0.4-SNAPSHOT/sirius-6.0.4-SNAPSHOT-linux64.zip)
<!--end download-->

All (including previous) releases can be found [here](https://github.com/sirius-ms/comet/releases).

A typical install time should not exceed 10min which is mostly dependend on the speed of the internet connection for downloading the installation files. The installation and functionality of COMET was successfully tested on Windows 10 (x64) and on Ubuntu 24.04.2 (x64).

#### [Installation Instructions](https://v6.docs.sirius-ms.io/install)
For  Windows and MacOS, the installer version of COMET (msi/pkg) should be preferred but might require administrator permissions.
Since we do not pay Microsoft/Apple for certification, you might have to confirm that you want to trust "software from
an unknown source" on Windows/MacOS when using the .msi/.pkg installers.

If you choose to download the .zip file corresponding to your operating system, you have to extract that .zip file into a directory where you have writing permissions, e.g. ```C:\COMET```. To start COMET, you have to execute ```sirius.exe``` in that folder.
As COMET is currently using the SIRIUS platform, you can also chek out the [documentation](https://v6.docs.sirius-ms.io/install) for more details about the installation procedure.

#### Dependencies
All installation versions of COMET include the Java Runtime Environment (JRE). Therefore, there is no need to install Java separately. 
In case you have screened your own combinatorial molecule library, you have to create a .csv file containing the building blocks of this library in order to use the COMET filter. A description of how such a .csv file looks like is given below.
As we recommend to use our own scripts to create such a .csv file, you need to have ```python```,```jupyter notebook```, and the python packages ```rdkit``` and ```pandas``` installed.

#### [Creating a user account](https://v6.docs.sirius-ms.io/account-and-license/)
User accounts can be created directly via the COMET/SIRIUS GUI. Please, use your **institutional email address**. SIRIUS
web services are free for academic/non-commercial use. Usually academic institutions are identified by their
email domain and access will be granted automatically. In some cases, further validation of your academic/non-commercial
may be required.
[See also SIRIUS Documentation – Account and License](https://v6.docs.sirius-ms.io/account-and-license/).

#### [Sources on GitHub](https://github.com/sirius-ms)
- [SIRIUS](https://github.com/sirius-ms/sirius)
- [SIRIUS-API Java SDK](sirius_nightsky_sdk/sirius_nightsky_sdk.openapi/README.md)
- [SIRIUS-API SDKs](https://github.com/sirius-ms/sirius-client-openAPI)

#### [Changelog](https://v6.docs.sirius-ms.io/changelog/)

### How to use COMET

#### Demo Data
This [Zenodo repository](https://doi.org/10.5281/zenodo.14070388) contains the obtained LC-MS/MS data, the scripts, and other supplementary files used in [(van der Nol et al., 2025)](https://doi.org/10.26434/chemrxiv-2025-s3v2z). Here, you can also find a file called ```demo_data.zip```. This archive contains three files belonging to a 500 membered combinatorial molecule library where each molecule consists of a benzimidazole scaffold decorated with one amino acid building block, one amine building block, and one aldehyde building block:
- ```ENL161_50uM_100fmol_SCE15-25_27112023.mzML``` is the LC-MS/MS data obtained by measuring the whole synthesized libary via nanoLC-MS/MS
- ```ENL161_CustomDB.tsv``` contains all structures of that library in the form of SMILES strings. Each structure has a unique ```id``` and ```name``` which represents its composition of building blocks.
- ```ENL161_CustomDB_BBs.csv``` contains the building blocks for each position.

Since each class of building block (e.g. amino acids) only occurs at a predefined position in the scaffold and this position never changes in the entire library, unique indices can be assigned to these positions. In case of this molecule library, this could be the index ```0``` for the amino acids, ```1``` for the amines and ```2``` for the aldehydes. In ```ENL161_CustomDB_BBs.csv```, you will find that each building block is assigned with such a position which is called ```bb_pos```. Additionally, this file contains for each building block its SMILES string (```smiles```), its corresponding molecular formula (```formula```), the formula of the loss when incorporated into the final molecule (```reaction_loss```), and an ```id``` specifying the exact building block in its class.

This is how ```ENL161_CustomDB_BBs.csv``` looks like: 
```
bb_pos,smiles,formula,reaction_loss,id
0,OC([C@H](NC(OCC1c2c(c3c1cccc3)cccc2)=O)CC4CCCCC4)=O,C24H27NO4NH,C15O3H11,1
0,O=C([C@@H](NC(OCC1C2=CC=CC=C2C3=CC=CC=C13)=O)C)O,C18H17NO4NH,C15O3H11,2
0,OC([C@@H]1CSCN1C(OCC2c3c(c4c2cccc4)cccc3)=O)=O,C19H17NO4SNH,C15O3H11,3
...
```
For example, the building block ```2``` at position ```0``` has the SMILES string ```O=C([C@@H](NC(OCC1C2=CC=CC=C2C3=CC=CC=C13)=O)C)O``` and corresponding formula ```C18H17NO4NH``` before synthesis (not incorporated into the final molecule). When incorporated into the final molecule, its molecular formula changes to ```C3H7N2O```. Therefore, ```C15O3H11``` is the molecular formula which describes this loss.

To create such a CSV file for your own libraries, you can use the Jupyter Notebook ```COMET_Building blocks_input.ipynb``` which is also part of that [Zenodo repository](https://doi.org/10.5281/zenodo.14070388).

#### LC-MS/MS Data Import
Once you have started the COMET GUI (by executing ```sirius.exe```), you can import your measured MS/MS data via the "Import" button or Drag and Drop to the left most panel. COMET/SIRIUS supports multiple MS data formats:
- .ms, .mgf, .mat/.msp, and Agilent’s .cef: These formats contain pre-processed peak lists for each feature.
- .mzml, and .mzxml: For these formats, feature detection and alignment will be performed.

 **Note that, all data must be centroided and that raw file formats are not supported.** For more information about the import of MS data, see [here](https://v6.docs.sirius-ms.io/io/#input).

**Example:**
Let's import ```ENL161_50uM_100fmol_SCE15-25_27112023.mzML``` via the "Import" button or drag-and-drop. Here, a small window with a progess bar opens up showing you the current steps of the preprocessing workflow. This takes about 5min on a laptop with a 11th Gen Intel(R) Core(TM) i7-11859H processor and 48Gb of RAM. After the import, the detected features are shown in the left most panel. At the bottom of this panel, you will find three numbers: "0 of 4865 (24619) selected". The number in paranthesis is the total number of detected features which is 24,619. The number 4,865 refers to the number of features which are obtained using the default filter settings (each feature has at least one MS/MS scan and has at least a decent feature quality). The first number (here 0) is the number of features in the selection. You can select a custom number of features which you might want to analyze.  

#### Feature Filtering with COMET
To use COMET's feature filtering, you have to click on the button with the three dots "..." next to the search text field ("Type and hit enter to search"). A panel called "Filter configuration" opens. Go to the "COMET" tab. Here, you have to specify the architecture of your combinatorial molecule library first. You do this by selecting the CSV file containing all the building blocks (described in section "Demo Data") and the molecular formula of the scaffold when incorporated in the final molecule. In case of the 500 membered benzimidazole library, this would be the file path to ```ENL161_CustomDB_BBs.csv``` and ```C8H3N2O``` as the molecular formula of the scaffold. Additionally, you have the option to set up the MS1 mass accurarcy and fallback adducts. Let's keep the default parameters of ```10ppm``` and ```[M + H]+```. When pressing "Apply" the first filter will be applied; i.e. all features will be filtered out those precursor mass doesn't match any compound of the library. This takes about 30sec and results in 806 filtered features. Note: when pressing "Apply", it's normal that the window freezes for a moment because there is currently no progress bar showing that something is computed.

To filter the features according to their fragmentation pattern and retain only those with a fragmentation pattern characteristic to the library molecules (i.e. fragmentation by cleaving off building blocks), you have to open the COMET filtering panel again. Here, you have to check "Enable peak matching filter". 
Now, you can specify which fragments you would expect; i.e. the fragments which are more likely to form and result by cleaving of single building blocks. If you leave this text field empty, all such fragments will be considered.\
Let's specify ```0,1,S[0;2],S[1;2]``` for the benzimidazole library from the example. Here, ```0``` specifies the single amino acid building block, ```1``` specifies the single amine, ```S[0;2]``` specifies the scaffold plus amino acid and aldehyde, and ```S[1;2]``` is the scaffold plus amine and aldehyde. See [(van der Nol et al., 2025)](https://doi.org/10.26434/chemrxiv-2025-s3v2z) for a more detailed description of how to specify fragments in the COMET filtering panel.\
In addition to the fragment specification, you can chose the minimum number of peaks which should be explainable by any of the specified fragments. You also have to provide how many of the highest intensity peaks should be considered per spectrum. Let's say that the ```5``` highest intensity peaks should be considered and at least ```2``` of those fragment peaks should be explainable by at least one candidate's building block fragment.\
As hydrogen rearrangments can occur during collision-induced dissociation (CID), you can also specify the maximum allowed number of hydrogen atom masses a theoretical fragment can deviate from its fragment peak. Let's say ```2``` allowed hydrogen shifts and let's change the MS2 mass accuracy to ```5ppm``` as the MS data was measured on an Orbitrap.\
If you specify an output location, a CSV file with the information about the filtering will be stored. This CSV file will contain all filtered features and it will tell you which fragment peaks can be explained with which building block fragments of which candidate structure. If you do not specifiy an output location, these information will be printed into the console.\
Applying the filter with these settings results in 285 filtered features. The filtering takes about 5min on a laptop with a 11th Gen Intel(R) Core(TM) i7-11859H processor and 48Gb of RAM. **Note again, it will look like that the window freezes. But in reality, the filtering is computed and will take some minutes.**\
For more information, see [(van der Nol et al., 2025)](https://doi.org/10.26434/chemrxiv-2025-s3v2z).

#### Structure Annotation

##### Creating a custom database with all library molecules
As we already know that the measured structures are a subset of all library compounds, we only want to annotate MS/MS spectra or features with structures of that library; i.e. we want to search in that library for potential candidates. [Here](https://v6.docs.sirius-ms.io/gui/#custom-database-import), you will find a description on how to create a custom database. Regarding the example library from above, we would open the ```Databases``` dialog, click on the ```Create custom Database``` button, enter the name ```ENL161``` as the database name in COMET, and specify the location of that custom database file (called ```enl161.siriusdb```). Then, we would add the TSV file ```ENL161_CustomDB.tsv``` and press the ```Import structures and spectra``` button. After creation, this new custom database ```enl161``` is shown together with the information of having 500 compounds with 455 different molecular formulas.

##### Compute Panel
To perform the structure annoation, you have to open the ```Compute``` dialog. Here, you have to activate the molecular formula identification with SIRIUS by clicking on the ```SIRIUS``` button. It's important to use the option ```Database search``` in the drop-down list for ```Molecular formula generation``` and define your custom database as the only database to search in, e.g. ```ENL161``` from the example above. 

Furthermore, you have to click on the ```Predict``` button for the fingerprint prediction and the ```Sarch DBs``` button for the structure annotation. Regarding the latter, we recommend you to use the option ```Rank with EPIMETHEUS``` in case of combinatorial molecule library. Again, as search database you should only use your created custom database, e.g. ```ENL161``` from the example above.\
Now press the ```Compute``` button and the MS/MS spectra or features will be annotated with their corresponding candidate strucutes. Those candidate structures are ranked according to their assigned score.

The computation of those remaining 285 features takes about 10min on a laptop with a 11th Gen Intel(R) Core(TM) i7-11859H processor and 48Gb of RAM.

### Integration of CSI:FingerID, CANOPUS and MSNovelist

Fragmentation trees and spectra can be directly uploaded from SIRIUS to the CSI:FingerID, CANOPUS and MSNovelist web services.
Results are retrieved from the web service and can be displayed in the SIRIUS graphical user interface. This functionality is
also available for the SIRIUS command-line tool. Training structures for CSI:FingerID's predictors are available through the CSI:FingerID web API:
<!--begin training-->

- https://www.csi-fingerid.uni-jena.de/v3.0/api/fingerid/trainingstructures?predictor=1 (training structures for positive ion mode)
- https://www.csi-fingerid.uni-jena.de/v3.0/api/fingerid/trainingstructures?predictor=2 (training structures for negative ion mode)

<!--end training-->

### Fragmentation Tree Computation

The manual interpretation of tandem mass spectra is time-consuming and
non-trivial. SIRIUS analyses the fragmentation pattern resulting in
a hypothetical fragmentation tree, in which nodes are annotated with
molecular formulas of the fragments and arcs (edges) represent fragmentation
events (losses). SIRIUS allows for the automated and high-throughput analysis of
small-compound MS data beyond elemental composition without requiring
compound structures or a mass spectral database.

### Isotope Pattern Analysis

SIRIUS deduces molecular formulas of small compounds by ranking isotope
patterns from mass spectra of high resolution. After preprocessing, the
output of a mass spectrometer is a list of peaks which corresponds to
the masses of the sample molecules and their abundance. In principle,
elemental compositions of small molecules can be identified using only
accurate masses. However, even with very high mass accuracy, many
formulas are obtained in higher mass regions. High resolution mass
spectrometry allows us to determine the isotope pattern of sample
molecule with outstanding accuracy and apply this information to
identify the elemental composition of the sample molecule. SIRIUS can be
downloaded either as graphical user interface (see Sirius GUI) or as
command-line tool.

<!--begin cite-->
## Main citations

### Main citations for COMET related features

Edith van der Nol, Nils Alexander Haupt, Qing Qing Gao, Benthe A.M. Smit, Martin Andre Hoffmann, Martin Engler-Lukajewski
Marcus Ludwig, Sean McKenna, J. Miguel Mata, Olivier Bequignon, Gerard van Westen, Tiemen J. Wendel, Sylvie M. Noordermeer, Sebastian Böcker, and Sebastian Pomplun. 
[Barcode-free hit discovery from massive libraries enabled by automated small molecule structure annotation](https://doi.org/10.26434/chemrxiv-2025-s3v2z). 
*ChemRxiv*, 2025.

### Main citations for SIRIUS related features

Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Alexander A. Aksenov, Alexey V. Melnik, Marvin Meusel, Pieter C. Dorrestein, Juho Rousu, and Sebastian Böcker,
[SIRIUS 4: Turning tandem mass spectra into metabolite structure information.](https://doi.org/10.1038/s41592-019-0344-8)
*Nature Methods* 16, 299–302, 2019.

---
Stravs, Michael A. and Dührkop, Kai and Böcker, Sebastian and Zamboni, Nicola
[MSNovelist: de novo structure generation from mass spectra](https://doi.org/10.1038/s41592-022-01486-3)
*Nature Methods* 19, 865–870, 2022. (Cite if you are using: MSNovelist)

Martin A. Hoffmann and Louis-Félix Nothias and Marcus Ludwig and Markus Fleischauer and Emily C. Gentry and Michael Witting and Pieter C. Dorrestein and Kai Dührkop and Sebastian Böcker
[High-confidence structural annotation of metabolites absent from spectral libraries](https://doi.org/10.1101/2021.03.18.435634)
*Nature Biotechnology* 40, 411–421, 2022. (Cite if you are using: *CSI:FingerID*, *COSMIC*)

Kai Dührkop, Louis-Félix Nothias, Markus Fleischauer, Raphael Reher, Marcus Ludwig, Martin A. Hoffmann, Daniel Petras, William H. Gerwick, Juho Rousu, Pieter C. Dorrestein and Sebastian Böcker.
[Systematic classification of unknown metabolites using high-resolution fragmentation mass spectra.](https://doi.org/10.1038/s41587-020-0740-8)
*Nature Biotechnology*, 2020. (Cite if you are using *CANOPUS*)

Yannick Djoumbou Feunang, Roman Eisner, Craig Knox, Leonid Chepelev, Janna Hastings, Gareth Owen, Eoin Fahy, Christoph Steinbeck, Shankar Subramanian, Evan Bolton, Russell Greiner, David S. Wishart.
[ClassyFire: automated chemical classification with a comprehensive, computable taxonomy.](https://doi.org/10.1186/s13321-016-0174-y)
*Journal of Cheminformatics* 8, 61, 2016. (*ClassyFire* publication; cite this if you are using *CANOPUS*)

Marcus Ludwig, Louis-Félix Nothias, Kai Dührkop, Irina Koester, Markus Fleischauer, Martin A. Hoffmann, Daniel Petras, Fernando Vargas, Mustafa Morsy, Lihini Aluwihare, Pieter C. Dorrestein, Sebastian Böcker.
[Database-independent molecular formula annotation using Gibbs sampling through ZODIAC.](https://doi.org/10.1038/s42256-020-00234-6)
*Nature Machine Intelligence* 2, 629–641, 2020. (Cite if you are using *ZODIAC*)

Kai Dührkop and Sebastian Böcker.
[Fragmentation trees reloaded.](http://dx.doi.org/10.1007/978-3-319-16706-0_10)
*Journal of Cheminformatics* 8, 5, 2016. (Cite this for *fragmentation pattern analysis and fragmentation tree computation*)

Kai Dührkop, Huibin Shen, Marvin Meusel, Juho Rousu, and Sebastian Böcker.
[Searching molecular structure databases with tandem mass spectra using CSI:FingerID](http://dx.doi.org/10.1073/pnas.1509788112).
*Proceedings of the National Academy of Sciences U S A* 112(41), 12580-12585, 2015. (cite this when *using CSI:FingerID*)

Sebastian Böcker, Matthias C. Letzel, Zsuzsanna Lipták and Anton Pervukhin.
[SIRIUS: decomposing isotope patterns for metabolite identification.](http://bioinformatics.oxfordjournals.org/content/25/2/218.full)
*Bioinformatics* 25(2), 218-224, 2009. (Cite this for *isotope pattern analysis*)

### Additional citations

Marcus Ludwig, Kai Dührkop and Sebastian and Böcker.
[Bayesian networks for mass spectrometric metabolite identification via molecular fingerprints.](http://doi.org/10.1093/bioinformatics/bty245)
*Bioinformatics*, 34(13): i333-i340. 2018. Proc. of Intelligent Systems for Molecular Biology (ISMB 2018). (Cite for CSI:FingerID Scoring)

W. Timothy J. White, Stephan Beyer, Kai Dührkop, Markus Chimani and
Sebastian Böcker. [Speedy Colorful
Subtrees.](http://dx.doi.org/10.1007/978-3-319-16706-0_10) In *Proc. of
Computing and Combinatorics Conference (COCOON 2015)*, volume 9198 of
*Lect Notes Comput Sci*, pages 310-322. Springer, Berlin, 2015. (cite
this on *why computations are swift*, even on a laptop computer)

Huibin Shen, Kai Dührkop, Sebastian Böcker and Juho Rousu. [Metabolite
Identification through Multiple Kernel Learning on Fragmentation
Trees.](http://dx.doi.org/10.1093/bioinformatics/btu275)
*Bioinformatics*, 30(12):i157-i164, 2014. Proc. of *Intelligent Systems
for Molecular Biology* (ISMB 2014). (Introduces *the machinery behind
CSI:FingerID*)

Imran Rauf, Florian Rasche, François Nicolas and
Sebastian Böcker. [Finding Maximum Colorful Subtrees in
practice.](http://dx.doi.org/10.1089/cmb.2012.0083) *J Comput Biol*,
20(4):1-11, 2013. (More, earlier work on *why computations are swift*
today)

Heinonen, M.; Shen, H.; Zamboni, N.; Rousu, J. [Metabolite
identification and molecular fingerprint prediction through machine
learning](http://dx.doi.org/10.1093/bioinformatics/bts437).
*Bioinformatics*, 2012. Vol. 28, nro 18, pp. 2333-2341. (Introduces the
*idea of predicting molecular fingerprints* from tandem MS data)

Florian Rasche, Aleš Svatoš, Ravi Kumar Maddula, Christoph Böttcher, and
Sebastian Böcker. [Computing Fragmentation Trees from Tandem Mass
Spectrometry
Data](http://pubs.acs.org/doi/abs/10.1021/ac101825k). *Analytical
Chemistry* (2011) 83 (4): 1243–1251. (Cite this for *introduction of
fragmentation trees* as used by SIRIUS)

Sebastian Böcker and Florian Rasche. [Towards de novo identification of metabolites by analyzing
tandem mass
spectra](http://bioinformatics.oxfordjournals.org/content/24/16/i49.abstract).
*Bioinformatics* (2008) 24 (16): i49-i55. (The very *first paper to
mention fragmentation trees* as used by SIRIUS)

<!--end cite-->

## License

Starting with version 4.4.27, SIRIUS is licensed under the [GNU Affero General
Public License (GPL)](https://www.gnu.org/licenses/agpl-3.0.txt). If you integrate SIRIUS into other software, we
strongly encourage you to make the usage of SIRIUS as well as the literature to cite transparent to the user.

## Acknowledgements
#### Thanks for supporting the development of SIRIUS!
[![MSCJ Logo](https://www.mscj.uni-jena.de/wp-content/uploads/2015/05/logo-svg-text-horizontal.svg)](https://www.mscj.uni-jena.de)
