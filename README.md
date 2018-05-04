# Creating Consensus Tree Figures with Ancestral State Estimations in BayesTraits

These R scripts are used to generate figures of consensus trees with nodes that represent the ancestral state estimations and phylogenetic uncertainty. 

These are two steps to this process, with two separate R scripts. The first step is to generate the commands that you need to enter into BayesTraits before running analyses. This step requires a consensus tree with nodes ordered as you want them in your final images (see BEAST’s [TreeAnnotator](http://beast.community/treeannotator) for how to generate a consensus tree). An efficient way of entering commands into BayesTraits is to use a command file, rather than having to enter all commands line by line. For more information on using command files, check out the “Running BayesTraits with a command file” section the BayesTraits manual. One point of caution here though is that this can hide any warnings, so make sure you are happy with how it is running first. 
The second step involves reading in the output from [BayesTraits](http://www.evolution.rdg.ac.uk/BayesTraits.html) in [R](http://www.cran.r-project.org/), then using this output and the same consensus tree that you used in step one to create two pdf images representing the ancestral state estimation for each trait. These figures will need to be manually added together and edited. I recommend saving the figures in pdf format or some other vector based format to preserve image quality when rescaling and editing. 

### Example files
The example used in this script is based on the “Primates” dataset from BayesTraits. The tree files supplied here have been modified slightly to work with the scripts provided (modifications noted in the tree files). When using the fields supplied in this repository, the figures created should look like Primates_Trait_A.pdf and Primates_Trait_B.pdf.

### Notes:
* This script only works for the Discrete Dependent analyses in v3 of [BayesTraits](http://www.evolution.rdg.ac.uk/BayesTraits.html). The node orders and commands are different in earlier versions of BayesTraits. It shouldn’t be too hard for others to modify it to work with older version of BayesTraits or for use in the Discrete Independent or Multistate analyses.
* This is designed to be run with a sample of trees in BayesTraits, not a single tree. The only time a single tree is used is for the purposes of figure diagrams in R. 
* By “phylogenetic uncertainty” I mean the proportion of trees in the sample of trees for which that a corresponding node with the same descendants in found on the consensus tree.
* If you use this code, please cite the authors of [R](http://www.cran.r-project.org/), [BayesTraits](http://www.evolution.rdg.ac.uk/BayesTraits.html), and [ape](https://cran.r-project.org/web/packages/ape/index.html).
* I am not a programmer and I know that my code is not elegant or efficient. I appreciate suggested improvements, though it may take some time for me to work through them. If you spot any mistakes please log them as issues on the issues page so that others are aware of them too. 

### Examples of some papers with figures based on this code:
* [Sheehan, O., Watts, J., Gray, R. D., & Atkinson, Q. D. (2018). Coevolution of landesque capital intensive agriculture and sociopolitical hierarchy. Proceedings of the National Academy of Sciences, 115(14), 3628-3633.](https://doi.org/10.1073/pnas.1714558115)
* [Watts, J., Sheehan, O., Atkinson, Q. D., Bulbulia, J., & Gray, R. D. (2016). Ritual human sacrifice promoted and sustained the evolution of stratified societies. Nature, 532(7598), 228.](http://dx.doi.org/10.1038/nature17159)

### Other GitHub pages worth taking a looking at for BayesTraits and R users:
* [Sam Passmore’s excdr R package for useful R functions for BayesTraits.](https://github.com/SamPassmore/excdr)
