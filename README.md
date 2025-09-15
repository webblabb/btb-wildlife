# Code associated with Stochastic modeling of bovine tuberculosis dynamics in white-tailed deer

\textit{Mycobacterium bovis}, the causative agent of the disease bovine tuberculosis (bTB), poses potential economic threats to the cattle industry.  Although currently prevalent at low levels in cattle herds in the United States, the presence of the pathogen in wildlife reservoirs provides the seed for continued outbreaks in livestock. The northern lower peninsula of Michigan is a known habitat for the white-tailed deer (\textit{Odocoileus virginianus}), which have been observed to serve as a reservoir host for \textit{M. bovis} in this region. It is suspected that this reservoir was established by spillover from livestock into wildlife populations, and it has since contributed to repeated spillback into livestock, hindering eradication of the disease in the United States. 

Here, we use a stochastic simulation to predict bTB dynamics within a wildlife population to illustrate the mechanistic drivers of bTB outbreaks and disease prevalence in wildlife populations. We account for seasonal variation in deer populations with a seasonal birth pulse and mortality due to hunter-harvest. Heterogeneity in behavior and immunology across individuals is also considered, as increased farm contact can result in increased infection risk and individual level variation in immune function can significantly affect the progression of a bTB infection.

We find that the probability of developing an endemic disease state depends on the number of infected animals rather than the size of the herd. Fadeout of bTB infections appears to be unlikely overall. We find that outbreaks that do fade out are likely driven by stochastic effects, as fade out was only observed in simulations where the outbreak failed to take off. It is also apparent that increased farm contact does not contribute to larger outbreaks in wildlife, but it may still have implications in spillback of the disease into livestock populations. Assuming perfect diagnostics and hunter compliance, we also find that prevalence estimates using postmortem diagnostics on hunter-harvested animals provides a systematic underestimate of true prevalence, especially when disease prevalence is low within the population. 

Preventive measures to reduce spillover events should be implemented given the low probability of fadeout in wildlife. Physical boundaries such as increased fencing and covered storage of feed and water troughs may reduce transmission potential via livestock-wildlife contact, possibly lowering bTB prevalence in both species. Although prevalence estimates from hunter harvest may be imperfect, change between years could be used as a tool to identify outbreaks or evaluate intervention effectiveness.

# How to use

pipeline.R runs all simulations and pub_plots.R generates all figures in the manuscript.

## Contact
Please contact [Lindsay Beck-Johnson](mailto:L.Beck-Johnson@colostate.edu) with any questions.
