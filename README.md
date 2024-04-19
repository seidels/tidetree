TiDeTree
========

TiDeTree is a BEAST 2 package that enables inference of time-scaled single-cell phylogenies and population dynamic parameters such as cell division, death, and differentiation rates from genetic lineage tracing data.

For further information please check out our [preprint](https://doi.org/10.1101/2022.02.14.480422) and find the code to reproduce the analyses in [this GitHub repo](https://github.com/seidels/tidetree-material).


To apply TiDeTree to your data, adapt the example XML file in the examples directory. For general guidance on setting prior distributions on parameters, look at [this BEAST2 tutorial](https://taming-the-beast.org/tutorials/Prior-selection/).


### Mini Tutorial

We migrated TiDeTree to BEAST 2.7. So you need to first [install JAVA 17](https://www.azul.com/downloads/?package=jdk#zulu) [as recommended here.](https://www.beast2.org/2022/08/22/what-will-change-in-v2-7-0-for-developers.html).

To run TiDeTree with the example.xml within the ./examples folder, use the following command:

```bash
java -jar bin/tidetree.jar examples/example.xml
```

### Parameter choices

A key component of TiDeTree is the editing and silencing model. Let's talk about editing first. TiDeTree can model any editing process, where an initially unedited site is edited once. Experimentally, this could be achieved by CRISPR-Cas9 editing, where a random edit is introduced. Or by recombinase editing, where an unedited site is either inverted or deleted.


We model this editing process as a 2-step process. The rate at which any edit occurs per unit of time is the clock rate r. 

![tidetree_editing_clock](https://github.com/seidels/tidetree/assets/13159214/5d0c9c47-f012-477a-903c-2eba60e35844)

Then, a particular edit $i$ is introduced with a relative editing rate $s_i$. You can also easily get the probability of introducing edit $i$  by taking $s_i / \sum_j s_j$   (This parameterisation is in TiDeTree's supplementary material page 5.  An alternative parameterisation is in the main text page 3, matrix Q 2.2).

![tidetree_editing](https://github.com/seidels/tidetree/assets/13159214/03e5f486-733c-4dbf-91d7-b48fddc39b12)


So how do we put this into our BEAST XML? We specify starting values for the clock and the editing rates. In this example, we have a clock rate of 1 per time unit and 3 editing rates for 3 different edited states.

```XML
<!-- the clock rate is the rate of acquiring edited state 2  -->
<parameter id="clockRate.c" spec="parameter.RealParameter" name="stateNode">
 1.0
</parameter>


<!--the editing rates -->
<parameter id="editRate" spec="parameter.RealParameter" lower="0.0" name="stateNode">
 0.8 0.1 0.1
</parameter>
        
```
These parameters are then used by the editing model (see the "@editRate" notation referencing the parameter with id="editRate"). 

```XML
<substModel id="substModel"
               spec="tidetree.substitutionmodel.EditAndSilencingModel"
               editRates="@editRate" silencingRate="@silencingRate"
               editHeight="54" editDuration="36">
    <frequencies spec="beast.base.evolution.substitutionmodel.Frequencies" frequencies="1 0 0 0" estimate="false"/>
  </substModel>
```
In the editing model we additinally have to specify the  editHeight and the editDuration. The editHeight is the amount of time between the start of editing ($t_1$) and the end of the experiment ($t_{end}$). The editDuration is the amount of time from the start of editing ($t_1$) until the end of editing ($t_2$). 

![editingSilencingModel](https://github.com/seidels/tidetree/assets/13159214/2c13a6fc-48e1-4349-9c8c-ccbdcb6ff4d7)

For example, let the entire experiment take 100 hours. Editing is induced after 10 hours and lasts for 10 hours. Then, the editHeight would be 90 hours and the editDuration 10 hours.

### TODO:
- [ ] integrate with [BEAUti](https://www.beast2.org/beauti/) to allow easy package installation
- [ ] Write a tutorial to be published on taming-the-beast.org 
