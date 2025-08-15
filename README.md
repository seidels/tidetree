TiDeTree
========

TiDeTree is a [BEAST 2[(https://www.beast2.org/) package that enables inference of time-scaled single-cell phylogenies and population dynamic parameters such as cell division, death, and differentiation rates from genetic lineage tracing data.

For further information please check out our [publication](https://doi.org/10.1098/rspb.2022.1844) or [preprint](https://doi.org/10.1101/2022.02.14.480422) and find the code to reproduce the analyses in [this GitHub repo](https://github.com/seidels/tidetree-material).


If you want to apply TiDeTree to your data, check out the tutorial [here](https://taming-the-beast.org/tutorials/TiDeTree-Tutorial/). If you are comfortable with BEAST xml hacking, take a look at the example XML file in the examples directory. For general guidance on setting prior distributions on parameters, look at [this BEAST2 tutorial](https://taming-the-beast.org/tutorials/Prior-selection/).


## Installation

1. **Install BEAST 2**  
   Download and install BEAST 2 by following the instructions on the [official BEAST 2 website](https://www.beast2.org/).

2. **Add TiDeTree to the Package Manager**  
   To install TiDeTree via BEAUTi:

   - Open **BEAUTi**.
   - Go to **File > Manage Packages**.
   - Under **Package Repositories**, click **Add URL**, and enter:  
     ```
     https://raw.githubusercontent.com/seidels/tidetree/main/package.xml
     ```

   This will add TiDeTree to the list of available packages.

3. **Install TiDeTree**  
   - Return to the list of available packages.
   - Select **TiDeTree** and click **Install**.
   - Close and reopen **BEAUTi** to complete the installation.

**You're all set to run analyses with TiDeTree!**


Mini Tutorial
------------

You can either set up your analyses using the graphical user interface BEAUTI, or, if you are more familiar with BEAST xml already, by editing the examle.xml provided within the ./examples folder.


### Parameter choices

A key component of TiDeTree is the editing and silencing model. Let's talk about editing first. TiDeTree can model any editing process, where an initially unedited site is edited once. Experimentally, this could be achieved by CRISPR-Cas9 editing, where a random edit is introduced. Or by recombinase editing, where an unedited site is either inverted or deleted.


We model this editing process as a 2-step process. The rate at which any edit occurs per unit of time is the clock rate r. 

![tidetree_editing_clock](https://github.com/seidels/tidetree/assets/13159214/5d0c9c47-f012-477a-903c-2eba60e35844)

Then, a particular edit $i$ is introduced with a relative editing rate $s_i$. You can also easily get the probability of introducing edit $i$  by taking $s_i / \sum_j s_j$   (This parameterisation is in TiDeTree's supplementary material page 5.  An alternative parameterisation is in the main text page 3, matrix Q 2.2).

![tidetree_editing](https://github.com/seidels/tidetree/assets/13159214/03e5f486-733c-4dbf-91d7-b48fddc39b12)

In the editing model we additinally have to specify the  editHeight and the editDuration. The editHeight is the amount of time between the start of editing ($t_1$) and the end of the experiment ($t_{end}$). The editDuration is the amount of time from the start of editing ($t_1$) until the end of editing ($t_2$). 

![editingSilencingModel](https://github.com/seidels/tidetree/assets/13159214/2c13a6fc-48e1-4349-9c8c-ccbdcb6ff4d7)

For example, let the entire experiment take 100 hours. Editing is induced after 10 hours and lasts for 10 hours. Then, the editHeight would be 90 hours and the editDuration 10 hours.


### Setting Up the Editing Model in BEAUTi

1. **Start a TiDeTree Analysis**  
   Open **BEAUTi** and go to **File > Template > tidetree**. This sets the analysis format for TiDeTree analyses.

2. **Load Your Alignment**  
   Import your data via **File > Import Alignment**, and select your `.tidetree` file (e.g., `./examples/example.tidetree`).

3. **Specify Tip Dates**  
   Navigate to the **Tip Dates** panel and check **Use Tip dates**.  
   For now, you can leave all entries set to `0`. After setting all the panels, go back to the XML and change the entries to the length of your experiment (in units of time).

4. **Configure Site Model**  
   In the **Site Model** panel, you will see that the **Edit Rates** and **Silencing Rate** have already been included and are set to be estimated during inference. Modify the **Edit Height** and **Edit Duration** to reflect your specific experimental setup.

5. **Other Panels**  
   Most other panels can be configured as described in the [standard BEAST 2 tutorials](https://taming-the-beast.org/tutorials/Introduction-to-BEAST2/). Below we highlight a few TiDeTree-specific considerations:

6. **Clock Model**  
   In the **Clock Model** panel, choose a clock model appropriate for your analysis.  
   If unsure, start with the **Strict clock**.

7. **Tree Initialisation**  
   Go to the **Initialisation** panel and expand the **Tree** section.  
   Enter values for the **Root Height**, **Edit Duration**, and **Edit Height** to ensure a reasonable starting tree is used.

8. **Set Priors**  
   In the **Priors** panel, define prior distributions for your parameters based on biological knowledge.  
   If you're uncertain about choosing priors, refer to the [Prior Selection tutorial](https://taming-the-beast.org/tutorials/Prior-selection/).

9. **Configure MCMC Settings**  
    Set the **MCMC chain length**, sampling frequency, and other relevant settings to match your analysis goals.

10. **Save and Run**  
    Save your analysis via **File > Save As**, which will generate an `.xml` file.  
    Run this file with **BEAST 2**.


### Setting up the editing model via modifying the XML

So how do we put this into our BEAST XML? We specify starting values for the clock and the editing rates. In this example, we have a clock rate of 1 per time unit and 3 editing rates for 3 different edited states.

```XML
<!-- the clock rate  -->
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
