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

### TODO:
- [ ] integrate with [BEAUti](https://www.beast2.org/beauti/) to allow easy package installation
- [ ] Write a tutorial to be published on taming-the-beast.org 
