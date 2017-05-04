
Example
~~~~~~~~~~~

.. code:: python

    import fwdpy
    import numpy as np
    rng = fwdpy.GSLrng(100)
    ##Some basic parameters
    N=1000
    theta=100.0
    rho=100.0
    ##All neutral muts are [0,1)
    nregions = [ fwdpy.Region(0,1,1) ]
    #Selected mutations.  All are additive, to keep this example simple.
    ##Strongly deleterious mutations to the "left"
    ##Weaker mutations (2Ns = 10 on average) to the "right"
    ## 1% of selected mutations will be positively selected
    ## and uniform throughout the region.  The distribution
    ## of s will be exponential with mean 1e-3
    smodels = [fwdpy.ConstantS(-1,0,0.99/2,-0.1),fwdpy.ExpS(1,2,0.99/2,-10),fwdpy.ExpS(-1,2,0.01,0.001)]
    ##Recombination models--10x hotspot in the middl
    rregions = [fwdpy.Region(-1,1,1),fwdpy.Region(0.45,0.55,10)]
    #set up list of population sizes,
    #which are NumPy arrays of ints
    popsizes = np.array([N],dtype=np.uint32) 
    popsizes = np.tile(popsizes,10*N)
    pops = fwdpy.evolve_regions(rng,1,N,popsizes[0:],theta/(4*N),0.1*theta/(4*N),rho/(4*N),nregions,smodels,rregions)
    #Take a sample of size n = 10 from the population via list comprehension
    popSample = [fwdpy.get_samples(rng,i,100) for i in pops]

If we pass *the same neutral and selected mutation rates to evolve.regions*, then the above model satisfies:

* The total number of mutations occurring in Exon 2 is 2x the number occuring in Exon 1.
* Within an exon, 3/4 of all new mutations are deleterious.

.. code:: python

    #Let's convert from tuples to pandas DataFrames.
    #Ideally, one would further split each tuple element into a list,
    #but this example let's us get the point...
    import pandas
    neutralMuts = pandas.DataFrame.from_records(popSample[0][0],columns=['pos','genotypes'])
    selectedMuts = pandas.DataFrame.from_records(popSample[0][1],columns=['pos','genotypes'])

.. code:: python

    print(neutralMuts.head())


.. parsed-literal::

            pos                                          genotypes
    0  0.006588  0000000000000000000000001000010000000000000000...
    1  0.012971  1000000000010101000001000010101000000000000010...
    2  0.014566  0110000010000010000010010000000111001000010000...
    3  0.014729  0000000000000000000000000000000000000000100000...
    4  0.020637  0000010000000000000000000000000000000000000000...


.. code:: python

    print(selectedMuts.head())


.. parsed-literal::

            pos                                          genotypes
    0  1.016680  0101110011001110000011010011011001010000011100...
    1  1.192438  1010001100110001111100101100100110101111100011...
    2  1.197182  0101110011001110000011010011011001010000011100...
    3  1.365021  0101110011001110000011010011011001010000011100...
    4  1.402332  1010001100110001111100101100100110101111100011...


