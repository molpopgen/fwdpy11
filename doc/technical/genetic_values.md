(gvalue-technical-details)=

# Technical overview of genetic value calculations

## A mock overview in Python

This section illustrates the technical details of genetic value objects
by presenting a mock up of their `API` in Python. To help document
key concepts, we make use of Python 3's type hinting features from
{mod}`typing` and decorators from {mod}`abc`:

* A class inheriting from {class}`abc.ABC` is an abstract base class.
  It is not possible to create instances of ABCs.  Rather, you must
  subclass them.
* A function decorated with {func}`abc.abstractmethod` is
  a part of an ABC's interface that must be defined by a subclass.
  In other words, to be a valid `MockDiploidGeneticValue`, subclasses
  must define the abstract methods.  In C++ terminology,
  these are "pure virtual" functions.
* A function decorated with {func}`typing.final` must not be redefined
  in a subclass. Doing so is considered an error, although Python
  will allow it.
* In Python, any undecorated functions are viewed as default implementations
  of the class' behavior.  A subclass may or may not redefine
  them as needed.

These concepts do not map perfectly from Python to C++, but it is close
enough for now. Simply put, C++ is much stricter than Python.

The following code block shows a Python mock-up corresponding to
{class}`fwdpy11.DiploidGeneticValue`. (The next section shows the actual
definition of {class}`fwdpy11.DiploidGeneticValue`.)

This class is an abstract base class.  The initialization semantics are:

* The object is initialized with the dimensionality of the model, a
  genetic value to fitness map, and a noise function.
* If the genetic value to fitness map is `None`, then the genetic
  value is assumed to be fitness itself (a "standard population
  genetic" simulation).  Thus, we supply {class}`fwdpy11.GeneticValueIsFitness`.
* If the noise object is `None`, then the model has no random effects and
  we supply {class}`fwdpy11.NoNoise`.

The remaining elements of the class are its public interface, consisting
of a mixture of final functions, abstract methods, and default implementations.
Each function takes a single argument, most of which are declared as some
sort of mocked data type.  These data types hold data that are relevant
to carrying out our calculations.  Below, the C++ versions get defined and
the Python versions are listed at the end of this section.

```{code-block} python

import abc
import typing
import attr
import numpy as np
import fwdpy11


@attr.s(auto_attribs=True)
class MockDiploidGeneticValue(abc.ABC):
    ndim: int
    genetic_value_to_fitness_map: typing.Union[None, fwdpy11.GeneticValueToFitnessMap]
    noise_function: typing.Union[None, fwdpy11.GeneticValueNoise]

    def __attrs_post_init__(self):
        self.gvalues = np.zeros(self.ndim)
        if self.genetic_value_to_fitness_map is None:
            self.genetic_value_to_fitness_map = fwdpy11.GeneticValueIsFitness(self.ndim)
        if self.noise_function is None:
            self.noise_function = fwdpy11.NoNoise()

    @typing.final
    def __call__(self, data: MockDiploidGeneticValueData) -> None:
        # Set the genetic value metadata for offspring
        data.offspring_metadata.g = self.calculate_gvalue(data)
        # Set the random effects metadata for offspring
        data.offspring_metadata.e = self.noise(MockDiploidGeneticValueNoiseData(data))
        # Set the fitness metadata for offspring
        data.offspring_metadata.w = self.genetic_value_to_fitness(
            MockDiploidGeneticValueToFitnessData(data)
        )

    @abc.abstractmethod
    def calculate_gvalue(self, data: MockDiploidGeneticValueData) -> float:
        """
        This function must populate self.gvalues
        and return a single value for the offspring's DiploidMetadata.g
        """
        raise NotImplementedError("calculate_gvalue not implemented")

    def genetic_value_to_fitness(
        self, data: MockDiploidGeneticValueToFitnessData
    ) -> float:
        return self.genetic_value_to_fitness_map(data)

    def noise(self, data: MockDiploidGeneticValueNoiseData) -> float:
        return self.noise_function(data)

    @abc.abstractmethod
    def update(self, pop: fwdpy11.DiploidPopulation) -> None:
        raise NotImplementedError("update not implemented")

```

To make things complete, a Python analogue of {class}`fwdpy11.GeneticValueToFitnessMap`
would look like this:

```{code-block} python

@attr.s(auto_attribs=True)
class MockGeneticValueToFitnessMap(abc.ABC):
    ndim: int

    @abc.abstractmethod
    def __call__(self, data: MockDiploidGeneticValueToFitnessData) -> float:
        raise NotImplementedError("__call__ not implemented")

    @abc.abstractmethod
    def update(self, pop: fwdpy11.DiploidPopulation) -> None:
        raise NotImplementedError("update not implemented")

```

Similarly, this is how we'd implement {class}`fwdpy11.GeneticValueNoise` in Python:

```{code-block} python

class MockGeneticValueNoise(abc.ABC):
    @abc.abstractmethod
    def __call__(self, data: MockDiploidGeneticValueNoiseData) -> float:
        raise NotImplementedError("__call__ not implemented")

    @abc.abstractmethod
    def update(self, pop: fwdpy11.DiploidPopulation) -> None:
        raise NotImplementedError("update not implemented")

```

(gvalue-simulation-behavior)=

### Behavior during a simulation

During a simulation, the following operations happen:

```{code-block} python

gv.update(pop)
gv.genetic_value_to_fitness_map.update(pop)
gv.noise_function.update(pop)

```

#offspring_metadata is a list of fwdpy11.DiploidMetadata 

: for o in offspring_metadata:

#Execute the calculations for each offspring 

: gv(MockDiploidGeneticValueData(o, pop))

The definition of our classes and the way that `update` functions are applied
implies that all three objects are updated *independently from one another*.
This independence covers a large number of use cases.  For example, a moving
trait optimum ({class}`fwdpy11.GaussianStabilizingSelection`) only needs acces to
{attr}`fwdpy11.DiploidPopulation.generation` to know if it needs to update its
internal state.

For some models, this independent updating doesn't do the job.  For example, if
fitness depends on an individual's genetic value *in comparison* to some or all
other genetic values in the population (as may be the case for models of
social interactions), then we don't have an easy way to
implement that using a class derived from {class}`fwdpy11.GeneticValueIsTrait`.

To handle these more complex models, subclasses may simply override the
functions `genetic_value_to_fitness` and/or `noise` as needed.
"Snowdrift" models of social interactions are one example where this method
is necessary (see {ref}`here <snowdrift-cpp-example>` for a C++ implementation
and {ref}`here <pysnowdrift>` for a Python version).

### The data objects

The member functions of our mock classes take different argument types.
These types hold things like the offspring metadata, parental metadata, etc.,
that may be needed to calculate return values.

For a genetic value object written in Python:

* `MockDiploidGeneticValueData` is {class}`fwdpy11.PyDiploidGeneticValueData`.
* `MockDiploidGeneticValueToFitnessData` is {class}`fwdpy11.DiploidGeneticValueToFitnessData`.
* `MockDiploidGeneticValueNoiseData` is {class}`fwdpy11.DiploidGeneticValueNoiseData`.

The data types closely match those used by the C++ API, which are
described in the next section.  The difference is that
{class}`fwdpy11.PyDiploidGeneticValueData` supports the buffer protocol
in order to have write access to {attr}`fwdpy11.DiploidGeneticValue.genetic_values`.
For the other two classes, their Python versions are exactly like their C++
back ends.

See {ref}`here <gvalues_python>` for details on writing Python genetic
value classes.

## The C++ side

This section shows the exact C++ implementations behind `fwdpy11` `ABCs`.
For those users implementing custom models, subclassing these types is
the way to maximize performance.  Concrete examples of subclasses are shown
below.

:::{note}

The C++ classes described here have functions that interact
with Python.  A future version of `fwdpy11` will likely
remove these functions from the C++ side and move them to
the Python-only side in order to promote thread safety.

:::

The C++ back-end for {class}`fwdpy11.DiploidGeneticValue` is found in the header file
`<fwdpy11/genetic_values/DiploidGeneticValue.hpp>`:

```{literalinclude} ../../fwdpy11/headers/fwdpy11/genetic_values/DiploidGeneticValue.hpp
:language: cpp
:lines: 19-

```

Some notes:

* Only functions marked `virtual` may be overridden by a subclass.
* A function with the following form is a "pure virtual function", which
  is the C++ analogue of {func}`abc.abstractmethod`:

  ```
  virtual void member_function() = 0;
  ```

  Attempting to allocate an instance of a subclass that has not defined
  all pure virtual functions will fail to compile.
* Classes like {class}`fwdpy11.Additive` take instances of {class}`fwdpy11.GeneticValueIsTrait`
  instead of {class}`fwdpy11.GeneticValueIsFitness`.  The reason is that most use cases
  requiring a non-`None` argument here are models of traits.  For other cases,
  the default ({class}`fwdpy11.GeneticValueIsFitness`) suffices.  However, the C++
  API allows you to subclass the base class however you wish.

The interface class for genetic value to fitness maps is `<fwdpy11/headers/fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp>`:

```{literalinclude} ../../fwdpy11/headers/fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp
:language: cpp
:lines: 19-

```

The C++ definition of {class}`fwdpy11.GeneticValueIsTrait` inherits from the above, and is another interface class (because it does not define any of the pure virtual functions):

```{literalinclude} ../../fwdpy11/headers/fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp
:language: cpp
:lines: 19-

```

The three data types `MockDiploidGeneticValueData`, `MockDiploidGeneticValueNoiseData`,
and `MockDiploidGeneticValueToFitnessData` map to the following three classes,
respectively, which hold references to objects from the simulation:

```{literalinclude} ../../fwdpy11/headers/fwdpy11/genetic_value_data/genetic_value_data.hpp
:language: cpp
:lines: 19-

```

### Example: strict additive effects

This example comes from the test suite:

```{literalinclude} ../../tests_with_cpp/custom_additive.cc
:language: cpp

```

(snowdrift-cpp-example)=

### Example: a snowdrift model

This example is based on {cite}`Doebeli2004-ny` and comes from the `fwdpy11` test suite.

```{eval-rst}
.. todo:: 

    Describe in words the implementation details.
```

The low-level details in C++ are:

```{literalinclude} ../../tests_with_cpp/ll_snowdrift.cc
:language: cpp

```

The user-facing Python class is implemented using `attrs`, which
is handy because we don't have to write C++ code to pickle/unpickle:

```{literalinclude} ../../tests_with_cpp/snowdrift.py
:language: python

```


