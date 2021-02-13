(genetics_vignettes_intro)=

# Introduction

Simulation parameters are stored in instances of {class}`fwdpy11.ModelParams`.
Once initialized, instances of this class are immutable, or "frozen".

This class can only be initialized via explicit use of `kwargs`.
In many of the vignettes, you will see that a {class}`dict` is filled first, and then used to construct a `ModelParams` object with the `**` syntax.
The reason to use a `dict` is because of the immutability of the `ModelParams` instances.
Sometimes, it is convenient to add to existing lists of things in the `dict`.
Such modifications are not allowed on an initialized `ModelParams`.

Many of the vignettes that follow have a `note` section at the top.
These notes tell you what `kwarg` for the `ModelParams` initializer is being discussed.
