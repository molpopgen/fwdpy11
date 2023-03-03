---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(demography_debugger_vignette)=

# Debugging demographic models

```{note}
In fwdpy11 0.20.0 the DemographyDebugger got a lot less useful.
The reason is that you cannot have an invalide model because
the demes back end will catch it.

In the future, we will add more features back to this type or
remove it altogether if we can't think of anything useful to 
do.
```
