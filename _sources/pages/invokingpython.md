# Invoking Python

fwdpy11 is under rapid development.  As such, API features may change from version to version. Where possible, we will
issue deprecation warnings and point you to the new API. However, Python suppresses deprecation warnings by default. To
see deprecation warnings, invoke python like this:

```{code-block} bash

python -Wd

```

To treat all warning as errors:

```{code-block} bash

python -We

```
Our opinion is that is useful to use each of these options at different times.
Treating warnings as errors forces you to update your code to what upstream developers think are best practices for their packages.
It is also good to be able to see the deprecation warnings (`-Wd`) because they warn you that your current code will eventually break.
