def setup_mako(cfg):
    # import fwdpy11 so we can find its C++ headers
    import fwdpy11 as fp11
    import subprocess
    import sys

    gsl_prefix = (
        subprocess.check_output(["gsl-config", "--prefix"]).strip().decode("utf-8")
    )
    # add fwdpy11 header locations to the include path
    cfg["include_dirs"].extend(
        [fp11.get_includes(), fp11.get_fwdpp_includes(), gsl_prefix + str("/include")]
    )
    gsl_libs = subprocess.check_output(["gsl-config", "--libs"]).strip().decode("utf-8")
    gsl_libs = gsl_libs.split(" ")
    gsl_runtime_libs = [i.replace("-l", "") for i in gsl_libs[1:]]
    cfg["libraries"].extend(gsl_runtime_libs)
    cfg["extra_link_args"].extend([gsl_libs[0]])
    # On OS X using clang, there is more work to do.  Using gcc on OS X
    # gets rid of these requirements. The specifics sadly depend on how
    # you initially built fwdpy11, and what is below assumes you used
    # the provided setup.py + OS X + clang:
    if sys.platform == "darwin":
        cfg["compiler_args"].extend(["-stdlib=libc++", "-mmacosx-version-min=10.7"])
        cfg["linker_args"] = ["-stdlib=libc++", "-mmacosx-version-min=10.7"]
