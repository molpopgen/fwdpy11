(sec_deployment)=

# Deployment tools

:::{note}
`Docker` commands require `sudo` on most systems.
On trusted machines, you may wish to add yourself to the `docker` group as described [here](https://docs.docker.com/engine/install/linux-postinstall/).

In general, [podman](https://podman.io) will prove easier to use than `Docker`.
`podman` is standard on systems like Fedora Linux.
All commands shown below work with `podman` because that tool can alias the `Docker` commands.
On Fedora, this aliasing requires `podman-docker` to be installed. (`sudo dnf install podman-docker`).

:::

## Docker tags and fwdpy11 versions

The policy is:

* The `:latest` tag corresponds to the latest commit to the `main` branch.
  In general, this will be the same as the latest release as most development happens in a separate branch.
* Beginning with the 0.19.x series, we also have tags for specific versions/releases.
  For example, `:v0.19.0` will correspond to that `fwdpy11` version 0.19.0.
  In general, this will be the same as `:latest` unless we are prepping for a release and/or releasing
  bug fixes for folks to test.
* The versioned tags *do not* distinguish stable releases from pre-releases!
  For example, the latest release could be `v0.19.0` and the latest pre-release
  is `v0.20.0a0` (alpha release 0 of `v0.20.0`).
  To guarantee a stable release, peruse the [tags](https://hub.docker.com/r/molpopgen/fwdpy11/tags) at `Docker` hub.

## Singularity images

HPC systems are moving towards containerization, with more institutional clusters supporting [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) images.

These images are simple to build.

To bootstrap a `singularity` image from a `fwdpy11` `Docker` image and install `msprime`, the following file will suffice:

```
Bootstrap: docker
From: molpopgen/fwdpy11:latest

%post
    . /venv/bin/activate
    python -m pip install msprime
```

If we call the above file `example.def`, then we can build an image:

```{code-block} bash
singularity build --fakeroot example.sif example.def
```

You would execute the previous command on a Linux machine of your own.
For users of systems like macOS, you could run the Linux installation in a virtual machine.

Once built, you can transfer the image to your cluster to execute jobs.

To run a script in the image, imagine we have a simulation defined in `example.py`:

```{code-block} bash
singularity exec example.sif bash -c ". /venv/bin/activate; python example.py"
```

## Docker images

### Images from Docker Hub

We distribute minimal environments based on Ubuntu LTS and Python 3 with `fwdpy11` installed.

To obtain an image based on the latest commits merged into `fwdpy11/main`:

```sh
docker pull molpopgen/fwdpy11
```

The above is equivalent to:

```sh
docker pull molpopgen/fwdpy11:latest
```

This tag corresponds to the latest commit to the main branch.
See above for details about `fwdpy11` versions and Docker tags.

### The virtual environment

In the image, `fwdpy11` is installed into a virtual environment located at `/venv`.

#### Extending existing images

Arguably the best use for the existing images is to create your own reproducible environments from them.
For example, the following `Dockerfile` creates an image with `matplotlib` added to the Python environment:

```
FROM molpopgen/fwdpy11:latest

WORKDIR /app

RUN . /venv/bin/activate
  && python -m pip install matplotlib
```

For maximal reproducibility, you may want to "pin" your own images to specific versions of the `fwdpy11` image:


```
FROM molpopgen/fwdpy11@sha256:106b32cb879f6f42bd5fc34a4a44ac40371dc66f51a89ee04dbf33b906dbcf69

WORKDIR /app

RUN . /venv/bin/activate
  && python -m pip install matplotlib
```
    
To sleuth out the `sha256` value that you need, simply pull the image:

```sh
docker pull molpopgen/fwdpy11:latest
```

The output will look something like:

```
latest: Pulling from molpopgen/fwdpy11
83ee3a23efb7: Already exists 
db98fc6f11f0: Already exists 
f611acd52c6c: Already exists 
1902ff943d3d: Pull complete 
2f1378e5b22a: Pull complete 
bfbf9abf860e: Pull complete 
55f3583009c4: Pull complete 
Digest: sha256:106b32cb879f6f42bd5fc34a4a44ac40371dc66f51a89ee04dbf33b906dbcf69
Status: Downloaded newer image for molpopgen/fwdpy11:latest
docker.io/molpopgen/fwdpy11:latest
```

The `sha256` info that you need is on the line starting with `Digest:`.

(sec_building_docker_images)=

### Building images from scratch

The directory `deployment/docker` contains a `Dockerfile` that will install a minimal environment for running `fwpdy11`.

From the `root` of the `fwdpy11` repository:

```sh
docker build --rm . -f deployment/docker/Dockerfile -t IMAGE_NAME
```

The result with be a Docker image named `IMAGE_NAME`.

Once this image is built, you may use it to run scripts on your local machine.
If you have a script `$HOME/tmp/foo.py`, you can run it via:

```sh
docker run --rm -v $HOME/tmp:/app IMAGE_NAME /bin/bash -c ". /venv/bin/activate; python /app/foo.py"
```

## Linux wheel building

The directory `deployment/linux_wheels` defines a `Docker` work flow to build binary wheels on Linux.
These scripts are what we use to build wheels for each release.
They may also be used to make wheels for development versions on a local machine.

From the `root` of the `fwdpy11` repository.

```sh
bash deployment/linux_wheels/run_wheel_workflow.sh VERSION
```

`VERSION` must be a valid Python version present in the docker image.
For example, `python3.10`.

The above command will require `sudo` on some systems.
On trusted machines, you may wish to add yourself to the `docker` group as described [here](https://docs.docker.com/engine/install/linux-postinstall/).

This work flow builds binary wheels for Python [manylinux](https://github.com/pypa/manylinux).

The wheels will be found in `dist/wheelhouse` upon completion.

**Do not** run this work flow in parallel on the same filesystem.

### Extracting the wheels

```sh
bash deployment/linux_wheels/extract_wheels.sh
```

The output will be a file called `wheels.tar`.
