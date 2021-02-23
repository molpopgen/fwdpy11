(sec_deployment)=

# Deployment tools

:::{note}
`Docker` commands require `sudo` on most systems.
On trusted machines, you may wish to add yourself to the `docker` group as described [here](https://docs.docker.com/engine/install/linux-postinstall/).
:::

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

To get an image with the latest stable release installed:

```sh
docker pull molpopgen/fwdpy11:stable
```

Most people, most of the time, will want to use the `:stable` tag.
The `:latest` tag is closer to the bleeding edge, containing features added since the latest release and probably a bug or two that we haven't sorted out yet!
If you need features in the `dev` branch that are not yet merged into `main`, you can build your own images based on the `dev` branch.
See {ref}`here <sec_building_docker_images>`.

#### Extending existing images

Arguably the best use for the existing images is to create your own reproducible environments from them.
For example, the following `Dockerfile` creates an image with `matplotlib` added to the Python environment:

```
FROM molpopgen/fwdpy11:latest

WORKDIR /app

RUN python3 -m pip install matplotlib
```

For maximal reproducibility, you may want to "pin" your own images to specific versions of the `fwdpy11` image:


```
FROM molpopgen/fwdpy11@sha256:106b32cb879f6f42bd5fc34a4a44ac40371dc66f51a89ee04dbf33b906dbcf69

WORKDIR /app

RUN python3 -m pip install matplotlib
```
    
To sleuth out the `sha256` value that you need, navigate [here](https://hub.docker.com/r/molpopgen/fwdpy11/tags), click on the desired `tag`, then copy the `sha256` information from up at the top near the word `DIGEST`.
This `sha256` will give you the unique identifier corresponding to the latest build of the `tag` you are looking at.

Some notes:

* The `sha256` will differ for the most recent `:latest` and `:stable` `tags`, even if they correspond to the exact same version of `fwdpy11`.
  This difference is because the `sha256` correspond to a build of an image, and not to the version of `fwdpy11`.

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
docker run --rm -v $HOME/tmp:/app IMAGE_NAME python3 /app/foo.py
```

## Linux wheel building

The directory `deployment/linux_wheels` defines a `Docker` work flow to build binary wheels on Linux.
Although this work flow *could* be used to make wheels for each new release, its main purposes are for testing and easily getting development versions to collaborators.

From the `root` of the `fwdpy11` repository:

```sh
bash deployment/linux_wheels/build_docker_image.sh
```

The above command will require `sudo` on some systems.
On trusted machines, you may wish to add yourself to the `docker` group as described [here](https://docs.docker.com/engine/install/linux-postinstall/).

This work flow builds binary wheels for Python 3.6, 3.7 and 3.8 using [manylinux](https://github.com/pypa/manylinux) 2014.
After building and auditing each wheel, the wheels are installed into virtual environments and the `fwdpy11` test suite is run.

### Extracting the wheels

```sh
bash deployment/linux_wheels/extract_wheels.sh
```

The output will be a file called `wheels.tar`.

### Caveats and limitations

The methods used to build wheels may not be totally ready for prime time.
Ideally, we'd use [oldest-supported-numpy](https://pypi.org/project/oldest-supported-numpy/), but we haven't tested it yet.


