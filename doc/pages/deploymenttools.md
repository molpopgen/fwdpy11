# Deployment tools

:::{note}
`Docker` commands require `sudo` on most systems.
On trusted machines, you may wish to add yourself to the `docker` group as described [here](https://docs.docker.com/engine/install/linux-postinstall/).
:::

## Docker images

The directory `deployment/docker` contains a `Dockerfile` that will install a minimal environment for running `fwpdy11`.
In the long run, this file will be the basis for auto-deployed Docker images.

From the `root` of the `fwdpy11` repository:

```sh
docker build --rm IMAGE_NAME -f deployment/docker/Dockerfile .
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


