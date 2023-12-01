PYTHON=$1

docker run --rm -v `pwd`:/project -w /project quay.io/pypa/manylinux_2_28_x86_64:2023-11-29-1ba608e /bin/sh ./deployment/linux_wheels/wheel_workflow.sh $PYTHON
