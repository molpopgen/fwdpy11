PYTHON=$1

docker run --rm -v `pwd`:/project -w /project quay.io/pypa/manylinux_2_28_x86_64:2022-10-02-69a0972 /bin/sh ./deployment/linux_wheels/wheel_workflow.sh $PYTHON
