ls dist/wheelhouse
pwd
for py in cp36-cp36m cp37-cp37m cp38-cp38
do
    PYPATH=/opt/python/${py}
    PYBIN=${PYPATH}/bin/python 
    ${PYBIN} -m venv test-${py}
    source test-${py}/bin/activate
    which pip
    pip install --upgrade --no-cache pip wheel
    pip install --no-cache $(ls dist/wheelhouse/*-${py}-manylinux2014_x86_64.whl)
    pip install --no-cache msprime pytest pytest-xdist
    ${PYBIN} -m pytest tests -n 4
    deactivate
    rm -rf test-${py}
done
