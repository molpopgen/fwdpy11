ls dist/wheelhouse
rm -rf test_wheels
# Have to do this from a separate dir
# so that the source repo doesn't interfere.
# The repo gets seen as the fwdpy11 package.
mkdir test_wheels
cd test_wheels
# Copy tests over b/c some of the test rely on paths
cp -r ../tests .
pwd
for py in cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
    PYPATH=/opt/python/${py}
    PYBIN=${PYPATH}/bin/python 
    ${PYBIN} -m venv test-${py}
    source test-${py}/bin/activate
    which pip
    pip install --upgrade --no-cache-dir pip wheel
    pip install fwdpy11 --pre --only-binary fwdpy11 -f ../dist/wheelhouse
    pip install --no-cache-dir msprime pytest pytest-xdist
    python -m pytest tests -n 4
    deactivate
    rm -rf test-${py}
done
cd ..
rm -rf test_wheels
