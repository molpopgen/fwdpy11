rm -rf dist/

for py in cp37-cp37m cp38-cp38 cp39-cp39
do
    rm -rf build
    PYPATH=/opt/python/${py}
    PYBIN=${PYPATH}/bin/python 
    ${PYBIN} -m pip install --no-cache-dir -r requirements.txt
    PATH=${PYPATH}:$PATH ${PYBIN} setup.py build_ext -i
    PATH=${PYPATH}:$PATH ${PYBIN} setup.py bdist_wheel
done

cd dist
for whl in *.whl; do
    auditwheel -v repair "$whl"
    rm "$whl"
done
cd ..
