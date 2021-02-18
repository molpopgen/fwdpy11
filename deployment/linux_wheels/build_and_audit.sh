for py in cp36-cp36m cp37-cp37m cp38-cp38
do
    PYPATH=/opt/python/${py}
    PYBIN=${PYPATH}/bin/python 
    ${PYBIN} -m pip install --no-cache -r requirements.txt
    PATH=${PYPATH}:$PATH ${PYBIN} setup.py build_ext -i
    PATH=${PYPATH}:$PATH ${PYBIN} setup.py bdist_wheel
done

cd dist
for whl in *.whl; do
    auditwheel -v repair "$whl"
    rm "$whl"
done
cd ..
