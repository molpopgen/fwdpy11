docker create -ti --name extract_fwdpy11_wheels fwdpy11_build_wheels
docker start extract_fwdpy11_wheels
docker exec -ti extract_fwdpy11_wheels tar cvf wheels.tar dist
docker cp extract_fwdpy11_wheels:/app/wheels.tar .
docker rm -f extract_fwdpy11_wheels
