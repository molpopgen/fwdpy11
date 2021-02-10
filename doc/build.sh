jupyter-book build -n .
RV=$?
if [ $RV -ne 0 ];
then
	echo "Error!"
fi
exit $RV

