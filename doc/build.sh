jupyter-book build -W .
RV=$?
if [ $RV -ne 0 ];
then
	echo "Error!"
fi
exit $RV

