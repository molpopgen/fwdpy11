REPORTDIR=_build/html/reports
jupyter-book build -n -W .
RV=$?
if [ $RV -ne 0 ];
then
	echo "Error!"
    if [ -e $REPORTDIR ]; then
        cat $REPORTDIR/*
    fi
else
    rm -rf $REPORTDIR/*
fi
exit $RV

