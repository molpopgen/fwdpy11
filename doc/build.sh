REPORTDIR=_build/html/reports
jupyter-book build -W -n .
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

