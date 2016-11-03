#!/bin/bash

#threelL=/home/grining/AA/geminal_slater/program/other/ThreeElectron/threel.00.x
#threelP=/home/grining/AA/geminal_slater/program/other/Puch/threel
threelQ=/home/grining/Robota/Grining/Quch/threel

case $1 in
"L")
    $threelL $2 > /dev/null
    case $3 in
    "J")
	mv fileJ0.F $4
	;;
    "K")
	mv fileK0.F $4
	mv fileK1.F $5
	;;
    *) 
	echo "Unknown integral class!!!"
	;;
    esac
    ;;
"P")
    $threelP $2 
    ;;
"Q")
    $threelQ $2 
    ;;
*)
    echo "Unknown three-electron program!!!"
    ;;
esac
