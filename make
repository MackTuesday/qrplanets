#!/bin/zsh

compile () {
    if [[ $1.cpp -nt $1.o || ! -a $1.o ]]
    then
        echo "Compiling $1.cpp"
        g++ $2 $3 $4 $5 $6 $7 $8 $9 -c $1.cpp -o $1.o
    fi
}

compile BrentsRand -std=c++14 -g
compile ChooseMass -std=c++14 -g
compile Star -std=c++14 -g
compile GlobalConstants -std=c++14 -g
compile Planet -std=c++14 -g
compile StarSystem -std=c++14 -g
ar rvs qrplanets.a BrentsRand.o ChooseMass.o Star.o GlobalConstants.o Planet.o StarSystem.o
