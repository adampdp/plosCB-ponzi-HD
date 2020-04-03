#!/bin/bash

# Script to manage local parameter searches.

##########################
# ### Parameters ###
##########################
numMPI=1
numTHR=25
seed=88782
changeFile="MahonNet.gsl"

variables=`echo $(seq 0.001 0.004 1)`

cb=0.2; 
synb=0.08;
deltatsim=0.2;
deltatinp=1;
size=50;
tscale=1.0;
donoise=0;
noiselev=0.00;
inpbaselow=0.04381;
inpbasehi=0.12381;
inprangelow=0.00;
inprangehi=0.00;
tlen=1060000;
var1off=0.001;
inpbufsz=65536;
inpbufoff=262144;
vk=-90.0;
vkrp=-77.5;
vkir=-90.0;


for var1 in $variables; do
cp $changeFile $changeFile.run;
sed -i 's/VKIRXX/'"$vkir"'/g' ./$changeFile.run
sed -i 's/VKRPXX/'"$vkrp"'/g' ./$changeFile.run
sed -i 's/VKXX/'"$vk"'/g' ./$changeFile.run
sed -i 's/INPBUFSZXX/'"$inpbufsz"'/g' ./$changeFile.run
sed -i 's/INPBUFOFFXX/'"$inpbufoff"'/g' ./$changeFile.run
sed -i 's/DELTATINPXX/'"$deltatinp"'/g' ./$changeFile.run
sed -i 's/DELTATSIMXX/'"$deltatsim"'/g' ./$changeFile.run
sed -i 's/CONPROBXX/'"$cb"'/g' ./$changeFile.run
sed -i 's/SIZEXX/'"$size"'/g' ./$changeFile.run
sed -i 's/TSCALEXX/'"$tscale"'/g' ./$changeFile.run
sed -i 's/NLEVXX/'"$noiselev"'/g' ./$changeFile.run
sed -i 's/TLENXX/'"$tlen"'/g' ./$changeFile.run
sed -i 's/DONOISEXX/'"$donoise"'/g' ./$changeFile.run
sed -i 's/INPBASELOWXX/'"$inpbaselow"'/g' ./$changeFile.run
sed -i 's/INPBASEHIXX/'"$inpbasehi"'/g' ./$changeFile.run
sed -i 's/INPRANGELOWXX/'"$inprangelow"'/g' ./$changeFile.run
sed -i 's/INPRANGEHIXX/'"$inprangehi"'/g' ./$changeFile.run
f2=$(expr "scale = 5; $var1" | bc | sed 's/^\./0./' );
f3=$(expr "scale = 5; $f2+$var1off" | bc | sed 's/^\./0./' );
sed -i 's/WEILOXX/'"$f2"'/g' ./$changeFile.run
sed -i 's/WEIHIXX/'"$f3"'/g' ./$changeFile.run 
sl=$(expr "scale = 5; $synb" | bc | sed 's/^\./0./' );
sh=$(expr "scale = 5; $synb+0.01" | bc | sed 's/^\./0./' );
sed -i 's/SYNBLOXX/'"$sl"'/g' ./$changeFile.run
sed -i 's/SYNBHIXX/'"$sh"'/g' ./$changeFile.run 
ulimit -Sn 4096
mpiexec -n $numMPI ../../../gsl/bin/gslparser -t $numTHR -f $changeFile.run -s $seed
rm $changeFile.run;
seed=$((seed+1))
done
