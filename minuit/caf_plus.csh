#!/bin/tcsh -f

setenv PATH /usr/bin:/bin
ln -s /usr/lib/libreadline.so.4 libreadline.so.3 
ln -s /usr/lib/libhistory.so.4 libhistory.so.3 
ln -s /usr/lib/libncurses.so.5 libncurses.so.4 
ln -s /usr/lib/libreadline.so.4 libreadline.so.4.1 
ln -s /usr/lib/libhistory.so.4 libhistory.so.4.1 
ln -s /usr/lib/libncurses.so.5 libncurses.so.4.1 

source ~cdfsoft/cdf2.cshrc
setup cdfsoft2

./toy_plus 

rm lib* toy_plus
exit 0
 
