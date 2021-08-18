for (( j=1 ; j<=9 ; j++))
do
root -l -b <<-EOF
.L DataFrame.cc
.L PadMap.cc
.x grawToTree.cc( $1, $j)
.q
EOF
done
#examle :
#   . run.sh 100  : 100 events 
#   . run.sh -1   : all events 
