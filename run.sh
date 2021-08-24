#for (( j=1 ; j<=1 ; j++))
for (( j=1 ; j<=9 ; j++))
do
root -l -b <<-EOF
.L GETAnalyzer.cc
.L GETDecoder.cc
.L GETPad.cc 
.x grawToTree.cc($j)
.q
EOF
done
#examle :
#   . run.sh 100  : 100 events 
#   . run.sh -1   : all events 
