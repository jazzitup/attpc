root -l -b <<-EOF
.L DataFrame.cc
.L PadMap.cc
.x grawToTree.cc( $1)
.q
EOF

#examle :
#   . run.sh 100  : 100 events 
#   . run.sh -1   : all events 
