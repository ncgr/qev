

cat $1 | awk -F"\t" '{if($0~ />/){a=$0} if(a~/^>'$2'$/){print $0}}'
