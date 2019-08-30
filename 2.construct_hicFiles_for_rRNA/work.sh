perl sam_to_loops.pl pairTag.sam rRNA 0 20000 > rRNAloops
perl interaction_from_loop_to_JuiceBox.pl rRNA.loops > rRNA.list
cat rRNA.list | sort -k3,3 -k7,7 > rRNA.sort.list
java -jar juicebox_tools.jar pre -r 1,2,5,10,20,50 rRNA.sort.list rRNA.hic hs.size