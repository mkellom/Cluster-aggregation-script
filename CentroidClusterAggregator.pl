#!usr/bin/perl
#All directory paths in the script will need to be customized for other users, whether specified in line or not.
use List::Util qw(shuffle);
$infile=$ARGV[0]; #Input FASTA sequence file; Important that every sequence has a unique header.
$iters=$ARGV[1]-1; #Number of clustering iterations.
$scale=$ARGV[2]; #Name of dataset.
open(FILE,"<$infile"); #Reads input FASTA file into memory
@lines=<FILE>;
close FILE;
$single=join('',@lines);
@fastas=split(/>/,$single); #Splits FASTA file into individual header/sequence strings (array values).
shift @fastas;
@mxn=grep{/>SOLEXA/} @lines; #This line is meant to grab all sequence headers in order. The string between "{/" and "/}" should be a string that is in every sequence header and not the sequence code. Will most likely need to be customized for other users.
$total=scalar @mxn; #Counts total headers (and in effect sequences).
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Creation of log file that will keep track of steps and outputs.
print LOG "Mapping inputs...$timestamp\n";
close LOG;
print "Mapping inputs...\n";
for($i=0;$i<=$total-1;$i++){ #Creates Index Hash of headers as keys and unique numerical identifiers as values.
	($sol,$junk)=split(/ /,$mxn[$i],2); #This line was used to shorten headers for our particular dataset. May not be needed or may need to be customized for other users.
	chomp $sol;
	$hash{$sol}="$i";
}
system("mkdir -m 777 /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale"); #Creates linux directory. Will need to be customized for other users.
for($n=0;$n<=$iters;$n++){ #Runs clustering for specified number of iterations.
	$timestamp=localtime(time);
	open(LOG,">>log$scale.txt"); #Prints current clustering iteration into log.
	print LOG "Clustering Run $n...$timestamp\n";
	close LOG;
	print "Clustering Run $n...\n";
	system("mkdir -m 777 /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n"); #Creates linux directory. Will need to be customized for other users.
	@fastas=shuffle(@fastas); #Shuffles (randomizes) sequence input order.
	$fasta=join('',@fastas);
	$fasta=~s/SOLEXA/>SOLEXA/g; #This line is meant to replace the ">" symbol in the sequence headers since we used it to split the FASTA file into individual header/sequence array values. Will most likely need to be customized for other users.
	open(OUTFILE,">>temp$scale.fasta"); #This is a temporary file that contains the randomized sequence input order for this particular clustering iteration.
	print OUTFILE "$fasta";
	close OUTFILE;
	system("/mnt/datA/Applications/usearch8.0.1517_i86linux64 -cluster_fast temp$scale.fasta -id 0.95 -msaout /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/cluster_ -threads 60"); #Runs USEARCH for the randomized sequence input order. Will need to be customized for other users.
	unlink("temp$scale.fasta"); #Deletes the temporary file containing the randomized sequence order for this particular clustering iteration.
	$timestamp=localtime(time);
	open(LOG,">>log$scale.txt"); #Prints into log that the concatenation step is beginning. Everything in the rest of this 'for loop' is specific to the use of USEARCH output file, possible specific to the version marked here.
	print LOG "Concatenating Run $n...$timestamp\n";
	close LOG;
	print "Concatenating Run $n...\n";
	system("sed -i '1i > ' /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/*"); #linux command to add "> " as the first line of every clustering output file for this clustering iteration. This is used to separate clusters in the master output file for this iteration when concatenated. Will need to be customized for other users.
	system("cat /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/* > /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/cluster_master"); #Linux command to concatenate all clustering output files into a master output file. Will need to be customized for other users.
	system("echo '> ' >> /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/cluster_master"); #Linux command to add "> " to the end of the master output file. Will need to be customized for other users.
	open(FILE,"</mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/cluster_master"); #Reads master output file into memory
	@lines=<FILE>;
	close FILE;
	unlink("/mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/cluster_master"); #Deletes master output file. The next regex lines alter each line of the master output; Might need to be customized for other users.
	$single=join('',@lines);
	$single=~s/[\-]\n/-/g;
	$single=~s/A\n/A/g;
	$single=~s/T\n/T/g;
	$single=~s/C\n/C/g;
	$single=~s/G\n/G/g;
	$single=~s/[\-]//g;
	$single=~s/>/\n>/g;
	$single=~s/\n//;
	$single=~s/> \n/> /g;
	open(OUTFILE,">>/mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/cluster_master"); #Writes new master output file with the regex altered lines.
	print OUTFILE "$single";
	close OUTFILE;
}
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Writes into log that the compiling of edges process has started.
print LOG "Compiling edges...$timestamp\n";
close LOG;
print "Compiling edges...\n";
system("mkdir -m 777 /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/edges"); #Creates linux directory. Will need to be customized for other users.
for($n=0;$n<=$iters;$n++){ #For the number of iterations specified, cluster output files from each iteration are opened and edges stored
	@cluster_files=</mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/$n/*>; #Grabs all output files.
	foreach $cluster_file (@cluster_files){
		next if ($cluster_file=~/master/); #This line ensures that master output files are skipped while the script iterates through all the output files.
		open(FILE,"<$cluster_file");
		@lines=<FILE>;
		close FILE;
		@lines=grep{/>/} @lines; #Grabs all the sequence headers and cluster separators in the order that they appear in the output file.
		next if (scalar @lines<=2); #This line tells the script to skip clusters if it contains a number of sequences that is less than or equal to the number given in this line. This can be customized by the user. The larger the number, the faster the script will complete but at the cost of stringency. Remember that the output files have a first line of "> " so a value of 2 in this line equates to a single sequence, or singleton cluster.
		shift @lines;
		foreach $node (@lines){
			chomp $node;
			$code=$hash{"$node"}; #For each sequence header in this cluster, the numerical identifier is looked up in th index hash.
			push(@codes,$code); #The numerical identifier from the index hash is pushed into an array for this cluster.
		}
		@codes=sort{$a<=>$b}@codes; #The numerical identifiers in this cluster are sorted from smallest to largest.
		$total_codes=scalar @codes; #Counts the size of this cluster.
		for($i=0;$i<=$total_codes-2;$i++){ #In this loop, and the one nested within it, every cluster edge combination is written to master edge lists. The smaller numerical idetifier in the edge pairing is used as the name of the edge list file (AKA Hub) and the larger numerical identifier is written into that file as a line (AKA Node).
			open(EDGE_FILE,">>/mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/edges/$codes[$i].txt");
			for($j=$i+1;$j<=$total_codes-1;$j++){
				print EDGE_FILE "$codes[$j]\n"
			}
			close EDGE_FILE;
		}
		undef @codes; #The cluster array is undefined so that it can be used for the next cluster.
	}
}
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Writes into log that the counting of edges process has started.
print LOG "Counting edges...$timestamp\n";
close LOG;
print "Counting edges...\n";
system("mkdir -m 777 /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/counts");
@codes=</mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/edges/*>; #Grabs all edge files.
$max_level=0;
$n=0;
foreach $file (@codes){ #For each edge file:
	($node1)=$file=~/(\d+)/; #Grab the Hub numerical idetifier from the edge file name.
	open(FILE,"<$file"); #Open the edge file.
	@edges=<FILE>; #Read the Node numerical identifiers into an array.
	close FILE;
	chomp @edges;
	map($counts{$_}++, @edges); #Convert the array of Node numerical identifiers into a hash. As it does this, it counts the number of times a Node numerical identifer is added to the hash and updates this count as the value for this Node key. This is counting the number of times this Hub-Node edge was stored in this edge list.
	foreach $node2 (keys %counts){ #For each Node key:
		$level=$counts{$node2}; #Lookup the count of each edge appearance
		if($level>=$max_level){ #This small loop is used to keep track of the most times any edge was found.
			$max_level=$level;
		}
		open(OUTFILE,">>/mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/counts/$level.txt"); #Each Hub-Node edge is written into a file with their count as the filename. This keeps track of the number of times each edge occured throught all clustering iterations.
		print OUTFILE "$node1 $node2\n";
		close OUTFILE;
		$n++; #Used to keep track of the number of edges total.
	}
	undef %counts
}
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Writes into log the number of edges processed.
print LOG "Processed $n edges total...$timestamp\n";
close LOG;
print "Processed $n edges total...\n";
open(LOG,">>log$scale.txt"); #Writes into log that the cluster aggregation process has started.
print LOG "Final clustering...$timestamp\n";
close LOG;
print "Final clustering...\n";
%rhash=reverse %hash; #Created a Reverse Index Hash (inverse of the index hash).
@levels=</mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/counts/*>; #Grabs all the files containing the edges and counts as filenames.
$c=0;
for($level=$max_level;$level>=2;$level--){ #Starting with the most frequently found edges (in the file with the max_level filename:
	open(FILE,"</mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/counts/$level.txt") or next; #If the file exists, open it
	@edges=<FILE>;
	close FILE;
	$timestamp=localtime(time);
	open(LOG,">>log$scale.txt"); #Writes into log that the aggregating process has started for this counts file.
	print LOG "Level: $level $timestamp\n";
	close LOG;
	print "Level: $level\n";
	chomp @edges;
	@edges=shuffle(@edges); #Randomized the edges in the current file so that the largest Node numerical identifier isn't always the first edge to be built into the aggregate clusters.
	foreach $edge (@edges){ #For each edge in the file:
		($hub,$node)=split(/ /,$edge,2); #Split the Hub and Node into separate strings.
		next if ((exists $antifinals{$hub}) and (exists $antifinals{$node})); #1)	Skip to the next edge if both the Hub and Node have already been assigned to clusters in the Tracking Hash.
		if (exists $antifinals{$hub}){ #2)	If the Hub has already been assigned to a cluster in the Tracking Hash (implying with step 1 that the Node has not been assigned yet):
			$avalue=$antifinals{$hub}; #a)	Get the Cluster Numerical Identifier value that the Hub Numerical Identifier key has been assigned to in the Tracking Hash.
			$finals{$avalue}.=$rhash{$node}; #b)	Get the sequence header value for the Node Numerical Identifier key from the Inverse Index Hash and append it to the value for the Cluster Numerical Identifier (from step 2a) key in the Aggregate Cluster Hash.
			$antifinals{$node}=$avalue; #c)	Append this Node Numerical Identifier key - Cluster Numerical Identifier value pair to the Tracking Hash
		}elsif (exists $antifinals{$node}){ #3)	If the Node has already been assigned to a cluster in the Tracking Hash (implying with step 1 that the Hub has not been assigned yet):
			$avalue=$antifinals{$node}; #a)	Get the Cluster Numerical Identifier value that the Node Numerical Identifier key has been assigned to in the Tracking Hash.
			$finals{$avalue}.=$rhash{$hub}; #b)	Get the Sequence Header Value for the Hub Numerical Identifier key from the Inverse Index Hash and append it to the value for the Cluster Numerical Identifier (from step 3a) key in the Aggregate Cluster Hash.
			$antifinals{$hub}=$avalue; #c)	Append this Hub Numerical Identifier key - Cluster Numerical Identifier value pair to the Tracking Hash.
		}else{ #4)	If neither the Hub nor Node have been previously assigned to a cluster in the Tracking Hash:
			$finals{$c}="$rhash{$hub}$rhash{$node}"; #a)	Create an Aggregate Cluster Hash pair with a Cluster Numerical Identifier as the key and the sequence headers for the Hub and Node Numerical Identifiers from the Inverse Index Hash as the value.
			$antifinals{$hub}=$c; #b)	Append the Hub Numerical Identifier key - Cluster Numerical Identifier value to the Tracking Hash.
			$antifinals{$node}=$c; #c)	Append the Node Numerical Identifier key - Cluster Numerical Identifier value to the Tracking Hash.
			$c++; #d)	Assign the next Cluster Numerical Identifier to be +1 greater than the current one (to create a new cluster).
		}
	}
}
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Writes into log that the cluster aggregation process has completed.
print LOG "Writing final clusters...$timestamp\n";
close LOG;
print "Writing final clusters...\n";
system("mkdir -m 777 /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/final");
for($key=0;$key<=$c-1;$key++){ #This loop prints out the aggragate cluster files.
	$value=$finals{$key};
	$value=~s/>/\n>/g;
	open(OUTFILE,">>/mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/final/cluster_$key");
	print OUTFILE "$value\n";
	close OUTFILE;
	$n=$key;
}
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Writes into log that the aggrageted cluster files are being concatenated into a master file. This is a similar process to that of the master files for the clusering outputs for the individual iterations.
print LOG "Concatenating final clusters...$timestamp\n";
close LOG;
print "Concatenating final clusters...\n";
system("sed -i '1i > ' /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/final/*");
system("cat /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/final/* > /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/final/cluster_master");
system("echo '> ' >> /mnt/datA/matt/landsort_deep/ftp_jcvi_org/clusters$scale/final/cluster_master");
$timestamp=localtime(time);
open(LOG,">>log$scale.txt"); #Writes into log that the script is finished.
print LOG "Done...$timestamp\n";
close LOG;
print "Done\n";
