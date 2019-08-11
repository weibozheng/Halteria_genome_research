#!/usr/bin/perl
use Bio::SeqIO;
use Bio::SearchIO;
$output_folder=$ARGV[2];
open(out1,'>',$output_folder . 'denovo_new.fa');#De novo predicted genes that is not involved in ncbi database
open(out2,'>',$output_folder . 'denovo_same.fa');#Genes predicted by both database-based method and augustus
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => $ARGV[0]);#genes predicted by database-based method 
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$seq_database{$idthis}=$result_Trans->seq();

}
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => $ARGV[1]);#genes predicted by de novo method
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$idthis=~/(.*?)\|gi\|.*?\|sp.*\|.*?\|.*?\|(\d+)\.\.(\d+)\|/;#This regulation expression extract 1. chromosome, 2. start, 3. end
	$tag_denovo=$1;
	$start_denovo=$2;
	$end_denovo=$3;
	$seq_denovo{$idthis}=$result_Trans->seq();
	$ok=0;
	foreach$keys(keys%seq_database){
		$keys=~/(.*?)\|gi\|.*?\|sp.*\|.*?\|.*?\|(\d+)\.\.(\d+)\|/;#ditto
		$tag_database=$1;
		$start_database=$2;
		$end_database=$3;
		if($tag_denovo eq $tag_database){
			if($start_database<=$start_denovo and $end_database>=$start_denovo){
				$ok=1;
				last;
			}elsif($start_database>=$start_denovo and $end_database<=$end_denovo){
				$ok=1;
				last;
			}elsif($start_database>$start_denovo and $start_database<$end_denovo and $end_database>=$end_denovo){
				$overlap=$end_denovo-$start_database+1;
				if($overlap/($end_denovo-$start_denovo+1)>0.4){
					$ok=1;
					last;
				}
			}elsif($start_denovo>$start_database and $start_denovo<$end_database and $end_denovo>=$end_database){
				$overlap=$end_database-$start_denovo+1;
				if($overlap/($end_denovo-$start_denovo+1)>0.4){
					$ok=1;
					last;
				}
			}
		}
	}
	if($ok==0){
		print out1 ">$idthis\n$seq_denovo{$idthis}\n";
	}else{
		print out2 ">$idthis\n$seq_denovo{$idthis}\n";
	}
}
