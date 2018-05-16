#!/usr/bin/perl
use Bio::SeqIO;
use Bio::SearchIO;
open(out1,'>C:\Users\wbzhe\Stentor\denovo_new.fa');
open(out2,'>C:\Users\wbzhe\Stentor\denovo_same.fa');
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => 'C:\Users\wbzhe\Stentor\Stentor_intronfree_db_CDS_v1.fas');
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();
	$seq_database{$idthis}=$result_Trans->seq();

}
my $stream_Trans = Bio::SeqIO->new(-format => 'fasta',-file   => 'C:\Users\wbzhe\Stentor\Stentor_intronfree_denovo_CDS_v1.fa');
while(my $result_Trans=$stream_Trans->next_seq()){
	$idthis=$result_Trans->id();#>NODE_50304_length_663_cov_165.749|gi|75017804|sp|Q8T277.3|PRKAG_DICDI|300..524|-
	$idthis=~/(NODE.*?)\|gi\|.*?\|sp.*\|.*?\|.*?\|(\d+)\.\.(\d+)\|/;
	$tag_denovo=$1;
	$start_denovo=$2;
	$end_denovo=$3;
	$seq_denovo{$idthis}=$result_Trans->seq();
	$ok=0;
	foreach$keys(keys%seq_database){
		$keys=~/(NODE.*?)\|gi\|.*?\|sp.*\|.*?\|.*?\|(\d+)\.\.(\d+)\|/;
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
