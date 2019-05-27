use Bio::SeqIO;
use List::Util qw /sum/;
open(FH,'/N/dc2/projects/ciliate/Halteria_grandinella/Halteria_trans/Hal_teloremoved.sort.sam'); #sorted samtools output file
open(OUTFILE1,'>/N/dc2/projects/ciliate/Halteria_grandinella/Halteria_trans/Hal_teloremoved_RNA.coverage');#Temporary file 1
open(OUTFILE2,'>/N/dc2/projects/ciliate/Halteria_grandinella/Halteria_trans/Hal_teloremoved_RNA.coverage.intergrated'); #Temporary file 2
@data=<FH>;
close(FH);
$stream_querys = Bio::SeqIO->new(-file => '/N/dc2/projects/ciliate/Halteria_grandinella/Halteria_Nuohe/Hal_SPCA/Hal_SPCA_v3_2telo_teloremoved_strand_modified.txt', -format => 'fasta');#all 2-telo chromosomes
while ($seqa=$stream_querys->next_seq){
	$seqname=$seqa->id();
	$length=$seqa->length();
	$length_s{$seqname}=$length;
}
foreach my $key1 (keys %length_s){
	for ($i=1;$i<=$length_s{$key1};$i++){
		$hash_all{$key1}{$i}=0;
	}
}
foreach$keys(@data){ #
	@this=undef;
	$keys=~/(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t/;
	$tag=$3;
	$start=$4;
	$sim=$6;
	@this=split(/\D/,$sim);
	$sum_D=sum @this;
	for($i=$start;$i<=$start+$sum_D;$i++){  
		$hash_all{$tag}{$i}++;
	}

}
foreach my $key1 (keys %hash_all){
	my $hash2 = $hash_all{$key1};
	foreach my $key2 (sort{$a<=>$b} keys %$hash2){
		$spe_hash{$key2}=$spe_hash{$key2}+$hash2->{$key2};
		print OUTFILE1 "$key1\t$key2\t$hash2->{$key2}\n";
	}
}
foreach$keyss(sort{$a<=>$b}keys%spe_hash){
	print OUTFILE2 "$keyss\t$spe_hash{$keyss}\n";
}

open(INF,'/N/dc2/projects/ciliate/Halteria_grandinella/Halteria_trans/Hal_teloremoved_RNA.coverage');#path of Temporary file 1, tab-separated file, format: tag_name\tnucleotide_position\tsequencing_depth
open (OUTFILE3,'>/N/dc2/projects/ciliate/Halteria_grandinella/Halteria_trans/Hal_SPCA_v3_2telo_teloremoved_strand_modified_RNA_reads_start_end_basedon_23_1.ext'); #Final Output
@data=<INF>;
foreach$keys(@data){
	$keys=~/(.*?)\|.*?\t(.*?)\t(.*)/;
	$hash{$1}{$2}=$3;
}
$stream_querys = Bio::SeqIO->new(-file => 'C:\Users\wbzheng\Halteria_Nuohe\Hal_SPCA\subtelomere\based_on_all_CDS\Hal_SPCA_v3_2telo_teloremoved_modified_allCDS.txt', -format => 'fasta');#all single gene chromosomes - only sense strand
while ($seqa=$stream_querys->next_seq){
	$seqname=$seqa->id();
	$seqname=~/(NODE.*?)\|/;
	$lengthx=$seqa->length();
	$length_s{$1}=$lengthx;
}
foreach$keyx(keys%length_s){
	chomp($keyx);
	if($keyx=~/(NODE.*)/){
		$tag=$1;
		$hash2=$hash{$tag};
		$out_5_pos=0;
		$out_3_pos=0;
		foreach my $j (sort{$a<=>$b}keys%$hash2){
			if($hash2->{$j}==0){
				$out_5_pos=$j;
			}else{
				last;
			}
		}

		foreach my $k (sort{$b<=>$a}keys%$hash2){
			if($hash2->{$k}!=0){
				$out_3_pos=$length_s{$tag}-$k;
				if($out_3_pos<0){
					$out_3_pos=0;
				}
				last;
			}
		}
		if($out_5_pos == $length_s{$tag}){
			$out_5_pos='NA';
			$out_3_pos='NA';
		}
		if($out_5_pos != $length_s{$tag} and $out_5_pos > 400){
			$out_5_pos='NA';
		}
		if($out_5_pos != $length_s{$tag} and $out_3_pos > 400){
			$out_3_pos='NA';
		}
	}
	print OUTFILE3 "$keyx\t$out_5_pos\t$out_3_pos\n";
}
