import sys
open(INFX,sys.argv[1]);#Genome annotation file
@data_gff=<INFX>;
foreach$keys(@data_gff){
	$keys=~/(.*?)\t.*?\t(.*?)\t(.*?)\t(.*?)\t.*?\t.*?\t.*?\t/;#chr1	.	exon	1065	1313	.	-	.	ID=NODE_206_length_95276_cov_226.532|gi|NOHIT|sp|NOHIT|NOHIT|533..1591|-
	$tag=$1;
	$type=$2;
	$start=$3;
	$end=$4;
	$hash{$tag}{$start}=$keys;
}
open(INF1,sys.argv[2]);
open(OUT,sys.argv[3]);
print OUT "chromosome\tposition\tminor_allele\tMAF\tmajor_allele\tMAAF\tposition_type\tpoly_type\treads_num\tgene\tPi\n";
@data=<INF1>;
$cutoff=0.02;
foreach$keys(@data){
	chomp($keys);
	if($keys=~/>(N.*)/){
		$current_tag=$1;
		#print OUT1 "$keys\n";
	}else{
		$keys=~/(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)/;
		$pos=$1;
		$A_num=$2;
		$C_num=$3;
		$G_num=$4;
		$T_num=$5;
		$sum=$A_num+$C_num+$G_num+$T_num;
		$freq_A=$A_num/$sum;
		$freq_C=$C_num/$sum;
		$freq_G=$G_num/$sum;
		$freq_T=$T_num/$sum;
		$good=0;
		$ef=0;
		$MAF=1;
		$MMAF=0;
		$type_this='inter_genic';
		$gene_id='NULL';
		$hash2=$hash{$current_tag};
		foreach$key2(keys%$hash2){
			$hash{$current_tag}{$key2}=~/(.*?)\t.*?\t(.*?)\t(.*?)\t(.*?)\t.*?\t.*?\t.*?\t/;#NODE_206_length_95276_cov_226.532	.	exon	1065	1313	.	-	.	ID=NODE_206_length_95276_cov_226.532|gi|NOHIT|sp|NOHIT|NOHIT|533..1591|-
			$tag=$1;
			$type=$2;
			$start=$3;
			$end=$4;
			if($current_tag eq $tag and $start<=$pos and $end>=$pos){
				if($type eq 'exon' or $type eq 'CDS'){
					$type_this='exon';
					$hash{$current_tag}{$key2}=~/ID=(.*)/;
					$gene_id=$1;
				}elsif($type eq 'intron'){
					$type_this='intron';
					$gene_id="NULL"
				}else{
					$gene_id="NULL"
				}
				last;
			}
		}
		if($freq_A*$sum>3){
			$ef++;
			if($freq_A<$MAF){
				$MAF=$freq_A;
				$minor='A';
			}
			if($freq_A>$MMAF){
				$MMAF=$freq_A;
				$major='A';
			}
		}
		if($freq_C*$sum>3){
			$ef++;
			if($freq_C<$MAF){
				$MAF=$freq_C;
				$minor='C';
			}
			if($freq_C>$MMAF){
				$MMAF=$freq_C;
				$major='C';
			}
		}
		if($freq_G*$sum>3){
			$ef++;
			if($freq_G<$MAF){
				$MAF=$freq_G;
				$minor='G';
			}
			if($freq_G>$MMAF){
				$MMAF=$freq_G;
				$major='G';
			}
			
		}
		if($freq_T*$sum>3){
			$ef++;
			if($freq_T<$MAF){
				$MAF=$freq_T;
				$minor='T';
			}
			if($freq_T>$MMAF){
				$MMAF=$freq_T;
				$major='T';
			}
		}
		
		$count++;
		if($ef==2){
			$type_sites='bi';
		}
		if($ef==3){
			$type_sites='tri';
		}
		if($ef==4){
			$type_sites='tetra';
		}
		if($ef>1){
			$pi_v=2*$MAF*(1-$MAF);
			print OUT "$current_tag\t$pos\t$minor\t$MAF\t$major\t$MMAF\t$type_this\t$type_sites\t$sum\t$gene_id\t$pi_v\n";
		}
	}
}
