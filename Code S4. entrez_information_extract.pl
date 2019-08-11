use LWP::Simple;
open(FH1,$ARGV[0]);#Gi and expression level of all genes
open(out,'>',$ARGV[1]);#Tab-like output file
@datassssss=<FH1>;


$db = 'protein';
print out "Protein_GI\tLogFC\tTitle\tCatalytic activity\tFunction\tSubcellular location\tMiscellaneous\tSimilarity\tTaxname\n";
$count=0;
foreach$keys(@datassssss){
$keys=~/(.*?)\t(.*?)\t(.*?)\n/;
$tag=$1;
$id_list = $2;

$exp=$3;
$count++;
print out "$tag\t$id_list\t$exp\t";

#print"$count\t$id_list\t$exp\t";
# Download protein records corresponding to a list of GI numbers.

#assemble the epost URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "epost.fcgi?db=$db&id=$id_list";

#post the epost URL
$output = get($url);

#parse WebEnv and QueryKey
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);



### include this code for EPost-EFetch
#assemble the efetch URL
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
$url .= "&rettype=COMMENT&retmode=text";

#post the efetch URL
$data = get($url);
if($data=~/title "RecName: (.*?)"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\t";
if($data=~/comment "\[CATALYTIC ACTIVITY\](.*?)"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\t";
if($data=~/comment "\[FUNCTION\](.*?)"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\t";
if($data=~/comment "\[SUBCELLULAR LOCATION\](.*?)"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\t";
if($data=~/comment "\[MISCELLANEOUS\](.*?)"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\t";
if($data=~/comment "\[SIMILARITY\](.*?)"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\t";
if($data=~/taxname \"(.*?)\"/ms){
	$thisout=$1;
	$thisout=~s/\n//g;
	print out "$thisout";
}
print out "\n";

}
