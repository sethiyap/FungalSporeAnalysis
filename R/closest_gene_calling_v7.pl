use strict;
use warnings;
use LWP::Simple;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use URI::Escape;

use Pod::Usage;

###############  Start Time  ########################################################
my $Time_Start="";
my $TS=time();
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";

##################  Main     ######################################################
use threads;
our $g_opts;
sub parse_opts{
    my $result = GetOptions(
                    "gff=s" => \$g_opts->{'gff'},#string
                    

                    "peak_file|pf=s" => \$g_opts->{'peak_file'},#string
                    "summit_column|sc=i" => \$g_opts->{'summit_column'},#string
                    
                    "calling_gene_by_closest|cgbc" => \$g_opts->{'calling_gene_by_closest'},#string
                    #"calling_gene_by_distance|cgbd" => \$g_opts->{'calling_gene_by_distance'},#string
                    #"calling_gene_by_both|cgbb" => \$g_opts->{'calling_gene_by_both'},#string
                    
                    "abcd" => \$g_opts->{'abcd'},#string
                    
                    "chromosome_name_file|cnf=s" => \$g_opts->{'chromosome_name_file'},#string
                    
                    "thread_limit|tl=s" => \$g_opts->{'thread_limit'},#string
                    "distance|d=s"  => \$g_opts->{'distance'},#string
                    "common_name|cn" => \$g_opts->{'common_name'},#string
                    
                    "cds" => \$g_opts->{'cds'},#string
                    
					"output_folder|of=s"=> \$g_opts->{'output_folder'},
					"outout_name|on=s"=> \$g_opts->{'outout_name'},
					
                    "quiet"     => sub { $g_opts->{'verbose'} = 0 },
                    "help|?"    => \$g_opts->{'help'}
                  );
    if(!($g_opts->{'gff'}) or !($g_opts->{'peak_file'})){
        &pod2usage( -verbose => 1);#exit status will be 1
    }
    if($g_opts->{'help'}){
        &pod2usage( -verbose => 1);#exit status will be 1
    }
}

&parse_opts();

my $gff=$g_opts->{'gff'};
my $chromosome_name_file=0;
$chromosome_name_file=$g_opts->{'chromosome_name_file'} if($g_opts->{'chromosome_name_file'});
my $peak_file=$g_opts->{'peak_file'};




my $summit_column=-1;
$summit_column=$g_opts->{'summit_column'}-1 if($g_opts->{'summit_column'});

my $cds=0;
$cds=1 if($g_opts->{'cds'});

my $calling_gene_by_closest=0;
$calling_gene_by_closest=1 if($g_opts->{'calling_gene_by_closest'});

my $calling_gene_by_both=0;
$calling_gene_by_both=1 if($g_opts->{'calling_gene_by_both'});

my $calling_gene_by_distance=0;
$calling_gene_by_distance=1 if($g_opts->{'calling_gene_by_distance'});


my $distance=2000;
$distance=$g_opts->{'distance'} if($g_opts->{'distance'});

my $common_name=0;
$common_name=1 if($g_opts->{'common_name'});

my $abcd=0;
$abcd=1 if($g_opts->{'abcd'});

my $thread_limit=2; #default thread_limit
$thread_limit=$g_opts->{'thread_limit'} if($g_opts->{'thread_limit'});


my $output_folder=dirname($peak_file);
my $outout_name=basename($peak_file);

$output_folder=$g_opts->{'output_folder'} if(($g_opts->{'output_folder'}));
$output_folder=~s/\/$//;

$outout_name=$g_opts->{'outout_name'} if(($g_opts->{'outout_name'}));

print "Parameter:\n\n";
print "gff: $gff\n";
print "peak_file: $peak_file\n";
print "summit_column: ".($summit_column+1)."\n" if(($g_opts->{'summit_column'}));
print "cds: True\n" if(($g_opts->{'cds'}));
print "chromosome_name_file: $chromosome_name_file\n" if(($g_opts->{'chromosome_name_file'}));
print "calling_gene_by_closest: $calling_gene_by_closest\n";
print "abcd: True\n" if(($g_opts->{'abcd'}));
print "calling_gene_by_distance: $calling_gene_by_distance\n";
print "calling_gene_by_both: $calling_gene_by_both\n";
print "distance: $distance\n" if(($g_opts->{'distance'}));
print "thread_limit: $thread_limit\n";
print "output_folder: $output_folder\n";
print "outout_name: $outout_name\n";
print "\n";



my %chr_name=();
if(($g_opts->{'chromosome_name_file'})){
    open(IN,"$chromosome_name_file")||die;
    while(my $t=<IN>){
        chomp $t;
        my @a=split(/\t/,$t);
        $chr_name{$a[0]}=$a[1];
    }
}

my $genome=read_gff($gff,\%chr_name);





open(IN,"$peak_file")||die;
my @lines=();
my $title="";
while(my $t=<IN>){
	next if($t=~/^\#/);
	next if($t=~/^"\#/);
	next if($t=~/^[\t\s]*$/);
	next if($t=~/^$/);
	next if($t=~/^\t+$/);
	chomp $t;
	$t=~s/\r//;
	$t=~s/\n//;
	if($t=~/^chr\tstart/i){
		$title=$t;
		next;
	}
	if($t=~/^chrom.*?\tstart/i){
		$title=$t;
		next;
	}
	
	
	push(@lines,$t);
}

my $splitarr=split_arr(\@lines,$thread_limit);
my $thread_num=0;
if($calling_gene_by_closest==1){
	my @thread_res=();
	my %res_finial=();
	print "calling_gene_by_closest\n";
	foreach(my $i=0;$i<$thread_limit;$i++){
		$thread_res[$i]=threads->create(\&call_closest_gene,\@{$$splitarr[$i]},$genome,$summit_column,$cds,$abcd);
		$thread_num++;
	}
	
	for my $t (@thread_res) {
		my $res=$t->join();
		#print "$res\n";
		foreach my $l (keys %$res){
			$res_finial{$l}=$$res{$l};
			#print "$l\n";
		}
	}
	
	open(OUT,">$output_folder/$outout_name.calling_gene_by_closest")||die;
	print OUT "$title\tupstream_closest_gene\tupstream_distance\tdownstream_closest_gene\tdownstream_distance\tgene_body\tfirst_closest_gene\tfirst_closest_distance\tsecond_closest_gene\tsecond_closest_distance\n";
	foreach my $l (@lines){
		print OUT "$l\t$res_finial{$l}\n";
	}
}

# if($calling_gene_by_distance==1){
# 	my @thread_res=();
# 	my %res_finial=();
# 	print "calling_gene_by_distance\n";
# 	foreach(my $i=0;$i<$thread_limit;$i++){
# 		$thread_res[$i]=threads->create(\&call_gene_by_dist,\@{$$splitarr[$i]},$genome,$distance,$summit_column,$cds);
# 		$thread_num++;
# 	}
# 	
# 	for my $t (@thread_res) {
# 		my $res=$t->join();
# 		#print "$res\n";
# 		foreach my $l (keys %$res){
# 			$res_finial{$l}=$$res{$l};
# 		}
# 	}
# 	
# 	open(OUT,">$output_folder/$outout_name.calling_gene_by_$distance")||die;
# 	print OUT "$title\tupstream_gene,distance\tdownstream_gene,distance\n";
# 	foreach my $l (@lines){
# 		print OUT "$l\t$res_finial{$l}\n";
# 	}
# }

# if($calling_gene_by_both==1){
# 	my @thread_res=();
# 	my %res_finial=();
# 	print "calling_gene_by_both\n";
# 	foreach(my $i=0;$i<$thread_limit;$i++){
# 		$thread_res[$i]=threads->create(\&call_gene_by_both,\@{$$splitarr[$i]},$genome,$distance,$summit_column,$cds);
# 		$thread_num++;
# 	}
# 	
# 	for my $t (@thread_res) {
# 		my $res=$t->join();
# 		#print "$res\n";
# 		foreach my $l (keys %$res){
# 			$res_finial{$l}=$$res{$l};
# 		}
# 	}
# 	
# 	open(OUT,">$output_folder/$outout_name.calling_gene_by_both_$distance")||die;
# 	print OUT "$title\tupstream_closest_gene\tupstream_distance\tdownstream_closest_gene\tdownstream_distance\tfirst_closest_gene\tfirst_closest_distance\tsecond_closest_gene\tsecond_closest_distance\n";
# 	foreach my $l (@lines){
# 		print OUT "$l\t$res_finial{$l}\n";
# 	}
# }










sub call_gene_by_dist{
	my ($line,$genome,$distance,$summit_column,$cds)=@_;
	my %res=();
	if($summit_column==-1){
		if($cds==0){
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $mid=int(($a[1]+$a[2])/2);
				my $downstream_gene="";
				for(my $i=$mid;$i<=$mid+$distance;$i++){
					if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
						my $dist=$i-$mid;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,$dist;";
						}
						
					}
				}
				my $upstream_gene="";
				for(my $i=$mid-1;$i>=$mid-$distance;$i--){
					if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
						my $dist=$mid-$i;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,$dist;";
						}
						
					}
				}
				if($upstream_gene eq ""){
					$upstream_gene="-";
				}
				if($downstream_gene eq ""){
					$downstream_gene="-";
				}
				$res{$l}="$upstream_gene\t$downstream_gene";
			}
		}
		else{
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $mid=int(($a[1]+$a[2])/2);
				my $downstream_gene="";
				for(my $i=$mid;$i<=$mid+$distance;$i++){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
						my $dist=$i-$mid;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,$dist;";
						}
						
					}
				}
				my $upstream_gene="";
				for(my $i=$mid-1;$i>=$mid-$distance;$i--){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
						my $dist=$mid-$i;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,$dist;";
						}
						
					}
				}
				if($upstream_gene eq ""){
					$upstream_gene="-";
				}
				if($downstream_gene eq ""){
					$downstream_gene="-";
				}
				$res{$l}="$upstream_gene\t$downstream_gene";
			}
		}
		
	}
	else{
		if($cds==0){
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $summit=$a[$summit_column];
				my $downstream_gene="";
				for(my $i=$summit;$i<=$summit+$distance;$i++){
					if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
						my $dist=$i-$summit;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,$dist;";
						}
						
					}
				}
				my $upstream_gene="";
				for(my $i=$summit-1;$i>=$summit-$distance;$i--){
					if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
						my $dist=$summit-$i;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,$dist;";
						}
						
					}
				}
				if($upstream_gene eq ""){
					$upstream_gene="-";
				}
				if($downstream_gene eq ""){
					$downstream_gene="-";
				}
				$res{$l}="$upstream_gene\t$downstream_gene";
			}
		}
		else{
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $summit=$a[$summit_column];
				my $downstream_gene="";
				for(my $i=$summit;$i<=$summit+$distance;$i++){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
						my $dist=$i-$summit;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,$dist;";
						}
						
					}
				}
				my $upstream_gene="";
				for(my $i=$summit-1;$i>=$summit-$distance;$i--){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
						my $dist=$i-$summit;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,$dist;";
						}
						
					}
				}
				if($upstream_gene eq ""){
					$upstream_gene="-";
				}
				if($downstream_gene eq ""){
					$downstream_gene="-";
				}
				
				$res{$l}="$upstream_gene\t$downstream_gene";
			}
		}
	}
	return \%res;
}

sub call_gene_by_both{
	my ($line,$genome,$distance,$summit_column,$cds)=@_;
	my %res=();
	if($summit_column==-1){
		if($cds==0){
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $mid=int(($a[1]+$a[2])/2);
				my $downstream_gene="";
				my $down_dist=0;
				for(my $i=$mid;$i<=$mid+$distance;$i++){
					if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
						$down_dist=$i-$mid;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,";
						}
						last;
					}
				}
				my $upstream_gene="";
				my $up_dist=0;
				for(my $i=$mid-1;$i>=$mid-$distance;$i--){
					if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
						$up_dist=$mid-$i;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,";
						}
						last;
					}
				}
				$upstream_gene=~s/,$//;
				$downstream_gene=~s/,$//;
				my $first="";
				my $second="";
				
				
				if($upstream_gene ne "" and $downstream_gene ne ""){
					if($up_dist>=$down_dist){
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					else{
						$second="$downstream_gene\t$down_dist";
						$first="$upstream_gene\t$up_dist";
					}
				}
				elsif($upstream_gene eq "" and $downstream_gene ne ""){
					$upstream_gene="-";
					$up_dist="-";
					$first="$downstream_gene\t$down_dist";
					$second="-\t-";
				}
				elsif($upstream_gene ne "" and $downstream_gene eq ""){
					$downstream_gene="-";
					$down_dist="-";
					$first="$upstream_gene\t$up_dist";
					$second="-\t-";
				}
				else{
					$upstream_gene="-";
					$up_dist="-";
					$downstream_gene="-";
					$down_dist="-";
					$first="-\t-";
					$second="-\t-";
				}
				
				
				$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$first\t$second";
			}
		}
		else{
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $mid=int(($a[1]+$a[2])/2);
				my $downstream_gene="";
				my $down_dist=0;
				for(my $i=$mid;$i<=$mid+$distance;$i++){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
						$down_dist=$i-$mid;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,";
						}
						last;
					}
				}
				my $upstream_gene="";
				my $up_dist=0;
				for(my $i=$mid-1;$i>=$mid-$distance;$i--){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
						$up_dist=$mid-$i;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,";
						}
						last;
					}
				}
				$upstream_gene=~s/,$//;
				$downstream_gene=~s/,$//;
				my $first="";
				my $second="";
				
				
				if($upstream_gene ne "" and $downstream_gene ne ""){
					if($up_dist>=$down_dist){
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					else{
						$second="$downstream_gene\t$down_dist";
						$first="$upstream_gene\t$up_dist";
					}
				}
				elsif($upstream_gene eq "" and $downstream_gene ne ""){
					$upstream_gene="-";
					$up_dist="-";
					$first="$downstream_gene\t$down_dist";
					$second="-\t-";
				}
				elsif($upstream_gene ne "" and $downstream_gene eq ""){
					$downstream_gene="-";
					$down_dist="-";
					$first="$upstream_gene\t$up_dist";
					$second="-\t-";
				}
				else{
					$upstream_gene="-";
					$up_dist="-";
					$downstream_gene="-";
					$down_dist="-";
					$first="-\t-";
					$second="-\t-";
				}
				
				
				$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$first\t$second";
			}
		}
		
	}
	else{
		if($cds==0){
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $summit=$a[$summit_column];
				my $downstream_gene="";
				my $down_dist=0;
				for(my $i=$summit;$i<=$summit+$distance;$i++){
					if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
						$down_dist=$i-$summit;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,";
						}
						last;
					}
				}
				my $upstream_gene="";
				my $up_dist=0;
				for(my $i=$summit-1;$i>=$summit-$distance;$i--){
					if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
						$up_dist=$summit-$i;
						foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,";
						}
						last;
					}
				}
				$upstream_gene=~s/,$//;
				$downstream_gene=~s/,$//;
				my $first="";
				my $second="";
				
				
				if($upstream_gene ne "" and $downstream_gene ne ""){
					if($up_dist>=$down_dist){
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					else{
						$second="$downstream_gene\t$down_dist";
						$first="$upstream_gene\t$up_dist";
					}
				}
				elsif($upstream_gene eq "" and $downstream_gene ne ""){
					$upstream_gene="-";
					$up_dist="-";
					$first="$downstream_gene\t$down_dist";
					$second="-\t-";
				}
				elsif($upstream_gene ne "" and $downstream_gene eq ""){
					$downstream_gene="-";
					$down_dist="-";
					$first="$upstream_gene\t$up_dist";
					$second="-\t-";
				}
				else{
					$upstream_gene="-";
					$up_dist="-";
					$downstream_gene="-";
					$down_dist="-";
					$first="-\t-";
					$second="-\t-";
				}
				
				
				$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$first\t$second";
			}
		}
		else{
			foreach my $l (@$line){
				my @a=split(/\t/,$l);
				my $summit=$a[$summit_column];
				my $downstream_gene="";
				my $down_dist=0;
				for(my $i=$summit;$i<=$summit+$distance;$i++){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
						$down_dist=$i-$summit;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
							$downstream_gene.="$g,";
						}
						last;
					}
				}
				my $upstream_gene="";
				my $up_dist=0;
				for(my $i=$summit-1;$i>=$summit-$distance;$i--){
					if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
						$up_dist=$summit-$i;
						foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
							$upstream_gene.="$g,";
						}
						last;
					}
				}
				$upstream_gene=~s/,$//;
				$downstream_gene=~s/,$//;
				my $first="";
				my $second="";
				if($upstream_gene ne "" and $downstream_gene ne ""){
					if($up_dist>=$down_dist){
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					else{
						$second="$downstream_gene\t$down_dist";
						$first="$upstream_gene\t$up_dist";
					}
				}
				elsif($upstream_gene eq "" and $downstream_gene ne ""){
					$upstream_gene="-";
					$up_dist="-";
					$first="$downstream_gene\t$down_dist";
					$second="-\t-";
				}
				elsif($upstream_gene ne "" and $downstream_gene eq ""){
					$downstream_gene="-";
					$down_dist="-";
					$first="$upstream_gene\t$up_dist";
					$second="-\t-";
				}
				else{
					$upstream_gene="-";
					$up_dist="-";
					$downstream_gene="-";
					$down_dist="-";
					$first="-\t-";
					$second="-\t-";
				}
				
				$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$first\t$second";
			}
		}
	}
	return \%res;
}

sub call_closest_gene{
	my ($line,$genome,$summit_column,$cds,$abcd)=@_;
	my %res=();
	
	if($abcd==1){
		if($summit_column==-1){
			if($cds==0){
				print "A:calling_closest_gene\t$abcd\t$summit_column\t$cds\t1\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $mid=int(($a[1]+$a[2])/2);
					my $downstream_gene="";
					my $down_dist=0;
					
					
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,+,";
							}
							last;
						}
						
						if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
							$down_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
								$downstream_gene.="$g,-,";
							}
							last;
						}
						
						
					
						
					}
					
					
					my $upstream_gene="";
					my $up_dist=0;
					
					
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,-,";
							}
							last;
						}
						
						if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
							$up_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
								$upstream_gene.="$g,+,";
							}
							last;
						}
					}
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					
					if($upstream_gene=~/,\+/ and $downstream_gene=~/,\+/){
						$upstream_gene="";
					}
					if($upstream_gene=~/,-/ and $downstream_gene=~/,-/){
						$downstream_gene="";
					}
					
					
					my $first="";
					my $second="";
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
			else{
				print "B:calling_closest_gene\t$abcd\t$summit_column\t$cds\t2\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $mid=int(($a[1]+$a[2])/2);
					my $downstream_gene="";
					my $down_dist=0;
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,+,";
							}
							last;
						}
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
							$down_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
								$downstream_gene.="$g,-,";
							}
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,-,";
							}
							last;
						}
						
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
							$up_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
								$upstream_gene.="$g,+,";
							}
							last;
						}
					}
					
					
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					
					if($upstream_gene=~/,\+/ and $downstream_gene=~/,\+/){
						$upstream_gene="";
					}
					if($upstream_gene=~/,-/ and $downstream_gene=~/,-/){
						$downstream_gene="";
					}
					
					
					
					
					my $first="";
					my $second="";
					if($up_dist>=$down_dist){
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					else{
						$second="$downstream_gene\t$down_dist";
						$first="$upstream_gene\t$up_dist";
					}
				
				
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
				
				
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
		
		}
		else{
			if($cds==0){
				print "C:calling_closest_gene\t$abcd\t$summit_column\t$cds\t3\n";
				foreach my $l (@$line){
					#print "$l\n";
					my @a=split(/\t/,$l);
					my $summit=$a[$summit_column];
					my $downstream_gene="";
					my $down_dist=0;
					if(!defined $$genome{"chr_length"}{$a[0]}){
						print "$l\n";
					}
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,+,";
							}
							last;
						}
						
						if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
							$down_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
								$downstream_gene.="$g,-,";
							}
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,-,";
							}
							last;
						}
						
						if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
							$up_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
								$upstream_gene.="$g,+,";
							}
							last;
						}
					}
					
					
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$downstream_gene_end";
						
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					print "$summit||$upstream_gene_end||$downstream_gene||$downstream_gene_end||$upstream_gene||||$gene_body\n";
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					my $first="";
					my $second="";
					
					if($upstream_gene=~/,\+/ and $downstream_gene=~/,\+/){
						$upstream_gene="";
					}
					if($upstream_gene=~/,-/ and $downstream_gene=~/,-/){
						$downstream_gene="";
					}
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
				
				
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
			else{
				print "D:calling_closest_gene\t$abcd\t$summit_column\t$cds\t4--\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $summit=$a[$summit_column];
					my $downstream_gene="";
					my $down_dist=0;
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,+," if(!exists $$genome{"exclude"}{$g});
							}
							
						}
						
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
							$down_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
								$downstream_gene.="$g,-," if(!exists $$genome{"exclude"}{$g});
							}
							
						}
						if($downstream_gene eq ""){
							next; 
						}
						else{
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$summit-$i;
							
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,-," if(!exists $$genome{"exclude"}{$g});
							}
							
						}
						
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
							$up_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
								$upstream_gene.="$g,+," if(!exists $$genome{"exclude"}{$g});
							}
							
						}
						if($upstream_gene eq ""){
							next; 
						}
						else{
							last;
						}
					}
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					my $first="";
					my $second="";
					
					if(($upstream_gene=~/,-/ and $downstream_gene=~/,+/) or ($upstream_gene=~/,+/ and $downstream_gene=~/,-/)){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene=~/,\+/ and $downstream_gene=~/,\+/){
						
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					elsif($upstream_gene=~/,-/ and $downstream_gene=~/,-/){

						$first="$upstream_gene\t$up_dist";
						$second="$downstream_gene\t$down_dist";
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}

					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
		}
	}
	else{
	
		if($summit_column==-1){
			if($cds==0){
				print "E:calling_closest_gene\t$abcd\t$summit_column\t$cds\t5\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $mid=int(($a[1]+$a[2])/2);
					my $downstream_gene="";
					my $down_dist=0;
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,";
							}
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,";
							}
							last;
						}
					}
					
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					my $first="";
					my $second="";
					
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
			else{
				print "F:calling_closest_gene\t$abcd\t$summit_column\t$cds\t6\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $mid=int(($a[1]+$a[2])/2);
					my $downstream_gene="";
					my $down_dist=0;
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,";
							}
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,";
							}
							last;
						}
					}
					
					
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$mid;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$mid;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$mid-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$mid-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					my $first="";
					my $second="";
					
					if($up_dist>=$down_dist){
						$first="$downstream_gene\t$down_dist";
						$second="$upstream_gene\t$up_dist";
					}
					else{
						$second="$downstream_gene\t$down_dist";
						$first="$upstream_gene\t$up_dist";
					}
					
					
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
				
				
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
		
		}
		else{
			if($cds==0){
				print "G:calling_closest_gene\t$abcd\t$summit_column\t$cds\t7\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $summit=$a[$summit_column];
					my $downstream_gene="";
					my $down_dist=0;
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,";
							}
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,";
							}
							last;
						}
					}
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
					
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					my $first="";
					my $second="";
				
					
					
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
				
				
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
			else{
				print "H:calling_closest_gene\t$abcd\t$summit_column\t$cds\t8\n";
				foreach my $l (@$line){
					my @a=split(/\t/,$l);
					my $summit=$a[$summit_column];
					my $downstream_gene="";
					my $down_dist=0;
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}){
							$down_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene.="$g,";
							}
							last;
						}
					}
					my $upstream_gene="";
					my $up_dist=0;
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}){
							$up_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"cds_start2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene.="$g,";
							}
							last;
						}
					}
					
					
					my $downstream_gene_end="";
					my $down_end_dist=0;
					my $upstream_gene_end="";
					my $up_end_dist=0;
					
					for(my $i=$summit;$i<$$genome{"chr_length"}{$a[0]};$i++){
						if(exists $$genome{"end2gene"}{$a[0]}{"+"}{$i}){
							$down_end_dist=$i-$summit;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"+"}{$i}}){
								$downstream_gene_end.="$g,+,";
							}
							last;
						}
					}
					
					
					for(my $i=$summit-1;$i>0;$i--){
						if(exists $$genome{"end2gene"}{$a[0]}{"-"}{$i}){
							$up_end_dist=$summit-$i;
							foreach my $g (keys %{$$genome{"end2gene"}{$a[0]}{"-"}{$i}}){
								$upstream_gene_end.="$g,-,";
							}
							last;
						}
					}
					
					my $gene_body="";
					
					
					if($upstream_gene_end eq $downstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($downstream_gene_end eq $upstream_gene){
						$gene_body.="$upstream_gene_end";
					}
					
					if($gene_body eq ""){
						$gene_body="-";
					}
				
				
					$upstream_gene=~s/,$//;
					$downstream_gene=~s/,$//;
					my $first="";
					my $second="";
					
					if($upstream_gene ne "" and $downstream_gene ne ""){
						if($up_dist>=$down_dist){
							$first="$downstream_gene\t$down_dist";
							$second="$upstream_gene\t$up_dist";
						}
						else{
							$second="$downstream_gene\t$down_dist";
							$first="$upstream_gene\t$up_dist";
						}
					}
					elsif($upstream_gene eq "" and $downstream_gene ne ""){
						$upstream_gene="-";
						$up_dist="-";
						$first="$downstream_gene\t$down_dist";
						$second="-\t-";
					}
					elsif($upstream_gene ne "" and $downstream_gene eq ""){
						$downstream_gene="-";
						$down_dist="-";
						$first="$upstream_gene\t$up_dist";
						$second="-\t-";
					}
					else{
						$upstream_gene="-";
						$up_dist="-";
						$downstream_gene="-";
						$down_dist="-";
						$first="-\t-";
						$second="-\t-";
					}
				
				
					$res{$l}="$upstream_gene\t$up_dist\t$downstream_gene\t$down_dist\t$gene_body\t$first\t$second";
				}
			}
		}
	}
	return \%res;
}



sub read_gff{
    my ($gff,$chr_name)=@_;
    open(IN,"$gff")||die "can't open $gff\n";
    my %chr_length=();

    #my %chr_name=();
    my %gene_cds_start=();
    my %gene_cds_end=();
    my %gene_info=();
    my %gene=();
    my %gene_name2gene_id=();
    my %mRNA=();
    my %gene2ID=();
    my %ID2gene=();
    my %gene2common=();
	my %exon=();
	my %CDS=();
	my %promoter_region=();
	my %tail_region=();
    my %all_gene_name=();
    my %start2gene=();
    my %end2gene=();
    my %cds_start2gene=();
    my %cds_end2gene=();
    my %gene_cds=();
    my %gene_desc=();
    my %gene_hash=();
    my %exclude=();
    $/="\n";
    while(my $t=<IN>){
            chomp $t;
            next if($t=~/^\#/);
            my @a=split(/\t/,$t);
            if($a[2] eq "region" or $a[2] eq "chromosome"){
				if(exists $$chr_name{$a[0]}){
					$a[0]=$$chr_name{$a[0]};
				}
                if(exists $chr_length{$a[0]}){
                    $chr_length{$a[0]}=$a[4] if($a[4]>$chr_length{$a[0]});
                }
                else{
                     $chr_length{$a[0]}=$a[4];
                }
                #print "$a[0]\t$a[4]\n";
            }
        	if($a[2] eq "ncRNA" or $a[2] eq "tRNA" or $a[2] eq "rRNA" or $a[2] eq "snRNA" or $a[2] eq "sRNA" or $a[2] eq "snoRNA" or $a[2] eq "siRNA" or $a[2] eq "scaRNA" or $a[2] eq "piRNA" or $a[2] eq "centromere" or $a[2] eq "repeat_region" or $a[2] eq "retrotransposon" or $a[2] eq "long_terminal_repeat" or $a[2] eq "blocked_reading_frame"){
        		$a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }    
                $exclude{$tmphash{"ID"}}=0 if(exists $tmphash{"ID"});
                $exclude{$tmphash{"Parent"}}=0 if(exists $tmphash{"Parent"});
                $exclude{$tmphash{"Name"}}=0 if(exists $tmphash{"Name"});
                $exclude{$tmphash{"name"}}=0 if(exists $tmphash{"name"});
        	}
        	
        	
        
            if($a[2] eq "CDS"){
                
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }      

                
                my $mrna_ID=$tmphash{"Parent"};

            
                if(exists $mRNA{$mrna_ID}){
                    next if(!exists $ID2gene{$mRNA{$mrna_ID}});
                    $CDS{$a[0]}{$ID2gene{$mRNA{$mrna_ID}}}{$a[3]."\t".$a[4]}=0;
                    if(exists $gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}}){ 
                        $gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}}=$a[3] if($a[3]<$gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}});
                        $gene_cds_end{$ID2gene{$mRNA{$mrna_ID}}}=$a[4] if($a[4]>$gene_cds_end{$ID2gene{$mRNA{$mrna_ID}}});
                        #print "$mRNA{$mrna_ID}\t$a[3]\n";
                        
                    }
                    else{
                    	
                        $gene_cds_start{$ID2gene{$mRNA{$mrna_ID}}}=$a[3];
                        $gene_cds_end{$ID2gene{$mRNA{$mrna_ID}}}=$a[4];
                        #print "$mrna_ID\t$a[3]\n";
                    }
                }
                else{
                	next if(!exists $ID2gene{$mrna_ID});
                	$CDS{$a[0]}{$ID2gene{$mrna_ID}}{$a[3]."\t".$a[4]}=0;
                    if(exists $gene_cds_start{$ID2gene{$mrna_ID}}){
                        
                        $gene_cds_start{$ID2gene{$mrna_ID}}=$a[3] if($a[3]<$gene_cds_start{$ID2gene{$mrna_ID}});
                        $gene_cds_end{$ID2gene{$mrna_ID}}=$a[4] if($a[4]>$gene_cds_end{$ID2gene{$mrna_ID}});
                        #print "$mRNA{$mrna_ID}\t$a[3]\n";
                    }
                    else{
                        
                        $gene_cds_start{$ID2gene{$mrna_ID}}=$a[3];
                        $gene_cds_end{$ID2gene{$mrna_ID}}=$a[4];
                        #print "$mrna_ID\t$a[3]\n";
                    }
                }
            
            
            
        
            }
            if($a[2] eq "mRNA"){
                
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }      
                
                my $ID=$tmphash{"Parent"};
                my $mrna_ID=$tmphash{"ID"};

                $mRNA{$mrna_ID}=$ID;

        
            }
            
            if($a[2] eq "exon"){
                
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }
                
                
                my $mrna_ID=$tmphash{"Parent"};

                
				if(exists $$chr_name{$a[0]}){
					$a[0]=$$chr_name{$a[0]};
				}

                
                if(exists $mRNA{$mrna_ID}){
                	next if(!exists $ID2gene{$mRNA{$mrna_ID}});
                    $exon{$a[0]}{$ID2gene{$mRNA{$mrna_ID}}}{$a[3]."\t".$a[4]}=0;
                }
                else{
                	next if(!exists $ID2gene{$mrna_ID});
                    $exon{$a[0]}{$ID2gene{$mrna_ID}}{$a[3]."\t".$a[4]}=0;
                }               
        
            }
            
            

            if($a[2] eq "gene" or $a[2] eq "pseudogene"){
                if($chr_name){
                    if(exists $$chr_name{$a[0]}){
                        $a[0]=$$chr_name{$a[0]};
                    }
                }
                $a[8]=uri_unescape($a[8]);
                my @aaa=split(/;/,$a[8]);
                my %tmphash=();
                foreach my $ttt (@aaa){
                	my ($tmpna,$tmpval)=split(/=/,$ttt);
                	next if(!defined $tmpna or !defined $tmpval);
                	$tmphash{$tmpna}=$tmpval;
                }
                if($a[8]=~/GeneID:(\d+)/){
                	$tmphash{"entrz_id"}=$1; 
                }
                if($a[8]=~/HGNC:(\d+)/){
                	$tmphash{"hgnc_id"}=$1; 
                }
                
                if($a[0] eq "ChrM" or $a[0] eq "NC_012920.1" or $a[0] eq "mito_A_nidulans_FGSC_A4" or $a[0] eq "Ca19-mtDNA" or $a[0] eq "NC_010339.1" or $a[0] eq "NC_005089.1"){
                	$exclude{$tmphash{"ID"}}=0 if(exists $tmphash{"ID"});
                	$exclude{$tmphash{"Parent"}}=0 if(exists $tmphash{"Parent"});
                	
                }
                if(exists $tmphash{"pseudo"} and $tmphash{"pseudo"} eq "true"){
                	$exclude{$tmphash{"ID"}}=0 if(exists $tmphash{"ID"});
                	$exclude{$tmphash{"Parent"}}=0 if(exists $tmphash{"Parent"});
                }
              	
              	
            	my $ID=$tmphash{"ID"};
                my $gene_name=$ID;
                $gene_name=$tmphash{"Name"} if(exists $tmphash{"Name"});
                $gene_name=$tmphash{"gene_name"} if(exists $tmphash{"gene_name"});

                
                my $gene_id=$tmphash{"entrz_id"} if(exists $tmphash{"entrz_id"});
                my $common="";
                $common=$tmphash{"gene"} if(exists $tmphash{"gene"});
                $common=$tmphash{"Gene"} if(exists $tmphash{"Gene"});

                $gene{$gene_name}=$gene_name;
                

                $gene_name2gene_id{$gene_name}=$gene_id if(defined $gene_id);
                
                $all_gene_name{$gene_name}=0;
                
                $gene2ID{$gene_id}=$ID if(defined $gene_id);
                
                $gene2ID{$gene_name}=$ID;
                
                $ID2gene{$ID}=$gene_name;
                
                $gene2common{$gene_id}=$common if($common ne "" and defined $gene_id);
                
                $gene{$gene_id}=$gene_name if(defined $gene_id);
                
                $gene2common{$gene_name}=$common if($common ne "");
                $gene{$common}=$gene_name if($common ne "");
                $gene{uc($common)}=$gene_name if($common ne "");

                %{$gene_hash{$gene_name}}=%tmphash;
                
                $gene_info{$gene_name}="$a[0]\t$a[3]\t$a[4]\t$a[6]";
                $gene_info{$gene_id}="$a[0]\t$a[3]\t$a[4]\t$a[6]" if(defined $gene_id);
                if($a[6] eq "+"){
                	$start2gene{$a[0]}{$a[6]}{$a[3]}{$gene_name}=0;
                	$end2gene{$a[0]}{$a[6]}{$a[4]}{$gene_name}=0;
                }
                else{
                	$start2gene{$a[0]}{$a[6]}{$a[4]}{$gene_name}=0;
                	$end2gene{$a[0]}{$a[6]}{$a[3]}{$gene_name}=0;
                }
                my $alias=$tmphash{"Alias"} if(exists $tmphash{"Alias"});
                my $note=$tmphash{"Note"} if(exists $tmphash{"Note"});
                my $desc=$tmphash{"description"} if(exists $tmphash{"description"});

                
                if(defined $alias and $alias ne ""){
                    my @aa=split(/,/,$alias);
                    foreach my $tt (@aa){
                            if(!exists $gene{$tt}){
                                    $gene{$tt}=$gene_name;
                            }
                    }
                }
                
                
                if(defined $note){
                    $gene_desc{$gene_name}="$note";
                    $gene_desc{$gene_id}="$note" if(defined $gene_id);
                }
                
                if(defined $desc){
                    $gene_desc{$gene_name}="$desc";
                    $gene_desc{$gene_id}="$desc" if(defined $gene_id);
                }
            }
        
    }
    
    
    foreach my $g (keys %all_gene_name){
    	
    	my @info=split(/\t/,$gene_info{$g});
    	if(exists $gene_cds_start{$g}){
    		$gene_cds{$info[0]}{$g}=$gene_cds_start{$g}."\t".$gene_cds_end{$g};
    		if($info[3] eq "+"){
    			$cds_start2gene{$info[0]}{$info[3]}{$gene_cds_start{$g}}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$gene_cds_end{$g}}{$g}=0;
    		}
    		else{
    			$cds_start2gene{$info[0]}{$info[3]}{$gene_cds_end{$g}}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$gene_cds_start{$g}}{$g}=0;
    		}
    		
    	}
    	else{
    		$gene_cds{$info[0]}{$g}=$info[1]."\t".$info[2];
			if($info[3] eq "+"){
    			$cds_start2gene{$info[0]}{$info[3]}{$info[1]}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$info[2]}{$g}=0;	
    		}
    		else{
    			$cds_start2gene{$info[0]}{$info[3]}{$info[2]}{$g}=0;
    			$cds_end2gene{$info[0]}{$info[3]}{$info[1]}{$g}=0;
    		}
    	}
    }
    
    

    my %genome=();
    %{$genome{"gene"}}=%gene;
    %{$genome{"exon"}}=%exon;
    %{$genome{"gene_cds"}}=%gene_cds;
    
    %{$genome{"promoter_region"}}=%promoter_region;
    %{$genome{"tail_region"}}=%tail_region;
    %{$genome{"exon"}}=%exon;
    %{$genome{"CDS"}}=%CDS;
    %{$genome{"ID2gene"}}=%ID2gene;
    %{$genome{"all_gene_name"}}=%all_gene_name;
    
    %{$genome{"chr_length"}}=%chr_length;
    %{$genome{"gene_info"}}=%gene_info;
    %{$genome{"gene_cds_start"}}=%gene_cds_start;
    %{$genome{"gene_cds_end"}}=%gene_cds_end;
    
    %{$genome{"cds_start2gene"}}=%cds_start2gene;
    %{$genome{"cds_end2gene"}}=%cds_end2gene;
    %{$genome{"start2gene"}}=%start2gene;
    %{$genome{"end2gene"}}=%end2gene;
    
    %{$genome{"gene2common"}}=%gene2common;
    %{$genome{"gene2ID"}}=%gene2ID;
    %{$genome{"gene_name2gene_id"}}=%gene_name2gene_id;
    %{$genome{"gene_hash"}}=%gene_hash;
    %{$genome{"gene_desc"}}=%gene_desc;
    %{$genome{"exclude"}}=%exclude;

    
    return \%genome;


}


sub split_arr{
	my ($a,$threads)=@_;
	my $n=@$a;
	print "$n\n";
	my $sub_n=int($n/$threads)+1;
	my @res=();
	my $j=0;
	for(my $i=0;$i<$n;$i++){
		push(@{$res[$j]},$$a[$i]);
		$j++ if(($i+1)%$sub_n==0);
	}
	return \@res;
}


###############  End Time  ########################################################
my $Time_End="";
my $TE=time();
my $UseTime=$TE-$TS;
my $danwei="Seconds";
my $UT="";
if($UseTime>120 and $UseTime<7200){
    my $m=int($UseTime/60);
    my $s=$UseTime%60;
    $UT= "Used Time: $m minutes $s seconds";
}
elsif($UseTime<=120){
    $UT= "Used Time: $UseTime seconds";
}
elsif($UseTime>=7200){
       my $h=int($UseTime/3600);
       my $m=int(($UseTime%3600)/60);
       my $s= ($UseTime%3600)%60 ;
       $UT= "Used Time: $h hours $m minutes $s seconds";
}
 $Time_End = sub_format_datetime(localtime(time()));
 print "\nEnd Time :[$Time_End]\nTime-consuming: $UT\n";

###############################################################################
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}








__END__

=head1 NAME

script - Using Getopt::Long and Pod::Usage

=head1 SYNOPSIS

script [options] [args ...]

Options: 

   -gff                                   * gff file 
   -peak_file                             * ur input peak bed files, "*chr_name  *start  *end  length  summit  .."
   -summit_column                         Int. which column is summit info in your peak_file, count from 1. Default is unset
                                          without setting summit_column, script will use middle site(mid) of peak to call gene

   -calling_gene_by_closest               Flag. call closest one gene from peak summit or mid(cis regulated).Default is TRUE
   -abcd                                  Flag. required by -calling_gene_by_closest,call closest genes cis and trans.

   -thread_limit                          Int. thread_limit, how many threads want to be used, Default is 2
   -distance                              Int. distance, required by -calling_gene_by_distance and -calling_gene_by_both, default is 2000 (bp)



   -chromosome_name_file                  chromosome_name_file,if chr_name in gff is diff with ur chr_name in peak_file,
                                          plz set it, "chr_name(in gff)  chr_name(in peak_file)"
                                          
   -common                                Flag. if want to output gene common name
   -cds                                   Flag. if want to use cds coordinate as gene coordinate

                                          * means required.
                                          
   -help            brief help message

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut