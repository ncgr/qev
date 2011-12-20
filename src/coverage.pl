#! /usr/bin/env perl
use warnings;
use strict;
use Bio::DB::Sam;
#use lib "/sw/perl/lib/site_perl/5.10.1/x86_64-linux/Math";
use Math::CDF;

my $x2 = 1;
my $dof = 20;
 
#print Math::CDF::pbinom(3,1000,.001) . "\n";



#unless ($ARGV[1]){
#	die "usage:$0 file.bam ref.fa\n";
#}

my $sam = Bio::DB::Sam->new(-autoindex=>1, -bam=>$ARGV[0]);
my @seq_ids = $sam->seq_ids;

foreach my $ref (@seq_ids){
	my @align = $sam->features(-type=>'match', -seq_id=>$ref);
#foreach my $a (@align) {
#my ($rid, $rseq, $rstart, $rstop) = ($a->query->seq_id, $a->query->dna, $a->start, $a->end);
#print "$rid, $rseq, $rstart, $rstop\n";
#}
	my $len = $sam->length($ref);
	my @pairs = $sam->features(-type=>'read_pair', -seq_id=>$ref);
	my @pcov = map{0}(1 .. $len);
	for (my $pos=1;$pos<=$len;$pos++){
		$pcov[$pos]=0;
	}
	my $read_length=0.0;
	my $pair_count;
	foreach my $p (@pairs){
		next unless $p;
		my @ss = $p->get_SeqFeatures;
		next unless $ss[1];
		my $start = $ss[0]->start;
		my $stop = $ss[1]->end;
		$read_length+=$stop -$start +1;
		foreach my $i ($start .. $stop){
			++$pcov[$i-1];
		}
		$pair_count++;
	}
	$read_length =  int($read_length/$pair_count);
	#print $read_length;
	my @hist=(0,0,0,0,0,0,0,0);
	for (my $pos=101;$pos<=$len-100;$pos++){
		#print "$pos $pcov[$pos]\n";
		$hist[$pcov[$pos]]++;
	}
#	my $med=median(@pcov);
	my $Pr=1.0;
	
	my $arg=2.0;
	my $ci=0;
	my $ci_pos=1;
	my $Gcum=0.0;
	my $chicum=0.0;
	my $chidegrees=0.0;
	my $G_degrees=0.0;
	$chidegrees=$G_degrees=$len-1;
	my @model = map{0}(1 .. $len);
        my @model_exp = map{0}(1 .. $len);
        my $summodel;
        my $pcov_sum;
        $pcov_sum=0.0;
	for (my $pos=1;$pos<=$len;$pos++){
                my $endl=$read_length;
                if($pos<$endl){
                        $endl=$pos;
                }
                if($len-$pos+1<$endl){
                        $endl = ($len-$pos+1);
                }
                
		$model[$pos]=$endl/$len;
 		$pcov_sum+=$pcov[$pos];
        }
	for (my $pos=1;$pos<=$len;$pos++){
		$model_exp[$pos]=$model[$pos]*$pair_count;
        }

  
	for (my $pos=1;$pos<=$len;$pos++){
		$Pr = Math::CDF::pbinom($pcov[$pos], $pair_count, $model[$pos]);
		#Accumulate G-statistic
		my $ent;
		if($pcov[$pos]==0){
			$ent =0;
		}else{
			$ent = 2*$pcov[$pos]*log($pcov[$pos]/$model_exp[$pos]);
		}
		$Gcum += $ent;
                print "pos $pos pcov $pcov[$pos] p $model[$pos] exp $model_exp[$pos] Pr $Pr ent $ent Gcum $Gcum\n";
		if($Pr < $arg){
			$arg = $Pr;
			$ci=$pcov[$pos];
			$ci_pos=$pos;
		}
        }

	my $stddev=sqrt(2*$G_degrees);
	my $StandardGstat=($Gcum-$G_degrees)/$stddev;
        my $Pr_model = Math::CDF::pnorm($StandardGstat);
	#print "Pr_model $Pr_model gcum $Gcum gdegrees $G_degrees stddev $stddev stand $StandardGstat\n";	
	print ">$ref model $Pr_model am $arg ci $ci mpos $ci_pos len $len pairs $pair_count hist $hist[0] $hist[1] $hist[2] $hist[3] $hist[4] $hist[5] $hist[6]\n";
        #print "@pcov\n";
        my ($cov) = $sam->features(-type=>'coverage', -seq_id=>$ref);
        my @bcov = $cov->coverage;
}

sub median {
	@ == 1 or die ('Sub usage: $median = median(\@array);');
	my ($array_ref) = @_;
	my $count = scalar @$array_ref;
	# Sort a COPY of the array, leaving the original untouched
	my @array = sort { $a <=> $b } @$array_ref;
	if ($count % 2) {
		return $array[int($count/2)];
	} else {
		return ($array[$count/2] + $array[$count/2 - 1]) / 2;
	}
} 
