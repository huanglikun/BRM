#!/usr/bin/perl
use strict;
use warnings;
use 5.014;
use File::Basename;

# usage
my $script = basename $0;
die "Usage: perl $script vcf2bsa_conf.txt markers.vcf markers.bsa\n" if(not @ARGV);

# input and output
my $file_cnf = shift @ARGV;
my $file_vcf = shift @ARGV;# both full or slim version
my $file_out = shift @ARGV;

# read global options from configure file
my %conf = &read_configure($file_cnf);

# check contents of configure
die "Error: pool1/pool2 is missing\n" if(not exists $conf{'Table2x2.pool1'} and not exists $conf{'Table2x2.pool2'});
my $s1_list = $conf{'Table2x2.pool1'} if(exists $conf{'Table2x2.pool1'});
my $s2_list = $conf{'Table2x2.pool2'} if(exists $conf{'Table2x2.pool2'});
my $p1_list = $conf{'Table2x2.parent1'} if(exists $conf{'Table2x2.parent1'});
my $p2_list = $conf{'Table2x2.parent2'} if(exists $conf{'Table2x2.parent2'});

# global array
my @RG = ();
my @S1_LIST = ();
my @S2_LIST = ();
my @P1_LIST = ();
my @P2_LIST = ();
#
@S1_LIST = split /\s*,\s*/, $s1_list if($s1_list);
@S2_LIST = split /\s*,\s*/, $s2_list if($s2_list);
@P1_LIST = split /\s*,\s*/, $p1_list if($p1_list);
@P2_LIST = split /\s*,\s*/, $p2_list if($p2_list);

# 
my $dir_out = dirname $file_out;
do {mkdir $dir_out or die ""} if(not -d $dir_out);

# read vcf file
open FILE, $file_vcf or die "";
open OUT, ">$file_out" or die "";
while(<FILE>){
	chomp ;
	# read header info
	if($_=~m/^##/){
		next;
	}
	# read title
	if($_=~m/^#(.*)/){
		my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
		@RG = @samples;
		# check
		my %samples = &assign_samples(\@samples);
		die "Error: pool1 is not match\n" if( @S1_LIST and &check_samples(\%samples, \@S1_LIST) );
		die "Error: pool2 is not match\n" if( @S2_LIST and &check_samples(\%samples, \@S2_LIST) );
		die "Error: parent1 is not match\n" if( @P1_LIST and &check_samples(\%samples, \@P1_LIST) );
		die "Error: parent2 is not match\n" if( @P2_LIST and &check_samples(\%samples, \@P2_LIST) );
		next;
	}
	my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
	my @smp1 = &new_samples(\@samples, \@S1_LIST) if(@S1_LIST);
	my @smp2 = &new_samples(\@samples, \@S2_LIST) if(@S2_LIST);
	my @par1 = &new_samples(\@samples, \@P1_LIST) if(@P1_LIST);
	my @par2 = &new_samples(\@samples, \@P2_LIST) if(@P2_LIST);

	#
	my $p1_allele = &which_parent_allele($tag, \@par1) if(@par1);
	my $p2_allele = &which_parent_allele($tag, \@par2) if(@par2);
	#
	my @mix = (@smp1, @smp2);
	$p2_allele = &which_another_allele($tag,\@mix,$p1_allele) if(@par1 and not @par2);
	$p1_allele = &which_another_allele($tag,\@mix,$p2_allele) if(@par2 and not @par1);
	#
	($p1_allele, $p2_allele) = (0, 1) if(not defined $p1_allele and not defined $p2_allele);
	#
	next if($p1_allele eq "NA");
	next if($p2_allele eq "NA");
	#
	# two column of samples in one file.
	my ($s1p1, $s1p2) = &sort_dp_by_allele($tag, \@smp1, $p1_allele, $p2_allele) if(@smp1);
	my ($s2p1, $s2p2) = &sort_dp_by_allele($tag, \@smp2, $p1_allele, $p2_allele) if(@smp2);
	my $str_s1p1 = &sum_array(@$s1p1) if(@smp1);
	my $str_s1p2 = &sum_array(@$s1p2) if(@smp1);
	my $str_s2p1 = &sum_array(@$s2p1) if(@smp2);
	my $str_s2p2 = &sum_array(@$s2p2) if(@smp2);
	#
	print OUT "$chr\t$loc";
	print OUT "\t$str_s1p1\t$str_s1p2" if(@smp1);
	print OUT "\t$str_s2p1\t$str_s2p2" if(@smp2);
	print OUT "\n";
}
close FILE;
close OUT;

########################################################################
# functions
########################################################################
sub read_configure{
	# read configure file
	my $file_conf = shift @_;
	my %CONF = ();
	open CONF, "<", $file_conf or die "";
	while(<CONF>){
		chomp ;
		next if($_=~m/^#/);
		next if($_=~m/^\s/);
		next if($_ eq "");
		$_=~s/\s*$//;
		$_=~s/\s*#.*$//;
		$_=~s/\s*;\s*$//;
		$_=~s/\s*([=\t])\s*/$1/;
		my ($id, $item) = (split /[=\t]/, $_);
		next  if(not $id or not $item);
		$CONF{$id}=$item;
	}
	close CONF;
	return %CONF;
}

sub assign_samples{
	my $samples = shift @_;
	#
	my %hash = ();
	for(my $i=0;$i<scalar @RG;$i++){
		$hash{$RG[$i]}=$$samples[$i];
	}
	return %hash;
}

sub check_samples{
	my $samples = shift @_;
	my $rg_arr  = shift @_;
	my $err = 0;
	foreach my $rg (@$rg_arr){
		if(not exists $$samples{$rg}){
			print STDERR "$rg not exists\n";
			$err++;
		}
	}
	return $err;	
}

sub new_samples{
	my $samples = shift @_;
	my $list    = shift @_;
	#
	my %samples=&assign_samples($samples);
	# seletect specific samples
	my @new=();
	foreach my $rg (@$list){
		push @new, $samples{$rg};
	}
	return @new;
}

sub which_parent_allele{
	my $tag_str = shift @_;
	my $samples = shift @_;
	#
	my @tag = split /:/, $tag_str;
	my %allele = ();
	foreach my $str (@$samples){
		next if($str =~ m/^\./);
		my @arr  = split /:/, $str;
		my %hash = &assign_by_tag(\@tag,\@arr);
		my ($a1,$a2) = split /\/|\|/, $hash{'GT'};
		$allele{$a1}++;
		$allele{$a2}++;
	}
	my @allele = keys %allele;
	my $num=scalar @allele;
	my $allele = "NA";
	$allele=$allele[0] if($num==1);
	return $allele;
}

sub which_another_allele{
	my $tag_str = shift @_;
	my $samples = shift @_;
	my $parent1 = shift @_;
	#
	my @tag = split /:/, $tag_str;
	my @depth = ();
	foreach my $str (@$samples){
		next if($str =~ m/^\./);
		my @arr = split /:/, $str;
		my @cov = &read_allele_depth(\@tag,$str);
		for(my $i=0;$i<scalar @cov;$i++){
			$depth[$i]+=$cov[$i];
		}
	}
	my $another="NA";
	return $another if ($parent1 eq 'NA');
	my $max=0;
	for(my $a=0;$a<scalar @depth;$a++){
		next if($a==$parent1);
		if($max<$depth[$a]){
			$max=$depth[$a];
			$another=$a;
		}
	}
	return $another;
}

sub assign_by_tag{
	my $tag = shift @_;
	my $arr = shift @_;
	#
	my %hash=();
	for(my $i=0;$i<scalar @$tag;$i++){
		$hash{$$tag[$i]}=$$arr[$i];
	}
	return %hash;
}

sub sort_dp_by_allele{
	my ($tag_str,$samples,$p1,$p2) = @_;
	my @tag = split /:/, $tag_str;
	my @dp1 = ();
	my @dp2 = ();
	foreach my $str (@$samples){
		if($str =~ m/^\./){
			push @dp1,"NA";
			push @dp2,"NA";
			next;
		}
		my @cov=&read_allele_depth(\@tag,$str);
		push @dp1, $cov[$p1];
		push @dp2, $cov[$p2];
	}
	return (\@dp1,\@dp2);
}

sub read_allele_depth{
	my $tag = shift @_;
	my $str = shift @_;
	#
	my @arr=split /:/, $str;
	my %hash=&assign_by_tag($tag,\@arr);
	#
	my @cov = ();
	if (exists $hash{'AD'}){
		@cov = split /,/,$hash{'AD'};
	}else{
		@cov = ($hash{'RO'});
		push @cov, (split /,/,$hash{'AO'});
	}
	return @cov;
}

sub sum_array{
	my $mk  = 0;
	my $sum = 0;
	foreach my $n (@_){
		next if($n eq "NA");
		$mk=1;
		$sum+=$n;
	}
	$sum = "NA" if(not $mk);
	return $sum;
}

