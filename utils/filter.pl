#!/usr/bin/perl -w
#!/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;

# option
my ($helpFlag,$iTmpDir,$iFirstFQ,$iSecondFQ,$iOut1,$iOut2,$iSum,$iLog);


GetOptions(
        "h|?|help"              => \$helpFlag,
        "d|tempdir=s"              => \$iTmpDir,
        "r1|read1=s"               => \$iFirstFQ,            ## Input paired-end  fastq file 1
        "r2|read2=s"               => \$iSecondFQ,            ## Input paired-end fastq file 2
        "o1|out1=s"               => \$iOut1,            ## output file 1
        "o2|out2=s"               => \$iOut2,            ## output file 2
        "s|sum=s"               => \$iSum,            ## output file 2
        "l|log=s"               => \$iLog,            ## output file 2
) || die "\n";

my $msg = "Arguments: [-d tmpdir] [-r1 first pe file .fq.gz or .fq or fastq or fastq.gz] [-r2 second pe file .fq.gz or .fq] [-o1 first output read file] [-o2 2nd output read file] [-s summary file]\t\n";

checkOptions();
my @delList;

my $PWD= dirname(__FILE__);

my $nfilter = "python $PWD/filterNCountAndQscoreFromFastqGZ.v02.py";


print_content("\nStart NFilter at ",$iLog);

`date >> $iLog`;
print_content("\n",$iLog);
print_content("command:perl filter.pl -d $iTmpDir -r1 $iFirstFQ -r2 $iSecondFQ -o1 $iOut1 -o2 $iOut2 -s $iSum -l $iLog\n",$iLog);

# Making gzip for paired end input files
#-----------------------------
my $firstfq = $iFirstFQ;
if ($firstfq !~/.gz$/){
	$firstfq=gzip_func($iFirstFQ,$iTmpDir."/".getFullFileName($iFirstFQ),$iLog);
	push(@delList,$firstfq);
}
my $secondfq = $iSecondFQ;
if ($secondfq !~/.gz$/){
	$secondfq=gzip_func($iSecondFQ,$iTmpDir."/".getFullFileName($iSecondFQ),$iLog);
	push(@delList,$secondfq);
}


# always gzip for nfilter output
my $out1 = $iOut1;
if ($out1 !~/.gz$/){
	$out1 .=".gz";
}
my $out2 = $iOut2;
if ($out2!~/.gz$/){
	$out2 .=".gz";
}

my $cmd_filter="$nfilter $firstfq $secondfq $out1 $out2 $iSum >> $iLog";
print_content("real command:$cmd_filter ....\n",$iLog);

if (system($cmd_filter) !=0){
	print_content("Error in $cmd_filter\n" ,$iLog);
	die ("$cmd_filter has been failed!");
}

# process trimmed output file and format
#----------------------------------------------
if ($iOut1 !~/.gz$/){
	my $zcatcmd = "zcat $out1 > $iOut1";
	print_content("$zcatcmd ....\n" ,$iLog);
	if (system($zcatcmd) !=0){
		print_content("Error in $zcatcmd\n" ,$iLog);
		die ("$zcatcmd has some problem\t\n");
	}
	push(@delList,$out1);
}
if ($iOut2 !~/.gz$/){
	my $zcatcmd = "zcat $out2 > $iOut2";
	print_content("$zcatcmd ....\n" ,$iLog);
	if (system($zcatcmd) !=0){
		print_content("Error in $zcatcmd\n" ,$iLog);
		die ("$zcatcmd has some problem\t\n");
	}
	push(@delList,$out2);
}
#----------------------------------------------

# delete delLists
#----------------------------------------------
foreach(@delList){
	print_content("rm -f $_...\n",$iLog);
	if (-f $_){
		`rm -f $_`;
	}
	else{
		print_content("Error in rm -f $_...\n",$iLog);
	}
}
#----------------------------------------------
print_content("Completed NFilter at",$iLog);
`date >> $iLog`;
print_content("\n----------------------\n",$iLog);

sub makeTar{
	my $res = shift();
	my @input= @{shift()};
	my $tarcmd= "tar cfvz $res ".join(" ",@input);
	if (system($tarcmd)!=0){
		print "Die:".$tarcmd."\n";
	} 
}

sub cp_res{
	my $input = shift();
	my $output = shift();
	my $logfile = shift();

	if ($input ne $output){
		my $cpcmd = "cp $input $output";
		print_content($cpcmd."...\n",$logfile);
		if (system ($cpcmd) !=0){
			print_content("Error in $cpcmd\n",$logfile);
			die ("cp $input $output has been failed!");
		}
	}
	return $output;
}
sub mv_res{
	my $input = shift();
	my $output = shift();
	my $logfile = shift();
	my $mvcmd = "mv $input $output";
	print_content($mvcmd."...\n",$logfile);
	if (system ($mvcmd) !=0){
		print_content("Error in $mvcmd\n",$logfile);
		die ("mv $input $output has been failed!");
	}
	return $output;
}



sub openAppend
{
        my ($fileName) = @_;

	return *STDOUT unless defined $fileName;

	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c >> $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c >> $fileName" : ">>$fileName") ||
	die("Open error: $fileName");
	return $fd;
}
sub openOutput
{
	my ($fileName) = @_;

	return *STDOUT unless defined $fileName;
	
	my ($fd);
	open($fd, $fileName =~ /.gz$/ ? "| gzip -c > $fileName" : $fileName =~ /.bz(ip)?2$/ ? "| bzip2 -z -c > $fileName" : ">$fileName") ||
	die("Open error: $fileName");
	return $fd;
}
sub print_new_content{
	my $msg = shift();
	my $fname = shift();
	my $fscr = openOutput($fname);
	print $fscr $msg;
	close($fscr);
}
sub print_content{
	my $msg = shift();
	my $fname = shift();
	my $fscr = openAppend($fname);
	print $fscr $msg;
	close($fscr);
}
sub gzip_func{
	my $fqfile = shift();
	my $move_file = shift();
	my $log_file = shift();
		
	
	cp_res($fqfile,$move_file,$log_file);
	
	my $gzipcmd = "gzip ".$move_file;
	print_content("$gzipcmd...\n",$iLog);
	if (system ($gzipcmd) !=0){
		print_content("Error in $gzipcmd...\n",$iLog);
		die ($gzipcmd." has been failed!");
	}
	return $move_file.".gz";
}

sub trim {
        my @result = @_;

        foreach (@result) {
                s/^\s+//;
                s/\s+$//;
        }

        return wantarray ? @result : $result[0];
}
sub getFullFileName{
	my $fpath = shift();
	my ($name,$dir,$suffix) = fileparse($fpath,qr"\..[^.]*$");
	return $name.$suffix;
}
sub getDirName{
	my $fpath = shift();
	my ($name,$dir,$suffix) = fileparse($fpath,qr"\..[^.]*$");
	return $dir;
}
sub getFileName{
	my $fpath = shift();
	my ($name,$dir,$suffix) = fileparse($fpath,qr"\..[^.]*$");
	return $name;
}
# compare file type 
# @param fpath file path
# @param ext file extension to be compared. Note that this parameter starts with . for file extension
sub checkFileType{
	my $fpath = shift();
	my $ext = shift();
	my ($name,$dir,$suffix) = fileparse($fpath,qr"\..[^.]*$");
	if ($suffix eq $ext){
		return 1;
	}
	return 0;
}

sub checkOptions
{

	$iTmpDir = shift(@ARGV) if !defined $iTmpDir && @ARGV>0;
	$iFirstFQ = shift(@ARGV) if !defined $iFirstFQ && @ARGV>0;
	$iSecondFQ = shift(@ARGV) if !defined $iSecondFQ && @ARGV > 0;
	$iOut1 = shift(@ARGV) if !defined $iOut1 && @ARGV > 0;
	$iOut2 = shift(@ARGV) if !defined $iOut2 && @ARGV > 0;
	$iSum = shift(@ARGV) if !defined $iSum && @ARGV > 0;

	if ($helpFlag || !$iTmpDir || !$iFirstFQ || !$iSecondFQ || !$iOut1 || !$iOut2 ||!$iSum ){
		die ($msg);
	}

}

