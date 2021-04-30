#! /usr/bin/perl -w
#
# annotateByMuSeqBox.pl
# Version of April 8, 2015.  Volker Brendel

use strict;
use Getopt::Std;


my $USAGE="\nUsage: $0 -m msbfile -q queryfile -o outfile\


** This script reads input from a MuSeqBox BLASTx output file specified\
   by the argument to -m and changes the FASTA headers in the corresponding\
   query file specified by the argument to -q.  Output will be written to\
   the mandatory argument to -o and is a copy of queryfile with modified\
   headers that show the best qualifying matches, as read from the MuSeqBox\
   file.

   \n";


my %args;
getopts('m:q:o:', \%args);

my ($FHMSBF,$FHQF,$FHOUTF);

if (!defined($args{m})) {
  print "\n!!! No MuSeqBox file specified.\n\n";
  die $USAGE;
}
my $msbfile = $args{m};
if (! -e $msbfile) {
  print "\n!!! MuSeqBox file $msbfile does not exist.\n\n";
  die $USAGE;
}
if (!defined($args{q})) {
  print "\n!!! No query file specified.\n\n";
  die $USAGE;
}
my $queryfile = $args{q};
if (! -e $queryfile) {
  print "\n!!! Query file $queryfile does not exist.\n\n";
  die $USAGE;
}
if (!defined($args{o})) {
  print "\n!!! No output file specified.\n\n";
  die $USAGE;
}
else {
  my $outputfile = $args{o};
  open FHOUTF,  ">$outputfile" || die ("Cannot open file: $outputfile"); 
}


my $line="";
my (@QueryID, @SubjectID, @QLen, @HSP, @HLen, @CovQ, @Qx, @Qy, @Sx, @Sy, @SLen, @CovS, @Pid, @Psi, @NGap, @Frame, @Score, @Eval, @Db, @Annotation);
my ($qid, $sid, $ann);
my %qid2sid;
my %qid2ann;

open FHMSBF,  "<$msbfile" || die ("Cannot open file: $msbfile"); 
if (defined($line=<FHMSBF>)){
  if ($line !~ /^BLASTx/) {
    print "\n\n ! MuSeqBox BLASTx output expected.  Please check your input file $msbfile.\n";
    die $USAGE;
  }
}

while(defined($line=<FHMSBF>)){
    if ($line =~ /^QueryID/) {
      $QueryID[0] = 0;
      $line =~ /SubjectID/; $QueryID[1] = $-[0]-2; $SubjectID[0] = $-[0];
      $line =~ /QLen/; $SubjectID[1] = $-[0]-3;
      $QLen[0] = $-[0]-1; $QLen[1] = $QLen[0] + 4;
      $HSP[0]  = $QLen[1]+2; $HSP[1] = $HSP[0] + 4;
      $HLen[0] = $HSP[1]+2; $HLen[1] = $HLen[0] + 6;
      $CovQ[0] = $HLen[1]+2; $CovQ[1] = $CovQ[0] + 4;
      $Qx[0]   = $CovQ[1]+2; $Qx[1] = $Qx[0] + 4;
      $Qy[0]   = $Qx[1]+2; $Qy[1] = $Qy[0] + 4;
      $Sx[0]   = $Qy[1]+2; $Sx[1] = $Sx[0] + 4;
      $Sy[0]   = $Sx[1]+2; $Sy[1] = $Sy[0] + 4;
      $SLen[0] = $Sy[1]+2; $SLen[1] = $SLen[0] + 4;
      $CovS[0] = $SLen[1]+2; $CovS[1] = $CovS[0] + 4;
      $Pid[0]  = $CovS[1]+2; $Pid[1] = $Pid[0] + 4;
      $Psi[0]  = $Pid[1]+2; $Psi[1] = $Psi[0] + 4;
      $NGap[0] = $Psi[1]+2; $NGap[1] = $NGap[0] + 3;
      $Frame[0]= $NGap[1]+2; $Frame[1] = $Frame[0] + 4;
      $Score[0]= $Frame[1]+2; $Score[1] = $Score[0] + 7;
      $Eval[0] = $Score[1]+2; $Eval[1] = $Eval[0] + 6;
      $Db[0] = $Eval[1]+2;
      $line =~ /Annotation/; $Db[1] = $-[0]-2; $Annotation[0] = $-[0];
    }
    elsif ($line !~ /BLAST/  &&  $line !~ /^ /  &&  $line !~ /^-/  &&  $line !~ /no_hit/  &&  $line !~ /^No query/  &&  $line !~ /^Database/  &&  $line !~ /^$/) {
      $qid = substr $line, $QueryID[0], $QueryID[1]-$QueryID[0]+1;
      ($qid) = $qid =~ /([^\s]+)/;
      $sid = substr $line, $SubjectID[0], $SubjectID[1]-$SubjectID[0]+1;
      ($sid) = $sid =~ /([^\s]+)/;
      $ann = substr $line, $Annotation[0];
      $ann =~ s/\s+$//;
      if (! exists $qid2sid{$qid}) {
        $qid2sid{$qid} = $sid;
        $qid2ann{$qid} = $ann;
      }
    }
}


open FHQF,  "<$queryfile" || die ("Cannot open file: $queryfile"); 

my $header    = "";
my $qseqname  = "";
while(defined($line=<FHQF>)){
  if ($line =~ /^>/) {
    $header = $line;
    if    ($header =~ /^>lcl/) {($qseqname) = $header =~ /^>lcl\|([^\s]+)\s*.*/;}
    elsif ($header =~ /^>gnl/) {($qseqname) = $header =~ /^>gnl\|[^\s]*\|([^\s\|]+).*/;}
    elsif ($header =~ /^>/)    {($qseqname) = $header =~ /^>([^\s]+)\s*.*/;}
    else {print "\nUnsupported FASTA header line: $header.  Please check.";}
    if (exists $qid2sid{$qseqname}) {
      $header =~ s/$/ matches $qid2sid{$qseqname} \($qid2ann{$qseqname}\)/;
    }
    print FHOUTF $header;
  }
  else {
    print FHOUTF $line;
  }
}
