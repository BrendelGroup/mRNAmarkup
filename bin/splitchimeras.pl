#! /usr/bin/perl -w
print "\n\n";
#
# splitchimeras.pl
# Version of April 8, 2015.  Volker Brendel

use strict;
use Getopt::Std;


my $USAGE="\nUsage: $0 -m msbfile -q queryfile -o outfile\


** This script reads input from a MuSeqBox BLASTx output file specified\
   by the argument to -m and searches for lines indicating potential chimeras.\
   The corresponding sequences are pulled out of the FASTA file supplied by\
   the argument to -q.  Output consists of a copy of the FASTA file, modified\
   to replace the chimeric entries by several shorter segments.  The shorter\
   segments are name according to their parent name, appended by 'a' for the\
   5'-segment, 'b' for the 3'-segment, and where applicable 'c' for the\
   central segment.  'p' is appended for potential chimeras not split by the\
   criteria of the progam.

   \n";


my %args;
getopts('m:q:o:', \%args);

my $splitgap = 700;
my $LINE_LENGTH = 80;
my $inputSEQ = "";

my ($FHMSBF,$FHOUTF);

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
  open FHOUTF,  ">>$outputfile" || die ("Cannot open file: $outputfile");
}


my $line="";
my ($ort1,$ort2,$query,$subject1,$subject2)=("","","","","",""); 
my ($qfrom1,$qto1,$qfrom2,$qto2)=(0,0,0,0);

open FHMSBF,  "<$msbfile" || die ("Cannot open file: $msbfile"); 


while(defined($line=<FHMSBF>)){
    if ($line =~ /  !Potential chimera/) {
print "\n\n$line";
      my ($ort1,$ort2,$query,$subject1,$qfrom1,$qto1,$subject2,$qfrom2,$qto2)= $line =~ /  !Potential chimera(.)(.): Query (.+) matches (.+) from (\d+) to (\d+) and (.+) from (\d+) to (\d+)./;
      $subject1 = (split (/ /,$subject1))[0];
      $subject2 = (split (/ /,$subject2))[0];
      &splitquery($query,$ort1,$qfrom1,$qto1,$ort2,$qfrom2,$qto2);
    }
}



sub splitquery {
  my $qname  = $_[0];
  my $ort1   = $_[1];
  my $from1  = $_[2];
  my $to1    = $_[3];
  my $ort2   = $_[4];
  my $from2  = $_[5];
  my $to2    = $_[6];

  my ($FHQF);
  open FHQF,  "<$queryfile" || die ("Cannot open file: $queryfile"); 

  my $fullheader = "";
  my $header = "";
  my $sequence = "";
  my $foundflag = 0;
  while(defined($line=<FHQF>)){
    if ($line =~ /$qname/) {
      ($fullheader) = $line;
      ($header) = $line =~ /(^>[^\s]*)/;
      while(defined($line=<FHQF>) && $line !~ /^>/) {
        chomp($sequence .= $line);
      }
    $foundflag = 1;
    last;
    }
  }

  if ($foundflag == 0) {
    print "\n!!! Query file $qname does not seem to be in file $queryfile.  Please check.\n";
    exit;
  }
  my $db = "";
  my $sgmntname = "";
  if    ($header =~ /^>lcl/) {($sgmntname) = $header =~ /^>lcl\|([^\s]+)\s*.*/;}
  elsif ($header =~ /^>gnl/) {($db,$sgmntname) = $header =~ /^>gnl\|([^\s]*)\|([^\s\|]+).*/;}
  elsif ($header =~ /^>/) {($sgmntname) = $header =~ /^>([^\s]+)\s*.*/;}
  else {print "\nUnsupported FASTA header line: $header.  Please check.";}

  $inputSEQ = $sequence;
  $inputSEQ =~ s/[\s\d]//g;
  $inputSEQ =~ tr/a-z/A-Z/;
  my $seq_length = $inputSEQ =~ tr/a-zA-Z/a-zA-Z/;

  if ($ort1 eq "+"  &&  $ort2 eq "-") {
print "\nPLUSMINUS $sgmntname";
    if ($to2 > $to1) {
# completely to the right
print "\nTO THE RIGHT GAP ", $to2-$to1, "\n";
if ($to2-$to1 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to1+200, " of $qname";
      &PULL_SEG(1,$to1+200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to2-200, " of $qname";
      &PULL_SEG($seq_length,$to2-200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $to1+100 , " to " , $to2-100, " of $qname";
      &PULL_SEG($to1+100,$to2-100,*FHOUTF);

}
else {
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to2-1, " of $qname";
      &PULL_SEG(1,$to2-1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to1+1, " of $qname";
      &PULL_SEG($seq_length,$to1+1,*FHOUTF);
}
    }
# now $to2 <= $to1
    elsif ($from2 > $to1) {
print "\nRIGHT OVERLAP\n";
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to1, " of $qname";
      &PULL_SEG(1,$to1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to2, " of $qname";
      &PULL_SEG($seq_length,$to2,*FHOUTF);
    }
# now $from2 <= $to1
    elsif ($to2 > $from1) {
# inclusion
print "\nINCLUSION\n";
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to1 , " of $qname";
      &PULL_SEG(1,$to1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to2, " of $qname";
      &PULL_SEG($seq_length,$to2,*FHOUTF);
    }
    elsif ($from2 > $from1) {
# partial overlap
print "\nLEFT OVERLAP\n";
      print FHOUTF ">" , $sgmntname , "a segment " , $from1, " to " , $seq_length, " of $qname";
      &PULL_SEG($from1,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from2 , " to " , "1" , " of $qname";
      &PULL_SEG($from2,1,*FHOUTF);
    }
    else {
# completely to the left
print "\nTO THE LEFT GAP ", $from1-$from2, "\n";
if ($from1-$from2 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , $from1-200 , " to " , $seq_length, " of $qname";
      &PULL_SEG($from1-200,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from2+200 , " to " , "1" , " of $qname";
      &PULL_SEG($from2+200,1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $from2+100 , " to " , $from1-100, " of $qname";
      &PULL_SEG($from2+100,$from1-100,*FHOUTF);
}
else {
      print FHOUTF ">" , $sgmntname , "a segment " , $from2+1 , " to " , $seq_length, " of $qname";
      &PULL_SEG($from2+1,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from2 , " to " , "1" , " of $qname";
      &PULL_SEG($from2,1,*FHOUTF);
}

    }
  }
  if ($ort1 eq "-"  &&  $ort2 eq "+") {
print "\nMINUSPLUS";
    if ($to1 > $to2) {
# completely to the right
print "\nTO THE RIGHT GAP ", $to1-$to2, "\n";
if ($to1-$to2 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to2+200, " of $qname";
      &PULL_SEG(1,$to2+200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to1-200, " of $qname";
      &PULL_SEG($seq_length,$to1-200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $to2+100 , " to " , $to1-100, " of $qname";
      &PULL_SEG($to2+100,$to1-100,*FHOUTF);

}
else {
      print FHOUTF ">" . $sgmntname . "a segment " . "1" . " to " , $to1-1, " of $qname";
      &PULL_SEG(1,$to1-1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to2+1, " of $qname";
      &PULL_SEG($seq_length,$to2+1,*FHOUTF);
}
    }
# now $to1 <= $to2
    elsif ($from1 > $to2) {
print "\nRIGHT OVERLAP\n";
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to2, " of $qname";
      &PULL_SEG(1,$to2,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to1, " of $qname";
      &PULL_SEG($seq_length,$to1,*FHOUTF);
    }
# now $from1 <= $to2
    elsif ($to1 > $from2) {
# inclusion
print "\nINCLUSION\n";
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to2, " of $qname";
      &PULL_SEG(1,$to2,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $to1, " of $qname";
      &PULL_SEG($seq_length,$to1,*FHOUTF);
    }
    elsif ($from1 > $from2) {
# partial overlap
print "\nLEFT OVERLAP\n";
      print FHOUTF ">" , $sgmntname , "a segment " , $from2 , " to " , $seq_length, " of $qname";
      &PULL_SEG($from2,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from1 , " to " , "1" , " of $qname";
      &PULL_SEG($from1,1,*FHOUTF);
    }
    else {
# completely to the left
print "\nTO THE LEFT GAP ", $from2-$from1, "\n";
if ($from2-$from1 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , $from2-200 , " to " , $seq_length, " of $qname";
      &PULL_SEG($from2-200,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from1+200 , " to " , "1" , " of $qname";
      &PULL_SEG($from1+200,1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $from1+100 , " to " , $from2-100, " of $qname";
      &PULL_SEG($from1+100,$from2-100,*FHOUTF);
}
else {
      print FHOUTF ">" , $sgmntname , "a segment " , $from1+1 , " to " , $seq_length, " of $qname";
      &PULL_SEG($from1+1,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from1 , " to " , "1" , " of $qname";
      &PULL_SEG($from1,1,*FHOUTF);
}
    }
  }

  if ($ort1 eq "+"  &&  $ort2 eq "+") {
print "\nPLUSPLUS";
    if ($from2 > $to1) {
# completely to the right
print "\nTO THE RIGHT GAP ", $from2-$to1, "\n";
if ($from2-$to1 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , 1 , " to " , $to1+200, " of $qname";
      &PULL_SEG(1,$to1+200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from2-200 , " to " , $seq_length , " of $qname";
      &PULL_SEG($from2-200,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $to1+100 , " to " , $from2-100, " of $qname";
      &PULL_SEG($to1+100,$from2-100,*FHOUTF);
}
else {
      print FHOUTF ">" . $sgmntname . "a segment " . "1" . " to " , $from2-1, " of $qname";
      &PULL_SEG(1,$from2-1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $to1+1 , " to " , $seq_length, " of $qname";
      &PULL_SEG($seq_length,$to1+1,*FHOUTF);
}
    }
    elsif ($to2 < $from1) {
# completely to the left
print "\nTO THE LEFT GAP ", $from1-$to2, "\n";
if ($from1-$to2 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , 1 , " to " , $to2+200, " of $qname";
      &PULL_SEG(1,$to2+200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $from1-200 , " to " , $seq_length , " of $qname";
      &PULL_SEG($from1-200,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $to2+100 , " to " , $from1-100, " of $qname";
      &PULL_SEG($to2+100,$from1-100,*FHOUTF);
}
else {
      print FHOUTF ">" , $sgmntname , "a segment " , "1" , " to " , $to2, " of $qname";
      &PULL_SEG(1,$to2,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $to2+1 , " to " , $seq_length , " of $qname";
      &PULL_SEG($to2+1,$seq_length,*FHOUTF);
}
    }
    else {
print "\nLEAVE\n";
      if ($fullheader =~ /potential fusion transcript assembly/) {
        print FHOUTF ">" , $sgmntname, " ... potential fusion transcript assembly based on distinct BLASTx matches";
      }
      else {
        print FHOUTF ">" , $sgmntname, "p ... potential fusion transcript assembly based on distinct BLASTx matches";
      }
      &PULL_SEG(1,$seq_length,*FHOUTF);
    }
  }

  if ($ort1 eq "-"  &&  $ort2 eq "-") {
print "\nMINUSMINUS";
    if ($to2 > $from1) {
# completely to the right
print "\nTO THE RIGHT GAP ", $to2-$from1, "\n";
if ($to2-$from1 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , 1 , " to " , $from1+200, " of $qname";
      &PULL_SEG(1,$from1+200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $to2-200 , " to " , $seq_length , " of $qname";
      &PULL_SEG($to2-200,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $from1+100 , " to " , $to2-100, " of $qname";
      &PULL_SEG($from1+100,$to2-100,*FHOUTF);
}
else {
      print FHOUTF ">" , $sgmntname , "a segment " , $to2-1 , " to " , "1" , " of $qname";
      &PULL_SEG($to2-1,1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $from1+1, " of $qname";
      &PULL_SEG($seq_length,$from1+1,*FHOUTF);
}
    }
    elsif ($from2 < $to1) {
# completely to the left
print "\nTO THE LEFT GAP ", $to1-$from2, "\n";
if ($to1-$from2 >= $splitgap) {
print "\nTHREE SEGMENTS for $qname\n";
      print FHOUTF ">" , $sgmntname , "a segment " , 1 , " to " , $from2+200, " of $qname";
      &PULL_SEG(1,$from2+200,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $to1-200 , " to " , $seq_length , " of $qname";
      &PULL_SEG($to1-200,$seq_length,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "c segment " , $from2+100 , " to " , $to1-100, " of $qname";
      &PULL_SEG($from2+100,$to1-100,*FHOUTF);
}
else {
      print FHOUTF ">" , $sgmntname , "a segment " , $to1-1 , " to " , "1" , " of $qname";
      &PULL_SEG($to1-1,1,*FHOUTF);
      print FHOUTF ">" , $sgmntname , "b segment " , $seq_length , " to " , $from2+1 , " of $qname";
      &PULL_SEG($seq_length,$from2+1,*FHOUTF);
}
    }
    else {
print "\nLEAVE\n";
      if ($fullheader =~ /potential fusion transcript assembly/) {
        print FHOUTF ">" , $sgmntname, " ... potential fusion transcript assembly based on distinct BLASTx matches";
      }
      else {
        print FHOUTF ">" , $sgmntname, "p ... potential fusion transcript assembly based on distinct BLASTx matches";
      }
      &PULL_SEG($seq_length,1,*FHOUTF);
    }
  }

}


sub REV_COMP {
  my $seq = $_[0];
  $seq =~ tr/aAgGcCtT/tTcCgGaA/;
  my @seq = split(//,$seq);
  
  return join("",reverse(@seq));
}

sub FASTA_SL {
  my $seq = $_[0];
  my $seq_length = $seq =~ tr/a-zA-Z/a-zA-Z/;
  if(defined($_[1])){
    $LINE_LENGTH = $_[1];
  }

  my $i=0;
  my $line;
  while($i < $seq_length){
    if(($i + $LINE_LENGTH) <= $seq_length){
      $line .= substr($seq,$i,$LINE_LENGTH) . "\n";
      $i += $LINE_LENGTH;
    }else{
      $line .= substr($seq,$i,($seq_length-$i)) . "\n";
      last;
    }
  }

  return $line;
}

sub PULL_SEG {
  my $START = $_[0];
  my $END = $_[1];
  my $FH = $_[2];

  my $Comp_Rev = 0;
    if($START > $END){
      $Comp_Rev = $START;
      $START = $END;
      $END = $Comp_Rev;
      $Comp_Rev = 1;
    }

  my $segment = substr($inputSEQ,($START-1),($END - $START + 1));

  if($Comp_Rev){
    print $FH " (reverse-complement)\n";
    print $FH FASTA_SL(REV_COMP($segment),$LINE_LENGTH);
  }else{
    print $FH "\n";
    print $FH &FASTA_SL($segment,$LINE_LENGTH);
  }
}
