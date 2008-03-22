use strict;
use File::Find;

# get path containing "..\VC\BIN"
my @path = split(';', $ENV{PATH});
my $pathvc = "";
my $i;
for $i (0..$#path) {
   if ($path[$i] =~ /VC\\BIN$/) {
      $pathvc = $path[$i];
      last;
   }#if
}#for
#print ("found path: $pathvc\n");

# input file names
my $file2find   = "cl.exe";
my @directories = ($pathvc);

# output file names
my $file_found = "";
my $dir_found  = "";

# check if VC++ is installed in fixed location
find(\&wanted, @directories);

sub wanted {
   if ($_ =~ /$file2find$/) {
      $file_found = $_;
      $dir_found  = $File::Find::dir;
      last;
   }#if
}#sub

# return result
if ($file_found =~ /^cl.exe$/) {
   print ("$file_found");
} else {
   print("ERR_CLEXE");
}#if





