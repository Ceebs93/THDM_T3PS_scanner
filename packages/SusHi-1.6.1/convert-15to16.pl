#!/usr/bin/perl -w
#
# convert sushi-1.5 input files to sushi-1.6
#
$version = "1.0";
@FILES = @ARGV;

foreach $file (@FILES) {
    system("/bin/cp $file $file.bak");
    open(OUT,">$file");
    print {OUT} ("# converted 1.5 to 1.6 using convert-15to16.pl, v$version\n");
    open(IN,$file.".bak");
    $block = "";
    while (<IN>) {
        $conv = 0;
	$line = $_;
	$lineorig = $line;
	if ($line =~ /^Block +(\S+)/) {
	    $block = $1;
	    print {OUT} ($line);
	    next;
	}
	if ($block =~ /PDFSPEC/) {
# PDF set number is now in entry 10:	    
	    if ($line =~ s/^\s+4\s/ 10/) {
                $conv=1;	                    
		print("$block (OLD): ",$lineorig);
		print("$block (NEW): ",$line);
	    }
    	} elsif ($block =~ /SUSHI/) {
# approximate NNLO stop is now SUSHI(5)=12:
	    if ($line =~ s/^\s+5\s+3/ 5 12/) {
                $conv=1;	                    
		print("$block (OLD): ",$lineorig);
		print("$block (NEW): ",$line);
	    }
	}
	if ($conv) {print {OUT} ("#v1.5: $lineorig")}
	print {OUT} ($line);
    }
}


