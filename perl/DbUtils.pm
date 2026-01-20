sub checkGeneInfo {
    my %args = (
        fileHandle => "",
        expData => {},
        @_
    );

    my $fhdb = $args{fileHandle};
    my %expData = %{$args{expData}};
    while (my $line = <$fhdb>) {
        chomp $line;
        next if ($line =~ /Synonyms/);
        my @line = split(/\t/, $line);
        
        my $ens;
        if ($line[5] =~ /Ensembl:(.+?)\|/) { #capture through regexp usin ()
            $ens = $1;
        }

        if ($expData{$line[2]} || ($ens && $expData{$ens})) {
            my $key;
            if ($expData{$line[2]}) {
                $key = $line[2];
            } elsif ($expData{$ens}) {
                $key = $ens;
            }
            if ($line[1]) {
                $expData{$key}{trueId} = $line[1];
            } elsif ($ens) {     
                ${$expData{$key}}{trueId} = $ens;
            }
        }
    }
    return(%expData);
}

sub checkmiRNInfo {
    my %args = (
        fileHandle => "",
        expData => {},
        @_
    );
    
    my $fhdb = $args{fileHandle};
    my %expData = %{$args{expData}};

    while (my $line = <$fhdb>) {
        chomp $line;
        my @line = split(/\t/, $line);

        if ($expData{$line[1]}) {
            $expData{$line[1]}{trueId} = $line[2];
        }
    }
    return(%expData);
}

sub checkProtInfo{
    my %args = (
        fileHandle => "",
        expData => {},
        @_
    );
    
    my $fhdb = $args{fileHandle};
    my %expData = %{$args{expData}};

    while (my $line = <$fhdb>) {
        chomp $line;
        my @line = split(/\t/, $line);
        $line[5] =~ s/;//;

        if ($expData{$line[3]}) {
            $expData{$line[3]}{trueId} = $line[5];
        }
    }
    return(%expData);
}







1;