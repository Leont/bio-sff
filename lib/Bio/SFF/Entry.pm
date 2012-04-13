package Bio::SFF::Entry;

use Moo;
use Sub::Name;

use Scalar::Util qw/looks_like_number/;

for my $attr (qw/clip_qual_left clip_qual_right clip_adaptor_left clip_adaptor_right/) {
	has $attr => (
		is => 'ro',
		required => 1,
		isa => sub {
			return looks_like_number($_[0]);
		},
	);
}

for my $attr (qw/name bases/){ 
	has $attr => (
		is => 'ro',
		required => 1,
		isa => sub {
			return defined and ref($_[0]) eq '';
		},
	);
}

my %unpack = (
	flowgram_values     => 'n*',
	flow_index_per_base => 'C*',
	quality_scores      => 'C*',
);

for my $attr(qw/flowgram_values flow_index_per_base quality_scores/) {
	my $raw = "${attr}_raw";
	has $raw => (
		is => 'ro',
		required => 1,
		init_arg => $attr,
		isa => sub {
			return defined and ref($_[0]) eq '';
		},
	);
	my $meth = "_$attr";
	has "_$attr" => (
		is => 'ro',
		init_arg => undef,
		lazy => 1,
		default => sub {
			return [ unpack $unpack{$attr}, $_[0]->$raw ];
		},
	);
	no strict 'refs';
	*{$attr} = subname($attr, sub {
		return @{ $_[0]->$meth };
	});
}

1;

#ABSTRACT: An SFF entry

__END__

=head1 SYNOPSIS

 while(my $entry = $reader->next_entry) {
     say '>', $entry->name;
     say $entry->bases;
 }

=head1 DESCRIPTION

This object represents an entry in an SFF file. It contains both processed data (the nucleotides and quality scores) and raw data (e.g. flowgram values).

=attr name

The name of this sequence

=attr bases

The nucleotides of this sequence

=attr flowgram_values

Returns an array containing all flowgram values. 

=attr flow_index_per_base

This contains the flow positions for each base in the called sequence (i.e., for each base, the position in the flowgram whose estimate resulted in that base being called). These values are "incremental" values, meaning that the stored position is the offset from the previous flow index in the field. All position values (prior to their incremental encoding) use 1-based indexing, so the first flow is flow 1.

=attr quality_scores

The quality scores for each of the bases in the sequence, where the values use the standard -log10 probability scale.

=attr clip_qual_left

The first base after the clipping point for quality, using 1-based indexing.

=attr clip_qual_right

The last base before the clipping point for quality, using 1-based indexing.

=attr clip_adaptor_left

The first base after the clipping point for quality, using 1-based indexing.

=attr clip_adaptor_right

The last base before the clipping point for quality, using 1-based indexing.

=begin Pod::Coverage

flow_index_per_base_raw
flowgram_values_raw
quality_scores_raw

=end Pod::Coverage
