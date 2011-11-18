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

=head1 DESCRIPTION

=attr name

The name of this sequence

=attr bases

The nucleotides of this sequence

=attr clip_qual_left

The first base after the clipping point for quality, using 1-based indexing.

=attr clip_qual_right

The last base before the clipping point for quality, using 1-based indexing.

=attr clip_adaptor_left

The first base after the clipping point for quality, using 1-based indexing.

=attr clip_adaptor_right

The last base before the clipping point for quality, using 1-based indexing.

=attr flowgram_values

Returns an array containing all flowgram values. 

=attr flow_index_per_base

=attr quality_scores

