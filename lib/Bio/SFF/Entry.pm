package Bio::SFF::Entry;

use Sub::Name;
use Moo;

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
