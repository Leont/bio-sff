package Bio::SFF::Header;

use Sub::Name;
use Moo;

use Scalar::Util qw/looks_like_number/;

for my $attr (qw/magic version index_offset index_length number_of_reads header_length key_length number_of_flows flowgram_format_code/) {
	has $attr => (
		is => 'ro',
		required => 1,
		isa => sub {
			return looks_like_number($_[0]);
		},
	);
}

for my $attr(qw/flow_chars key_sequences/) {
	has $attr => (
		is => 'ro',
		required => 1,
		isa => sub {
			return defined and ref($_[0]) eq '';
		},
	);
}

1;
