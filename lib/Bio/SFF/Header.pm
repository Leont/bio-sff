package Bio::SFF::Header;

use Moo;
use Sub::Name;

use Scalar::Util qw/looks_like_number/;

for my $attr (qw/magic version index_offset index_length number_of_reads header_length number_of_flows_per_read flowgram_format_code/) {
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

#ABSTRACT: An SFF header

__END__

=head1 DESCRIPTION

This object represents the header of an SFF file. You probably don't want to deal with this in any way.

=attr magic

The magic bytes at the start of any SFF file, this is always C<779314790>.

=attr version

The version of SFF that is used. This must currently be 1.

=attr index_offset

The offset of the index, or C<0> if no index is present.

=attr index_length

The length of the index, or C<0> if no index is present.

=attr number_of_reads

The number of reads in the SFF file.

=attr header_length

The length of the header in bytes.

=attr number_of_flows

The number of flows in the entries.

=attr flowgram_format_code

Currently, this must always be C<1>.

=attr flow_chars

The array of nucleotide bases ('A', 'C', 'G' or 'T') that correspond to the nucleotides used for each flow of each read.

=attr key_sequences

The nucleotide bases of the key sequence used for the reads.

=attr number_of_flows_per_read

The number of flowgram values in each entry.

=cut
