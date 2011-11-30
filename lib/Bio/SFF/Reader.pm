package Bio::SFF::Reader;

use Moo;

use Bio::SFF::Entry;
use Bio::SFF::Header;
use Carp qw/croak/;
use Config;
use Const::Fast;
use Scalar::Util qw/reftype/;

const my $bits64 => 8;
const my $header_size => 31;
const my $entry_header_size => 4;
const my $idx_off_type => ($] >= 5.010 && $Config{use64bitint} ? 'Q>' : 'x[N]N');
const my $size_of_flowgram_value => 2;
const my $uses_number_of_bases => 3;

sub _roundup {
	my $number = shift;
	my $remain = $number % $bits64;
	return $number + ($remain ? $bits64 - $remain : 0);
}

has _fh => (
	is       => 'ro',
	required => 1,
	init_arg => 'file',
	isa      => sub {
		reftype($_[0]) eq 'GLOB';
	},
	coerce   => sub {
		my $val = shift;
		return $val if ref($val) and reftype($val) eq 'GLOB';
		open my $fh, '<:raw', $val or croak "Could open file $val: $!";
		return $fh;
	}
);

has header => (
	is => 'ro',
	init_arg => undef,
	builder => '_build_header',
	lazy => 1,
);

## no critic (Subroutines::ProhibitUnusedPrivateSubroutines,Subroutines::ProhibitBuiltinHomonyms)
sub _build_header {
	my $self = shift;
	my $templ = "a4N $idx_off_type N2n3C";
	my ($magic, $version, $index_offset, $index_length, $number_of_reads, $header_length, $key_length, $number_of_flows_per_read, $flowgram_format_code) = unpack $templ, $self->_read_bytes($header_size);
	my ($flow_chars, $key_sequence) = unpack sprintf('A%dA%d', $number_of_flows_per_read, $key_length), $self->_read_bytes($header_length - $header_size);

	my $header = Bio::SFF::Header->new(
		magic => $magic,
		version => $version,
		index_offset => $index_offset,
		index_length => $index_length,
		number_of_reads => $number_of_reads,
		header_length => $header_length,
		key_length => $key_length,
		number_of_flows => $number_of_flows_per_read,
		flowgram_format_code => $flowgram_format_code,
		flow_chars => $flow_chars,
		key_sequences => $key_sequence,
	);

	return $header;
}

has _number_of_reads => (
	is => 'ro',
	init_arg => undef,
	default => sub {
		my $self = shift;
		return $self->header->number_of_reads;
	},
	lazy => 1,
);

has _current_read => (
	is => 'rw',
	init_arg => undef,
	default => sub { 0 },
);

has _number_of_flows_per_read => (
	is => 'ro',
	lazy => 1,
	default => sub {
		my $self = shift;
		return $self->header->number_of_flows;
	},
);

sub _read_bytes {
	my ($self, $num) = @_;
	my $buffer;
	croak "Could not read SFF file: $!" if not defined read $self->_fh, $buffer, $num;
	return $buffer;
}

my $read_template = 'Nnnnn A%d';
my @header_keys = qw/clip_qual_left clip_qual_right clip_adaptor_left clip_adaptor_right name/;

sub next_entry {
	my $self = shift;
	return if $self->_current_read >= $self->_number_of_reads;
	
	my %entry;
	@entry{qw/read_header_length name_length/} = unpack 'nn', $self->_read_bytes($entry_header_size);
	(my ($number_of_bases), @entry{@header_keys}) = unpack sprintf($read_template, $entry{name_length}), $self->_read_bytes($entry{read_header_length} - $entry_header_size);

	my $data_template = sprintf 'A%dA%dA%dA%d', $size_of_flowgram_value * $self->_number_of_flows_per_read, ($number_of_bases) x $uses_number_of_bases;
	my $data_length = _roundup($size_of_flowgram_value * $self->_number_of_flows_per_read + $uses_number_of_bases * $number_of_bases);
	@entry{qw/flowgram_values flow_index_per_base bases quality_scores/} = unpack $data_template, $self->_read_bytes($data_length);

	$self->_current_read($self->_current_read + 1);
	return Bio::SFF::Entry->new(\%entry);
}

sub reset {
	my $self = shift;
	$self->header;
	seek $self->_fh, $self->file->header_length, 0 or croak "Couldn't seek: $!";
	return;
}

1;

#ABSTRACT: An SFF reader

__END__

=head1 SYNOPSIS

 my $reader = Bio::SFF::Reader(file => $filename);
 while (my $entry = $reader->next_entry) {
     say '>', $entry->name;
     say $entry->bases;
 }

=head1 DESCRIPTION

=method new(...)

This method creates a new SFF object

=over 4

=item * file

The file that should be read. This can either be a filename or a filehandle.

=back

=method next_entry()

Read an entry and return it as a Bio::SFF:Entry object.

=method reset()

Reset the iterator to the start of the file.

=method header()

Returns the L<Bio::SFF::Header|Bio::SFF::Header> object associated with this reader.

=cut
