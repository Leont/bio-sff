package Bio::SFF::Reader;
use strict;
use warnings;

use Moo;

use Bio::SFF::Entry;
use Bio::SFF::Header;
use Carp qw/croak/;
use Config;
use Scalar::Util qw/reftype/;

sub _roundup {
	my $number = shift;
	my $remain = $number % 8;
	return $number + ($remain ? 8 - $remain : 0);
}

around BUILDARGS => sub {
	my ($orig, $class, @args) = @_;
	 
	unshift @args, 'file' if @args % 2 == 1 and ref($_[0]) ne 'HASH';

	return $class->$orig(@args);
};

has _fh => (
	is       => 'ro',
	required => 1,
	init_arg => 'file',
	isa      => sub {
		reftype($_[0]) eq 'GLOB';
	},
	coerce   => sub {
		my $val = shift;
		return $val if ref($val);
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

sub _build_header {
	my $self = shift;
	my $idx_off_type = $] >= 5.010 && $Config{use64bitint} ? 'Q>' : 'x[N]N';
	my $templ = "a4N $idx_off_type N2n3C";
	my ($magic, $version, $index_offset, $index_length, $number_of_reads, $header_length, $key_length, $number_of_flows_per_read, $flowgram_format_code) = unpack $templ, $self->_read_bytes(31);
	my ($flow_chars, $key_sequence) = unpack sprintf('A%dA%d', $number_of_flows_per_read, $key_length), $self->_read_bytes($header_length - 31);

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

has number_of_reads => (
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
my @header_keys = qw/number_of_bases clip_qual_left clip_qual_right clip_adaptor_left clip_adaptor_right name/;

sub read_entry {
	my $self = shift;
	return if $self->_current_read >= $self->number_of_reads;
	
	my %entry;
	@entry{qw/read_header_length name_length/} = unpack "nn", $self->_read_bytes(4);
	@entry{@header_keys} = unpack sprintf($read_template, $entry{name_length}), $self->_read_bytes($entry{read_header_length} - 4);

	my $data_template = sprintf 'A%dA%dA%dA%d', 2 * $self->_number_of_flows_per_read, ($entry{number_of_bases}) x 3;
	my $data_length = _roundup(2 * $self->_number_of_flows_per_read + 3 * $entry{number_of_bases});
	@entry{qw/flowgram_values flow_index_per_base bases quality_scores/} = unpack $data_template, $self->_read_bytes($data_length);

	$self->_current_read($self->_current_read + 1);
	return Bio::SFF::Entry->new(\%entry);
}

1;
