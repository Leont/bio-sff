package Bio::SFF::Reader;

use Moo;

use Bio::SFF::Entry;
use Bio::SFF::Header;
use Bio::SFF::Index;
use Carp qw/croak/;
use Config;
use Const::Fast;
use Fcntl qw/SEEK_SET SEEK_CUR/;
use Scalar::Util qw/reftype/;

const my $padding_to => 8;
const my $index_header => 8;
const my $roche_offset => 5;
const my $base255 => 255;
const my $header_size => 31;
const my $entry_header_size => 4;
const my $idx_off_type => ($] >= 5.010 && $Config{use64bitint} ? 'Q>' : 'x[N]N');
const my $size_of_flowgram_value => 2;
const my $uses_number_of_bases => 3;

sub _roundup {
	my $number = shift;
	my $remain = $number % $padding_to;
	return $number + ($remain ? $padding_to - $remain : 0);
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
		return $val if ref $val and reftype($val) eq 'GLOB';
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
		number_of_flows_per_read => $number_of_flows_per_read,
		flowgram_format_code => $flowgram_format_code,
		flow_chars => $flow_chars,
		key_sequences => $key_sequence,
	);

	return $header;
}

for my $method (qw/number_of_reads number_of_flows_per_read index_offset index_length/) {
	has "_$method" => (
		is => 'ro',
		init_arg => undef,
		default => sub {
			my $self = shift;
			return $self->header->$method;
		},
		lazy => 1,
	);
}

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
	$self->_fh->seek($self->_index_length, SEEK_CUR) if $self->_fh->tell == $self->_index_offset;
	return if $self->_fh->eof;
	
	my %entry;
	@entry{qw/read_header_length name_length/} = unpack 'nn', $self->_read_bytes($entry_header_size);
	(my ($number_of_bases), @entry{@header_keys}) = unpack sprintf($read_template, $entry{name_length}), $self->_read_bytes($entry{read_header_length} - $entry_header_size);

	my $data_template = sprintf 'A%dA%dA%dA%d', $size_of_flowgram_value * $self->_number_of_flows_per_read, ($number_of_bases) x $uses_number_of_bases;
	my $data_length = _roundup($size_of_flowgram_value * $self->_number_of_flows_per_read + $uses_number_of_bases * $number_of_bases);
	@entry{qw/flowgram_values flow_index_per_base bases quality_scores/} = unpack $data_template, $self->_read_bytes($data_length);

	return Bio::SFF::Entry->new(\%entry);
}

has _index_info => (
	is => 'ro',
	init_arg => undef,
	builder => '_build_index_info',
	lazy => 1,
);

sub _build_index_info {
	my $self = shift;
	my ($index_offset, $index_length) = ($self->header->index_offset, $self->header->index_length);
	return if not $index_offset or not $index_length;
	
	my $tell = $self->_fh->tell;
	$self->_fh->seek($index_offset, SEEK_SET);
	my ($magic_number) = unpack 'A8', $self->_read_bytes($index_header);
	$self->_fh->seek($tell, SEEK_SET);
	return $magic_number;
}

has manifest => (
	is => 'ro',
	init_arg => undef,
	builder => '_build_manifest',
	lazy => 1,
);

sub _build_manifest {
	my $self = shift;
	return $self->index->manifest if $self->_has_index;
	my $magic_number = $self->_index_info;
	if ($magic_number eq '.mft1.00') { 
		my ($index_offset, $index_length) = ($self->_index_offset, $self->_index_length);
		my $tell = $self->_fh->tell;
		$self->_fh->seek($index_offset + $index_header, SEEK_SET);
		my $xml = $self->_read_manifest($magic_number);
		$self->_fh->seek($tell, SEEK_SET);
		return $xml;
	}
	return;
}

sub _read_manifest {
	my ($self, $magic_number) = @_;
	my $xmldata_head = $self->_read_bytes($index_header);
	if ( $magic_number eq '.mft1.00') {
		my ($xml_size, $data_size) = unpack 'NN', $xmldata_head;
		return $self->_read_bytes($xml_size);
	}
	return;
}

has index => (
	is => 'ro',
	init_arg => undef,
	builder => '_build_index',
	lazy => 1,
	predicate => '_has_index'
);

sub _build_index {
	my $self = shift;
	my $magic_number = $self->_index_info;
	my $has_roche_index = defined $magic_number and $magic_number =~ / \A \.[sm]ft 1\.00 \z /xm;
	return $has_roche_index ? $self->_read_roche_index($magic_number) : $self->_read_slow_index;
}

sub _read_roche_index {
	my ($self, $magic_number) = @_;

	my ($index_offset, $index_length) = ($self->header->index_offset, $self->header->index_length);
	my $tell = $self->_fh->tell;
	$self->_fh->seek($index_offset + $index_header, SEEK_SET);

	my $xml = $self->_read_manifest($magic_number);
	my ($buffer, %offset_for) = ('');
	my $counter = 0;
	while ($counter < $self->_number_of_reads) {
		read $self->_fh, $buffer, 8192, length $buffer or croak "Couldn\'t read index($counter)";
		while ($buffer =~ m/ (.+?) \xFF /gcxs) {
			my $name = $1;
			my @offset = unpack 'C5', substr $name, -$roche_offset, $roche_offset, '';
			$offset_for{$name} = $offset[-1] + 255 * $offset[-2] + 65025 * $offset[-3] + 16581375 * $offset[-4];
			$counter++;
		}
		$buffer = substr $buffer, pos $buffer;
	}
	$self->_fh->seek($tell, SEEK_SET);
	return Bio::SFF::Index->new(offsets => \%offset_for, manifest => $xml);
}

sub _read_slow_index {
	my $self = shift;

	my $tell = $self->_fh->tell;
	$self->reset;

	my %offset_for;
	for my $counter (1 .. $self->_number_of_reads) {
		my $offset = $self->_fh->tell;
		my $entry = $self->next_entry;
		$offset_for{ $entry->name } = $offset;
	}
	$self->_fh->seek($tell, SEEK_SET);
	return Bio::SFF::Index->new(offsets => \%offset_for, manifest => undef);
}

sub lookup {
	my ($self, $name) = @_;
	my $offset = $self->index->offset_of($name);
	return if not defined $offset;
	warn "offset is $offset\n";
	$self->_fh->seek($offset, SEEK_SET);
	return $self->next_entry;
}

sub reset {
	my $self = shift;
	$self->header;
	$self->_fh->seek($self->header->header_length, SEEK_SET) or croak "Couldn't seek: $!";
	return;
}

1;

#ABSTRACT: An SFF reader

__END__

=head1 SYNOPSIS

 my $reader = Bio::SFF::Reader->new(file => $filename);
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

Read an entry and return it as a L<Bio::SFF:Entry|Bio::SFF::Entry> object.

=method lookup($name)

This will look up a named sequence, and return it. Note that this will affect the iterator.

=method reset()

Reset the iterator to the start of the file.

=method manifest

This returns the (XML) manifest as a string or undef if none is present.

=method header()

Returns the L<Bio::SFF::Header|Bio::SFF::Header> object associated with this reader.

=cut
