package Bio::SFF::Reader::Random;

use Moo;

use Bio::SFF::Index;
use Carp qw/croak/;
use Const::Fast;
use Fcntl qw/SEEK_SET/;

const my $index_header => 8;
const my $roche_offset => 5;
const my $base255 => 255;

has _index => (
	is => 'ro',
	init_arg => undef,
	builder => '_build_index',
	lazy => 1,
	predicate => '_has_index'
);

with 'Bio::SFF::Reader';

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
	my ($counter, $buffer, %offset_for) = (0, '');
	while ($counter < $self->_number_of_reads) {
		read $self->_fh, $buffer, 8192, length $buffer or croak "Couldn\'t read index($counter)";
		while ($buffer =~ m/ (.+?) \xFF /gcxs) {
			my $name = $1;
			my @offset = unpack 'C5', substr $name, -$roche_offset, $roche_offset, '';
			$offset_for{$name} = $offset[-1] + 255 * $offset[-2] + 255**2 * $offset[-3] + 255**3 * $offset[-4];
			$counter++;
		}
		$buffer = substr $buffer, pos $buffer;
	}
	$self->_fh->seek($tell, SEEK_SET);
	return Bio::SFF::Index->new(offsets => \%offset_for, manifest => $xml);
}

sub _read_slow_index {
	my $self = shift;
	my $position = $self->header->header_length;

	my $tell = $self->_fh->tell;
	$self->_fh->seek($position, SEEK_SET) or croak "Couldn't seek: $!";

	my %offset_for;
	for my $counter (1 .. $self->_number_of_reads) {
		my $offset = $self->_fh->tell;
		$offset_for{ $self->_read_entry->name } = $offset;
	}
	$self->_fh->seek($tell, SEEK_SET);
	return Bio::SFF::Index->new(offsets => \%offset_for, manifest => undef);
}

sub lookup {
	my ($self, $name) = @_;
	my $offset = $self->_index->offset_of($name);
	return if not defined $offset;
	$self->_fh->seek($offset, SEEK_SET);
	return $self->_read_entry;
}

1;

#ABSTRACT: Random-access SFF reader

__END__

=head1 SYNOPSIS

 my $reader = Bio::SFF::Reader::Random->new(file => $filename);
 if (my $entry = $reader->lookup('ABCDEF')) {
     say '>', $entry->name;
     say $entry->bases;
 }
 else {
     say 'No such entry "ABCDEF"';
 }

=head1 DESCRIPTION

This class implements L<Bio::SFF::Reader|Bio::SFF::Reader>. It provides random access to an SFF file.

=method lookup($name)

This will look up a named sequence, and return it. Note that this will affect the iterator.

