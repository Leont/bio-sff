package Bio::SFF::Reader::Sequential;

use Moo;

with 'Bio::SFF::Reader';

use Carp qw/croak/;
use Const::Fast;
use Fcntl qw/SEEK_SET SEEK_CUR/;

has _current_read => (
	is => 'rw',
	init_arg => undef,
	default => sub { 0 },
);

sub next_entry {
	my $self = shift;
	return if $self->_current_read >= $self->_number_of_reads;
	$self->_fh->seek($self->_index_length, SEEK_CUR), if $self->_fh->tell == $self->_index_offset;
	my $ret = $self->_read_entry;
	$self->_current_read($self->_current_read + 1) if defined $ret;
	return $ret;
}

sub _has_index {
	return 0;
}

sub reset {
	my $self = shift;
	$self->header;
	$self->_fh->seek($self->header->header_length, SEEK_SET) or croak "Couldn't seek: $!";
	$self->_current_read(0);
	return;
};

1;

#ABSTRACT: Sequential SFF reader

__END__

=head1 SYNOPSIS

 my $reader = Bio::SFF::Reader::Sequential->new(file => $filename);
 while (my $entry = $reader->next_entry) {
     say '>', $entry->name;
     say $entry->bases;
 }

=head1 DESCRIPTION

This class implements L<Bio::SFF::Reader|Bio::SFF::Reader>. It provides sequential access to an SFF file.

=method next_entry()

Read an entry and return it as a L<Bio::SFF:Entry|Bio::SFF::Entry> object.

=method reset()

Reset the iterator to the start of the file.

