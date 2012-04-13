package Bio::SFF::Index;

use Moo;

has manifest => (
	is => 'ro',
	required => 1,
);

has _offsets => (
	is => 'ro',
	isa => sub { ref($_[0]) eq 'HASH' },
	init_arg => 'offsets',
	required => 1,
);

sub offset_of {
	my ($self, $name) = @_;
	return $self->_offsets->{$name};
}

1;

#ABSTRACT: SFF index object

__END__

=head1 DESCRIPTION

This class represents the index of an SFF file.

=method manifest()

This returns the (XML) manifest as a bytestring.

=method offset_of($name)

This returns the offset of a specific entry in the SFF file.

=cut
